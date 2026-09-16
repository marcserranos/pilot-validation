#!/usr/bin/env python3
"""Fixture tests for 24_novelty_by_field.py (SCHEMA.md hard rule 6).

Unit tests cover:
- field-class mapping at every depth, including an undetermined consensus
- per-call artifact labels
- translation, frameshift and premature-stop detection
- cds_known / protein_known detection against a small fake refdata directory, with mixed header
  styles
- protein-difference strings (Hamming and Levenshtein)
- groove annotation
- greedy unrelated selection
- strict ancestry
- count suppression

The end-to-end test builds a small synthetic outroot (Table 1, cds.fa.gz, refdata,
alleles.csv.gz, cohort_membership, relatedness). It then runs main() and checks the output files
and the key counts.

Run (no pytest needed):
    cd scripts/hla_popgen && python3 tests/test_novelty_by_field.py
"""
import gzip
import importlib.util
import inspect
import json
import os
import pathlib
import random
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m24 = _load("24_novelty_by_field.py", "hla_popgen_novelty_by_field")

# ---------------------------------------------------------------------------
# Synthetic sequences
# ---------------------------------------------------------------------------
_RNG = random.Random(24)
SENSE_CODONS = sorted(c for c, a in m24.CODON_TABLE.items() if a != "*")


def rand_cds(n_codons, rng=_RNG):
    return "ATG" + "".join(rng.choice(SENSE_CODONS) for _ in range(n_codons - 1)) + "TAA"


def mutate_codon(seq, p, synonymous):
    """Replace codon p (1-based) with a synonymous or missense (never stop) codon."""
    i = 3 * (p - 1)
    old = seq[i:i + 3]
    aa = m24.CODON_TABLE[old]
    for c in SENSE_CODONS:
        if c == old:
            continue
        if (m24.CODON_TABLE[c] == aa) == synonymous:
            return seq[:i] + c + seq[i + 3:]
    raise ValueError(f"no {'synonymous' if synonymous else 'missense'} codon for {old}")


def _first_codon_with_synonym(seq, start):
    p = start
    while True:
        old = seq[3 * (p - 1):3 * p]
        aa = m24.CODON_TABLE[old]
        if any(m24.CODON_TABLE[c] == aa and c != old for c in SENSE_CODONS):
            return p
        p += 1


A1 = rand_cds(365)  # 1098 nt including stop = canonical HLA-A length
A0101_02 = mutate_codon(A1, _first_codon_with_synonym(A1, 150), synonymous=True)
A0201 = mutate_codon(A1, 40, synonymous=False)
DRB1_1 = rand_cds(21)
X_SEQ = mutate_codon(A1, 250, synonymous=False)   # exon 4, not groove
Y_SEQ = mutate_codon(A1, 30, synonymous=False)    # exon 2, groove
Z_SEQ = mutate_codon(A1, 260, synonymous=False)   # named depth 3 but protein is new
S_SEQ = mutate_codon(A1, _first_codon_with_synonym(A1, 100), synonymous=True)
PK_SEQ = mutate_codon(A1, _first_codon_with_synonym(A1, 120), synonymous=True)
FS_SEQ = A1[:500] + A1[501:]
STOP_SEQ = A1[:3 * 49] + "TAA" + A1[3 * 50:]
DQB1_SEQ = rand_cds(30)


def write_fasta_gz(path, records):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for k in range(0, len(seq), 60):
                f.write(seq[k:k + 60] + "\n")


def make_refdata(root):
    ref = os.path.join(root, "refdata", "Data-2024Feb02")
    # Mixed header styles: IMGT-like, HLA-prefixed, bare with trailing tokens, bare fields only.
    # Real header style (verified on the VM): '>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp',
    # plus tolerated variants (no HLA- prefix, bare fields, extra tokens, frame!=1).
    write_fasta_gz(os.path.join(ref, "CDSseq", "HLA-A.fa.gz"), [
        ("HLA-A*01:01:01:01 HLA00001 frame=1 1098bp", A1),
        ("HLA-A*01:01:01:02 HLA00002 frame=1 1098bp", A1),
        ("A*01:01:02 extra tokens", A0101_02),
        ("HLA-A*02:01:01:01 HLA00005 frame=1 1095bp", A0201[:-3]),  # no stop: must still match
        ("HLA-A*03:07:01:01 HLA00099 frame=2 546bp", "C" + A1[:545]),  # partial, not translatable
    ])
    write_fasta_gz(os.path.join(ref, "CDSseq", "DRB1.fa.gz"), [("01:01:01 frame=1", DRB1_1)])
    cds = "CDS=101..173,301..570,701..976,1101..1376,1501..1703"
    exon = "exon=1:1..173,2:301..570,3:701..976,4:1101..1376,5:1501..2000"
    with gzip.open(os.path.join(ref, "alleles.csv.gz"), "wt") as f:
        for allele, good in (("HLA-A*01:01:01:01", "true"), ("HLA-A*01:01:01:02", "true"),
                             ("HLA-A*02:01:01:01", "false")):
            f.write("\t".join(["gene=HLA-A", "version=IPD-IMGT/HLA-V3.55.0", "ID=HLA00001",
                               f"allele={allele}", "type=gene", "geneRange=1..3000",
                               f"goodTemplate={good}", "startCodon=ATG:101..103", cds,
                               "UTR=1..100,1704..3000", "nexon=5", exon]) + "\n")
    return os.path.dirname(ref)


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------
def test_field_class_mapping():
    fc = m24.field_class_of
    assert fc("HLA-A*01:01:01:01") == "known"
    assert fc("HLA-A*01:01:01:new") == "f4_noncoding"
    assert fc("HLA-A*01:01:new") == "f3_synonymous"
    assert fc("HLA-A*01:new") == "f2_protein"
    assert fc("HLA-A*new") == "f1_undetermined"
    assert fc("undetermined") == "uncalled"
    assert fc(None) == "uncalled"
    assert fc(float("nan")) == "uncalled"
    assert fc("") == "uncalled"
    # explicit novelty_depth wins over the string
    assert fc("HLA-A*01:new", 3) == "f3_synonymous"
    assert fc("HLA-A*01:01:new", float("nan")) == "f3_synonymous"
    # vectorized path agrees
    df = pd.DataFrame({
        "consensus": ["HLA-A*01:01:01:01", "HLA-A*01:01:01:new", "HLA-A*01:01:new",
                      "HLA-A*01:new", "HLA-A*new", "undetermined", None],
        "novelty_depth": [None, 4, None, 2, 1, None, None],
        "template_warning": ["NA"] * 7, "cds_mut": [None] * 7,
    })
    out = m24.add_call_labels(df)
    assert list(out["field_class"]) == ["known", "f4_noncoding", "f3_synonymous", "f2_protein",
                                        "f1_undetermined", "uncalled", "uncalled"]


def test_artifact_labels():
    al = m24.artifact_label_of
    hp = "HLA-A*01:01|:100-a:200|K(aaa)<N(aa)"
    assert al("NA", hp) == "homopolymer_indel"
    assert al("partial_CDS", hp) == "homopolymer_indel"  # rule B checked first
    assert al("partial_CDS", None) == "partial_cds"
    assert al("no-start_codon,inframe_stop", None) == "inframe_stop"
    assert al("NA", "HLA-A*01:01|:10*ag:5-a:4|x") == "clean"
    assert al(None, None) == "clean"
    assert al("no-stop_codon", None) == "clean"
    df = pd.DataFrame({"consensus": ["HLA-A*01:new"] * 4, "novelty_depth": [2] * 4,
                       "template_warning": ["NA", "partial_CDS", float("nan"), "inframe_stop"],
                       "cds_mut": [hp, None, None, None]})
    assert list(m24.add_call_labels(df)["artifact_label"]) == [
        "homopolymer_indel", "partial_cds", "clean", "inframe_stop"]


def test_translation_frameshift_stop():
    assert m24.translate("ATGAAATAA") == "MK*"
    pi = m24.protein_info("ATGAAATAA")
    assert pi == {"protein": "MK", "frameshift": False, "premature_stop": False}
    pi = m24.protein_info("ATGTAAAAATAA")
    assert pi["premature_stop"] and not pi["frameshift"]
    pi = m24.protein_info("ATGAAAATAA")
    assert pi["frameshift"]
    assert m24.cds_core("ATGAAATAG") == "ATGAAA"
    assert m24.cds_core("ATGAAA") == "ATGAAA"
    assert m24.cds_core("ATGAATAG") == "ATGAATAG"  # out of frame: untouched
    assert m24.translate("ATGNNN") == "MX"


def test_refdata_known_detection(tmp_path):
    refdir = make_refdata(str(tmp_path))
    ref = m24.load_refdata(refdir)
    assert ref.n_records == 6 and ref.n_unparsed_headers == 0
    assert ref.n_nonstandard_frame == 1  # frame=2 record: CDS index only, no protein
    assert "HLA-A*03:07:01:01" in ref.name_cds and "HLA-A*03:07:01:01" not in ref.name_prot
    assert m24.header_frame(">HLA-A*01:01 HLA00001 frame=2 546bp") == 2
    assert m24.header_frame(">HLA-A*01:01") == 1
    assert set(ref.genes()) == {"A", "DRB1"}
    assert "HLA-A*01:01:02" in ref.name_cds
    assert "HLA-DRB1*01:01:01" in ref.name_cds
    c = m24.classify_sequence
    r = c(A0201, "A", 2, ref)  # observed carries a stop, refdata does not
    assert r["seq_class"] == "cds_known" and r["mapped_prot2"] == "HLA-A*02:01"
    assert r["mapped_cds3"] == "HLA-A*02:01:01"
    r = c(PK_SEQ, "A", 2, ref)
    assert r["seq_class"] == "protein_known" and r["protein_known"] and not r["cds_known"]
    assert r["mapped_prot2"] == "HLA-A*01:01"
    assert c(PK_SEQ, "A", 3, ref)["seq_class"] == "novel_cds_synonymous"
    assert c(X_SEQ, "A", 2, ref)["seq_class"] == "novel_protein"
    assert c(Z_SEQ, "A", 3, ref)["seq_class"] == "novel_protein"
    assert c(FS_SEQ, "A", 2, ref)["seq_class"] == "frameshift_or_stop"
    assert c(STOP_SEQ, "A", 2, ref)["seq_class"] == "frameshift_or_stop"
    assert c(DQB1_SEQ, "DQB1", 2, ref)["seq_class"] == "no_refdata"
    # header parser tolerance
    assert m24.parse_ref_header(">HLA:HLA00001 A*01:01 1098 bp", "X") == ("A", "HLA-A*01:01")
    assert m24.parse_ref_header(">01:02", "DQB1") == ("DQB1", "HLA-DQB1*01:02")
    assert m24.parse_ref_header(">MICA*001", "MICA") == ("MICA", "MICA*001")


def test_protein_diff_positions(tmp_path):
    assert m24.format_diffs(m24.hamming_diffs("MKVLA", "MRVLG")) == "2:K>R,5:A>G"
    d, diffs = m24.levenshtein_align("MKVLA", "MKLA")
    assert d == 1 and m24.format_diffs(diffs) == "3:V>-"
    d, diffs = m24.levenshtein_align("MKLA", "MKWLA")
    assert d == 1 and m24.format_diffs(diffs) == "2:->W"
    ref = m24.load_refdata(make_refdata(str(tmp_path)))
    xp = m24.protein_info(X_SEQ)["protein"]
    nn = m24.nearest_known_protein(xp, "A", ref)
    ap = m24.protein_info(A1)["protein"]
    assert nn["metric"] == "hamming" and nn["distance"] == 1
    assert nn["name"] == "HLA-A*01:01:01:01"
    assert m24.format_diffs(nn["diffs"]) == f"250:{ap[249]}>{xp[249]}"
    # in-frame deletion -> Levenshtein path, with the hint used as a candidate
    dp = ap[:10] + ap[11:]
    nn = m24.nearest_known_protein(dp, "A", ref, ["HLA-A*01:01:01:01"])
    assert nn["metric"] == "levenshtein" and nn["distance"] == 1
    assert nn["diffs"][0][2] == "-"


def test_groove_exon_sets():
    ge = m24.groove_exons_for
    assert ge("A") == {2, 3} and ge("E") == {2, 3} and ge("G") == {2, 3}
    assert ge("DRB1") == {2} and ge("DQB1") == {2} and ge("DPA1") == {2} and ge("DRA") == {2}
    assert ge("DRB5") == {2}
    assert ge("MICA") is None and ge("TAP1") is None and ge("HFE") is None


def test_cds_segments_and_exon_maps(tmp_path):
    segs = m24.cds_segments("101..173,301..570,701..976",
                            "1:1..173,2:301..570,3:701..976")
    assert segs == ((1, 73), (2, 270), (3, 276))
    # no exon= field: the k-th CDS range is exon k
    assert m24.cds_segments("1..73,100..369", "") == ((1, 73), (2, 270))
    refdir = make_refdata(str(tmp_path))
    em = m24.load_exon_maps(os.path.join(refdir, "Data-2024Feb02", "alleles.csv.gz"),
                            {"HLA-A*01:01:01:01"}, {"A"})
    assert em.version == "IPD-IMGT/HLA-V3.55.0"
    assert em.by_name["HLA-A*01:01:01:01"][:3] == ((1, 73), (2, 270), (3, 276))
    assert em.gene_mode["A"] == em.by_name["HLA-A*01:01:01:01"]
    assert em.genes == {"A"}


def test_groove_annotation():
    diffs = [(30, "K", "R"), (250, "A", "G"), (10, "M", "V")]
    empty = m24.ExonMaps()
    # canonical fallback only when the gene is absent from alleles.csv.gz entirely
    flags, method, _ = m24.annotate_groove(diffs, "A", "HLA-A*x", 1095, 1095, empty)
    assert method == "approx_canonical" and flags == [True, False, False]
    flags, method, reason = m24.annotate_groove(diffs, "A", "HLA-A*x", 1095, 1092, empty)
    assert method == "NA" and flags == [None] * 3 and "canonical" in reason
    assert m24.annotate_groove(diffs, "DRB1", "x", 798, 798, empty)[1] == "NA"
    em = m24.ExonMaps()
    em.genes = {"A", "DRB1"}
    # exon 1 = 87 nt, so codon 30's middle base (index 88) is the first codon inside exon 2
    em.by_name["HLA-DRB1*y"] = ((1, 87), (2, 270), (3, 441))
    em.gene_mode["DRB1"] = ((1, 87), (2, 270), (3, 441))
    flags, method, _ = m24.annotate_groove(diffs, "DRB1", "HLA-DRB1*y", 798, 798, em)
    assert method == "refdata_exon" and flags == [True, False, False]
    flags, method, _ = m24.annotate_groove(diffs, "DRB1", "HLA-DRB1*other", 798, 798, em)
    assert method == "refdata_gene_mode" and flags == [True, False, False]
    # gene present in alleles.csv.gz but no usable map -> NA, never the canonical fallback
    assert m24.annotate_groove(diffs, "A", "HLA-A*x", 1095, 1095, em)[1] == "NA"
    assert m24.annotate_groove(diffs, "MICA", "x", 1, 1, empty)[1] == "NA"


def test_greedy_unrelated():
    ids = ["a", "b", "c", "d", "e", "f"]
    pairs = [("a", "b", 0.25), ("a", "c", 0.25), ("a", "d", 0.1),  # a has 3 relatives
             ("e", "f", 0.05),                                    # tie -> drop 'e'
             ("b", "c", 0.0441),                                  # below threshold
             ("c", "zz", 0.5)]                                    # zz not long-read
    kept, removed = m24.greedy_unrelated(ids, pairs)
    assert removed == ["a", "e"], removed
    assert kept == {"b", "c", "d", "f"}
    # deterministic tie-break among equal degrees: sorted id
    kept, removed = m24.greedy_unrelated(["2", "10"], [("2", "10", 0.3)])
    assert removed == ["10"] and kept == {"2"}
    # triangle: first removal leaves one pair; next tie by id
    kept, removed = m24.greedy_unrelated(["x", "y", "z"],
                                         [("x", "y", 0.2), ("y", "z", 0.2), ("x", "z", 0.2)])
    assert removed == ["x", "y"] and kept == {"z"}
    # exactly the threshold counts as related
    kept, _ = m24.greedy_unrelated(["p", "q"], [("p", "q", 0.0442)])
    assert kept == {"q"}


def test_strict_ancestry():
    people = pd.DataFrame({
        "ancestry_pred": ["AFR", "EUR", "EAS", None, "SAS"],
        "p_afr": [0.95, 0.1, 0.0, 0.5, 0.0], "p_amr": [0.0] * 5, "p_eas": [0.0, 0.0, 0.9, 0, 0],
        "p_eur": [0.05, 0.89, 0.1, 0.5, 0.0], "p_mid": [0.0] * 5,
        "p_sas": [0.0, 0.0, 0.0, 0.0, float("nan")],
    })
    assert list(m24.strict_ancestry(people)) == ["AFR", None, "EAS", None, None]


def test_suppression():
    s = m24.suppress
    assert s(0) == "0" and s(1) == "<20" and s(19) == "<20" and s(20) == "20"
    assert s(12345) == "12345" and s(None) == "NA" and s(float("nan")) == "NA"
    df = pd.DataFrame({"n": [0, 5, 25], "frac": [0.1, 0.2, 0.3]})
    out = m24.suppress_df(df, ["n"])
    assert list(out["n"]) == ["0", "<20", "25"] and list(out["frac"]) == [0.1, 0.2, 0.3]
    tree = m24.suppress_tree({"a": 3, "b": [1, 40], "limit": 5, "f": 0.5, "t": True})
    assert tree == {"a": "<20", "b": ["<20", "40"], "limit": 5, "f": 0.5, "t": True}


# ---------------------------------------------------------------------------
# End-to-end fixture
# ---------------------------------------------------------------------------
T1_COLS = ["person_id", "hap", "contig", "gene", "copy_index", "gene_class", "consensus",
           "n_fields", "is_novel", "novelty_depth", "novelty_class", "template_allele",
           "template_distance", "cds_distance", "cds_mut", "n_aa_changes", "template_warning",
           "alleles", "n_tied", "strand"]
NCLASS = {1: "undetermined", 2: "protein_altering", 3: "synonymous", 4: "beyond_cds"}


def build_outroot(root):
    """Returns dict of paths. 12 people x 2 haps. HLA-A carries the scenarios; DRB1 is known."""
    pid = lambda i: f"PID{i:04d}"  # noqa: E731
    known_a = ("HLA-A*01:01:01:01", A1, "NA", None)
    # (person, hap) -> HLA-A (consensus, seq, warning, cds_mut)
    scen = {
        (1, "hap1"): ("HLA-A*01:new", X_SEQ, "NA", "HLA-A*01:01:01:01|:748*ga:349|A<G"),
        (2, "hap1"): ("HLA-A*01:new", X_SEQ, "NA", "HLA-A*01:01:01:01|:748*ga:349|A<G"),
        (3, "hap1"): ("HLA-A*01:new", Y_SEQ, "NA", "HLA-A*01:01:01:01|:88*ga:1009|K<R"),
        (4, "hap2"): ("HLA-A*01:new", Y_SEQ, "NA", "HLA-A*01:01:01:01|:88*ga:1009|K<R"),
        (5, "hap1"): ("HLA-A*01:01:new", S_SEQ, "NA", "HLA-A*01:01:01:01|:300*ga:797|L<L"),
        (6, "hap1"): ("HLA-A*02:new", A0201, "NA", "HLA-A*02:01:01:01|:10*ga:1087|x"),
        (7, "hap1"): ("HLA-A*01:new", PK_SEQ, "NA", "HLA-A*01:01:01:01|:360*ga:737|x"),
        (8, "hap1"): ("HLA-A*01:new", FS_SEQ, "NA", "HLA-A*01:01:01:01|:500-a:597|x"),
        (9, "hap1"): ("HLA-A*01:new", STOP_SEQ, "partial_CDS",
                      "HLA-A*01:01:01:01|:147*at:948|x"),
        (10, "hap1"): ("HLA-A*01:01:new", Z_SEQ, "NA", "HLA-A*01:01:01:01|:778*ga:319|x"),
        (11, "hap2"): ("HLA-A*new", A1, "NA", None),
        (12, "hap1"): ("HLA-A*01:01:01:new", A1, "NA", None),
        (12, "hap2"): ("undetermined", A1, "NA", None),
    }
    ancestry = {1: ("afr", "AFR"), 2: ("eur", "EUR"), 3: ("afr", "AFR"), 4: ("afr", "AFR"),
                5: ("eas", "EAS"), 6: ("eur", "EUR"), 7: ("eur", "EUR"), 8: ("amr", "AMR"),
                9: ("sas", "SAS"), 10: ("mid", "MID"), 11: ("eur", "EUR"), 12: ("eur", "EUR")}
    people_root = os.path.join(root, "people")
    rows = []
    for i in range(1, 13):
        for hap in ("hap1", "hap2"):
            contig = f"{hap}_ctg_{pid(i)}_a"
            cons, seq, warn, cds_mut = scen.get((i, hap), known_a)
            entries = [("HLA-A", 1, cons, seq, warn, cds_mut),
                       ("HLA-DRB1", 1, "HLA-DRB1*01:01:01", DRB1_1, "NA", None)]
            if (i, hap) == (11, "hap1"):
                # copy number 2 on one contig, both novel -> ambiguous join, excluded
                entries = [("HLA-A", 1, *known_a),
                           ("HLA-DRB1", 1, "HLA-DRB1*01:new", DRB1_1, "NA", "x|:3*ga:59|x"),
                           ("HLA-DRB1", 2, "HLA-DRB1*01:new", DRB1_1, "NA", "x|:3*ga:59|x")]
            if (i, hap) == (5, "hap2"):
                entries.append(("HLA-DQB1", 1, "HLA-DQB1*05:new", DQB1_SEQ, "NA", "x|:3*ga:9|x"))
            recs = []
            for k, (gene, ci, c, s, w, cm) in enumerate(entries, 1):
                depth = m24.depth_from_consensus(c)
                gclass = "classical_I" if gene == "HLA-A" else "classical_II"
                rows.append({
                    "person_id": pid(i), "hap": hap, "contig": contig, "gene": gene,
                    "copy_index": ci, "gene_class": gclass, "consensus": c,
                    "n_fields": None, "is_novel": depth is not None,
                    "novelty_depth": depth, "novelty_class": NCLASS.get(depth),
                    "template_allele": "HLA-A*01:01:01:01" if gene == "HLA-A" else c,
                    "template_distance": 1 if depth else 0, "cds_distance": 1 if cm else None,
                    "cds_mut": cm, "n_aa_changes": None, "template_warning": w,
                    "alleles": "HLA-A*01:01:01:01" if gene == "HLA-A" else None,
                    "n_tied": 1, "strand": "+"})
                recs.append((f"{contig}_{gene}_{k} GT_AG:GT_AG", s))
            write_fasta_gz(os.path.join(people_root, pid(i), "immuannot_output", hap,
                                        "cds.fa.gz"), recs)
    t1 = os.path.join(root, "hla_calls_rich.tsv")
    pd.DataFrame(rows, columns=T1_COLS).to_csv(t1, sep="\t", index=False)
    coh = []
    for i, (low, up) in ancestry.items():
        probs = {f"p_{a.lower()}": 0.0 for a in m24.vc.ANCESTRY_ORDER}
        # person 2 is EUR but admixed -> strict UNASSIGNED
        probs[f"p_{up.lower()}"] = 0.6 if i == 2 else 0.95
        coh.append({"person_id": pid(i), "in_sr": True, "in_lr": True, "ancestry_pred": low,
                    **probs})
    cm = os.path.join(root, "cohort_membership.tsv")
    pd.DataFrame(coh).to_csv(cm, sep="\t", index=False)
    rel = os.path.join(root, "samples_relatedness.tsv")
    pd.DataFrame({"i.s": [pid(3), pid(1), pid(5)], "j.s": [pid(4), "PID9999", pid(6)],
                  "kin": [0.25, 0.3, 0.03]}).to_csv(rel, sep="\t", index=False)
    return {"table1": t1, "cohort": cm, "rel": rel, "outroot": people_root,
            "refdata": make_refdata(root)}


def test_end_to_end(tmp_path):
    root = str(tmp_path)
    p = build_outroot(root)
    out_dir = os.path.join(root, "committed")
    local_dir = os.path.join(root, "local")
    S = m24.main([
        "--table1", p["table1"], "--cohort-membership", p["cohort"],
        "--outroot", p["outroot"], "--refdata", p["refdata"],
        "--relatedness-table", p["rel"], "--out-dir", out_dir, "--local-dir", local_dir,
        "--threads", "3", "--min-denominator", "1",
    ])
    committed = ["novelty_field_counts.tsv", "sequence_class_counts.tsv",
                 "protein_level_novel_clusters.tsv", "novelty_by_field_report.md",
                 "summary.json", "fig_field_class_by_gene.png",
                 "fig_clean_novelty_rate_by_ancestry.png", "fig_novelty_funnel.png"]
    local = ["unrelated_set.tsv", "novelty_field_counts.full.tsv",
             "protein_level_novel_clusters.full.tsv", "allele_counts_by_resolution.tsv",
             "novel_protein_clusters.faa"]
    for f in committed:
        assert os.path.getsize(os.path.join(out_dir, f)) > 0, f
    for f in local:
        assert os.path.getsize(os.path.join(local_dir, f)) > 0, f
    assert set(os.listdir(out_dir)) == set(committed), os.listdir(out_dir)
    # No person ids in any committed text file.
    for f in committed:
        if f.endswith((".tsv", ".md", ".json")):
            assert "PID" not in open(os.path.join(out_dir, f)).read(), f

    # Unrelated: PID0003/PID0004 related -> drop PID0003 (tie, smallest id)
    un = pd.read_csv(os.path.join(local_dir, "unrelated_set.tsv"), sep="\t")
    assert set(un.loc[~un["unrelated"], "person_id"]) == {"PID0003"}
    assert un.set_index("person_id").loc["PID0002", "anc_strict"] == "UNASSIGNED"
    assert un.set_index("person_id").loc["PID0001", "anc_strict"] == "AFR"
    assert S["n_people_removed_for_relatedness"] == 1

    # Field-class counts (full): HLA-A f2 pooled/all = persons 1,2,3,4,6,7,8,9 = 8 calls
    fc = pd.read_csv(os.path.join(local_dir, "novelty_field_counts.full.tsv"), sep="\t")
    sel = fc[(fc["gene"] == "HLA-A") & (fc["ancestry_scheme"] == "pred")
             & (fc["ancestry"] == "POOLED") & (~fc["unrelated_only"])]
    by = sel.groupby(["field_class", "artifact_label"])["n_haplotypes"].sum()
    assert by[("f2_protein", "clean")] == 6
    assert by[("f2_protein", "homopolymer_indel")] == 1
    assert by[("f2_protein", "partial_cds")] == 1
    assert by[("f3_synonymous", "clean")] == 2
    assert by[("f4_noncoding", "clean")] == 1
    assert by[("f1_undetermined", "clean")] == 1
    assert by[("uncalled", "clean")] == 1
    assert by[("known", "clean")] == 24 - 8 - 2 - 1 - 1 - 1
    # unrelated view drops PID0003's two calls (one f2, one known)
    selu = fc[(fc["gene"] == "HLA-A") & (fc["ancestry_scheme"] == "pred")
              & (fc["ancestry"] == "POOLED") & (fc["unrelated_only"])]
    assert selu["n_haplotypes"].sum() == 22
    # strict scheme: PID0002 is UNASSIGNED
    sels = fc[(fc["gene"] == "HLA-A") & (fc["ancestry_scheme"] == "strict")
              & (fc["ancestry"] == "UNASSIGNED") & (~fc["unrelated_only"])]
    assert sels["n_haplotypes"].sum() == 2
    # committed version is suppressed
    fcc = pd.read_csv(os.path.join(out_dir, "novelty_field_counts.tsv"), sep="\t", dtype=str)
    assert set(fcc["n_haplotypes"]) == {"<20"} or "<20" in set(fcc["n_haplotypes"])
    assert not fcc["n_haplotypes"].isin([str(i) for i in range(1, 20)]).any()

    # Headline / sequence classes
    h = S["headline"]
    assert h["n_novel_protein_clusters"] == 3, h
    assert h["n_novel_protein_clusters_clean_recurrent_unrelated"] == 1, h
    assert h["n_cds_known_naming_artifacts"] == 1
    assert h["n_frameshift_or_stop_calls"] == 2
    assert h["n_novel_protein_calls"] == 5  # X x2, Y x2, Z
    assert h["n_novel_cds_synonymous_clusters"] == 2  # S and PK
    assert S["n_depth3_protein_not_catalogued"] == 1
    assert S["n_depth2_protein_catalogued"] == 1
    assert S["match_stats"]["n_ambiguous_copy"] >= 2
    seqc = pd.read_csv(os.path.join(out_dir, "sequence_class_counts.tsv"), sep="\t", dtype=str)
    assert "no_refdata" in set(seqc["seq_class"])

    # Clusters
    cl = pd.read_csv(os.path.join(local_dir, "protein_level_novel_clusters.full.tsv"), sep="\t")
    npc = cl[cl["cluster_type"] == "novel_protein"].set_index("n_aa_diffs", drop=False)
    x = cl[cl["aa_diffs"].fillna("").str.startswith("250:")].iloc[0]
    y = cl[cl["aa_diffs"].fillna("").str.startswith("30:")].iloc[0]
    assert len(npc) == 3
    assert x["n_persons"] == 2 and x["n_persons_unrelated_clean"] == 2
    assert bool(x["clean_recurrent_unrelated"])
    assert x["nearest_known_allele"] == "HLA-A*01:01:01:01" and x["nearest_distance"] == 1
    assert x["groove_method"] == "refdata_exon" and str(x["in_groove"]) == "False"
    assert S["imgt_version_alleles_csv"] == "IPD-IMGT/HLA-V3.55.0"
    assert x["n_persons_pred_AFR"] == 1 and x["n_persons_pred_EUR"] == 1
    assert x["n_persons_strict_UNASSIGNED"] == 1
    assert y["n_persons"] == 2 and y["n_persons_unrelated"] == 1
    assert not bool(y["clean_recurrent_unrelated"])
    assert str(y["in_groove"]) == "True"
    clc = pd.read_csv(os.path.join(out_dir, "protein_level_novel_clusters.tsv"), sep="\t",
                      dtype=str)
    assert set(clc["n_persons"]) == {"<20"}

    # Allele counts (WS3)
    ac = pd.read_csv(os.path.join(local_dir, "allele_counts_by_resolution.tsv"), sep="\t")
    base = ac[(ac["gene"] == "HLA-A") & (ac["ancestry_scheme"] == "pred")
              & (ac["ancestry"] == "POOLED") & (~ac["unrelated_only"])]
    prot_all = base[(base["resolution"] == "protein") & (base["exclude_flagged"] == "all")]
    pa = dict(zip(prot_all["allele_id"], prot_all["n_haplotypes"]))
    # 11 plain known + f4 (1) + protein_known (1) + synonymous S (1) = 14
    assert pa["HLA-A*01:01"] == 14, pa
    assert pa["HLA-A*02:01"] == 1  # cds_known naming artifact mapped back
    assert sum(1 for k in pa if "_prot_" in k) == 3
    assert sum(1 for k in pa if "_nonfunc_" in k) == 2
    prot_clean = base[(base["resolution"] == "protein") & (base["exclude_flagged"] == "clean_only")]
    assert not prot_clean["allele_id"].str.contains("_nonfunc_").any()
    cds_all = base[(base["resolution"] == "cds") & (base["exclude_flagged"] == "all")]
    ca = dict(zip(cds_all["allele_id"], cds_all["n_haplotypes"]))
    assert ca["HLA-A*01:01:01"] == 12 and ca["HLA-A*02:01:01"] == 1
    assert sum(1 for k in ca if "_cds_" in k) == 5  # X, Y, Z, S, PK
    # dropped: f1 + uncalled + no_refdata DQB1 are not in HLA-A table counts
    assert sum(pa.values()) == 24 - 2 - 0  # 24 A calls minus f1 and uncalled

    summ = json.load(open(os.path.join(out_dir, "summary.json")))
    assert summ["headline"]["n_novel_protein_clusters"] == "<20"
    report = open(os.path.join(out_dir, "novelty_by_field_report.md")).read()
    assert "Headline numbers" in report and "refdata_exon" in report


def test_limit_pilot(tmp_path):
    root = str(tmp_path)
    p = build_outroot(root)
    S = m24.main([
        "--table1", p["table1"], "--cohort-membership", p["cohort"],
        "--outroot", p["outroot"], "--refdata", p["refdata"], "--skip-relatedness",
        "--out-dir", os.path.join(root, "c"), "--local-dir", os.path.join(root, "l"),
        "--limit", "2", "--threads", "1",
    ])
    assert S["n_people"] == 2 and S["relatedness_skipped"]
    assert os.path.exists(os.path.join(root, "c", "pilot_limit2", "summary.json"))
    assert os.path.exists(os.path.join(root, "l", "pilot_limit2", "unrelated_set.tsv"))
    assert S["headline"]["n_novel_protein_clusters_clean_recurrent_unrelated"] == 1


if __name__ == "__main__":
    for name, fn in inspect.getmembers(sys.modules[__name__], inspect.isfunction):
        if name.startswith("test_"):
            kw = ({"tmp_path": pathlib.Path(tempfile.mkdtemp())}
                  if "tmp_path" in inspect.signature(fn).parameters else {})
            fn(**kw)
            print(f"  PASS  {name}")
    print("ok")
