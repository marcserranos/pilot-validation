#!/usr/bin/env python3
"""Tests for 25_noncoding_novelty_paf.py.

The synthetic cs strings are produced by `make_alignment` below, which models minimap2 directly:
an alignment is a list of (allele_base, contig_base) columns in the ALLELE's forward orientation;
for a minus-strand PAF row minimap2 aligns revcomp(allele) to the contig, so the columns are
reversed and complemented, and the contig itself is the reverse complement. cs is then written
with minimap2's definition (reference = target = contig): `*<contig><allele>`, `+<allele-only>`,
`-<contig-only>`. This is independent of the parser under test.

Run (no pytest needed):
  cd scripts/hla_popgen && python3 -c "import sys,tempfile,pathlib,inspect; sys.path.insert(0,'tests'); \
import test_noncoding_novelty_paf as t; [f(**({'tmp_path': pathlib.Path(tempfile.mkdtemp())} \
if 'tmp_path' in inspect.signature(f).parameters else {})) for n,f in \
inspect.getmembers(t, inspect.isfunction) if n.startswith('test_')]; print('ok')"
"""
import gzip
import hashlib
import importlib.util
import json
import os
import random
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)


def _load(filename, modname):
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load("25_noncoding_novelty_paf.py", "noncoding_novelty_paf")

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-"}


def rc(s):
    return "".join(COMP[b] for b in reversed(s))


def rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def apply_edits(allele, edits):
    """edits (0-based allele coords, non-overlapping):
       ('snv', pos, contig_base) | ('del', pos, n) | ('ins', pos_before, seq)
    -> alignment columns [(allele_base or '-', contig_base or '-')] in allele orientation."""
    by_pos = {}
    for e in edits:
        by_pos.setdefault(e[1], []).append(e)
    cols = []
    i = 0
    while i < len(allele):
        evs = by_pos.get(i, [])
        for e in evs:
            if e[0] == "ins":
                cols.extend(("-", b) for b in e[2])
        dele = [e for e in evs if e[0] == "del"]
        snv = [e for e in evs if e[0] == "snv"]
        if dele:
            for k in range(dele[0][2]):
                cols.append((allele[i + k], "-"))
            i += dele[0][2]
            continue
        cols.append((allele[i], snv[0][2] if snv else allele[i]))
        i += 1
    return cols


def cs_from_cols(cols):
    """cols in TARGET forward orientation: (query_base, target_base)."""
    out = []
    run = 0
    i = 0
    while i < len(cols):
        q, t = cols[i]
        if q != "-" and t != "-" and q == t:
            run += 1
            i += 1
            continue
        if run:
            out.append(f":{run}")
            run = 0
        if q != "-" and t != "-":
            out.append(f"*{t.lower()}{q.lower()}")
            i += 1
        elif t == "-":
            j = i
            while j < len(cols) and cols[j][1] == "-":
                j += 1
            out.append("+" + "".join(c[0] for c in cols[i:j]).lower())
            i = j
        else:
            j = i
            while j < len(cols) and cols[j][0] == "-":
                j += 1
            out.append("-" + "".join(c[1] for c in cols[i:j]).lower())
            i = j
    if run:
        out.append(f":{run}")
    return "".join(out)


def make_alignment(allele, edits, strand, qstart=0, qend=None):
    """-> (contig_seq_aligned_part, cs, qstart, qend) for an alignment of allele[qstart:qend]."""
    qend = len(allele) if qend is None else qend
    sub_edits = [(e[0], e[1] - qstart) + tuple(e[2:]) for e in edits]
    cols = apply_edits(allele[qstart:qend], sub_edits)
    contig = "".join(t for _q, t in cols if t != "-")
    if strand == "-":
        cols = [(COMP[q], COMP[t]) for q, t in reversed(cols)]
        contig = rc(contig)
    return contig, cs_from_cols(cols), qstart, qend


def paf_line(qname, qlen, qs, qe, strand, tname, tlen, ts, te, nm, cs):
    return "\t".join(map(str, [qname, qlen, qs, qe, strand, tname, tlen, ts, te, 100, 100, 60,
                               f"NM:i:{nm}", "tp:A:P", f"cs:Z:{cs}"])) + "\n"


# ---------------------------------------------------------------------------
# 1. cs orientation
# ---------------------------------------------------------------------------
def test_cs_orientation_plus():
    ev = m.parse_cs_events(":10*ga:5+tt:3-c:2", 0, 23, "+")
    # *ga: contig g, allele a -> SNV allele A > contig G at allele pos 10
    assert ev[0] == (10, "SNV", "A", "G"), ev
    # +tt: allele bases TT absent from contig, at allele 16..17
    assert ev[1] == (16, "CONTIG_DEL", "TT", ""), ev
    # -c: contig has extra C before allele pos 21
    assert ev[2] == (21, "CONTIG_INS", "", "C"), ev
    assert m.event_token(ev[0]) == "11A>G"
    assert m.event_token(ev[1]) == "17_18delTT"
    assert m.event_token(ev[2]) == "21_22insC"


def test_cs_matches_synthetic_generator():
    allele = rand_seq(60, 1)
    edits = [("snv", 7, "A" if allele[7] != "A" else "C"), ("del", 20, 2), ("ins", 40, "GT")]
    _, cs, qs, qe = make_alignment(allele, edits, "+")
    ev = m.parse_cs_events(cs, qs, qe, "+")
    assert [e[1] for e in ev] == ["SNV", "CONTIG_DEL", "CONTIG_INS"]
    assert ev[0] == (7, "SNV", allele[7], edits[0][2])
    assert ev[1] == (20, "CONTIG_DEL", allele[20:22], "")
    assert ev[2] == (40, "CONTIG_INS", "", "GT")


def test_cs_bad_input_raises():
    for cs, qs, qe in ((":10", 0, 11), (":5x3", 0, 8), (":5~gt10ag", 0, 5)):
        try:
            m.parse_cs_events(cs, qs, qe, "+")
        except m.CsError:
            continue
        raise AssertionError(f"expected CsError for {cs}")


# ---------------------------------------------------------------------------
# 2. minus-strand normalization
# ---------------------------------------------------------------------------
def _non_repeat_allele():
    # no homopolymer runs >= 2, so indel placement is unique
    return ("ACGTACAGTCATGCATCGATGCTAGCATGACTGATCGTACGATCAGTCGATGCAT"
            "GCATCAGTGCATGCTAGCTGACTGCA")


def test_minus_strand_same_signature():
    allele = _non_repeat_allele()
    edits = [("snv", 5, "G" if allele[5] != "G" else "T"), ("del", 18, 3), ("ins", 33, "TTG"),
             ("snv", 50, "A" if allele[50] != "A" else "C")]
    sigs = {}
    for strand in "+-":
        _, cs, qs, qe = make_alignment(allele, edits, strand)
        sigs[strand] = m.build_signature("HLA-A", "HLA-A*01:01:01:01", cs, qs, qe, strand)
    assert sigs["+"]["signature"] == sigs["-"]["signature"], (sigs["+"]["signature"],
                                                              sigs["-"]["signature"])
    assert sigs["+"]["signature_id"] == sigs["-"]["signature_id"]
    assert sigs["+"]["n_snv"] == 2 and sigs["+"]["n_indel"] == 2


def test_minus_strand_partial_alignment():
    allele = _non_repeat_allele()
    edits = [("snv", 20, "C" if allele[20] != "C" else "A"), ("del", 30, 1)]
    a = m.parse_cs_events(make_alignment(allele, edits, "+", 5, 70)[1], 5, 70, "+")
    b = m.parse_cs_events(make_alignment(allele, edits, "-", 5, 70)[1], 5, 70, "-")
    assert a == b, (a, b)
    assert a[0][0] == 20 and a[1] == (30, "CONTIG_DEL", allele[30], "")


def test_minus_strand_repeat_needs_sequence():
    # Homopolymer AAAAA at 20..24; the same 1-bp deletion placed at different offsets inside
    # the run on the two strands (as minimap2 does: left-aligned in TARGET orientation).
    allele = _non_repeat_allele()[:20] + "AAAAA" + _non_repeat_allele()[25:]
    _, cs_p, qs, qe = make_alignment(allele, [("del", 20, 1)], "+")
    _, cs_m, _, _ = make_alignment(allele, [("del", 24, 1)], "-")
    no_seq_p = m.build_signature("HLA-A", "X", cs_p, qs, qe, "+")
    no_seq_m = m.build_signature("HLA-A", "X", cs_m, qs, qe, "-")
    assert no_seq_p["signature"] != no_seq_m["signature"]       # the documented limitation
    assert no_seq_p["homopolymer_context_checked"] is False
    s_p = m.build_signature("HLA-A", "X", cs_p, qs, qe, "+", allele_seq=allele)
    s_m = m.build_signature("HLA-A", "X", cs_m, qs, qe, "-", allele_seq=allele)
    assert s_p["signature"] == s_m["signature"] == "X|21delA", (s_p["signature"], s_m["signature"])
    assert s_p["signature_class"] == "homopolymer_only"
    assert s_p["event_flank_runs"] == "4"


def test_insertion_left_normalization():
    allele = "ACGTCAGGGGTACA"
    # insert G at the end of the GGGG run (before pos 10) vs at its start (before pos 6)
    ev = [(10, "CONTIG_INS", "", "G")]
    assert m.left_normalize(ev, allele) == [(6, "CONTIG_INS", "", "G")]
    ev2 = [(10, "CONTIG_INS", "", "TG")]   # dinucleotide not rotating past T at pos 9? allele[9]=G
    assert m.left_normalize(ev2, allele) == [(9, "CONTIG_INS", "", "GT")]


# ---------------------------------------------------------------------------
# 3. homopolymer flag
# ---------------------------------------------------------------------------
def test_homopolymer_context_and_content_only():
    allele = "ACGTAAAAGTC"
    ins_run = (4, "CONTIG_INS", "", "A")          # next to AAAA
    ins_lonely = (2, "CONTIG_INS", "", "T")       # between C and G; flank T? allele[1]=C, [2]=G
    assert m.homopolymer_info(ins_run, allele) == (True, 4)
    assert m.homopolymer_info(ins_lonely, allele) == (False, 0)
    # content-only: every single-base indel is flagged, context not checked
    assert m.homopolymer_info(ins_lonely, None) == (True, None)
    assert m.homopolymer_info((0, "CONTIG_DEL", "AC", ""), None) == (False, None)
    # deletion of 1 A from the run of 4: flank excludes the deleted base -> 3
    assert m.homopolymer_info((4, "CONTIG_DEL", "A", ""), allele) == (True, 3)
    # min_run configurable
    assert m.homopolymer_info((4, "CONTIG_DEL", "A", ""), allele, min_run=4) == (False, 3)
    sig = m.build_signature("G", "X", ":2-t:9", 0, 11, "+", allele_seq=allele)
    assert sig["signature_class"] == "indel_nonhomopolymer" and sig["homopolymer_context_checked"]
    sig = m.build_signature("G", "X", ":2-t:9", 0, 11, "+")
    assert sig["signature_class"] == "homopolymer_only" and not sig["homopolymer_context_checked"]


def test_inconsistent_sequence_falls_back():
    allele = "ACGTAAAAGTC"
    sig = m.build_signature("G", "X", ":1*ag:9", 0, 11, "+", allele_seq=allele)  # allele[1]=C != G
    assert sig["seq_consistent"] is False and sig["context_mode"] == "content_only"


# ---------------------------------------------------------------------------
# 4. row selection
# ---------------------------------------------------------------------------
def _row(qname, qs, qe, tname, ts, te, nm, strand="+"):
    return m.parse_paf_line(paf_line(qname, 100, qs, qe, strand, tname, 10000, ts, te, nm, ":1"))


def test_row_selection_template_and_fallback():
    rows = [
        _row("HLA-A*01:01:01:01", 0, 80, "ctg1", 1000, 1080, 3),
        _row("HLA-A*01:01:01:01", 0, 100, "ctg1", 1000, 1100, 5),   # longest -> chosen
        _row("HLA-A*01:01:01:01", 0, 100, "ctg2", 1000, 1100, 0),   # wrong contig
        _row("HLA-A*01:01:01:01", 0, 100, "ctg1", 5000, 5100, 0),   # outside gene span
        _row("HLA-A*01:01:01:02", 0, 90, "ctg1", 1000, 1090, 1),
        _row("HLA-A*01:01:01:03", 0, 100, "ctg1", 1000, 1100, 1),
        _row("HLA-A*01:01:010:01", 0, 100, "ctg1", 1000, 1100, 0),  # not the same 3 fields
    ]
    r, path = m.select_paf_row(rows, "ctg1", 1001, 1100, m.norm_allele("A*01:01:01:01"),
                               m.three_field_prefix("HLA-A*01:01:01:new"))
    assert path == "template" and r["qend"] == 100 and r["nm"] == 5 and r["tname"] == "ctg1"
    # template absent -> fallback: lowest NM (tie 1), then longest -> *03
    r, path = m.select_paf_row(rows, "ctg1", 1001, 1100, m.norm_allele("HLA-A*01:01:01:99"),
                               m.three_field_prefix("HLA-A*01:01:01:new"))
    assert path == "fallback_3field" and r["qname"] == "HLA-A*01:01:01:03", r
    # nothing on this contig
    r, path = m.select_paf_row(rows, "ctg9", 1, 100, "A*01:01:01:01", "A*01:01:01:")
    assert r is None and path == "no_row"
    # 1-based inclusive gene_end == PAF tstart+... boundary: gene 1..1000 does not overlap 1000..
    r, path = m.select_paf_row(rows, "ctg1", 1, 1000, "A*01:01:01:01", None)
    assert r is None


def test_name_normalization():
    assert m.norm_allele("HLA-A*01:01") == m.norm_allele("A*01:01") == m.norm_allele("hla-a*01:01")
    assert m.three_field_prefix("HLA-DRB1*15:01:01:new") == "DRB1*15:01:01:"
    assert m.three_field_prefix("HLA-A*01:new") is None
    assert m.key_prefix3("A*01:01:01:01N") == "A*01:01:01:"


# ---------------------------------------------------------------------------
# 5. unrelated set / ancestry / suppression
# ---------------------------------------------------------------------------
def test_greedy_unrelated():
    pairs = [("a", "b", 0.25), ("a", "c", 0.25), ("b", "c", 0.05), ("d", "e", 0.10),
             ("f", "g", 0.01), ("x", "h", 0.3)]
    kept = m.greedy_unrelated(pairs, {"a", "b", "c", "d", "e", "f", "g", "h"})
    # a has 2 relatives, b and c 2 each (b-c 0.05 >= 0.0442) -> tie, 'a' removed first by id,
    # then b-c remain -> tie 1/1 -> 'b' removed. d-e tie -> 'd'. f-g below threshold.
    # x is not long-read -> h stays.
    assert kept == {"c", "e", "f", "g", "h"}, kept


def test_strict_ancestry():
    df = pd.DataFrame({"ancestry_pred": ["AFR", "EUR", "EUR", pd.NA],
                       "p_afr": [0.95, 0.0, 0.2, 0.9], "p_eur": [0.05, 0.89, 0.8, 0.1]})
    s = m.strict_ancestry(df)
    assert s.iloc[0] == "AFR" and pd.isna(s.iloc[1]) and pd.isna(s.iloc[2]) and pd.isna(s.iloc[3])


def test_suppress():
    assert m.suppress(0) == 0 and m.suppress(19) == "<20" and m.suppress(20) == 20
    assert m.suppress(1) == "<20" and m.suppress(float("nan")) == "NA" and m.suppress(pd.NA) == "NA"


# ---------------------------------------------------------------------------
# 6. end to end
# ---------------------------------------------------------------------------
def _build_outroot(root):
    """Synthetic cohort: 30 persons, 1 depth-4 HLA-A call per hap1.
    - persons p00..p23: signature S1 (SNV + non-repeat del), half plus / half minus strand
      -> one cluster with 24 persons; p00/p01 are related (one dropped from unrelated).
    - p24..p27: homopolymer-only signature (1-bp del inside AAAAA), on both strands.
    - p28: template absent from PAF -> fallback row.
    - p29: no PAF row on the contig -> no_row.
    All share the same observed CDS -> one former novel_id.
    """
    people = root / "people"
    base = _non_repeat_allele()
    allele = base[:20] + "AAAAA" + base[25:] + rand_seq(40, 7)
    other = allele[:60] + ("C" if allele[60] != "C" else "G") + allele[61:]
    cds = "ATGGCC" + rand_seq(30, 3) + "TGA"
    sha1 = hashlib.sha1(cds.encode()).hexdigest()
    s1 = [("snv", 2, "T" if allele[2] != "T" else "G"), ("del", 33, 2)]
    hp = [("del", 22, 1)]
    t1_rows, rel_rows, coh_rows = [], [], []
    anc_cycle = ["AFR", "EUR", "EAS"]
    for i in range(30):
        pid = f"p{i:02d}"
        d = people / pid / "immuannot_output" / "hap1"
        d.mkdir(parents=True)
        contig = f"{pid}_ctg"
        strand = "+" if i % 2 == 0 else "-"
        template = "HLA-A*01:01:01:01"
        lines = []
        edits = hp if 24 <= i <= 27 else s1
        if i != 29:
            seq, cs, qs, qe = make_alignment(allele, edits, strand)
            ts = 1000
            te = ts + len(seq)
            qname = "HLA-A*01:01:01:05" if i == 28 else template
            lines.append(paf_line(qname, len(allele), qs, qe, strand, contig, 50000, ts, te,
                                  len(edits), cs))
            # decoys: other allele (higher NM) on same contig; template on a far locus
            _, cs2, _, _ = make_alignment(other, [], strand)
            lines.append(paf_line("HLA-A*02:01:01:01", len(other), 0, len(other), strand,
                                  contig, 50000, ts, ts + len(other), 9, cs2))
            lines.append(paf_line(template, len(allele), 0, 30, "+", contig, 50000, 40000,
                                  40030, 0, ":30"))
        else:
            lines.append(paf_line(template, len(allele), 0, len(allele), "+", "elsewhere",
                                  50000, 1000, 1000 + len(allele), 0, f":{len(allele)}"))
        with gzip.open(d / "mm2.ipd.gen.paf.gz", "wt") as fh:
            fh.writelines(lines)
        with gzip.open(d / "cds.fa.gz", "wt") as fh:
            fh.write(f">{contig}_HLA-A_1 GT..AG\n{cds.lower()}\n")
        t1_rows.append({
            "person_id": pid, "hap": "hap1", "contig": contig, "gene": "HLA-A", "copy_index": 1,
            "gene_class": "classical_I", "consensus": "HLA-A*01:01:01:new", "n_fields": 4,
            "is_novel": True, "novelty_depth": 4, "novelty_class": "beyond_cds",
            "template_allele": "HLA-A*01:01:01:99" if i == 28 else template,
            "template_distance": 2, "gene_start": 1001, "gene_end": 1000 + len(allele),
            "strand": strand, "template_warning": "NA"})
        # a known (non-novel) HLA-B row, must be ignored
        t1_rows.append({**t1_rows[-1], "gene": "HLA-B", "consensus": "HLA-B*07:02:01:01",
                        "is_novel": False, "novelty_depth": pd.NA, "novelty_class": pd.NA})
        a = anc_cycle[i % 3]
        coh_rows.append({"person_id": pid, "in_sr": True, "in_lr": True,
                         "ancestry_pred": a.lower(),
                         **{f"p_{x.lower()}": (0.95 if x == a else 0.01)
                            for x in ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]}})
    rel_rows = [{"i.s": "p00", "j.s": "p01", "kin": 0.25},
                {"i.s": "p00", "j.s": "outsider", "kin": 0.25}]
    pd.DataFrame(t1_rows).to_csv(root / "hla_calls_rich.tsv", sep="\t", index=False)
    pd.DataFrame(coh_rows).to_csv(root / "cohort_membership.tsv", sep="\t", index=False)
    pd.DataFrame(rel_rows).to_csv(root / "samples_relatedness.tsv", sep="\t", index=False)
    nid = f"HLA-A_nov_{sha1[:8]}"
    pd.DataFrame([{"novel_id": nid, "gene": "HLA-A", "gene_class": "classical_I",
                   "nearest_allele": "HLA-A*01:01:01:new", "cds_distance": 0.0,
                   "n_aa_changes": "", "novelty_class": "beyond_cds", "cds_seq_sha1": sha1,
                   "cds_len": len(cds), "n_haplotypes": 30, "n_persons": 30, "n_ancestries": 3,
                   "ancestry_counts": "{}", "is_homopolymer_indel_only": False,
                   "both_haps_one_person": False, "confidence_tier": "high",
                   "passes_qc": True, "passes_qc_singleton_ok": True}]) \
        .to_csv(root / "novel_alleles.tsv", sep="\t", index=False)
    refdata = root / "refdata"
    (refdata / "CDSseq").mkdir(parents=True)
    with gzip.open(refdata / "CDSseq" / "HLA-A.fa.gz", "wt") as fh:
        fh.write(">HLA-A*01:01:01:01\nATG\n")
    with gzip.open(refdata / "hla_gen.fa.gz", "wt") as fh:
        fh.write(f">HLA-A*01:01:01:01 {len(allele)} bp\n{allele[:50]}\n{allele[50:]}\n")
        fh.write(f">HLA-A*01:01:01:05\n{allele}\n")
        fh.write(f">HLA-B*07:02:01:01\n{other}\n")
    L = len(allele)
    with gzip.open(refdata / "alleles.csv.gz", "wt") as fh:
        # *01: exon1 1..10 (UTR5 1..4), intron1 11..40, exon2 41..60, intron2, exon3 100..L
        fh.write(f"gene=HLA-A\tallele=HLA-A*01:01:01:01\ttype=gene\tgeneRange=1..{L}\t"
                 f"CDS=5..10,41..60,100..110\tUTR=1..4,111..{L}\tnexon=3\t"
                 f"exon=1:1..10,2:41..60,3:100..{L}\tintron=GT-AG\n")
        # *05: CDS starts at 1 -> an event at position 3 touches CDS
        fh.write(f"gene=HLA-A\tallele=HLA-A*01:01:01:05\ttype=gene\tgeneRange=1..{L}\t"
                 f"CDS=1..10,41..60,100..110\tUTR=111..{L}\tnexon=3\t"
                 f"exon=1:1..10,2:41..60,3:100..{L}\n")
    return people, refdata, nid, allele


def _run(tmp_path, with_seq):
    people, refdata, nid, allele = _build_outroot(tmp_path)
    out = tmp_path / ("out_seq" if with_seq else "out_noseq")
    local = tmp_path / ("local_seq" if with_seq else "local_noseq")
    argv = ["--table1", str(tmp_path / "hla_calls_rich.tsv"),
            "--cohort-membership", str(tmp_path / "cohort_membership.tsv"),
            "--outroot", str(people),
            "--relatedness-table", str(tmp_path / "samples_relatedness.tsv"),
            "--novel-alleles", str(tmp_path / "novel_alleles.tsv"),
            "--refdata", str(refdata if with_seq else tmp_path / "nope"),
            "--out-dir", str(out), "--local-dir", str(local), "--threads", "4"]
    summary = m.main(argv)
    return summary, out, local, nid, allele


def test_end_to_end_with_sequence(tmp_path):
    summary, out, local, nid, allele = _run(tmp_path, with_seq=True)
    calls = pd.read_csv(local / "25_noncoding_calls.tsv", sep="\t")
    assert len(calls) == 30, len(calls)                      # HLA-B known row ignored
    st = calls.set_index("person_id")
    assert st.loc["p29", "status"] == "no_row"
    assert st.loc["p28", "selection_path"] == "fallback_3field"
    assert st.loc["p28", "used_allele"] == "HLA-A*01:01:01:05"
    assert (calls["novel_id"] == nid).all()
    ok = calls[calls["status"] == "ok"]
    assert (ok["context_mode"] == "sequence").all()
    s1 = ok[ok["person_id"].isin([f"p{i:02d}" for i in range(24)])]
    assert s1["signature_id"].nunique() == 1, s1["signature"].unique()   # both strands agree
    assert set(s1["paf_strand"]) == {"+", "-"}
    hp = ok[ok["person_id"].isin(["p24", "p25", "p26", "p27"])]
    assert hp["signature_id"].nunique() == 1, hp["signature"].unique()
    assert (hp["signature_class"] == "homopolymer_only").all()
    assert hp["event_tokens"].iloc[0] == "21delA"
    assert s1["event_regions"].iloc[0] == "UTR5,intron1", s1["event_regions"].iloc[0]
    assert (s1["region_class"] == "involves_utr").all()
    assert (hp["region_class"] == "intron_only").all() and hp["event_regions"].iloc[0] == "intron1"
    assert st.loc["p28", "region_class"] == "touches_cds"
    assert st.loc["p28", "event_regions"] == "exon1,intron1"

    full = pd.read_csv(local / "25_noncoding_signature_clusters_full.tsv", sep="\t")
    big = full.sort_values("n_persons", ascending=False).iloc[0]
    assert big["n_persons"] == 24 and big["n_persons_unrelated"] == 23
    assert big["n_plus_strand"] == 12 and big["n_minus_strand"] == 12
    assert big["n_persons_strict_AFR"] == 8

    committed = pd.read_csv(out / "noncoding_signature_clusters.tsv", sep="\t", dtype=str)
    assert len(committed) == 1 and committed["n_persons"].iloc[0] == "24"
    assert committed["n_persons_strict_AFR"].iloc[0] == "<20"
    pool = pd.read_csv(out / "beyond_cds_pooling.tsv", sep="\t", dtype=str)
    assert len(pool) == 1 and pool["n_distinct_signatures"].iloc[0] == "3"
    assert pool["old_n_persons"].iloc[0] == "30"
    frac = pd.read_csv(out / "homopolymer_fraction.tsv", sep="\t", dtype=str)
    row = frac[(frac["gene"] == "ALL") & (frac["ancestry"] == "POOLED")
               & (frac["ancestry_scheme"] == "pred") & (frac["unrelated_only"] == "False")].iloc[0]
    assert row["n_haplotypes"] == "29" and abs(float(row["frac_homopolymer_only"]) - 4 / 29) < 1e-9
    sfs = pd.read_csv(local / "25_nonhomopolymer_sfs_full.tsv", sep="\t")
    allrow = sfs[(sfs["gene"] == "ALL") & (~sfs["unrelated_only"])].set_index(
        "carrier_haplotypes_bin")["n_signatures"]
    assert allrow["21-50"] == 1 and allrow["1"] == 1 and allrow.sum() == 2
    for f in ("noncoding_novelty_report.md", "summary.json", "gene_summary.tsv",
              "pooling_distribution.tsv", "nonhomopolymer_sfs.tsv",
              "fig1_signatures_per_former_cluster.png", "fig2_homopolymer_fraction.png",
              "fig3_nonhomopolymer_sfs.png"):
        assert (out / f).exists(), f
    # nothing person-level in committed outputs
    for f in out.iterdir():
        if f.suffix in (".tsv", ".md", ".json"):
            txt = f.read_text()
            assert "p0" not in txt and "p1" not in txt and "_ctg" not in txt, f
    js = json.loads((out / "summary.json").read_text())
    assert js["n_depth4_haplotypes"] == 30 and js["frac_context_sequence"] == 1.0
    assert abs(js["frac_region_class"]["intron_only"] - 4 / 29) < 1e-3
    assert js["touches_cds"]["n"] == "<20"
    assert big["region_class"] == "involves_utr" and big["n_events_utr"] == 1 \
        and big["n_events_intron"] == 1 and big["n_events_cds"] == 0
    # resume: second run reads the checkpoint and gives the same answer
    summary2 = m.main(["--table1", str(tmp_path / "hla_calls_rich.tsv"),
                       "--cohort-membership", str(tmp_path / "cohort_membership.tsv"),
                       "--outroot", str(tmp_path / "people"),
                       "--relatedness-table", str(tmp_path / "samples_relatedness.tsv"),
                       "--novel-alleles", str(tmp_path / "novel_alleles.tsv"),
                       "--refdata", str(tmp_path / "refdata"), "--no-figures",
                       "--out-dir", str(out), "--local-dir", str(local), "--threads", "1"])
    assert summary2["n_signatures"] == summary["n_signatures"]


def test_end_to_end_without_sequence(tmp_path):
    summary, out, local, nid, allele = _run(tmp_path, with_seq=False)
    calls = pd.read_csv(local / "25_noncoding_calls.tsv", sep="\t")
    ok = calls[calls["status"] == "ok"]
    assert (ok["context_mode"] == "content_only").all()
    assert (~ok["homopolymer_context_checked"].astype(bool)).all()
    # non-repeat signature still agrees across strands without the sequence
    s1 = ok[ok["person_id"].isin([f"p{i:02d}" for i in range(24)])]
    assert s1["signature_id"].nunique() == 1
    # the 2-bp non-homopolymer deletion keeps S1 out of homopolymer_only
    assert (s1["signature_class"] == "snv_plus_indel").all()
    assert summary["n_depth4_haplotypes"] == 30


def test_limit_and_novel_matches(tmp_path):
    people, refdata, nid, allele = _build_outroot(tmp_path)
    t1 = pd.read_csv(tmp_path / "hla_calls_rich.tsv", sep="\t")
    nm = t1[t1["gene"] == "HLA-A"][["person_id", "hap", "contig", "gene"]].copy()
    nm["cds_seq_sha1"] = "f" * 40
    nm.to_csv(tmp_path / "matches.tsv", sep="\t", index=False)
    summary = m.main(["--table1", str(tmp_path / "hla_calls_rich.tsv"),
                      "--cohort-membership", str(tmp_path / "cohort_membership.tsv"),
                      "--outroot", str(people),
                      "--relatedness-table", str(tmp_path / "samples_relatedness.tsv"),
                      "--novel-alleles", str(tmp_path / "novel_alleles.tsv"),
                      "--novel-matches", str(tmp_path / "matches.tsv"),
                      "--refdata", str(refdata), "--limit", "5", "--no-figures",
                      "--out-dir", str(tmp_path / "o"), "--local-dir", str(tmp_path / "l")])
    calls = pd.read_csv(tmp_path / "l" / "25_noncoding_calls.tsv", sep="\t")
    assert len(calls) == 5 and (calls["novel_id"] == "HLA-A_nov_ffffffff").all()
    assert summary["n_persons"] == "<20"


def test_refdata_discovery_prefers_gen_fa(tmp_path):
    (tmp_path / "CDSseq").mkdir()
    for name in ("gen.fa.gz", "other.fa.gz", "CDSseq/HLA-A.fa.gz"):
        with gzip.open(tmp_path / name, "wt") as fh:
            fh.write(">HLA-A*01:01:01:01 HLA00001 x 10bp\nACGTACGTAC\n")
    assert m.find_allele_fastas(str(tmp_path)) == [str(tmp_path / "gen.fa.gz")]
    (tmp_path / "gen.fa.gz").unlink()
    assert m.find_allele_fastas(str(tmp_path)) == [str(tmp_path / "other.fa.gz")]
    seqs = m.load_allele_sequences([str(tmp_path / "other.fa.gz")], {"A*01:01:01:01"}, set())
    assert seqs == {"A*01:01:01:01": "ACGTACGTAC"}


def test_region_classification():
    sig = m.build_signature("G", "X", ":2*ag:17", 0, 20, "+", ann=None)
    assert sig["event_regions"] == "unannotated" and sig["region_class"] == "unannotated"
    assert m.build_signature("G", "X", ":20", 0, 20, "+")["region_class"] == "no_difference"
    ann = {"cds": [(301, 373), (504, 773)], "utr": [(1, 300), (774, 900)],
           "exons": m.parse_exons("1:1..373,2:504..900"), "gene_range": [(1, 900)]}
    lab = [m.region_of(p - 1, ann) for p in (5, 305, 400, 600, 800, 950)]
    assert lab == ["UTR5", "exon1", "intron1", "exon2", "UTR3", "outside"], lab
    assert m.classify_regions(["intron1", "intron3"]) == "intron_only"
    assert m.classify_regions(["intron1", "UTR3"]) == "involves_utr"
    assert m.classify_regions(["UTR5", "exon2"]) == "touches_cds"
    # 1-based position 3 of a 20-bp allele with CDS 3..5
    ann = {"cds": [(3, 5)], "utr": [(1, 2), (6, 20)], "exons": [(1, 1, 20)],
           "gene_range": [(1, 20)]}
    sig = m.build_signature("G", "X", ":2*ag:17", 0, 20, "+", ann=ann)
    assert sig["event_regions"] == "exon1" and sig["n_events_cds"] == 1
