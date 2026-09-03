#!/usr/bin/env python3
"""Real assertions against 01_extract_rich.py and 02_build_cohorts.py, run against synthetic
fixtures (tests/make_fixtures.py) and hand-crafted GTF snippets for cases the fixture generator
doesn't happen to produce (e.g. C4 rows -- the fixture never emits C4 genes at all).

Why this file exists: SCHEMA.md hard rule #6 -- "every script here must run against the synthetic
fixtures before it is handed to Marc." This file is the automatable half of that; the two
extraction scripts' own --validate modes cover the other half (row-grain uniqueness on real/sample
output).

Run:
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures_test -n 200
    python3 scripts/hla_popgen/tests/test_extraction.py --fixtures /tmp/hla_fixtures_test
"""
import argparse
import importlib.util
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


extract = _load_module("01_extract_rich.py", "hla_popgen_extract_rich")
cohorts = _load_module("02_build_cohorts.py", "hla_popgen_build_cohorts")

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        msg = f"  FAIL  {name}" + (f" -- {detail}" if detail else "")
        print(msg)
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Unit-level tests on hand-crafted attribute strings / GTF snippets
# ---------------------------------------------------------------------------
def test_unquoted_template_distance():
    attrs = extract.parse_attrs(
        'gene_id "IAG123456"; template_allele "HLA-A*01:01:01:01"; '
        'template_distance 7; gene_name "HLA-A";'
    )
    check("unquoted template_distance parsed as '7'", attrs.get("template_distance") == "7",
          detail=repr(attrs))


def test_cds_mut_not_split_naively():
    raw = ('HLA-A*01:01:01:01|:301*ag:302+cc|K(AAG)<E(GAG):D(GAT)<L(CTT)')
    attrs = extract.parse_attrs(f'gene_id "X"; cds_mut "{raw}";')
    check("cds_mut with |,<,:,* preserved verbatim", attrs.get("cds_mut") == raw,
          detail=attrs.get("cds_mut"))
    check("n_aa_changes counts '<' occurrences (2)", extract.n_aa_changes(raw) == 2)


def test_conditional_fields_absent_not_zero():
    attrs = extract.parse_attrs('gene_id "X"; transcript_id "Y"; gene_name "HLA-A"; '
                                 'consensus "HLA-A*01:01:01:01"; alleles "HLA-A*01:01:01:01";')
    check("cds_distance absent (not 0) when not in attrs", "cds_distance" not in attrs)
    check("template_warning absent (not empty string) when not in attrs",
          "template_warning" not in attrs)
    # Full pipeline: build_table1_row must emit None (-> NA on disk), not 0/False-as-data.
    row = {
        "contig": "ctg1", "gene_id": "IAG1", "gene_name_raw": "HLA-A",
        "template_allele": "HLA-A*01:01:01:01", "template_distance": "0",
        "gene_start": 100, "gene_end": 200, "strand": "+",
        "consensus": "HLA-A*01:01:01:01", "alleles": "HLA-A*01:01:01:01",
        "template_warning": None, "cds_distance": None, "cds_mut": None,
    }
    t1, warn = extract.build_table1_row("p1", "hap1", row)
    check("cds_distance is None (NA) when gene-level match was perfect",
          t1["cds_distance"] is None, detail=str(t1))
    check("cds_mut is None (NA) when absent", t1["cds_mut"] is None)
    check("template_warning is None (NA) when absent", t1["template_warning"] is None)
    check("has_warning is False when template_warning absent", t1["has_warning"] is False)


def test_novelty_depth_classification():
    cases = [
        ("HLA-A*new", 1, "undetermined"),  # regression: Fix 1, depth-1 gene-level unresolved
        ("HLA-A*02:new", 2, "protein_altering"),
        ("HLA-C*07:01:new", 3, "synonymous"),
        ("HLA-DPB1*17:01:01:new", 4, "beyond_cds"),
        ("HLA-A*01:01:01:01", None, None),  # exact match, no novelty
        ("undetermined", None, None),
    ]
    for consensus, exp_depth, exp_class in cases:
        n_fields, is_novel, depth, cls, warn = extract.parse_consensus(consensus)
        check(f"novelty depth for {consensus!r} == {exp_depth}", depth == exp_depth,
              detail=f"got depth={depth}")
        check(f"novelty class for {consensus!r} == {exp_class}", cls == exp_class,
              detail=f"got class={cls}")
    check("is_novel True only when 'new' present",
          extract.parse_consensus("HLA-A*02:new")[1] is True and
          extract.parse_consensus("HLA-A*01:01:01:01")[1] is False)
    check("depth-1 novel call ('HLA-A*new') is_novel True (not confused with the literal string "
          "'undetermined' consensus, which has no fields at all)",
          extract.parse_consensus("HLA-A*new")[1] is True)
    check("depth-1 does not raise an 'unexpected novelty depth' warning (it's a real, documented "
          "case per SCHEMA.md, not a parser anomaly)",
          extract.parse_consensus("HLA-A*new")[4] is None,
          detail=str(extract.parse_consensus("HLA-A*new")[4]))


def test_novelty_depth1_full_row_no_crash_no_misbucket():
    """Regression for Fix 1: novelty_depth==1 ('HLA-A*new', even the first/gene-resolution field
    unresolved -- consensusCall()'s commonprefix truncation colliding on tied candidates that
    disagree at field 1) must not crash build_table1_row, must not be silently dropped, and must
    not be mis-bucketed into protein_altering/synonymous/beyond_cds via a dict .get(depth) that
    only expected keys 2/3/4."""
    check("NOVELTY_CLASS_BY_DEPTH has depth 1 mapped to 'undetermined' (would otherwise KeyError-"
          "free but silently return None via .get(), mis-classifying a real, documented case)",
          extract.NOVELTY_CLASS_BY_DEPTH.get(1) == "undetermined",
          detail=str(extract.NOVELTY_CLASS_BY_DEPTH))
    row = {
        "contig": "ctg1", "gene_id": "IAG999001", "gene_name_raw": "HLA-A",
        "template_allele": "HLA-A*01:01:01:01", "template_distance": "5",
        "gene_start": 1000, "gene_end": 5000, "strand": "+",
        "consensus": "HLA-A*new", "alleles": "HLA-A*01:01:01:01,HLA-A*02:01:01:01",
        "template_warning": "partial_CDS", "cds_distance": None, "cds_mut": None,
    }
    t1, warn = extract.build_table1_row("p_depth1", "hap1", row)
    check("depth-1 row does not crash build_table1_row", t1 is not None)
    check("depth-1 row: is_novel True", t1["is_novel"] is True)
    check("depth-1 row: novelty_depth == 1", t1["novelty_depth"] == 1,
          detail=str(t1["novelty_depth"]))
    check("depth-1 row: novelty_class == 'undetermined', not None/mis-bucketed",
          t1["novelty_class"] == "undetermined", detail=str(t1["novelty_class"]))
    check("depth-1 row: no spurious 'unexpected novelty depth' warning attached",
          warn is None, detail=str(warn))


def test_c4_gene_name_stripping():
    for raw, exp_gene, exp_size in [("C4AL", "C4A", "L"), ("C4AS", "C4A", "S"),
                                     ("C4BL", "C4B", "L"), ("C4S", "C4", "S")]:
        gene, size = extract.clean_gene_name(raw)
        check(f"C4 strip {raw!r} -> gene={exp_gene!r}, size={exp_size!r}",
              gene == exp_gene and size == exp_size, detail=f"got ({gene!r}, {size!r})")
    gene, size = extract.clean_gene_name("HLA-A")
    check("non-C4 gene name untouched", gene == "HLA-A" and size is None)


def test_kir_absent_and_flagged_if_present():
    kir_present = [g for g in extract.GENE_CLASS if extract.GENE_CLASS[g] == "kir"]
    check("17 KIR genes classified", len(kir_present) == 17, detail=str(len(kir_present)))
    row = {
        "contig": "ctg1", "gene_id": "IAG1", "gene_name_raw": "KIR2DL1",
        "template_allele": "KIR2DL1*001", "template_distance": "0",
        "gene_start": 100, "gene_end": 200, "strand": "+",
        "consensus": "KIR2DL1*001", "alleles": "KIR2DL1*001",
        "template_warning": None, "cds_distance": None, "cds_mut": None,
    }
    t1, warn = extract.build_table1_row("p1", "hap1", row)
    check("KIR gene presence is flagged loudly via a warning", warn is not None and "KIR" in warn,
          detail=str(warn))


def test_copy_index_suffix():
    check("no suffix -> copy_index 1", extract.gene_id_copy_index("IAG123456") == ("IAG123456", 1))
    check("suffix .2 -> copy_index 2", extract.gene_id_copy_index("IAG123456.2") == ("IAG123456", 2))


def test_ancestry_case_normalization():
    check("cohorts module lowercases/uppercases ancestry_pred consistently",
          hasattr(cohorts, "normalize_ancestry"),
          detail="expected a normalize_ancestry() helper in 02_build_cohorts.py")
    if hasattr(cohorts, "normalize_ancestry"):
        check("normalize_ancestry('afr') == normalize_ancestry('AFR')",
              cohorts.normalize_ancestry("afr") == cohorts.normalize_ancestry("AFR"))


# ---------------------------------------------------------------------------
# Integration-level tests against the real fixture tree
# ---------------------------------------------------------------------------
def test_against_fixtures(fixtures_root):
    persons = extract.discover_persons(fixtures_root, limit=200)
    check("fixture discovery finds >0 person dirs", len(persons) > 0)

    multi_contig_found = False
    kir_found = False
    all_t1_rows = []
    for pid in persons[:80]:  # a meaningful chunk, not the whole 200, to keep this test fast
        t1, t2, n_unpairable, warnings = extract.process_person(
            pid, os.path.join(fixtures_root, pid, "immuannot_output"), strict_header_check=False)
        all_t1_rows.extend(t1)
        contigs_per_hap = {}
        for row in t1:
            contigs_per_hap.setdefault(row["hap"], set()).add(row["contig"])
        if any(len(c) > 1 for c in contigs_per_hap.values()):
            multi_contig_found = True
        if any(row["gene_class"] == "kir" for row in t1):
            kir_found = True
        for w in warnings:
            if "UNKNOWN GENE" in w:
                FAILURES.append(f"unexpected UNKNOWN GENE warning: {w}")

    check("multi-contig haplotypes produce multiple distinct contigs per hap (cis trap fixture)",
          multi_contig_found)
    check("KIR genes absent from fixture output (spec: must never appear)", not kir_found)
    check("row grain (person_id, hap, contig, gene, copy_index) unique across sampled rows",
          _grain_unique(all_t1_rows))

    # Regression (Fix 1): depth-1 ('undetermined') novelty calls are rare-but-real in the fixture
    # (~3% of novel variants per make_fixtures.py) and must survive the full parse -> classification
    # path without crashing, being dropped, or landing in the wrong novelty_class bucket.
    depth1_rows = [r for r in all_t1_rows if r["novelty_depth"] == 1]
    check("fixture produces at least one novelty_depth==1 row across the sampled persons "
          "(if this ever goes to 0, the fixture's own depth-1 generation may have changed shape)",
          len(depth1_rows) > 0, detail=f"sampled {len(all_t1_rows)} rows, 80 persons")
    check("every novelty_depth==1 row classifies as novelty_class=='undetermined', none None/"
          "mis-bucketed", all(r["novelty_class"] == "undetermined" for r in depth1_rows),
          detail=str({r["novelty_class"] for r in depth1_rows}))
    check("every novelty_depth==1 row still has is_novel True (not silently dropped from the "
          "novel-call set)", all(r["is_novel"] for r in depth1_rows))

    # Regression (Fix 2 fixture realism check): template_warning must be near-ubiquitous (~95%),
    # matching 00_recon_vm.py's real production measurement -- if this collapses back toward the
    # fixture's old ~4% rate, the token-aware-gate regression test downstream would stop being a
    # real test of anything.
    transcript_rows_with_warning = sum(1 for r in all_t1_rows if r["template_warning"])
    warn_rate = transcript_rows_with_warning / len(all_t1_rows) if all_t1_rows else 0
    check("fixture's template_warning presence rate is near the real ~95.3% production rate "
          "(checked >=0.85 to allow sampling noise over 80 persons)", warn_rate >= 0.85,
          detail=f"rate={warn_rate:.3f} ({transcript_rows_with_warning}/{len(all_t1_rows)})")

    # Table 2: pairing must respect contig, not just hap-file membership.
    total_t2 = 0
    total_unpairable = 0
    for pid in persons[:80]:
        t1, t2, n_unpairable, _ = extract.process_person(
            pid, os.path.join(fixtures_root, pid, "immuannot_output"), strict_header_check=False)
        total_t2 += len(t2)
        total_unpairable += n_unpairable
    check("Table 2 produced some cis pairs", total_t2 > 0)
    check("unpairable (cross-contig) cases were detected and counted, not silently dropped",
          total_unpairable > 0, detail=f"got {total_unpairable}")
    for pid in persons[:20]:
        t1, t2, _, _ = extract.process_person(
            pid, os.path.join(fixtures_root, pid, "immuannot_output"), strict_header_check=False)
        contig_of = {(r["hap"], r["gene"]): r["contig"] for r in t1}
        for pair_row in t2:
            gene_a, gene_b = pair_row["pair"].split("~")
            ga, gb = f"HLA-{gene_a}", f"HLA-{gene_b}"
            ca = contig_of.get((pair_row["hap"], ga))
            cb = contig_of.get((pair_row["hap"], gb))
            if ca is not None and cb is not None:
                check(f"cis pair {pair_row['person_id']}/{pair_row['hap']}/{pair_row['pair']} "
                      f"shares one contig", ca == cb == pair_row["contig"])


def _grain_unique(rows):
    seen = set()
    for r in rows:
        key = (r["person_id"], r["hap"], r["contig"], r["gene"], r["copy_index"])
        if key in seen:
            return False
        seen.add(key)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", default="/tmp/hla_fixtures")
    args = ap.parse_args()

    print("=== Unit tests ===")
    for fn in [test_unquoted_template_distance, test_cds_mut_not_split_naively,
               test_conditional_fields_absent_not_zero, test_novelty_depth_classification,
               test_novelty_depth1_full_row_no_crash_no_misbucket,
               test_c4_gene_name_stripping, test_kir_absent_and_flagged_if_present,
               test_copy_index_suffix, test_ancestry_case_normalization]:
        try:
            fn()
        except Exception:
            FAILURES.append(fn.__name__)
            traceback.print_exc()

    print(f"\n=== Integration tests (fixtures at {args.fixtures!r}) ===")
    if not os.path.isdir(args.fixtures):
        print(f"SKIPPED: {args.fixtures!r} not found -- run make_fixtures.py first.")
    else:
        try:
            test_against_fixtures(args.fixtures)
        except Exception:
            FAILURES.append("test_against_fixtures")
            traceback.print_exc()

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
