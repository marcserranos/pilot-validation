#!/usr/bin/env python3
"""Unit tests for 30_hla_structural_variation.py.

The whole scientific claim of that script is "we can tell a real deletion from a broken assembly".
So the tests are built around exactly that distinction:

  1. A gene missing from a contig that carries genes on BOTH sides is a deletion.
  2. The same gene missing from a contig that ends before it is NOT a deletion -- it is an
     unbridged absence, and must be counted separately. Getting this backwards would manufacture
     a deletion signal out of assembly fragmentation, and the resulting figure would look fine.
  3. Presence anywhere on the haplotype beats a bridged absence elsewhere (a gene split onto a
     second contig must not be reported as deleted on the first).
  4. The DR51/52/53 positive control recovers the textbook expectation, including the "no second
     DRB locus" groups, which are the easy ones to get wrong (absence is the correct answer).
  5. The KIR audit's verdict flips on a single KIR row -- it must not silently pass.

Run: python3 scripts/hla_popgen/tests/test_hla_structural_variation.py
"""
import importlib.util
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load_module("30_hla_structural_variation.py", "hla_sv")

FAILURES = []
# Canonical order stand-in: DRB1 < DRB3 < DQA1 < DQB1 < DPA1.
ORDER = {"DRB1": 0.0, "DRB3": 0.25, "DQA1": 0.5, "DQB1": 0.75, "DPA1": 1.0}
GOI = sorted(ORDER)


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def _t1(rows):
    df = pd.DataFrame(rows, columns=["person_id", "hap", "contig", "gene_bare", "copy_index",
                                     "consensus"])
    df["gene"] = "HLA-" + df["gene_bare"]
    return df


def _status(calls, pid, hap, gene):
    r = calls[(calls["person_id"] == pid) & (calls["hap"] == hap) & (calls["gene"] == gene)]
    return None if r.empty else r.iloc[0]["status"]


# ---------------------------------------------------------------------------
# bridged vs unbridged
# ---------------------------------------------------------------------------
def test_bridged_absence_is_called_a_deletion():
    """DRB1 and DQA1 are on one contig; DRB3 sits between them in canonical order and is absent.
    The assembly was built straight through DRB3's position, so DRB3 is deleted."""
    t1 = _t1([
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01"),
        ("p1", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*05:01"),
    ])
    calls, diag = m.bridged_absences(t1, ORDER, GOI)
    check("a gene absent between two present genes on one contig is called deleted",
          _status(calls, "p1", "hap1", "DRB3") == "deleted_bridged",
          str(_status(calls, "p1", "hap1", "DRB3")))
    check("...and is counted as a bridged absence", diag.get("bridged_absence", 0) >= 1, str(diag))


def test_unbridged_absence_is_not_a_deletion():
    """THE test. The contig carries only genes at the START of the region and stops. DQB1 and DPA1
    are past its end, so their absence says nothing about whether they exist on this haplotype.
    Calling them deleted would turn assembly fragmentation into a biological finding."""
    t1 = _t1([
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01"),
        ("p1", "hap1", "ctg1", "DRB3", 1, "HLA-DRB3*01:01"),
    ])
    calls, diag = m.bridged_absences(t1, ORDER, GOI)
    check("a gene past the end of the contig is NOT called deleted",
          _status(calls, "p1", "hap1", "DPA1") == "absent_unbridged",
          str(_status(calls, "p1", "hap1", "DPA1")))
    check("a gene between the contig's end and the region's end is also unbridged",
          _status(calls, "p1", "hap1", "DQB1") == "absent_unbridged",
          str(_status(calls, "p1", "hap1", "DQB1")))
    check("no bridged absence was recorded at all", diag.get("bridged_absence", 0) == 0, str(diag))


def test_single_gene_contig_bridges_nothing():
    t1 = _t1([("p1", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*05:01")])
    calls, diag = m.bridged_absences(t1, ORDER, GOI)
    check("one gene on a contig cannot bridge anything",
          diag.get("bridged_absence", 0) == 0, str(diag))
    check("everything else on that haplotype is unbridged",
          all(_status(calls, "p1", "hap1", g) == "absent_unbridged"
              for g in ["DRB1", "DRB3", "DQB1", "DPA1"]),
          str(calls.to_dict("records")))


def test_presence_on_a_second_contig_beats_bridged_absence():
    """A fragmented assembly can put DRB3 on its own contig while another contig spans DRB3's
    position. The gene is present on the haplotype; reporting it deleted would be wrong."""
    t1 = _t1([
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01"),
        ("p1", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*05:01"),
        ("p1", "hap1", "ctg2", "DRB3", 1, "HLA-DRB3*01:01"),
    ])
    calls, _ = m.bridged_absences(t1, ORDER, GOI)
    check("a gene present on any contig of the haplotype is 'present', not deleted",
          _status(calls, "p1", "hap1", "DRB3") == "present",
          str(_status(calls, "p1", "hap1", "DRB3")))


def test_haplotypes_are_independent():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01"),
        ("p1", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*05:01"),
        ("p1", "hap2", "ctg5", "DRB1", 1, "HLA-DRB1*15:01"),
        ("p1", "hap2", "ctg5", "DRB3", 1, "HLA-DRB3*01:01"),
        ("p1", "hap2", "ctg5", "DQA1", 1, "HLA-DQA1*01:02"),
    ])
    calls, _ = m.bridged_absences(t1, ORDER, GOI)
    check("hap1 shows the deletion", _status(calls, "p1", "hap1", "DRB3") == "deleted_bridged")
    check("hap2, which carries the gene, does not",
          _status(calls, "p1", "hap2", "DRB3") == "present")


# ---------------------------------------------------------------------------
# duplications
# ---------------------------------------------------------------------------
def test_duplication_counts_extra_copies_on_one_contig():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01"),
        ("p1", "hap1", "ctg1", "DRB1", 2, "HLA-DRB1*11:01"),
        ("p2", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*04:01"),
    ])
    dup = m.duplications(t1, ["DRB1"])
    r = dup[dup["gene"] == "DRB1"].iloc[0]
    check("a second copy on one contig is counted", int(r["n_multi_copy"]) == 1, str(r.to_dict()))
    check("the denominator is haplotype-contigs, not people",
          int(r["n_haplotype_contigs"]) == 2, str(r.to_dict()))
    check("the maximum copy index seen is reported", int(r["max_copy_seen"]) == 2)


def test_same_gene_on_two_contigs_is_not_a_duplication():
    """Two contigs each carrying one copy is a fragmented assembly, not a duplication."""
    t1 = _t1([
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01"),
        ("p1", "hap1", "ctg2", "DRB1", 1, "HLA-DRB1*03:01"),
    ])
    dup = m.duplications(t1, ["DRB1"])
    check("copies on separate contigs are not counted as multi-copy",
          int(dup.iloc[0]["n_multi_copy"]) == 0, str(dup.to_dict("records")))


# ---------------------------------------------------------------------------
# DR51/52/53 positive control
# ---------------------------------------------------------------------------
def test_drb_positive_control_recovers_expectation():
    t1 = _t1([
        # DR52 group: DRB1*03 expects DRB3, and has it -> concordant
        ("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01:01"),
        ("p1", "hap1", "ctg1", "DRB3", 1, "HLA-DRB3*01:01:02"),
        ("p1", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*05:01"),
        # DR51 group: DRB1*15 expects DRB5, but DRB5 is absent (bridged) -> discordant
        ("p2", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*15:01:01"),
        ("p2", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*01:02"),
        # DR1 group: DRB1*01 expects NOTHING, and has nothing -> concordant
        ("p3", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*01:01:01"),
        ("p3", "hap1", "ctg1", "DQA1", 1, "HLA-DQA1*01:01"),
    ])
    order = dict(ORDER, DRB5=0.3)
    calls, _ = m.bridged_absences(t1, order, sorted(order))
    detail, summ = m.drb_expectation_check(t1, calls)
    got = {r["drb1_group"]: r["concordant"] for r in detail.to_dict("records")}
    check("DRB1*03 with DRB3 present is concordant", got.get("03") is True, str(got))
    check("DRB1*15 without DRB5 is discordant", got.get("15") is False, str(got))
    check("DRB1*01 with no second DRB locus is concordant (absence is the right answer)",
          got.get("01") is True, str(got))
    check("the summary reports a percentage per group",
          "pct_concordant" in summ.columns, str(summ.columns.tolist()))


def test_first_field_parsing():
    check("first_field zero-pads a single digit", m.first_field("HLA-DRB1*3:01") == "03",
          str(m.first_field("HLA-DRB1*3:01")))
    check("first_field reads a two-digit group", m.first_field("HLA-DRB1*15:01:01") == "15")
    check("a field-1-novel call has no group", m.first_field("HLA-DRB1*new") is None)
    check("junk yields None", m.first_field("undetermined") is None)


# ---------------------------------------------------------------------------
# KIR audit (ask A5)
# ---------------------------------------------------------------------------
def test_kir_audit_expects_zero():
    t1 = _t1([("p1", "hap1", "ctg1", "DRB1", 1, "HLA-DRB1*03:01")])
    t1["gene_class"] = "classical_II"
    out = m.kir_audit(t1)
    check("no KIR rows gives the expected verdict", out["n_kir_rows"] == 0, str(out))
    check("...and the verdict explains why (chr19, outside the trim window)",
          "chr19" in out["verdict"], out["verdict"])


def test_kir_audit_flips_loudly_on_a_single_row():
    t1 = _t1([("p1", "hap1", "ctg1", "KIR2DL1", 1, "KIR2DL1*001")])
    t1["gene_class"] = "kir"
    out = m.kir_audit(t1)
    check("one KIR row is enough to flip the verdict", out["n_kir_rows"] == 1, str(out))
    check("...and the verdict says to stop, not to continue",
          out["verdict"].startswith("UNEXPECTED"), out["verdict"])
    check("...and names the gene that appeared", out["kir_genes_seen"] == ["KIR2DL1"], str(out))


# ---------------------------------------------------------------------------
# C4
# ---------------------------------------------------------------------------
def test_c4_copy_number_per_haplotype():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "C4A", 1, "C4A"),
        ("p1", "hap1", "ctg1", "C4B", 1, "C4B"),
        ("p2", "hap1", "ctg1", "C4A", 1, "C4A"),
        ("p2", "hap1", "ctg1", "C4A", 2, "C4A"),
    ])
    t1["c4_size"] = ["L", "S", "L", "L"]
    dist, size = m.c4_copy_number(t1)
    rows = {(int(r["C4A"]), int(r["C4B"])): r for r in dist.to_dict("records")}
    check("a 1A/1B haplotype is counted", (1, 1) in rows, str(rows.keys()))
    check("a 2A/0B haplotype is counted", (2, 0) in rows, str(rows.keys()))
    check("long/short composition is reported per gene",
          not size.empty and set(size["c4_size"]) <= {"L", "S"}, str(size.to_dict("records")))


def test_suppression_rule():
    check("zero is disclosable", m.suppress(0) == "0")
    check("1-19 is suppressed", m.suppress(19) == "<20")
    check("20 is disclosed", m.suppress(20) == "20")


def main():
    print("test_hla_structural_variation.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
