#!/usr/bin/env python3
"""62 -- Phased HLA class II haplotype units for fine-mapping (S05 WS-G).

One row per (person_id, hap, unit), unit in {DQA1~DQB1, DRB1~DQA1~DQB1}, with the allele strings
(2-field and full), the physical-evidence level of the phasing, the DQ G1/G2 class, a DQ
structural-anomaly flag, a regional coverage flag, platform, ancestry and the unrelated flag.

This per-person table is INTERNAL: it is written only under --internal-dir (default
~/s05/internal/62), refused anywhere inside the repo or any reports/ path
(_disclosure.assert_internal_path_allowed), and is meant for the in-perimeter collaborator bucket.
This script does NOT upload anything. The PUBLIC companion (--out-dir, default ~/s05/results/62)
holds only aggregates (masked) and the data dictionary.

Rules (defaults = the proposed S05 rules in context/DECISIONS.md, each behind a flag):
  Physical evidence (57b trio truth). A unit is physical only when every gene of the unit has a
  call on this haplotype, all of them (copy 1) sit on ONE assembly contig, and every gene pair
  inside the unit is "allowed" by --pair-rules. Pair statuses:
      allowed      same-contig pairs whose genes came from different parents in <7.7% of trio
                   children (class II pairs, DP, B~C)
      flagged      same contig, but a known elevated trio mismatch (DQB1~DPB1, ~21%)
      not_allowed  same contig is not enough (A~B, A~C, ~40%); needs statistical phasing / trios
  `evidence_level` of a unit:
      physical            one contig, all pairs allowed
      physical_flagged    one contig, >= 1 pair flagged, none not_allowed
      not_allowed_by_rule one contig, >= 1 pair not_allowed
      split_contigs       genes on >= 2 contigs: hap1/hap2 label carried across a contig break,
                          which 50e showed is not reliable (12.6% vs 42.4% expected cross-group)
      incomplete          a gene of the unit has no call on this haplotype
  Pairs not listed in --pair-rules get --unlisted-pair-status (default not_allowed).

  DQ G1/G2 (Petersdorf 2022, as 37/50e): G1 = DQA1*02-06 / DQB1*02-04, G2 = DQA1*01 / DQB1*05-06,
  field 1 only. dq_g1g2 = G1 | G2 | predicted_incompatible | unclassified | not_called (copy 1).

  DQ structural anomaly (50e follow-up): on same-contig haplotypes 98% of DQA1-DQB1 spacings lie in
  16.6-24.2 kb; the only same-contig cross-group haplotypes had spacing >= 100 kb and a second DQB1
  copy. dq_structure =
      anomalous_spacing   same contig and |gene_start(DQB1) - gene_start(DQA1)| >= --dq-spacing-max-kb
      anomalous_multicopy DQA1 or DQB1 has copy_index > 1 on this haplotype
      canonical           same contig, spacing below the cut, single copies
      not_assessable      DQA1/DQB1 missing or on different contigs
  (both anomalies -> anomalous_spacing;anomalous_multicopy)

  Regional coverage, two sources:
    (1) --exclude-people = 59's INTERNAL exclude_people.tsv (preferred). Its window_MHC rows give,
        per person, the depth value on 59's scale and 59's proposed cutoff for the MHC window.
        coverage_flag = low_coverage if value < cutoff, else pass (pass for everyone when 59
        proposed no MHC cutoff); coverage_cutoff_status = set_by_59. People absent from the file:
        not_measured. Only the MHC rows are used (59's `exclude` column also reflects the KIR
        window, which says nothing about class II units).
    (2) fallback, when --exclude-people is not given: 58's INTERNAL per_person_depth.tsv,
        --coverage-region (default window_MHC) at MAPQ >= --coverage-mapq (default 0), normalised
        like 59: rel_depth / the cohort median rel_depth of that region on the person's platform
        (all people in the depth table; platforms with < 20 people use the all-platform median).
        coverage_flag = low_coverage if normalised depth < --coverage-cutoff (default 0.6), pass
        otherwise, not_measured when the file is absent or the person is not in it;
        coverage_cutoff_status = provisional_until_59. --coverage-scale abs uses raw rel_depth.

  Unrelated: greedy removal of the highest-degree node at kin >= 0.0442 (_s05_common), over the
  people in Table 1. `not_assessed` when the relatedness table is missing (never silently True).

Inputs (VM, read-only): ~/pipeline_outputs/hla_calls_rich.tsv (Table 1: person_id, hap, contig,
gene, copy_index, consensus, gene_start, gene_end), ~/pipeline_outputs/cohort_membership.tsv,
the AoU relatedness table, ~/s05/internal/59/exclude_people.tsv (preferred) or
~/s05/internal/58/per_person_depth.tsv (fallback; both optional).

Outputs
  INTERNAL (--internal-dir): haplotype_units.tsv.gz (first line "# INTERNAL ..."), README.md (the
      data dictionary, identical to the repo copy), MANIFEST.tsv (bytes, rows, md5, sha256).
  PUBLIC (--out-dir): unit_freq_by_ancestry.tsv (2-field unit frequency per ancestry; unrelated,
      evidence `physical`, coverage not low by default; units named only with >= 20 haplotypes),
      unit_qc_summary.tsv (haplotype counts per evidence level / G1G2 / structure / coverage,
      masked), README.md (data dictionary).
  --dictionary-only PATH writes just the dictionary (used for reports/hla_popgen/62_*/README.md).

Upload (text only, a human runs it after review and Marc's consent; same bucket and pattern as 49):
  gsutil -m rsync -n -r ~/s05/internal/62/ \\
      gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/
  then the same without -n, then verify with -n -c.

Usage (VM, from the repo clone; see sprints/S05_qc_phasing_kir_trb/vm_parts/62.sh):
  python3 scripts/hla_popgen/62_cole_haplotype_units.py
  python3 scripts/hla_popgen/62_cole_haplotype_units.py --limit 200 \\
      --internal-dir ~/s05/internal/62_smoke --out-dir ~/s05/results/62_smoke   # smoke test

Tests: scripts/hla_popgen/tests/test_62_cole_haplotype_units.py (synthetic).
"""
import argparse
import gzip
import hashlib
import io
import os
import sys
import time
from collections import Counter

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _disclosure as _disc  # noqa: E402
import _s05_common as C  # noqa: E402

_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s05/internal/62")
DEFAULT_OUT_DIR = os.path.expanduser("~/s05/results/62")
DEFAULT_DEPTH = os.path.expanduser("~/s05/internal/58/per_person_depth.tsv")
REPO_DICTIONARY = os.path.join(_REPO_ROOT, "reports", "hla_popgen", "62_cole_haplotype_units",
                               "README.md")
SUPPRESS = 20

UNITS = (("DQA1~DQB1", ("DQA1", "DQB1")),
         ("DRB1~DQA1~DQB1", ("DRB1", "DQA1", "DQB1")))
UNIT_GENES = ("DRB1", "DQA1", "DQB1")
DEFAULT_PAIR_RULES = ("DQA1~DQB1=allowed,DRB1~DQA1=allowed,DRB1~DQB1=allowed,"
                      "DPA1~DPB1=allowed,B~C=allowed,DQB1~DPB1=flagged,DQA1~DPB1=flagged,"
                      "A~B=not_allowed,A~C=not_allowed")
PAIR_STATUSES = ("allowed", "flagged", "not_allowed")
EVIDENCE_LEVELS = ("physical", "physical_flagged", "not_allowed_by_rule", "split_contigs",
                   "incomplete")
COVERAGE_STATUS = "provisional_until_59"
COVERAGE_STATUS_59 = "set_by_59"
DEFAULT_EXCLUDE = os.path.expanduser("~/s05/internal/59/exclude_people.tsv")
MIN_PLATFORM_N_FOR_MEDIAN = 20  # same as 59

COLUMNS = [
    "person_id", "hap", "unit",
    "allele_DRB1", "allele_DQA1", "allele_DQB1",
    "two_field_DRB1", "two_field_DQA1", "two_field_DQB1",
    "unit_full", "unit_two_field",
    "n_contigs", "same_contig", "pair_evidence", "evidence_level",
    "dq_g1g2", "dq_spacing_kb", "dq_structure",
    "any_novel", "multi_copy",
    "mhc_rel_depth", "coverage_scale", "coverage_value", "coverage_flag", "coverage_cutoff",
    "coverage_cutoff_status",
    "platform", "ancestry_pred", "unrelated",
]

DICTIONARY = [
    ("person_id", "str", "AoU research id. INTERNAL table only; never in PUBLIC outputs."),
    ("hap", "str", "hap1 / hap2: which assembled haplotype (hifiasm, no trio/Hi-C)."),
    ("unit", "str", "DQA1~DQB1 or DRB1~DQA1~DQB1. Every (person, hap) present in Table 1 gets both rows."),
    ("allele_<GENE>", "str", "Full Immuannot consensus call of DRB1 / DQA1 / DQB1 on this hap (all copies, '/'-joined, copy order). Empty if no call or gene not in the unit. May contain 'new' (novel at that field)."),
    ("two_field_<GENE>", "str", "Same call truncated to 2 fields (e.g. DQB1*06:02). Empty when field 1 or 2 is novel/undetermined."),
    ("unit_full", "str", "'~'-joined full calls of the unit's genes (copy 1). Empty if incomplete."),
    ("unit_two_field", "str", "'~'-joined 2-field calls (copy 1). Empty if incomplete or any gene unresolved at 2 fields. Typing claims are supported at 2 fields (57b trio: ~4% Mendelian error, 12.6% at full resolution)."),
    ("n_contigs", "int", "Distinct assembly contigs carrying the unit's genes (copy 1) on this hap."),
    ("same_contig", "bool", "All unit genes (copy 1) on one contig."),
    ("pair_evidence", "str", "Per gene pair of the unit: 'GENE~GENE:<same_contig|split|missing>/<allowed|flagged|not_allowed>', ';'-joined."),
    ("evidence_level", "str", "physical | physical_flagged | not_allowed_by_rule | split_contigs | incomplete (rule in the script docstring). Use `physical` for haplotype-level claims."),
    ("dq_g1g2", "str", "Petersdorf 2022 DQ group of the copy-1 DQA1/DQB1 pair: G1 | G2 | predicted_incompatible | unclassified | not_called."),
    ("dq_spacing_kb", "float", "|gene_start(DQB1) - gene_start(DQA1)| in kb when both (copy 1) are on one contig; empty otherwise."),
    ("dq_structure", "str", "canonical | anomalous_spacing | anomalous_multicopy | anomalous_spacing;anomalous_multicopy | not_assessable."),
    ("any_novel", "bool", "Any call of the unit's genes contains 'new'."),
    ("multi_copy", "bool", "Any unit gene has copy_index > 1 on this hap."),
    ("mhc_rel_depth", "float", "Long-read depth over the MHC window / the person's genome-wide control depth (58, raw). Empty if 58's table was not read or the person is not in it."),
    ("coverage_scale", "str", "Scale of coverage_value and coverage_cutoff: norm (depth / cohort median of the MHC window on the person's platform; 1.0 = typical) or abs (raw rel_depth)."),
    ("coverage_value", "float", "The person's MHC-window depth on coverage_scale, as compared with the cutoff. Empty if not measured."),
    ("coverage_flag", "str", "pass | low_coverage | not_measured."),
    ("coverage_cutoff", "float", "Cutoff used for coverage_flag (on coverage_scale). Empty when 59 proposed no MHC cutoff."),
    ("coverage_cutoff_status", "str", "'set_by_59': cutoff taken from 59's exclude_people.tsv (its MHC-window proposal); 'provisional_until_59': fallback cutoff on 58's depth, a placeholder."),
    ("platform", "str", "revio / sequel2e (sequel2 people have no HLA calls yet). Platform affects novelty rates (50c); use as a covariate."),
    ("ancestry_pred", "str", "Predicted genetic ancestry (AoU), upper-cased."),
    ("unrelated", "str", "True / False from the greedy kin >= 0.0442 rule over HLA-called people; not_assessed when no relatedness table."),
]


def log(msg):
    print(f"[62 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------
def parse_pair_rules(text):
    """'A~B=status,...' -> {frozenset({A,B}): status}. Order-free pair keys."""
    rules = {}
    for item in filter(None, (x.strip() for x in str(text).split(","))):
        pair, _, status = item.partition("=")
        genes = pair.split("~")
        if len(genes) != 2 or status not in PAIR_STATUSES:
            raise ValueError(f"bad --pair-rules item {item!r}; expected GENE~GENE=" +
                             "|".join(PAIR_STATUSES))
        rules[frozenset(genes)] = status
    return rules


def pair_status(a, b, rules, unlisted="not_allowed"):
    return rules.get(frozenset((a, b)), unlisted)


def evidence_for_unit(genes, contig_of, rules, unlisted="not_allowed"):
    """genes: tuple of unit genes; contig_of: {gene: contig of copy 1}. Returns
    (evidence_level, pair_evidence str, n_contigs, same_contig)."""
    pairs = [(genes[i], genes[j]) for i in range(len(genes)) for j in range(i + 1, len(genes))]
    parts, statuses = [], []
    for a, b in pairs:
        st = pair_status(a, b, rules, unlisted)
        if a not in contig_of or b not in contig_of:
            where = "missing"
        elif contig_of[a] == contig_of[b]:
            where = "same_contig"
        else:
            where = "split"
        parts.append(f"{a}~{b}:{where}/{st}")
        statuses.append(st)
    present = [g for g in genes if g in contig_of]
    n_ctg = len({contig_of[g] for g in present})
    if len(present) < len(genes):
        return "incomplete", ";".join(parts), n_ctg, False
    if n_ctg > 1:
        return "split_contigs", ";".join(parts), n_ctg, False
    if "not_allowed" in statuses:
        lvl = "not_allowed_by_rule"
    elif "flagged" in statuses:
        lvl = "physical_flagged"
    else:
        lvl = "physical"
    return lvl, ";".join(parts), n_ctg, True


def dq_group(gene, consensus):
    """Petersdorf 2022, field 1 only (same rule as 37/50)."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    f1 = C.norm_gene(consensus).split("*", 1)[1].split(":")[0]
    if not f1.isdigit():
        return None
    f1 = int(f1)
    if gene == "DQA1":
        return "G2" if f1 == 1 else ("G1" if 2 <= f1 <= 6 else None)
    if gene == "DQB1":
        return "G1" if 2 <= f1 <= 4 else ("G2" if f1 in (5, 6) else None)
    return None


def classify_dq(dqa1, dqb1):
    """37's classify_pair vocabulary."""
    if dqa1 is None or dqb1 is None:
        return "not_called"
    ga, gb = dq_group("DQA1", dqa1), dq_group("DQB1", dqb1)
    if ga is None or gb is None:
        return "unclassified"
    return ga if ga == gb else "predicted_incompatible"


def dq_structure(rec_a, rec_b, multicopy, spacing_max_kb):
    """rec_*: dict(contig, start) of copy 1 or None. Returns (label, spacing_kb or None)."""
    flags = []
    spacing = None
    if rec_a is None or rec_b is None or rec_a["contig"] != rec_b["contig"]:
        if multicopy:
            return "anomalous_multicopy", None
        return "not_assessable", None
    if rec_a["start"] is not None and rec_b["start"] is not None:
        spacing = abs(rec_b["start"] - rec_a["start"]) / 1000.0
        if spacing >= spacing_max_kb:
            flags.append("anomalous_spacing")
    if multicopy:
        flags.append("anomalous_multicopy")
    return (";".join(flags) if flags else "canonical"), spacing


def coverage_flag(rel, cutoff):
    if rel is None or (isinstance(rel, float) and np.isnan(rel)):
        return "not_measured"
    return "low_coverage" if rel < cutoff else "pass"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def normalise_depth(depth, platform_of, min_n=MIN_PLATFORM_N_FOR_MEDIAN):
    """{pid: rel} -> {pid: rel / median rel of the person's platform}, 59's rule: medians over all
    people in the depth table; a platform with < min_n people (or unknown) uses the all-platform
    median."""
    if not depth:
        return {}
    s = pd.Series(depth, dtype=float)
    plat = s.index.map(lambda p: platform_of.get(p, "unknown"))
    allmed = float(s.median())
    out = {}
    for pl, g in s.groupby(plat):
        med = float(g.median()) if (len(g) >= min_n and pl != "unknown") else allmed
        for pid, v in g.items():
            out[pid] = v / med if med > 0 else float("nan")
    return out


def load_exclude_mhc(path, region="window_MHC"):
    """59's INTERNAL exclude_people.tsv -> ({pid: (value, cutoff or None)}, scale). Only the
    `region` rows are used."""
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str}, comment="#")
    need = {"person_id", "region", "scale", "cutoff", "value"}
    if not need <= set(df.columns):
        raise SystemExit(f"FATAL: {path} lacks 59's columns {sorted(need - set(df.columns))}")
    sub = df[df["region"] == region]
    if sub.empty:
        raise SystemExit(f"FATAL: {path} has no {region} rows")
    scale = str(sub["scale"].iloc[0])
    cut = pd.to_numeric(sub["cutoff"], errors="coerce")
    val = pd.to_numeric(sub["value"], errors="coerce")
    out = {p: (None if pd.isna(v) else float(v), None if pd.isna(c) else float(c))
           for p, v, c in zip(sub["person_id"], val, cut)}
    return out, scale


def coverage_info(pids, rel_raw, depth_available, platform_of, exclude_path=None,
                  cutoff=0.6, scale="norm", region="window_MHC"):
    """Per person: dict(mhc_rel_depth, coverage_scale, coverage_value, coverage_flag,
    coverage_cutoff, coverage_cutoff_status). 59's file wins when given; else the normalised
    58 fallback."""
    info = {}
    if exclude_path:
        ex, ex_scale = load_exclude_mhc(exclude_path, region)
        for p in pids:
            rel = rel_raw.get(p)
            base = {"mhc_rel_depth": "" if rel is None else round(rel, 4),
                    "coverage_scale": ex_scale, "coverage_cutoff_status": COVERAGE_STATUS_59}
            if p not in ex or ex[p][0] is None:
                info[p] = dict(base, coverage_value="", coverage_flag="not_measured",
                               coverage_cutoff="" if p not in ex or ex[p][1] is None else ex[p][1])
                continue
            v, c = ex[p]
            flag = "pass" if c is None else coverage_flag(v, c)
            info[p] = dict(base, coverage_value=round(v, 4), coverage_flag=flag,
                           coverage_cutoff="" if c is None else c)
        return info
    vals = rel_raw if scale == "abs" else normalise_depth(rel_raw, platform_of)
    for p in pids:
        rel = rel_raw.get(p) if depth_available else None
        v = vals.get(p) if depth_available else None
        info[p] = {"mhc_rel_depth": "" if rel is None else round(rel, 4),
                   "coverage_scale": scale,
                   "coverage_value": "" if v is None else round(v, 4),
                   "coverage_flag": coverage_flag(v, cutoff) if depth_available else "not_measured",
                   "coverage_cutoff": cutoff, "coverage_cutoff_status": COVERAGE_STATUS}
    return info


def load_depth(path, region="window_MHC", mapq=0):
    """58's INTERNAL per_person_depth.tsv -> ({person_id: rel_depth}, available)."""
    if not path or not os.path.exists(path):
        return {}, False
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str}, comment="#")
    if not {"person_id", "region", "mapq", "rel_depth"} <= set(df.columns):
        log("WARNING: depth table lacks expected columns; coverage_flag = not_measured")
        return {}, False
    sub = df[(df["region"] == region) & (pd.to_numeric(df["mapq"], errors="coerce") == mapq)]
    rel = pd.to_numeric(sub["rel_depth"], errors="coerce")
    return {p: float(v) for p, v in zip(sub["person_id"], rel) if pd.notna(v)}, True


def load_table1(path, limit=None):
    C.TABLE1_USECOLS = list(dict.fromkeys(C.TABLE1_USECOLS + ["gene_start", "gene_end"]))
    t1 = C.load_table1(path, limit_people=limit)
    t1["person_id"] = t1["person_id"].astype(str)
    return t1


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------
def build_units(t1, people_info, coverage, rules, unlisted="not_allowed", spacing_max_kb=100.0):
    """t1: Table 1 rows (gene names bare). people_info: {pid: dict(platform, ancestry_pred,
    unrelated)}. coverage: {pid: coverage_info() dict}. Returns the INTERNAL DataFrame (COLUMNS)."""
    all_haps = t1[["person_id", "hap"]].drop_duplicates().sort_values(["person_id", "hap"])
    t = t1[t1["gene"].isin(UNIT_GENES)].copy()
    t = t.sort_values(["person_id", "hap", "gene", "copy_index"])
    calls = {}
    for r in t.itertuples(index=False):
        start = getattr(r, "gene_start", None)
        start = None if start is None or pd.isna(start) else int(start)
        calls.setdefault((r.person_id, r.hap), {}).setdefault(r.gene, []).append(
            {"consensus": r.consensus, "contig": r.contig, "copy": int(r.copy_index),
             "start": start})
    rows = []
    for pid, hap in all_haps.itertuples(index=False):
        g = calls.get((pid, hap), {})
        first = {gene: recs[0] for gene, recs in g.items() if recs}
        contig_of = {gene: rec["contig"] for gene, rec in first.items()}
        info = people_info.get(pid, {})
        cov = coverage.get(pid) or {"mhc_rel_depth": "", "coverage_scale": "", "coverage_value": "",
                                    "coverage_flag": "not_measured", "coverage_cutoff": "",
                                    "coverage_cutoff_status": ""}
        dq_multi = any(len(g.get(x, [])) > 1 or any(c["copy"] > 1 for c in g.get(x, []))
                       for x in ("DQA1", "DQB1"))
        dq_struct, spacing = dq_structure(first.get("DQA1"), first.get("DQB1"), dq_multi,
                                          spacing_max_kb)
        g1g2 = classify_dq(first["DQA1"]["consensus"] if "DQA1" in first else None,
                           first["DQB1"]["consensus"] if "DQB1" in first else None)
        for unit, genes in UNITS:
            row = {"person_id": pid, "hap": hap, "unit": unit}
            for x in UNIT_GENES:
                recs = g.get(x, []) if x in genes else []
                row[f"allele_{x}"] = "/".join(str(c["consensus"]) for c in recs) or ""
                tf = [C.two_field(c["consensus"]) for c in recs]
                row[f"two_field_{x}"] = "/".join(tf) if recs and all(tf) else ""
            lvl, pe, n_ctg, same = evidence_for_unit(genes, contig_of, rules, unlisted)
            complete = all(x in first for x in genes)
            row["unit_full"] = "~".join(str(first[x]["consensus"]) for x in genes) if complete else ""
            tfs = [C.two_field(first[x]["consensus"]) for x in genes] if complete else []
            row["unit_two_field"] = "~".join(tfs) if complete and all(tfs) else ""
            row.update({
                "n_contigs": n_ctg, "same_contig": same, "pair_evidence": pe,
                "evidence_level": lvl, "dq_g1g2": g1g2,
                "dq_spacing_kb": "" if spacing is None else round(spacing, 2),
                "dq_structure": dq_struct,
                "any_novel": any("new" in str(c["consensus"]).lower()
                                 for x in genes for c in g.get(x, [])),
                "multi_copy": any(len(g.get(x, [])) > 1 or any(c["copy"] > 1 for c in g.get(x, []))
                                  for x in genes),
                **cov,
                "platform": info.get("platform", ""), "ancestry_pred": info.get("ancestry_pred", ""),
                "unrelated": info.get("unrelated", "not_assessed"),
            })
            rows.append(row)
    return pd.DataFrame(rows, columns=COLUMNS)


# ---------------------------------------------------------------------------
# PUBLIC aggregates
# ---------------------------------------------------------------------------
def suppress_partition(counts, floor=SUPPRESS):
    """Complementary suppression (same rule as 61): mask 1-19 cells, then the smallest shown
    non-zero cells until >= 2 cells are masked and they sum to >= floor."""
    masked = {c for c, v in counts.items() if 0 < v < floor}
    if not masked:
        return masked
    while len(masked) < 2 or sum(counts[c] for c in masked) < floor:
        rest = [(v, str(c), c) for c, v in counts.items() if c not in masked and v > 0]
        if not rest:
            break
        masked.add(min(rest)[2])
    return masked


def public_unit_freq(units, evidence=("physical",), exclude_low_coverage=True, unrelated_only=True):
    """Per unit x ancestry (+ALL): 2-field unit frequency among eligible haplotypes. A unit string
    is named only when >= 20 eligible haplotypes carry it; the rest pooled (masked if 1-19)."""
    d = units[units["evidence_level"].isin(evidence)]
    if unrelated_only:
        d = d[d["unrelated"].astype(str) == "True"]
    if exclude_low_coverage:
        d = d[d["coverage_flag"] != "low_coverage"]
    rows = []
    for unit, _ in UNITS:
        du = d[d["unit"] == unit]
        strata = [(a, du[du["ancestry_pred"] == a]) for a in C.ANCESTRY_ORDER] + [("ALL", du)]
        for anc, s in strata:
            n_elig = len(s)
            resolved = s[s["unit_two_field"] != ""]
            cnt = Counter(resolved["unit_two_field"])
            total = len(resolved)
            named = [(u, k) for u, k in cnt.most_common() if k >= SUPPRESS]
            rest = total - sum(k for _, k in named)
            while 0 < rest < SUPPRESS and named:  # else total - shown recovers a 1-19 count
                rest += named.pop()[1]
            for rank, (u, k) in enumerate(named, 1):
                r = C.rate_cells(k, total, "public", digits=2)
                rows.append({"unit": unit, "ancestry": anc, "rank": rank, "haplotype": u,
                             "n_haps": r["n"], "n_resolved_haps": r["d"], "freq_pct": r["pct"],
                             "freq_lo": r["lo"], "freq_hi": r["hi"]})
            r = C.rate_cells(rest, total, "public", digits=2)
            rows.append({"unit": unit, "ancestry": anc, "rank": "", "haplotype": "other",
                         "n_haps": r["n"], "n_resolved_haps": r["d"], "freq_pct": r["pct"],
                         "freq_lo": r["lo"], "freq_hi": r["hi"]})
            # share of eligible haplotypes not resolved at 2 fields; the eligible total is not
            # printed (with n_resolved_haps it would recover a masked unresolved count)
            r = C.rate_cells(n_elig - total, n_elig, "public", digits=2)
            rows.append({"unit": unit, "ancestry": anc, "rank": "",
                         "haplotype": "unresolved_at_2_fields", "n_haps": r["n"],
                         "n_resolved_haps": "", "freq_pct": r["pct"], "freq_lo": r["lo"],
                         "freq_hi": r["hi"]})
    out = pd.DataFrame(rows)
    out.attrs["filters"] = dict(evidence=list(evidence), exclude_low_coverage=exclude_low_coverage,
                                unrelated_only=unrelated_only)
    return out


def public_qc_summary(units):
    """Haplotype counts per (unit, field, value, ancestry), masked; plus % of the stratum.

    Complementary suppression, two passes: (1) within each (unit, field, ancestry), the values sum to
    the stratum (suppress_partition, as 61); (2) across rows, cells a reader can tie together are
    checked with _disclosure.suppress_equations: per (unit, field, value) the ancestries sum to ALL;
    dq_g1g2 `not_called` equals evidence `incomplete` in the DQA1~DQB1 unit; and fields that do not depend on the unit
    (dq_g1g2, dq_structure, coverage_flag, platform, unrelated) repeat the same counts in both units."""
    import _disclosure as _disc
    rows = []
    for unit, _ in UNITS:
        du = units[units["unit"] == unit]
        for field in ("evidence_level", "dq_g1g2", "dq_structure", "coverage_flag", "platform",
                      "unrelated"):
            for anc in C.ANCESTRY_ORDER + ["ALL"]:
                s = du if anc == "ALL" else du[du["ancestry_pred"] == anc]
                n = len(s)
                counts = {v: int(k) for v, k in s[field].astype(str).value_counts().sort_index().items()}
                sup = suppress_partition(counts)
                for val, k in counts.items():
                    r = C.rate_cells(k, n, "public", digits=2)
                    if val in sup and not (0 < k < SUPPRESS):
                        r = dict(r, n="<20*", pct="")
                    rows.append({"unit": unit, "field": field, "value": val, "ancestry": anc,
                                 "n_haps": r["n"], "of_haps": r["d"], "pct": r["pct"], "_raw_n_haps": k})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    eqs = _disc.partition_equations(df, ["unit", "field", "ancestry"], "value", None, ["n_haps"])
    eqs += _disc.partition_equations(df, ["unit", "field", "value"], "ancestry", "ALL", ["n_haps"],
                                     members=list(C.ANCESTRY_ORDER))
    idx = {(r.unit, r.field, r.value, r.ancestry): i for i, r in df.iterrows()}

    def cell(*k):
        return (idx[k], "n_haps") if k in idx else None
    u1 = UNITS[0][0]
    for anc in C.ANCESTRY_ORDER + ["ALL"]:
        a, b = cell(u1, "dq_g1g2", "not_called", anc), cell(u1, "evidence_level", "incomplete", anc)
        if a and b:
            eqs.append([a, b])
        for field in ("dq_g1g2", "dq_structure", "coverage_flag", "platform", "unrelated"):
            vals = set(df.loc[df["field"] == field, "value"])
            for v in vals:
                same = [c for c in (cell(u, field, v, anc) for u, _ in UNITS) if c]
                if len(same) > 1:
                    eqs.append(same)
    df = _disc.suppress_equations(df, eqs, {"n_haps": ["pct"]})
    return _disc.drop_raw(df)


# ---------------------------------------------------------------------------
# Dictionary / writing
# ---------------------------------------------------------------------------
def dictionary_markdown(params=None):
    p = params or {}
    lines = [
        "# 62 — Phased class II haplotype units: data dictionary",
        "",
        "Script: `scripts/hla_popgen/62_cole_haplotype_units.py`. This file is the data dictionary",
        "only; it contains no data. The per-person table (`haplotype_units.tsv.gz`) is INTERNAL:",
        "it lives on the VM under `~/s05/internal/62/` and in the in-perimeter collaborator bucket,",
        "never in this repo.",
        "",
        "## Grain",
        "One row per (`person_id`, `hap`, `unit`); `unit` is `DQA1~DQB1` or `DRB1~DQA1~DQB1`.",
        "Every haplotype present in HLA Table 1 gets both units, including incomplete ones.",
        "The file's first line is `# INTERNAL ...`; read with `comment.char = \"#\"` (R) or",
        "`comment=\"#\"` (pandas).",
        "",
        "## Columns",
        "| column | type | meaning |",
        "|---|---|---|",
    ]
    lines += [f"| `{c}` | {t} | {m} |" for c, t, m in DICTIONARY]
    lines += [
        "",
        "## Rules and parameters",
        f"- Pair rules (`--pair-rules`): `{p.get('pair_rules', DEFAULT_PAIR_RULES)}`; unlisted pairs:"
        f" `{p.get('unlisted', 'not_allowed')}`. Source: 57b trio truth (same-contig pairs whose",
        "  genes came from different parents: class II pairs and B~C <7.7%, DQB1~DPB1 21%,",
        "  A~B / A~C ~40%). Proposed S05 rule, not yet ratified (context/DECISIONS.md).",
        f"- DQ spacing cut (`--dq-spacing-max-kb`): {p.get('spacing', 100.0)} kb, measured between",
        "  `gene_start` coordinates on the shared contig (50e follow-up: 98% of same-contig",
        "  haplotypes lie in 16.6–24.2 kb).",
        "- Coverage, preferred source: 59's INTERNAL `exclude_people.tsv` (`--exclude-people`). Its",
        "  `window_MHC` rows give each person's MHC depth on 59's scale and 59's proposed MHC cutoff;",
        "  `coverage_flag = low_coverage` below it, status `set_by_59`. 59's `exclude` column (which",
        "  also reflects the KIR window) is not used here.",
        f"- Coverage, fallback (no 59 file): 58's depth over `{p.get('region', 'window_MHC')}` at"
        f" MAPQ ≥ {p.get('mapq', 0)}, scale `{p.get('scale', 'norm')}` (norm = rel_depth / the cohort",
        "  median on the person's platform, platforms with < 20 people use the all-platform median),",
        f"  cutoff {p.get('cutoff', 0.6)}, status `provisional_until_59`.",
        "- The status of the source actually used is in `coverage_cutoff_status` on every row.",
        "- Unrelated: greedy removal at kin ≥ 0.0442 among HLA-called people.",
        "",
        "## Recommended use",
        "- Haplotype-level association: `evidence_level == \"physical\"`,",
        "  `dq_structure == \"canonical\"`, `coverage_flag != \"low_coverage\"`, `unrelated == True`.",
        "- Use `unit_two_field` as the default allele resolution; full resolution needs read support.",
        "- Keep `platform` (and ancestry) as covariates.",
        "",
        "## Known limits",
        "- The 991 `sequel2` people (mostly AFR) have no HLA calls yet (Tier-3 trim; S05 WS-K).",
        "- Without 59's file the coverage cutoff is a placeholder (`provisional_until_59`).",
        "- `split_contigs` haplotypes keep their hap label, but the label is not reliable across a",
        "  contig break (50e).",
        "",
        "## Companion PUBLIC files (`~/s05/results/62/`)",
        "- `unit_freq_by_ancestry.tsv`: 2-field unit frequency by ancestry among unrelated, physical,",
        "  not-low-coverage haplotypes; units named only with ≥ 20 haplotypes (smaller ones pooled as",
        "  `other`), counts 1–19 shown `<20`.",
        "- `unit_qc_summary.tsv`: haplotype counts per evidence level, DQ class, DQ structure, coverage",
        "  flag, platform and unrelated flag, per ancestry. Counts 1–19 shown `<20`; a cell masked only",
        "  so that it cannot be recovered by subtraction from the total is shown `<20*`.",
        "",
        "## Upload (run by a person, after review; not done by the script)",
        "```",
        "gsutil -m rsync -n -r ~/s05/internal/62/ \\",
        "    gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/",
        "gsutil -m rsync -r ~/s05/internal/62/ \\",
        "    gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/",
        "gsutil -m rsync -n -r -c ~/s05/internal/62/ \\",
        "    gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/",
        "```",
        "",
    ]
    return "\n".join(lines)


def write_internal_gz(df, path):
    """gzip TSV with the INTERNAL header line; refuses repo / reports paths."""
    _disc.assert_internal_path_allowed(path)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    buf = io.StringIO()
    buf.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
    df.to_csv(buf, sep="\t", index=False)
    with gzip.open(path, "wt") as fh:
        fh.write(buf.getvalue())


def write_manifest(internal_dir, files):
    rows = []
    for name, nrows in files:
        full = os.path.join(internal_dir, name)
        data = open(full, "rb").read()
        rows.append({"file": name, "bytes": len(data), "rows": nrows,
                     "md5": hashlib.md5(data).hexdigest(),
                     "sha256": hashlib.sha256(data).hexdigest()})
    path = os.path.join(internal_dir, "MANIFEST.tsv")
    _disc.assert_internal_path_allowed(path)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def people_table(cohort, t1_people, relatedness_path):
    info = {}
    coh = cohort.set_index("person_id") if len(cohort) else cohort
    have_pairs = any(x and os.path.exists(x) for x in (relatedness_path, C.FALLBACK_RELATEDNESS))
    if have_pairs:  # load_relatedness_pairs falls back to the cached LR pair list itself
        unrelated = C.greedy_unrelated(t1_people, C.load_relatedness_pairs(relatedness_path))
        unrel = {p: str(p in unrelated) for p in t1_people}
    else:
        log("WARNING: relatedness table and cached pair list missing; unrelated = not_assessed")
        unrel = {p: "not_assessed" for p in t1_people}
    for p in t1_people:
        r = coh.loc[p] if p in coh.index else None
        info[p] = {"platform": "" if r is None else str(r.get("platform", "")),
                   "ancestry_pred": "" if r is None else str(r.get("ancestry_pred", "")),
                   "unrelated": unrel[p]}
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--depth-table", default=DEFAULT_DEPTH,
                    help="58 INTERNAL per_person_depth.tsv (absent -> coverage not_measured)")
    ap.add_argument("--coverage-region", default="window_MHC")
    ap.add_argument("--coverage-mapq", type=int, default=0)
    ap.add_argument("--exclude-people", default=None,
                    help="59's INTERNAL exclude_people.tsv (e.g. " + DEFAULT_EXCLUDE + "); when "
                         "given, coverage_flag comes from its window_MHC rows (set_by_59)")
    ap.add_argument("--coverage-scale", choices=["norm", "abs"], default="norm",
                    help="fallback only: norm = rel_depth / cohort median per platform")
    ap.add_argument("--coverage-cutoff", type=float, default=0.6,
                    help="fallback only: depth below this (on --coverage-scale) -> low_coverage "
                         "(PROVISIONAL until 59)")
    ap.add_argument("--pair-rules", default=DEFAULT_PAIR_RULES)
    ap.add_argument("--unlisted-pair-status", choices=PAIR_STATUSES, default="not_allowed")
    ap.add_argument("--dq-spacing-max-kb", type=float, default=100.0)
    ap.add_argument("--freq-evidence", default="physical",
                    help="comma list of evidence levels counted in unit_freq_by_ancestry")
    ap.add_argument("--freq-keep-low-coverage", action="store_true")
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--dictionary-only", metavar="PATH", nargs="?", const=REPO_DICTIONARY,
                    help="write only the data dictionary (default: the repo report README)")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)
    params = dict(pair_rules=args.pair_rules, unlisted=args.unlisted_pair_status,
                  spacing=args.dq_spacing_max_kb, region=args.coverage_region,
                  mapq=args.coverage_mapq, cutoff=args.coverage_cutoff, scale=args.coverage_scale)
    if args.dictionary_only:
        os.makedirs(os.path.dirname(os.path.abspath(args.dictionary_only)), exist_ok=True)
        with open(args.dictionary_only, "w") as fh:
            fh.write(dictionary_markdown(params))
        log(f"wrote {args.dictionary_only}")
        return 0

    _disc.assert_internal_path_allowed(args.internal_dir)  # fail before any heavy work
    rules = parse_pair_rules(args.pair_rules)
    t1 = load_table1(args.table1, args.limit)
    log(f"Table 1: {len(t1)} rows, {t1['person_id'].nunique()} people")
    cohort = C.load_cohort(args.cohort_membership)
    cohort["person_id"] = cohort["person_id"].astype(str)
    info = people_table(cohort, set(t1["person_id"]), args.relatedness_table)
    depth, avail = load_depth(args.depth_table, args.coverage_region, args.coverage_mapq)
    log(f"depth table available: {avail} ({len(depth)} people)")
    if args.exclude_people and not os.path.exists(args.exclude_people):
        raise SystemExit(f"FATAL: --exclude-people {args.exclude_people} not found")
    platform_of = dict(zip(cohort["person_id"], cohort["platform"]))
    cov = coverage_info(set(t1["person_id"]), depth, avail, platform_of, args.exclude_people,
                        args.coverage_cutoff, args.coverage_scale, args.coverage_region)
    log("coverage source: " + ("59 exclude_people.tsv" if args.exclude_people else
                               f"58 depth, {args.coverage_scale} < {args.coverage_cutoff} (provisional)"))
    units = build_units(t1, info, cov, rules, args.unlisted_pair_status, args.dq_spacing_max_kb)

    gz = os.path.join(args.internal_dir, "haplotype_units.tsv.gz")
    write_internal_gz(units, gz)
    readme_int = os.path.join(args.internal_dir, "README.md")
    _disc.assert_internal_path_allowed(readme_int)
    open(readme_int, "w").write(dictionary_markdown(params))
    write_manifest(args.internal_dir, [("haplotype_units.tsv.gz", len(units)), ("README.md", "")])
    log(f"INTERNAL: {len(units)} rows -> {gz}")

    ev = tuple(x for x in args.freq_evidence.split(",") if x)
    freq = public_unit_freq(units, ev, exclude_low_coverage=not args.freq_keep_low_coverage)
    os.makedirs(args.out_dir, exist_ok=True)
    _disc.write_tsv(freq, os.path.join(args.out_dir, "unit_freq_by_ancestry.tsv"), _disc.PUBLIC)
    _disc.write_tsv(public_qc_summary(units), os.path.join(args.out_dir, "unit_qc_summary.tsv"),
                    _disc.PUBLIC)
    open(os.path.join(args.out_dir, "README.md"), "w").write(dictionary_markdown(params))
    log(f"PUBLIC tables -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
