#!/usr/bin/env python3
"""WS5 follow-on (S03 call #8 sect 7) -- aggregate the full-cohort KIR Immuannot run.

Reads the per-person hap1/hap2.gtf.gz outputs a production-orchestrator run wrote under
--outroot (chr19 KIR window, --out-suffix .kir, Tier-3/self_align_needed people excluded
upstream -- see reports/hla_popgen/41_kir_scoping/README.md sect 4) and writes ONLY aggregate,
disclosure-safe TSVs to --out-dir. Reuses 41_kir_pilot.py's parsing/novelty/unrelatedness
functions via importlib (loaded from THIS script's own directory -- the orchestrator deploys
both files side by side under ~/s03/) instead of re-deriving them, so the novelty-tier and
framework-miss definitions never drift from the already-validated pilot logic.

Person discovery is directory-based: every subdirectory of --outroot IS a person the full run
attempted (whatever tier-filtering happened upstream), so this script does not need to
re-implement trim_tier filtering -- it just reports what it finds, including missing/corrupt
per-person output as its own QC row.

Disclosure (hard rule, matches 40b_novelty_rate_export.py's convention):
  - Any raw count 1-19 is written as the literal string "<20" (never as 0, never omitted).
  - Any rate whose numerator OR denominator is 1-19 is left blank (a disclosed denominator next
    to a rate can back-reveal a masked numerator) -- both the point estimate and its Wilson CI.
  - No person_id is ever written to an output file or printed to stdout/stderr. Progress lines
    are counts/percentages/timings only -- grep this file for `print(`/`log(` to audit.

Usage (VM, from ~/s03; imports 41_kir_pilot.py functions from this same directory):
  cd ~/s03 && PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    pixi run -e specimmune -- python3 43_kir_full_aggregate.py \\
      --outroot ~/pipeline_outputs_kir --out-dir ~/s03/results/43 --workers 32

Smoke test: --limit 50 restricts to the first 50 discovered person directories (sorted order,
so it is deterministic across re-runs, not a random sample -- fine for a code-path smoke test,
not for trusting the numbers).
"""
import argparse
import gzip
import importlib.util
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _disclosure as _disc

SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
# Set once, at the top of main(), from --disclosure. Module-level (not threaded through every
# call) so suppressed()/rate_row() keep their existing zero-argument call sites unchanged for
# every caller (44_kir_recurrence_saturation.py's m43().suppressed()/m43().rate_row() included) --
# see _disclosure.py / context/DECISIONS.md "two disclosure versions" for the policy this routes
# through. Defaults to PUBLIC, matching this script's behaviour before the disclosure layer
# existed.
_DISCLOSURE_MODE = _disc.PUBLIC

# Populated once per process (main process at import time; each worker via _worker_init) --
# module-level so worker functions (which multiprocessing pickles by reference, not closure) can
# reach 41_kir_pilot.py's functions without re-importing on every single task.
_kir = None


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def _load_kir_module(script_path=None):
    path = script_path or os.path.join(_THIS_DIR, "41_kir_pilot.py")
    spec = importlib.util.spec_from_file_location("kir_pilot_41", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _worker_init(script_path):
    global _kir
    _kir = _load_kir_module(script_path)


# ---------------------------------------------------------------------------
# Disclosure helpers -- one place, applied uniformly to every output table.
# ---------------------------------------------------------------------------
def suppressed(n, count_type=_disc.PARTICIPANT):
    """A true zero is not a disclosure risk and must stay '0' (this project's own
    feedback_suppressed_counts_are_not_zero.md: a masked <20 count silently collapsing to 0 is
    the bug to avoid -- the converse, a genuine 0 getting masked to '<20', would be the same class
    of information loss in the other direction and is equally wrong). Only 1-19 is masked --
    UNDER PUBLIC MODE and only for a disclosive count_type (default PARTICIPANT, this function's
    historical behaviour). Routes through _disclosure.mask() so the two-disclosure-versions
    policy (context/DECISIONS.md) lives in exactly one place; `count_type` lets a call site opt
    into the loosened exact-reporting rule (e.g. ALLELE_DISTINCT/RICHNESS/QC_TALLY) without
    duplicating the policy locally. `_DISCLOSURE_MODE` is set once in main() from --disclosure."""
    return _disc.mask(n, count_type, _DISCLOSURE_MODE).display


def wilson_ci(count, nobs, z=1.96):
    """Wilson score interval, scalar-friendly. Same formula as _viz_common.wilson_ci (05/16/17/
    24/38/40b's shared helper) reproduced here rather than imported, so this script's only
    dependency beyond stdlib is pandas/numpy -- _viz_common.py pulls in matplotlib, unneeded for
    an aggregate-only export and a heavier ask of the VM's pixi env than this script needs."""
    if nobs <= 0:
        return float("nan"), float("nan"), float("nan")
    p = count / nobs
    denom = 1 + z**2 / nobs
    center = (p + z**2 / (2 * nobs)) / denom
    half = (z * np.sqrt(p * (1 - p) / nobs + z**2 / (4 * nobs**2))) / denom
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return p, lo, hi


def rate_row(n, d, pct_decimals=1):
    """Returns (n_str, d_str, pct_str, lo_str, hi_str) with the hard disclosure rule applied
    UNDER PUBLIC mode: a rate is blanked when either n or d falls in the disclosive 1-19 band,
    OR when d is 0 (undefined). A genuine n==0 with d>=20 is NOT disclosive (see suppressed()'s
    docstring) and is reported as a real 0.0% rather than blanked. Under INTERNAL mode (set via
    --disclosure), returns exact n/d/pct/CI always (never blanked) -- see _disclosure.mask_rate().
    Routes through _disclosure.mask_rate() (count_type=RATE_NUM_DENOM); `_DISCLOSURE_MODE` is set
    once in main() from --disclosure."""
    r = _disc.mask_rate(n, d, _DISCLOSURE_MODE, pct_decimals=pct_decimals)
    return r["n"], r["d"], r["pct"], r["lo"], r["hi"]


# ---------------------------------------------------------------------------
# Step 1: discover people (directory-based -- ground truth of who the full run attempted).
# ---------------------------------------------------------------------------
def discover_people(outroot, limit=None):
    if not os.path.isdir(outroot):
        return []
    pids = sorted(d for d in os.listdir(outroot)
                  if os.path.isdir(os.path.join(outroot, d)))
    return pids[:limit] if limit else pids


# ---------------------------------------------------------------------------
# Step 2: per-person parse (runs in a worker process). Returns a compact, person_id-tagged dict
# that stays in-process memory only (never printed/written) -- merged into aggregate counters by
# the main process in Step 3, after which the person_id is discarded.
# ---------------------------------------------------------------------------
def parse_person(pid, outroot):
    """Standalone (no multiprocessing needed) -- also the function tests call directly against a
    synthetic --outroot. Requires _kir to already be loaded (module-level or via _worker_init)."""
    haps = []
    n_invalid = 0
    for hap in ("hap1", "hap2"):
        gtf_path = os.path.join(outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
        if not os.path.exists(gtf_path):
            n_invalid += 1
            continue
        try:
            rows = _kir.parse_hap_gtf(gtf_path)
        except (OSError, gzip.BadGzipFile, EOFError):
            n_invalid += 1
            continue
        kir_rows = [r for r in rows if r["gene"] in _kir.KIR_GENES]
        genes_called = sorted({r["gene"] for r in kir_rows})
        haps.append({"genes_called": genes_called, "kir_rows": kir_rows})
    return {"person_id": str(pid), "haps": haps, "n_invalid_haps": n_invalid}


def _parse_person_task(args):
    pid, outroot = args
    return parse_person(pid, outroot)


def parse_all(pids, outroot, workers, kir_script_path):
    """Multiprocessing map over parse_person. Falls back to serial (workers<=1) so a --limit
    smoke test or a single-core VM doesn't pay pool-startup overhead."""
    if workers <= 1:
        global _kir
        if _kir is None:
            _kir = _load_kir_module(kir_script_path)
        return [parse_person(pid, outroot) for pid in pids]
    results = []
    t0 = time.perf_counter()
    with mp.Pool(processes=workers, initializer=_worker_init,
                 initargs=(kir_script_path,)) as pool:
        for i, r in enumerate(pool.imap_unordered(
                _parse_person_task, [(p, outroot) for p in pids], chunksize=32), 1):
            results.append(r)
            if i % 1000 == 0 or i == len(pids):
                log(f"[43] parsed {i}/{len(pids)} people ({100.0*i/len(pids):.1f}%) "
                    f"in {time.perf_counter()-t0:.0f}s")
    return results


# ---------------------------------------------------------------------------
# Step 3: aggregation. One generic per-subset "bundle" pass, reused for the all/unrelated split
# in kir_run_summary and for each ancestry slice in kir_gene_by_ancestry, so the counting logic
# lives in exactly one place.
# ---------------------------------------------------------------------------
def compute_bundle(pids_subset, results_by_pid, kir):
    b = {
        "n_hap_attempted": 0, "n_hap_valid": 0, "n_hap_with_calls": 0,
        "genes_per_hap": [], "content": defaultdict(int),
        "gene_hap_presence": defaultdict(int),
        "gene_tier_counts": defaultdict(lambda: defaultdict(int)),
        "gene_known_allele_counts": defaultdict(lambda: defaultdict(int)),
        "gene_novel_protein_seqs": defaultdict(set),
        "framework_present": defaultdict(int),
        "framework_miss_class": defaultdict(int),
        "n_2dl2_and_2dl3": 0, "n_haps_2dl2_or_2dl3": 0,
        "n_3dl1_and_3ds1": 0, "n_haps_3dl1_or_3ds1": 0,
    }
    for pid in pids_subset:
        r = results_by_pid.get(pid)
        if not r:
            continue
        b["n_hap_attempted"] += 2
        for hap in r["haps"]:
            b["n_hap_valid"] += 1
            genes_called = hap["genes_called"]
            kir_rows = hap["kir_rows"]
            if not kir_rows:
                continue
            b["n_hap_with_calls"] += 1
            b["genes_per_hap"].append(len(genes_called))
            for g in genes_called:
                b["gene_hap_presence"][g] += 1
            for g in kir.FRAMEWORK_GENES:
                if g in genes_called:
                    b["framework_present"][g] += 1
                else:
                    b["framework_miss_class"][kir.classify_framework_miss(g, genes_called)] += 1
            content = "cB" if kir.B_CONTENT_GENES.intersection(genes_called) else "cA"
            b["content"][content] += 1
            has_2dl2, has_2dl3 = "KIR2DL2" in genes_called, "KIR2DL3" in genes_called
            has_3dl1, has_3ds1 = "KIR3DL1" in genes_called, "KIR3DS1" in genes_called
            if has_2dl2 or has_2dl3:
                b["n_haps_2dl2_or_2dl3"] += 1
            if has_2dl2 and has_2dl3:
                b["n_2dl2_and_2dl3"] += 1
            if has_3dl1 or has_3ds1:
                b["n_haps_3dl1_or_3ds1"] += 1
            if has_3dl1 and has_3ds1:
                b["n_3dl1_and_3ds1"] += 1
            for row in kir_rows:
                tier = row["novelty_tier"]
                b["gene_tier_counts"][row["gene"]][tier] += 1
                if tier == "known":
                    b["gene_known_allele_counts"][row["gene"]][row["consensus"]] += 1
                elif tier == "novel_protein":
                    b["gene_novel_protein_seqs"][row["gene"]].add(row["consensus"])
    return b


# ---------------------------------------------------------------------------
# Step 4: build each output table from bundle(s).
# ---------------------------------------------------------------------------
def build_run_summary(bundles):
    """bundles: {"all": bundle, "unrelated": bundle}."""
    rows = []
    for subset, b in bundles.items():
        n_novel = sum(sum(c for t, c in tiers.items() if t not in ("known", "undetermined"))
                      for tiers in b["gene_tier_counts"].values())
        n_undet = sum(tiers.get("undetermined", 0) for tiers in b["gene_tier_counts"].values())
        n_calls_total = sum(sum(tiers.values()) for tiers in b["gene_tier_counts"].values())
        n_ca, n_cb = suppressed(b["content"].get("cA", 0)), suppressed(b["content"].get("cB", 0))
        n_hap_str, d_hap_str, hap_pct, hap_lo, hap_hi = rate_row(
            b["n_hap_with_calls"], b["n_hap_valid"])
        n_nov_str, d_nov_str, nov_pct, nov_lo, nov_hi = rate_row(n_novel, n_calls_total)
        n_und_str, d_und_str, und_pct, und_lo, und_hi = rate_row(n_undet, n_calls_total)
        rows.append({
            "subset": subset,
            "n_haplotypes_attempted": b["n_hap_attempted"],
            "n_haplotypes_parsed_valid": b["n_hap_valid"],
            "n_haplotypes_with_ge1_kir_call": n_hap_str,
            "pct_haplotypes_with_ge1_kir_call": hap_pct,
            "pct_ci_lo": hap_lo, "pct_ci_hi": hap_hi,
            "mean_kir_genes_per_hap": (round(sum(b["genes_per_hap"]) / len(b["genes_per_hap"]), 2)
                                        if b["genes_per_hap"] else ""),
            "total_kir_calls": n_calls_total,
            "n_novel_calls": n_nov_str, "pct_novel": nov_pct,
            "novel_ci_lo": nov_lo, "novel_ci_hi": nov_hi,
            "n_undetermined_calls": n_und_str, "pct_undetermined": und_pct,
            "undetermined_ci_lo": und_lo, "undetermined_ci_hi": und_hi,
            "n_haplotypes_cA": n_ca, "n_haplotypes_cB": n_cb,
        })
    return pd.DataFrame(rows)


def build_gene_summary(bundle, kir):
    rows = []
    for gene in kir.KIR_GENES:
        tiers = bundle["gene_tier_counts"].get(gene, {})
        n_calls = sum(tiers.values())
        n_pres_str, d_pres_str, pres_pct, pres_lo, pres_hi = rate_row(
            bundle["gene_hap_presence"].get(gene, 0), bundle["n_hap_with_calls"])
        n_prot_str, d_prot_str, prot_pct, prot_lo, prot_hi = rate_row(
            tiers.get("novel_protein", 0), n_calls)
        known_pct = (round(100.0 * tiers.get("known", 0) / n_calls, 1) if n_calls else "")
        genomic_pct = (round(100.0 * tiers.get("novel_genomic_known_cds", 0) / n_calls, 1)
                       if n_calls else "")
        syn_pct = (round(100.0 * tiers.get("novel_cds_synonymous", 0) / n_calls, 1)
                   if n_calls else "")
        n_distinct_known = len(bundle["gene_known_allele_counts"].get(gene, {}))
        n_distinct_novel_prot = len(bundle["gene_novel_protein_seqs"].get(gene, set()))
        rows.append({
            "gene": gene,
            "n_haplotypes_carrying": n_pres_str, "n_haplotypes_total": d_pres_str,
            "presence_pct": pres_pct, "presence_ci_lo": pres_lo, "presence_ci_hi": pres_hi,
            "n_calls": n_calls,
            "pct_known": known_pct,
            "pct_novel_genomic_only": genomic_pct,
            "pct_novel_cds_synonymous": syn_pct,
            "pct_novel_protein": prot_pct,
            "novel_protein_ci_lo": prot_lo, "novel_protein_ci_hi": prot_hi,
            # allele/richness counts, not participant counts -- exact in BOTH disclosure modes
            # per AOU_SMALL_CELL_POLICY.md's loosened rule (reference/AOU_SMALL_CELL_POLICY.md
            # table: "Per-gene allele richness (distinct-allele counts)" -> allowed exact).
            "n_distinct_known_alleles": _disc.mask(n_distinct_known, _disc.RICHNESS,
                                                    _DISCLOSURE_MODE).display,
            "n_distinct_novel_protein_seqs": _disc.mask(n_distinct_novel_prot, _disc.RICHNESS,
                                                         _DISCLOSURE_MODE).display,
        })
    return pd.DataFrame(rows)


def build_gene_by_ancestry(bundles_by_ancestry, kir):
    rows = []
    for gene in kir.KIR_GENES:
        for anc in ANCESTRY_ORDER:
            b = bundles_by_ancestry.get(anc)
            if b is None:
                continue
            tiers = b["gene_tier_counts"].get(gene, {})
            n_calls = sum(tiers.values())
            n_pres_str, d_pres_str, pres_pct, pres_lo, pres_hi = rate_row(
                b["gene_hap_presence"].get(gene, 0), b["n_hap_with_calls"])
            n_prot_str, d_prot_str, prot_pct, prot_lo, prot_hi = rate_row(
                tiers.get("novel_protein", 0), n_calls)
            rows.append({
                "gene": gene, "ancestry": anc,
                "n_haplotypes_carrying": n_pres_str, "n_haplotypes_total": d_pres_str,
                "presence_pct": pres_pct, "presence_ci_lo": pres_lo, "presence_ci_hi": pres_hi,
                "n_novel_protein": n_prot_str, "n_calls": d_prot_str,
                "novel_protein_pct": prot_pct,
                "novel_protein_ci_lo": prot_lo, "novel_protein_ci_hi": prot_hi,
            })
    return pd.DataFrame(rows)


def build_allele_freq(bundle, kir):
    """count_type=CARRIER_NAMED_ALLELE (a named allele next to its carrier count) -- masked in
    PUBLIC mode regardless of _DISCLOSURE_MODE's general looseness elsewhere
    (AOU_SMALL_CELL_POLICY.md answer (3): "no, without an [RAB] exception"), so PUBLIC keeps the
    original pool-under-20-into-'other' behaviour (never prints a named allele next to a <20
    count, not even as the masked string). INTERNAL mode instead lists every allele individually,
    exact, with an lt20 flag column -- that is the whole point of the VM-only INTERNAL view."""
    rows = []
    for gene in kir.KIR_GENES:
        counts = bundle["gene_known_allele_counts"].get(gene, {})
        total_known = sum(counts.values())
        if _DISCLOSURE_MODE == _disc.INTERNAL:
            for allele, c in sorted(counts.items(), key=lambda kv: -kv[1]):
                r = _disc.mask_rate(c, total_known, _DISCLOSURE_MODE)
                rows.append({"gene": gene, "allele": allele, "n_haplotypes": r["n"],
                            "n_known_calls_total_gene": r["d"], "freq_pct": r["pct"],
                            "freq_ci_lo": r["lo"], "freq_ci_hi": r["hi"], "lt20": r["lt20_n"]})
            continue
        kept = {a: c for a, c in counts.items() if c >= SUPPRESS_BELOW}
        pooled = sum(c for a, c in counts.items() if c < SUPPRESS_BELOW)
        for allele, c in sorted(kept.items(), key=lambda kv: -kv[1]):
            n_str, d_str, pct, lo, hi = rate_row(c, total_known)
            rows.append({"gene": gene, "allele": allele, "n_haplotypes": n_str,
                        "n_known_calls_total_gene": d_str, "freq_pct": pct,
                        "freq_ci_lo": lo, "freq_ci_hi": hi})
        if pooled > 0:
            n_str, d_str, pct, lo, hi = rate_row(pooled, total_known)
            rows.append({"gene": gene, "allele": "other (<20 each)", "n_haplotypes": n_str,
                        "n_known_calls_total_gene": d_str, "freq_pct": pct,
                        "freq_ci_lo": lo, "freq_ci_hi": hi})
    return pd.DataFrame(rows)


def build_qc(bundle, kir, n_missing_or_corrupt, n_people_total):
    rows = []
    for g in kir.FRAMEWORK_GENES:
        n_str, d_str, pct, lo, hi = rate_row(bundle["framework_present"].get(g, 0),
                                             bundle["n_hap_with_calls"])
        rows.append({"metric": "framework_gene_presence", "item": g, "n": n_str, "d": d_str,
                    "pct": pct, "ci_lo": lo, "ci_hi": hi})
    for cls in ("flanked", "ambiguous_partial_flank", "edge_fragmented"):
        rows.append({"metric": "framework_miss_classification", "item": cls,
                    "n": suppressed(bundle["framework_miss_class"].get(cls, 0)),
                    "d": "", "pct": "", "ci_lo": "", "ci_hi": ""})
    n_str, d_str, pct, lo, hi = rate_row(bundle["n_2dl2_and_2dl3"], bundle["n_haps_2dl2_or_2dl3"])
    rows.append({"metric": "cooccurrence", "item": "KIR2DL2_and_KIR2DL3", "n": n_str, "d": d_str,
                "pct": pct, "ci_lo": lo, "ci_hi": hi})
    n_str, d_str, pct, lo, hi = rate_row(bundle["n_3dl1_and_3ds1"], bundle["n_haps_3dl1_or_3ds1"])
    rows.append({"metric": "cooccurrence", "item": "KIR3DL1_and_KIR3DS1", "n": n_str, "d": d_str,
                "pct": pct, "ci_lo": lo, "ci_hi": hi})
    n_zero = bundle["n_hap_valid"] - bundle["n_hap_with_calls"]
    n_str, d_str, pct, lo, hi = rate_row(n_zero, bundle["n_hap_valid"])
    rows.append({"metric": "haplotypes_zero_kir_calls", "item": "", "n": n_str, "d": d_str,
                "pct": pct, "ci_lo": lo, "ci_hi": hi})
    # QC tally (pipeline/assay outcome, no link to any phenotype/trait) -- exact allowed under
    # the loosened rule (AOU_SMALL_CELL_POLICY.md: "QC call tallies... no link to
    # phenotype/trait -> allowed exact"). Report n/d exact and derive pct/CI directly rather than
    # going through rate_row() (which is scoped to RATE_NUM_DENOM/participant-style rates).
    n_str = _disc.mask(n_missing_or_corrupt, _disc.QC_TALLY, _DISCLOSURE_MODE).display
    d_str = _disc.mask(n_people_total, _disc.QC_TALLY, _DISCLOSURE_MODE).display
    if n_people_total > 0:
        p, lo_f, hi_f = _disc.wilson_ci(n_missing_or_corrupt, n_people_total)
        pct, lo, hi = "%.1f" % (100 * p), "%.1f" % (100 * lo_f), "%.1f" % (100 * hi_f)
    else:
        pct, lo, hi = "", "", ""
    rows.append({"metric": "people_missing_or_corrupt_gtf", "item": "", "n": n_str, "d": d_str,
                "pct": pct, "ci_lo": lo, "ci_hi": hi})
    return pd.DataFrame(rows)


def build_content_by_ancestry(bundles_by_ancestry):
    rows = []
    for anc in ANCESTRY_ORDER:
        b = bundles_by_ancestry.get(anc)
        if b is None:
            continue
        n_ca, n_cb = b["content"].get("cA", 0), b["content"].get("cB", 0)
        n_str, d_str, pct, lo, hi = rate_row(n_ca, n_ca + n_cb)
        rows.append({"ancestry": anc, "n_cA": suppressed(n_ca), "n_cB": suppressed(n_cb),
                    "n_total": d_str, "pct_cA": pct, "pct_cA_ci_lo": lo, "pct_cA_ci_hi": hi})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5: ancestry/unrelatedness lookups -- reuses 41_kir_pilot.py's own logic verbatim.
# ---------------------------------------------------------------------------
def load_ancestry_map(cohort_path):
    cohort = pd.read_csv(cohort_path, sep="\t", dtype=str, low_memory=False)
    cohort["ancestry_pred"] = cohort["ancestry_pred"].fillna("UNASSIGNED").str.upper()
    return dict(zip(cohort["person_id"].astype(str), cohort["ancestry_pred"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs_kir"),
                    help="Where the full KIR orchestrator run wrote <person_id>/immuannot_output/"
                         "hap{1,2}.gtf.gz (one subdir per person attempted).")
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s03/results/43"))
    ap.add_argument("--cohort-membership",
                    default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table",
                    default=os.path.expanduser(
                        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/"
                        "aux/relatedness/samples_relatedness.tsv"))
    ap.add_argument("--kir-pilot-script", default=os.path.join(_THIS_DIR, "41_kir_pilot.py"),
                    help="Path to 41_kir_pilot.py -- functions are imported from here via "
                         "importlib, not re-derived.")
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--limit", type=int, default=None,
                    help="Smoke test: only process the first N discovered person directories "
                         "(sorted, deterministic).")
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args()

    global _DISCLOSURE_MODE
    _DISCLOSURE_MODE = args.disclosure
    if args.disclosure == _disc.INTERNAL and args.out_dir == ap.get_default("out_dir"):
        # caller didn't override --out-dir -- route to the INTERNAL default rather than the
        # PUBLIC default (~/s03/results/43), per _disclosure.default_out_dir().
        args.out_dir = _disc.default_out_dir(_disc.INTERNAL, "43", args.out_dir)
    if args.disclosure == _disc.INTERNAL:
        _disc.assert_internal_path_allowed(args.out_dir)

    os.makedirs(args.out_dir, exist_ok=True)

    global _kir
    _kir = _load_kir_module(args.kir_pilot_script)

    t0 = time.perf_counter()
    pids = discover_people(args.outroot, args.limit)
    n_people_total = len(pids)
    log(f"[43] discovered {n_people_total} person directories under {args.outroot!r}"
        + (f" (--limit {args.limit} applied)" if args.limit else ""))
    if n_people_total == 0:
        log("[43] nothing to do -- exiting.")
        return

    results = parse_all(pids, args.outroot, args.workers, args.kir_pilot_script)
    results_by_pid = {r["person_id"]: r for r in results}
    n_missing_or_corrupt = sum(1 for r in results if len(r["haps"]) == 0)
    n_people_with_outputs = n_people_total - n_missing_or_corrupt
    log(f"[43] {n_people_with_outputs}/{n_people_total} people have >=1 valid parsed haplotype "
        f"({100.0*n_people_with_outputs/n_people_total:.1f}%)")

    log("[43] loading ancestry map ...")
    ancestry_of = load_ancestry_map(args.cohort_membership)

    log("[43] computing unrelated subset ...")
    pairs = _kir.load_relatedness_pairs(args.relatedness_table)
    kept, removed = _kir.greedy_unrelated(pids, pairs, kin_min=_kir.KIN_MIN)
    unrelated_pids = [p for p in pids if p in kept]
    log(f"[43] unrelated subset: {len(unrelated_pids)}/{n_people_total} "
        f"({len(removed)} relatives dropped)")

    log("[43] aggregating (all vs unrelated) ...")
    bundle_all = compute_bundle(pids, results_by_pid, _kir)
    bundle_unrelated = compute_bundle(unrelated_pids, results_by_pid, _kir)

    log("[43] aggregating by ancestry (unrelated subset) ...")
    pids_by_ancestry = defaultdict(list)
    for p in unrelated_pids:
        anc = ancestry_of.get(p, "UNASSIGNED")
        if anc in ANCESTRY_ORDER:
            pids_by_ancestry[anc].append(p)
    bundles_by_ancestry = {a: compute_bundle(ps, results_by_pid, _kir)
                           for a, ps in pids_by_ancestry.items()}
    for a in ANCESTRY_ORDER:
        n = len(pids_by_ancestry.get(a, []))
        log(f"[43]   {a}: {n} people (unrelated)")

    log(f"[43] writing output tables (--disclosure {args.disclosure}) ...")
    _disc.write_tsv(build_run_summary({"all": bundle_all, "unrelated": bundle_unrelated}),
                     os.path.join(args.out_dir, "kir_run_summary.tsv"), _DISCLOSURE_MODE)
    _disc.write_tsv(build_gene_summary(bundle_unrelated, _kir),
                     os.path.join(args.out_dir, "kir_gene_summary.tsv"), _DISCLOSURE_MODE)
    _disc.write_tsv(build_gene_by_ancestry(bundles_by_ancestry, _kir),
                     os.path.join(args.out_dir, "kir_gene_by_ancestry.tsv"), _DISCLOSURE_MODE)
    _disc.write_tsv(build_allele_freq(bundle_unrelated, _kir),
                     os.path.join(args.out_dir, "kir_allele_freq.tsv"), _DISCLOSURE_MODE)
    _disc.write_tsv(build_qc(bundle_unrelated, _kir, n_missing_or_corrupt, n_people_total),
                     os.path.join(args.out_dir, "kir_qc.tsv"), _DISCLOSURE_MODE)
    _disc.write_tsv(build_content_by_ancestry(bundles_by_ancestry),
                     os.path.join(args.out_dir, "kir_content_by_ancestry.tsv"), _DISCLOSURE_MODE)

    log(f"[43] done in {time.perf_counter()-t0:.0f}s -- 6 tables written to {args.out_dir!r}")


if __name__ == "__main__":
    main()
