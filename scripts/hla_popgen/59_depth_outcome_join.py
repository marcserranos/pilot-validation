#!/usr/bin/env python3
"""Does low regional long-read depth explain the HLA/KIR calls we distrust? Where should a cutoff go?

Sprint S05, WS-C (SPRINT.md board). 58 measures each person's long-read depth over the MHC and the
KIR region relative to their own genome-wide depth (rel_depth; 1.0 = covered like the rest of the
genome). This script joins that number, per person, to four kinds of call that earlier QC already
suspects are technical, and asks whether thin coverage explains them:

  (a) KIR (rules from 51_kir_structure_qc.py, imported, never re-derived here)
      kir_incomplete           a haplotype that is not 51-"complete" (KIR calls on >1 contig, or one
                               contig without both ends KIR3DL3 + KIR3DL2, or no KIR call at all);
      kir_fw_missing           a framework gene (51.FRAMEWORK: 3DL3, 2DL4, 3DP1, 3DL2) missing on a
                               haplotype;
      kir_fw_missing_complete  the same, restricted to 51-complete haplotypes. 51 found most of these
                               are contracted haplotypes (biology), so depth should NOT explain them;
                               it is reported as an internal negative control, not used for the cutoff.
  (b) DQ (rules from 50_qc_foundation.py: 50.dq_group = the Petersdorf G1/G2 classifier on field 1;
      50e: a structurally anomalous same-contig DQ haplotype has DQA1-DQB1 spacing >= 100 kb)
      dq_split                 DQA1 and DQB1 both called on a haplotype but on different contigs;
      dq_anomalous             DQA1 and DQB1 on one contig with start-to-start spacing >= --dq-anomalous-kb;
      dq_cross_group_any       a haplotype pairing a G1 DQA1 with a G2 DQB1 (or the reverse), whatever
                               the contig layout (includes split haplotypes, whose cross-contig hap
                               label is arbitrary: 50e);
      dq_cross_group_physical  the same, but only when both genes sit on one contig: the 57b physical-
                               evidence rule (a DQA1~DQB1 pair counts as a physical haplotype only when
                               both genes are on one contig; trio mismatch < 7.7% for such pairs).
  (c) hla_absent               one of the always-present genes (--expected-genes, default A, B, C, DRB1,
                               DQA1, DQB1, DPA1, DPB1) absent from a haplotype. DRB3/4/5 are excluded:
                               their presence is real copy-number variation set by the DRB1 lineage.
      hla_hap_empty            a haplotype with no HLA call at all (the extreme case of the above).
  (d) MHC contig breaks (50b's mhc_contigs_per_hap: distinct contigs carrying copy-1 calls of
      50.MHC_GENES on one haplotype, built with 50.hap_gene_contigs)
      mhc_contigs_mean         continuous: the person's mean over their haplotypes;
      mhc_fragmented           binary: a haplotype with >= --frag-contigs contigs.
  Every outcome is 1/0 per person ("any haplotype"), blank when the person cannot be assessed (no
  KIR GTF pair; no Table 1 row; no haplotype with both DQ genes; no classifiable DQ pair).

Each outcome is judged against the depth of its own region (KIR outcomes: window_KIR; the rest:
window_MHC), at MAPQ >= 0 and MAPQ >= 20 (58's two floors). Windows only: single-gene depth is too
noisy for a person-level filter (58 smoke test: gene-level p5 0.4-0.6 in ordinary people).

Two depth scales, both reported everywhere (column `scale`):
  abs   rel_depth as 58 defines it (region depth / the person's own genome-wide depth);
  norm  rel_depth / the cohort median rel_depth of that region x MAPQ x platform (all people with
        depth, before the kinship filter; a platform with < 20 people uses the all-platform median).
        Regions run ~15% under genome depth even in ordinary people (58 smoke test: MHC median 0.825,
        p5-p95 0.76-0.91; KIR 0.86, 0.74-0.97), so 1.0 here = typical for the region. Per-platform
        medians keep a platform-wide depth offset from looking like a per-person defect.
  The cutoff rule is applied on `norm` by default (`--proposal-scale`), grid 0.50..0.95 step 0.05;
  `abs` keeps the grid 0.30..0.90 step 0.05.

Cohort: people in 58's INTERNAL depth table, joined to cohort_membership.tsv for ancestry and
platform (_s05_common.load_cohort). Unrelated only by default (_s05_common.greedy_unrelated, kin >=
0.0442, the repo-wide rule); `--include-related` skips it.

PUBLIC outputs (`--out-dir`, aggregate only; every participant count through _disclosure, so 1-19
shows as "<20" and every rate built on such a cell is blank):
    outcome_by_depth_bin.tsv   outcome rate by rel_depth bin, per ALL / ancestry / platform stratum,
                               for every outcome x MAPQ; mean/median for mhc_contigs_mean
    hla_gene_absent_by_bin.tsv per expected gene: share of haplotypes lacking it, by MHC depth bin
    cutoff_scan.tsv            per cutoff 0.30..0.90 step 0.05, ALL and per ancestry: people removed,
                               outcome rate kept vs removed, share of events removed, J
    logistic_or.tsv            OR of the outcome per 0.1 LOWER rel_depth, adjusted for platform and
                               ancestry (statsmodels Logit; Wald 95% CI)
    contigs_linear.tsv         change in mean MHC contigs per haplotype per 0.1 lower rel_depth (OLS)
    proposed_cutoff.tsv        the cutoff the stated rule picks, per outcome x MAPQ and per region
                               (rule text in every row)
    run_summary.tsv            inputs, people per stage, parameters, timing
INTERNAL (`--internal-dir`, never in the repo, never pulled):
    per_person_depth_outcome.tsv one row per person: rel_* and norm_* depth columns, ancestry,
                                 platform, every outcome
    exclude_people.tsv           the coverage filter, ready for 61/62: one row per person per region
                                 (window_MHC, window_KIR): person_id, region, mapq, scale, cutoff,
                                 value, below_cutoff, exclude. cutoff = that region's proposal on the
                                 primary scale (--proposal-scale) at --exclude-mapq (default 0);
                                 blank when no cutoff was proposed (nobody below). value = the
                                 person's depth on that scale. exclude (same on both of a person's
                                 rows) = below the cutoff in ANY region named by --exclude-region
                                 (default both). Covers every person with depth, relatives included
                                 (the cutoff is derived on the unrelated set, but 61/62 need a flag
                                 for everyone).
    kir_outcomes.tsv             per-person KIR parse cache (appended + fsynced; a rerun resumes)

Cutoff rule (CLI `--min-events-removed` X, `--max-people-removed` Y, `--min-ancestry-n` M). A
person is removed when their window depth (on the chosen scale) < c. Scan c upward over the
grid; propose the SMALLEST c that removes >= X of the outcome's events while removing <= Y of the people overall AND in every ancestry
with >= M people. If removing X of the events always costs more than Y of some ancestry, no cutoff
is proposed ("depth_does_not_explain"), and the row reports how much of the events the largest
admissible cutoff removes. Per region and MAPQ, the region proposal is the largest proposed cutoff
over that region's outcomes (excluding kir_fw_missing_complete, a negative control); since removal
only grows with c, that cutoff still satisfies the people limit in every ancestry. Smallest-c-first
keeps the filter as light as the rule allows.

Usage (VM, from the repo; pandas/numpy/statsmodels, all in the spechla pixi env):
    cd ~/repos/pilot-validation
    setsid nohup pixi run -e spechla -- python3 -u scripts/hla_popgen/59_depth_outcome_join.py \\
        --out-dir ~/s05/results/59 --internal-dir ~/s05/internal/59 --workers 8 \\
        > ~/s05/logs/59.log 2>&1 < /dev/null & disown
    Smoke test first: --people-limit 300 --out-dir ~/s05/results/59_smoke
        --internal-dir ~/s05/internal/59_smoke (own paths, ENVIRONMENT.md quirk #22).

Tests: scripts/hla_popgen/tests/test_59_depth_outcome_join.py (synthetic; no AoU data).
"""
import argparse
import importlib.util
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _disclosure as _disc  # noqa: E402
import _s05_common as C  # noqa: E402

DEFAULT_DEPTH = os.path.expanduser("~/s05/internal/58/per_person_depth.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/s05/results/59")
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s05/internal/59")

ANCESTRY_ORDER = C.ANCESTRY_ORDER
MAPQS = (0, 20)
REGIONS = ("window_MHC", "window_KIR")
EXPECTED_GENES = ("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1")
DEFAULT_BIN_EDGES = (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
DEFAULT_CUTOFFS = tuple(round(x, 2) for x in np.arange(0.30, 0.9001, 0.05))
# Cohort-normalised scale: rel_depth / median rel_depth of the same region x MAPQ x platform. The 58
# smoke test showed windows run ~15% under genome depth in ordinary people (MHC median 0.825, KIR
# 0.86), so 1.0 on this scale = "typical for this region on this platform".
DEFAULT_NORM_BIN_EDGES = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)
DEFAULT_NORM_CUTOFFS = tuple(round(x, 2) for x in np.arange(0.50, 0.9501, 0.05))
SCALES = ("norm", "abs")
MIN_PLATFORM_N_FOR_MEDIAN = 20  # smaller platform groups are normalised by the all-platform median
REL_CLIP = 1.5          # regression only: rel_depth above this is capped (duplication-like outliers
                        # say nothing about thin coverage and would dominate the slope)
OR_STEP = 0.1           # effects are per 0.1 LOWER rel_depth
MIN_EVENTS_FOR_MODEL = 5
MIN_LEVEL_EVENTS = 5    # covariate levels with fewer events or non-events are pooled into "other"

# outcome -> region whose depth it is judged against
OUTCOMES = [
    ("kir_incomplete", "window_KIR"),
    ("kir_fw_missing", "window_KIR"),
    ("kir_fw_missing_complete", "window_KIR"),
    ("hla_absent", "window_MHC"),
    ("hla_hap_empty", "window_MHC"),
    ("dq_split", "window_MHC"),
    ("dq_anomalous", "window_MHC"),
    ("dq_cross_group_any", "window_MHC"),
    ("dq_cross_group_physical", "window_MHC"),
    ("mhc_fragmented", "window_MHC"),
]
NEGATIVE_CONTROLS = {"kir_fw_missing_complete"}
CONTINUOUS = [("mhc_contigs_mean", "window_MHC")]


def log(msg):
    print(f"[59 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def depth_col(region, mapq, scale="abs"):
    return f"{'rel' if scale == 'abs' else 'norm'}_{region}_q{mapq}"


def add_normalised(pp):
    """norm_<region>_q<mapq> = rel / median(rel) over people of the same platform (all people with
    depth, before the kinship filter). Platforms with < MIN_PLATFORM_N_FOR_MEDIAN people use the
    all-platform median. Returns (pp, medians list for run_summary)."""
    meds = []
    for region in REGIONS:
        for q in MAPQS:
            col = depth_col(region, q)
            allmed = float(pp[col].median())
            norm = pd.Series(np.nan, index=pp.index)
            for plat, g in pp.groupby("platform"):
                v = g[col].dropna()
                med = float(v.median()) if len(v) >= MIN_PLATFORM_N_FOR_MEDIAN else allmed
                norm[g.index] = g[col] / med if med and np.isfinite(med) else np.nan
                meds.append((region, q, plat, round(med, 4), len(v) >= MIN_PLATFORM_N_FOR_MEDIAN))
            pp[depth_col(region, q, "norm")] = norm.astype(float)
    return pp, meds


# ---------------------------------------------------------------------------
# Sibling modules (importlib): each rule lives in exactly one script.
# ---------------------------------------------------------------------------
_mods = {}


def _load_module(filename, modname):
    if modname in _mods:
        return _mods[modname]
    spec = importlib.util.spec_from_file_location(modname, os.path.join(_THIS_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    _mods[modname] = mod
    return mod


def m50():
    return _load_module("50_qc_foundation.py", "qc_foundation_for_59")


def m51():
    return _load_module("51_kir_structure_qc.py", "kir_structure_qc_for_59")


def m58():
    return _load_module("58_regional_depth.py", "regional_depth_for_59")


# ---------------------------------------------------------------------------
# Lock (a second launch against the same internal dir fails loudly: quirk #23)
# ---------------------------------------------------------------------------
def acquire_lock(internal_dir, name="59.lock"):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, name)
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"FATAL: another run (pid {pid}) holds {lock}")
        except (ValueError, ProcessLookupError):
            log(f"removing stale lock {lock}")
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def load_depth(path, people_limit=None):
    """58's INTERNAL per_person_depth.tsv (columns = 58.PER_PERSON_COLS) -> one row per person with
    rel_window_MHC_q0/_q20 and rel_window_KIR_q0/_q20. `people_limit` keeps the first N person_ids
    in sorted order (deterministic smoke test)."""
    cols = m58().PER_PERSON_COLS
    d = pd.read_csv(path, sep="\t", comment="#", dtype={"person_id": str})
    missing = set(cols) - set(d.columns)
    if missing:
        raise SystemExit(f"FATAL: {path} lacks 58's columns {sorted(missing)}")
    d = d[d["region"].isin(REGIONS) & d["mapq"].isin(MAPQS)]
    d = d.drop_duplicates(["person_id", "region", "mapq"], keep="last")
    pids = sorted(d["person_id"].unique())
    if people_limit:
        pids = pids[:people_limit]
        d = d[d["person_id"].isin(set(pids))]
    w = d.pivot_table(index="person_id", columns=["region", "mapq"], values="rel_depth",
                      aggfunc="first")
    out = pd.DataFrame({"person_id": pids})
    for region in REGIONS:
        for q in MAPQS:
            key = (region, q)
            vals = w[key] if key in w.columns else pd.Series(dtype=float)
            out[depth_col(region, q)] = out["person_id"].map(vals).astype(float)
    return out


def attach_cohort(depth, cohort_path):
    coh = C.load_cohort(cohort_path)[["person_id", "ancestry_pred", "platform"]]
    coh = coh.drop_duplicates("person_id").rename(columns={"ancestry_pred": "ancestry"})
    pp = depth.merge(coh, on="person_id", how="left")
    pp["ancestry"] = pp["ancestry"].where(pp["ancestry"].isin(ANCESTRY_ORDER), "UNK")
    pp["platform"] = pp["platform"].fillna("unknown").replace({"nan": "unknown"})
    return pp


# ---------------------------------------------------------------------------
# Outcomes
# ---------------------------------------------------------------------------
def hla_outcomes(t1, pids, expected):
    """hla_absent / hla_hap_empty per person + per-gene absent-haplotype counts (absent_<gene>).
    A person with no Table 1 row at all is not assessable (NaN)."""
    in_t1 = set(t1["person_id"])
    genes_on = {k: set(g["gene"]) for k, g in t1.groupby(["person_id", "hap"], sort=False)}
    rows = []
    for pid in pids:
        if pid not in in_t1:
            rows.append({"person_id": pid})
            continue
        r = {"person_id": pid}
        n_abs = empty = 0
        per_gene = dict.fromkeys(expected, 0)
        for hap in ("hap1", "hap2"):
            genes = genes_on.get((pid, hap), set())
            empty += not genes
            for g in expected:
                if g not in genes:
                    per_gene[g] += 1
                    n_abs += 1
        r.update(hla_absent=int(n_abs > 0), hla_hap_empty=int(empty > 0), hla_n_absent_gene_haps=n_abs)
        r.update({f"absent_{g}": v for g, v in per_gene.items()})
        rows.append(r)
    cols = (["person_id", "hla_absent", "hla_hap_empty", "hla_n_absent_gene_haps"]
            + [f"absent_{g}" for g in expected])
    return pd.DataFrame(rows, columns=cols)


def dq_outcomes(t1, pids, anomalous_kb):
    """Per haplotype with DQA1 and DQB1 both called (copy 1): split / same-contig spacing / G1-G2
    group (50.dq_group). Person-level 'any haplotype' flags."""
    dq_group = m50().dq_group
    sub = t1[t1["gene"].isin(["DQA1", "DQB1"])]
    second_dqb1 = set(map(tuple, sub[(sub["gene"] == "DQB1") & (sub["copy_index"] >= 2)]
                          [["person_id", "hap"]].drop_duplicates().to_numpy()))
    c1 = sub[sub["copy_index"] == 1].drop_duplicates(["person_id", "hap", "gene"])
    by = {k: g.set_index("gene") for k, g in c1.groupby(["person_id", "hap"], sort=False)}
    rows = []
    for pid in pids:
        haps = [(h, by.get((pid, h))) for h in ("hap1", "hap2")]
        haps = [(h, g) for h, g in haps if g is not None and {"DQA1", "DQB1"} <= set(g.index)]
        r = {"person_id": pid}
        if not haps:
            rows.append(r)
            continue
        split = anom = cross_any = cross_phys = n_class = 0
        max_kb = np.nan
        for h, g in haps:
            a, b = g.loc["DQA1"], g.loc["DQB1"]
            is_split = a["contig"] != b["contig"]
            split += is_split
            if not is_split:
                kb = abs(float(b["gene_start"]) - float(a["gene_start"])) / 1000.0
                if np.isfinite(kb):
                    max_kb = kb if not np.isfinite(max_kb) else max(max_kb, kb)
                    anom += kb >= anomalous_kb
            ga, gb = dq_group("DQA1", a["consensus"]), dq_group("DQB1", b["consensus"])
            if ga and gb:
                n_class += 1
                if ga != gb:
                    cross_any += 1
                    cross_phys += not is_split
        r.update(dq_split=int(split > 0), dq_anomalous=int(anom > 0),
                 dq_max_same_contig_kb=round(max_kb, 1) if np.isfinite(max_kb) else np.nan,
                 dq_second_dqb1=int(any((pid, h) in second_dqb1 for h, _ in haps)))
        if n_class:
            r.update(dq_cross_group_any=int(cross_any > 0), dq_cross_group_physical=int(cross_phys > 0))
        rows.append(r)
    cols = ["person_id", "dq_split", "dq_anomalous", "dq_cross_group_any", "dq_cross_group_physical",
            "dq_max_same_contig_kb", "dq_second_dqb1"]
    return pd.DataFrame(rows, columns=cols)


def contig_outcomes(t1, pids, frag_contigs):
    """50.hap_gene_contigs -> distinct contigs carrying copy-1 MHC calls per haplotype."""
    hgc = m50().hap_gene_contigs(t1)
    n_by = {}
    for (pid, hap), gc in hgc.items():
        n_by.setdefault(pid, []).append(len(set(gc.values())))
    rows = []
    for pid in pids:
        v = n_by.get(pid)
        if not v:
            rows.append({"person_id": pid})
            continue
        rows.append({"person_id": pid, "mhc_contigs_mean": float(np.mean(v)),
                     "mhc_contigs_max": int(max(v)), "mhc_fragmented": int(max(v) >= frag_contigs)})
    return pd.DataFrame(rows, columns=["person_id", "mhc_contigs_mean", "mhc_contigs_max",
                                       "mhc_fragmented"])


KIR_COLS = ["person_id", "kir_n_haps", "kir_class_hap1", "kir_class_hap2", "kir_incomplete",
            "kir_fw_missing", "kir_fw_missing_complete", "kir_n_fw_missing"]


def kir_person(job):
    """51.person_records + 51.classify_hap + 51.FRAMEWORK for one person. Assessable only when both
    haplotype GTFs exist (an existing GTF with no KIR call is an 'empty' haplotype)."""
    pid, root = job
    k = m51()
    recs = k.person_records((pid, root))
    r = dict.fromkeys(KIR_COLS, "")
    r.update(person_id=pid, kir_n_haps=len(recs))
    if len(recs) != 2:
        return r
    n_miss = 0
    miss_complete = 0
    n_complete = 0
    incomplete = 0
    for _pid, hap, calls in recs:
        cls = k.classify_hap(calls)
        r[f"kir_class_{hap}"] = cls
        genes = {g for _, g, _, _ in calls}
        miss = [g for g in k.FRAMEWORK if g not in genes]
        n_miss += len(miss)
        incomplete += cls != "complete"
        if cls == "complete":
            n_complete += 1
            miss_complete += bool(miss)
    r.update(kir_incomplete=int(incomplete > 0), kir_fw_missing=int(n_miss > 0),
             kir_n_fw_missing=n_miss,
             kir_fw_missing_complete=int(miss_complete > 0) if n_complete else "")
    return r


def kir_outcomes(kir_outroot, pids, internal_dir, workers=1):
    """Parses each person's hap GTFs once; results appended + fsynced to the INTERNAL cache, so a
    rerun (after a crash or a VM stop) resumes."""
    cache = os.path.join(internal_dir, "kir_outcomes.tsv")
    done = set()
    if os.path.exists(cache):
        done = set(pd.read_csv(cache, sep="\t", comment="#", dtype={"person_id": str})["person_id"])
    todo = [p for p in pids if p not in done]
    log(f"KIR: {len(done)} cached, {len(todo)} to parse ({workers} workers)")
    m51()  # load in the parent before forking
    jobs = [(p, kir_outroot) for p in todo]
    t0 = time.time()
    if workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(workers) as ex:
            it = ex.map(kir_person, jobs, chunksize=64)
            _append_kir(cache, it, len(jobs), t0)
    else:
        _append_kir(cache, map(kir_person, jobs), len(jobs), t0)
    d = pd.read_csv(cache, sep="\t", comment="#", dtype={"person_id": str, "kir_class_hap1": str,
                                                         "kir_class_hap2": str})
    d = d.drop_duplicates("person_id", keep="last")
    return d[d["person_id"].isin(set(pids))]


def _append_kir(cache, results, n, t0, batch=500):
    new = not os.path.exists(cache)
    if new:
        _disc.assert_internal_path_allowed(cache)
    buf = []
    with open(cache, "a") as fh:
        if new:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n" + "\t".join(KIR_COLS) + "\n")
        for i, r in enumerate(results, 1):
            buf.append("\t".join(str(r[c]) for c in KIR_COLS))
            if len(buf) >= batch or i == n:
                fh.write("\n".join(buf) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
                buf = []
                log(f"  KIR parsed {i}/{n} ({i / max(1e-9, time.time() - t0):.1f}/s)")


# ---------------------------------------------------------------------------
# Disclosure helpers
# ---------------------------------------------------------------------------
def _small(*vals):
    return any(0 < int(v) < _disc.SUPPRESS_BELOW for v in vals)


def _rate_cells(prefix, n, d, mode):
    r = _disc.mask_rate(n, d, mode, pct_decimals=2)
    return {f"n_{prefix}": r["n"], f"pct_{prefix}": r["pct"], f"pct_{prefix}_lo": r["lo"],
            f"pct_{prefix}_hi": r["hi"]}


def _cnt(n, mode):
    return _disc.mask(int(n), _disc.PARTICIPANT, mode).display


# ---------------------------------------------------------------------------
# Aggregates (PUBLIC)
# ---------------------------------------------------------------------------
def bin_labels(edges):
    return [f"<{edges[0]:g}"] + [f"{a:g}-{b:g}" for a, b in zip(edges[:-1], edges[1:])] + \
           [f">={edges[-1]:g}"]


def assign_bins(x, edges):
    labels = bin_labels(edges)
    x = np.asarray(x, dtype=float)
    idx = np.searchsorted(np.asarray(edges, dtype=float), x, side="right")
    return [labels[i] if np.isfinite(v) else None for i, v in zip(idx, x)]


def _strata(d, with_platform=False):
    yield "ALL", "ALL", d
    for a in ANCESTRY_ORDER:
        s = d[d["ancestry"] == a]
        if len(s):
            yield "ancestry", a, s
    if with_platform:
        for p in sorted(d["platform"].unique()):
            yield "platform", p, d[d["platform"] == p]


def outcome_by_bin(pp, edges_by_scale, mode):
    rows = []
    for (outcome, region), scale in [(o, sc) for o in OUTCOMES + CONTINUOUS for sc in SCALES]:
        cont = (outcome, region) in CONTINUOUS
        labels = bin_labels(edges_by_scale[scale])
        edges = edges_by_scale[scale]
        for q in MAPQS:
            col = depth_col(region, q, scale)
            d = pp[pp[outcome].notna() & pp[col].notna()]
            if d.empty:
                continue
            d = d.assign(bin=assign_bins(d[col], edges))
            for stype, sval, s in _strata(d, with_platform=True):
                for lab in labels:
                    b = s[s["bin"] == lab]
                    n = len(b)
                    row = dict(outcome=outcome, region=region, scale=scale, mapq=q,
                               stratum_type=stype, stratum=sval, depth_bin=lab,
                               n_people=_cnt(n, mode))
                    if cont:
                        hide = n == 0 or (mode == _disc.PUBLIC and n < _disc.SUPPRESS_BELOW)
                        v = b[outcome].astype(float)
                        row.update(mean_value="" if hide else round(float(v.mean()), 3),
                                   median_value="" if hide else round(float(v.median()), 3))
                    else:
                        k = int(b[outcome].sum())
                        row.update(_rate_cells("events", k, n, mode))
                        if mode == _disc.PUBLIC and _small(n - k) and not _small(k):
                            # events + non-events = n: a 1-19 non-event count is a subtraction away
                            row.update(n_events="", pct_events="", pct_events_lo="", pct_events_hi="")
                    rows.append(row)
    return pd.DataFrame(rows)


def gene_absent_by_bin(pp, genes, edges_by_scale, mode, mapq=0):
    rows = []
    for scale in SCALES:
        rows += _gene_absent_one(pp, genes, edges_by_scale[scale], mode, mapq, scale)
    return pd.DataFrame(rows)


def _gene_absent_one(pp, genes, edges, mode, mapq, scale):
    col = depth_col("window_MHC", mapq, scale)
    d = pp[pp["hla_absent"].notna() & pp[col].notna()]
    rows = []
    if d.empty:
        return rows
    d = d.assign(bin=assign_bins(d[col], edges))
    for gene in genes:
        for lab in bin_labels(edges):
            b = d[d["bin"] == lab]
            n_haps, k = 2 * len(b), int(b[f"absent_{gene}"].sum())
            r = _disc.mask_rate(k, n_haps, mode, pct_decimals=2)
            row = dict(gene=gene, scale=scale, mapq=mapq, depth_bin=lab, n_people=_cnt(len(b), mode),
                       n_haps=r["d"], n_haps_absent=r["n"], pct_haps_absent=r["pct"],
                       pct_lo=r["lo"], pct_hi=r["hi"])
            if mode == _disc.PUBLIC and _small(n_haps - k) and not _small(k):
                row.update(n_haps_absent="", pct_haps_absent="", pct_lo="", pct_hi="")
            rows.append(row)
    return rows


def scan_exact(y, x, cutoffs):
    """Exact (unmasked) scan for one outcome x depth measure."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    n, e = len(y), float(y.sum())
    out = []
    for c in cutoffs:
        rem = x < c
        n_rem, e_rem = int(rem.sum()), float(y[rem].sum())
        fp = n_rem / n if n else float("nan")
        fe = e_rem / e if e else float("nan")
        out.append(dict(cutoff=c, n=n, n_removed=n_rem, events=int(e), events_removed=int(e_rem),
                        frac_people_removed=fp, frac_events_removed=fe,
                        J=(fe - fp) if e else float("nan")))
    return out


_SCAN_RENAME = {"pct_events_kept": "outcome_rate_kept_pct", "pct_events_kept_lo": "outcome_rate_kept_lo",
                "pct_events_kept_hi": "outcome_rate_kept_hi",
                "pct_events_removed": "outcome_rate_removed_pct",
                "pct_events_removed_lo": "outcome_rate_removed_lo",
                "pct_events_removed_hi": "outcome_rate_removed_hi",
                "pct_removed": "pct_people_removed", "pct_removed_lo": "pct_people_removed_lo",
                "pct_removed_hi": "pct_people_removed_hi"}


def _complementary_suppress(row, n, nr, e, er):
    """PUBLIC only. A row publishes a total and one part of it, so the other part is a subtraction
    away (kept = n - removed). When either half of a split is 1-19, every cell from which that half
    can be recovered is blanked too (a count that is really >= 20 is blanked, never shown as <20)."""
    tags = []
    if _small(nr, n - nr):
        tags.append("people")
        if not _small(nr):
            row["n_removed"] = ""
        for k in ("pct_people_removed", "pct_people_removed_lo", "pct_people_removed_hi"):
            row[k] = ""
    if _small(er, e - er) or (tags and e > 0):
        tags.append("events")
        for k in ("n_events_removed", "n_events_kept"):
            if row.get(k) not in ("<20", "0"):
                row[k] = ""
        for k in ("outcome_rate_kept_pct", "outcome_rate_kept_lo", "outcome_rate_kept_hi",
                  "outcome_rate_removed_pct", "outcome_rate_removed_lo", "outcome_rate_removed_hi",
                  "pct_of_events_removed", "J"):
            row[k] = ""
    row["complementary_suppression"] = ";".join(tags)


def cutoff_scan(pp, cutoffs_by_scale, mode):
    rows = []
    for (outcome, region), scale in [(o, sc) for o in OUTCOMES for sc in SCALES]:
        cutoffs = cutoffs_by_scale[scale]
        for q in MAPQS:
            col = depth_col(region, q, scale)
            d = pp[pp[outcome].notna() & pp[col].notna()]
            if d.empty:
                continue
            for _st, anc, s in _strata(d):
                y = s[outcome].to_numpy(dtype=float)
                x = s[col].to_numpy(dtype=float)
                for ex in scan_exact(y, x, cutoffs):
                    n, nr, e, er = ex["n"], ex["n_removed"], ex["events"], ex["events_removed"]
                    row = dict(outcome=outcome, region=region, scale=scale, mapq=q, ancestry=anc,
                               cutoff=ex["cutoff"], n_people=_cnt(n, mode))
                    for prefix, a, b in (("removed", nr, n), ("events_kept", e - er, n - nr),
                                         ("events_removed", er, nr)):
                        for k, v in _rate_cells(prefix, a, b, mode).items():
                            row[_SCAN_RENAME.get(k, k)] = v
                    fe = _disc.mask_rate(er, e, mode, pct_decimals=2)
                    row.update(n_events=_cnt(e, mode), pct_of_events_removed=fe["pct"])
                    blank = mode == _disc.PUBLIC and _small(n, nr, e, er, n - nr)
                    row["J"] = "" if (blank or not e) else round(ex["J"], 4)
                    row["complementary_suppression"] = ""
                    if mode == _disc.PUBLIC:
                        _complementary_suppress(row, n, nr, e, er)
                    rows.append(row)
    return pd.DataFrame(rows)


def rule_text(min_events, max_people, min_anc_n, scale="norm"):
    what = ("rel_depth / cohort median rel_depth of the region (same platform, same MAPQ)"
            if scale == "norm" else "rel_depth (region depth / own genome-wide depth)")
    return (f"scale={scale}: remove people with {what} < c; scan c upward; propose the smallest c that removes "
            f">= {100 * min_events:g}% of the outcome's events while removing <= {100 * max_people:g}% "
            f"of people overall and in every ancestry with >= {min_anc_n} people; region proposal = "
            f"largest proposed c over the region's outcomes (negative controls excluded)")


def propose_cutoff(y, x, anc, cutoffs, min_events, max_people, min_anc_n):
    """Applies the documented rule to exact numbers. Returns a dict (cutoff None when not proposed)."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    anc = np.asarray(anc)
    e = y.sum()
    base = dict(cutoff=None, status="", frac_events_removed=None, frac_people_removed=None,
                max_frac_people_removed_ancestry=None, binding_ancestry="",
                largest_admissible_cutoff=None)
    if e == 0:
        base["status"] = "no_events"
        return base
    groups = [a for a in ANCESTRY_ORDER if (anc == a).sum() >= min_anc_n]
    last_ok = None
    for c in sorted(cutoffs):
        rem = x < c
        fp = rem.mean()
        fpa = {a: rem[anc == a].mean() for a in groups}
        worst = max(fpa, key=fpa.get) if fpa else ""
        if fp > max_people or (fpa and fpa[worst] > max_people):
            base["binding_ancestry"] = worst if (fpa and fpa[worst] > max_people) else "ALL"
            break
        fe = y[rem].sum() / e
        last_ok = dict(c=c, fe=fe, fp=fp, fpa=fpa[worst] if fpa else None)
        if fe >= min_events:
            base.update(cutoff=c, status="proposed", frac_events_removed=fe, frac_people_removed=fp,
                        max_frac_people_removed_ancestry=last_ok["fpa"], binding_ancestry="",
                        largest_admissible_cutoff=None)
            return base
    if last_ok is None:
        base["status"] = "every_cutoff_removes_too_many_people"
        return base
    base.update(status="depth_does_not_explain", frac_events_removed=last_ok["fe"],
                frac_people_removed=last_ok["fp"], max_frac_people_removed_ancestry=last_ok["fpa"],
                largest_admissible_cutoff=last_ok["c"])
    return base


def proposed_table(pp, cutoffs_by_scale, min_events, max_people, min_anc_n, mode,
                   primary_scale="norm"):
    """Both scales are reported; `primary` marks the scale the rule is applied on by default."""
    rows = []
    for scale in SCALES:
        rows += _proposed_one(pp, cutoffs_by_scale[scale], min_events, max_people, min_anc_n, mode,
                              scale, scale == primary_scale)
    return pd.DataFrame(rows)


def _proposed_one(pp, cutoffs, min_events, max_people, min_anc_n, mode, scale, primary):
    rows = []
    text = rule_text(min_events, max_people, min_anc_n, scale)
    chosen = {}
    for outcome, region in OUTCOMES:
        for q in MAPQS:
            col = depth_col(region, q, scale)
            d = pp[pp[outcome].notna() & pp[col].notna()]
            row = dict(outcome=outcome, region=region, scale=scale, primary=primary, mapq=q, role=(
                "negative_control" if outcome in NEGATIVE_CONTROLS else "outcome"))
            if d.empty:
                row.update(status="not_assessable", rule=text)
                rows.append(row)
                continue
            y, x = d[outcome].to_numpy(float), d[col].to_numpy(float)
            p = propose_cutoff(y, x, d["ancestry"].to_numpy(), cutoffs, min_events, max_people,
                               min_anc_n)
            n, e = len(d), int(y.sum())
            if p["status"] == "proposed" and outcome not in NEGATIVE_CONTROLS:
                chosen.setdefault((region, q), []).append((p["cutoff"], outcome))
            at = p["cutoff"] if p["cutoff"] is not None else p["largest_admissible_cutoff"]
            small = mode == _disc.PUBLIC and _small(n, e)
            if mode == _disc.PUBLIC and at is not None:
                ex = scan_exact(y, x, [at])[0]  # fractions at c re-expose its kept/removed splits
                small = small or _small(ex["n_removed"], n - ex["n_removed"], ex["events_removed"],
                                        e - ex["events_removed"])
            row.update(status=p["status"],
                       proposed_cutoff="" if p["cutoff"] is None else p["cutoff"],
                       largest_admissible_cutoff=("" if p["largest_admissible_cutoff"] is None
                                                  else p["largest_admissible_cutoff"]),
                       binding_ancestry=p["binding_ancestry"],
                       n_people=_cnt(n, mode), n_events=_cnt(e, mode))
            for k in ("frac_events_removed", "frac_people_removed", "max_frac_people_removed_ancestry"):
                v = p.get(k)
                row[k] = "" if (v is None or small) else round(float(v), 4)
            row["rule"] = text
            rows.append(row)
    for region in REGIONS:
        for q in MAPQS:
            c = chosen.get((region, q), [])
            row = dict(outcome=f"REGION_{region}", region=region, scale=scale, primary=primary,
                       mapq=q, role="region_proposal", rule=text)
            if c:
                best = max(c)
                row.update(status="proposed", proposed_cutoff=best[0],
                           binding_ancestry="", n_people="", n_events="",
                           frac_events_removed="", frac_people_removed="",
                           max_frac_people_removed_ancestry="",
                           largest_admissible_cutoff=f"driven_by={best[1]}")
            else:
                row.update(status="no_outcome_explained_by_depth", proposed_cutoff="")
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Regression: outcome ~ rel_depth + platform + ancestry
# ---------------------------------------------------------------------------
def _norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _design(d, xcol, rel_clip, y=None):
    """Intercept, clipped rel_depth, treatment-coded platform and ancestry. For a binary y, levels
    with < MIN_LEVEL_EVENTS events or non-events are pooled into 'other' (a level with no events is
    perfectly separated and would stop the fit)."""
    x = np.minimum(d[xcol].to_numpy(dtype=float), rel_clip)
    cols, names = [np.ones(len(d)), x], ["const", "rel_depth"]
    for factor in ("platform", "ancestry"):
        f = d[factor].astype(str).to_numpy()
        if y is not None:
            f = f.copy()
            for lev in np.unique(f):
                m = f == lev
                ev = y[m].sum()
                if ev < MIN_LEVEL_EVENTS or (m.sum() - ev) < MIN_LEVEL_EVENTS:
                    f[m] = "other"
        levels = pd.Series(f).value_counts()
        if len(levels) < 2:
            continue
        ref = levels.index[0]
        for lev in sorted(levels.index):
            if lev != ref:
                cols.append((f == lev).astype(float))
                names.append(f"{factor}[{lev}]")
    return np.column_stack(cols), names


def logit_irls(X, y, max_iter=100, tol=1e-8):
    """Fallback Newton-Raphson logistic MLE (only used when statsmodels is missing)."""
    beta = np.zeros(X.shape[1])
    converged = False
    for _ in range(max_iter):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30)))
        H = X.T @ (X * (p * (1 - p))[:, None])
        g = X.T @ (y - p)
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, g, rcond=None)[0]
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            converged = True
            break
    p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30)))
    try:
        se = np.sqrt(np.clip(np.diag(np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None]))), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(X.shape[1], np.nan)
    if np.any(np.abs(beta) > 25):
        converged = False
    return beta, se, converged


def fit_logit(X, y):
    try:
        import statsmodels.api as sm
    except ImportError:
        b, s, c = logit_irls(X, y)
        return b, s, c, "numpy_irls"
    try:
        res = sm.Logit(y, X).fit(disp=0, maxiter=200)
        conv = bool(res.mle_retvals.get("converged", True))
        b, s = np.asarray(res.params), np.asarray(res.bse)
        if np.any(np.abs(b) > 25):
            conv = False
        return b, s, conv, "statsmodels"
    except Exception:  # perfect separation, singular matrix
        return np.full(X.shape[1], np.nan), np.full(X.shape[1], np.nan), False, "statsmodels"


def logistic_table(pp, mode, rel_clip=REL_CLIP):
    rows = []
    for (outcome, region), scale in [(o, sc) for o in OUTCOMES for sc in SCALES]:
        for q in MAPQS:
            col = depth_col(region, q, scale)
            d = pp[pp[outcome].notna() & pp[col].notna()]
            n = len(d)
            e = int(d[outcome].sum()) if n else 0
            row = dict(outcome=outcome, region=region, scale=scale, mapq=q, n_people=_cnt(n, mode),
                       n_events=_cnt(e, mode), or_per_0p1_lower=None, or_lo=None, or_hi=None,
                       p_value=None, covariates="", engine="", status="")
            if e < MIN_EVENTS_FOR_MODEL or (n - e) < MIN_EVENTS_FOR_MODEL:
                row["status"] = "too_few_events" if e < MIN_EVENTS_FOR_MODEL else "too_few_non_events"
                rows.append(row)
                continue
            if mode == _disc.PUBLIC and _small(e, n - e):
                row["status"] = "suppressed_lt20_events"
                rows.append(row)
                continue
            y = d[outcome].to_numpy(dtype=float)
            X, names = _design(d, col, rel_clip, y)
            beta, se, conv, used = fit_logit(X, y)
            b, s = beta[1], se[1]
            row.update(engine=used, covariates=";".join(names[2:]))
            if not conv or not np.isfinite(s) or s == 0:
                row["status"] = "did_not_converge"
                rows.append(row)
                continue
            lor, lse = -OR_STEP * b, OR_STEP * s
            row.update(or_per_0p1_lower=round(math.exp(lor), 4),
                       or_lo=round(math.exp(lor - 1.96 * lse), 4),
                       or_hi=round(math.exp(lor + 1.96 * lse), 4),
                       p_value=float(f"{2 * _norm_sf(abs(b / s)):.3g}"), status="ok")
            rows.append(row)
    return pd.DataFrame(rows)


def contigs_linear(pp, mode, rel_clip=REL_CLIP):
    """OLS: mhc_contigs_mean ~ rel_depth + platform + ancestry; effect per 0.1 LOWER rel_depth."""
    rows = []
    for (outcome, region), scale in [(o, sc) for o in CONTINUOUS for sc in SCALES]:
        for q in MAPQS:
            col = depth_col(region, q, scale)
            d = pp[pp[outcome].notna() & pp[col].notna()]
            n = len(d)
            row = dict(outcome=outcome, region=region, scale=scale, mapq=q, n_people=_cnt(n, mode),
                       delta_per_0p1_lower=None, lo=None, hi=None, p_value=None, covariates="",
                       status="")
            if n < _disc.SUPPRESS_BELOW:
                row["status"] = "too_few_people"
                rows.append(row)
                continue
            X, names = _design(d, col, rel_clip)
            y = d[outcome].to_numpy(dtype=float)
            try:
                import statsmodels.api as sm
                res = sm.OLS(y, X).fit()
                b, s = float(res.params[1]), float(res.bse[1])
            except ImportError:
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                resid = y - X @ beta
                sigma2 = resid @ resid / max(1, n - X.shape[1])
                b = float(beta[1])
                s = float(np.sqrt(sigma2 * np.linalg.inv(X.T @ X)[1, 1]))
            row.update(delta_per_0p1_lower=round(-OR_STEP * b, 4),
                       lo=round(-OR_STEP * b - 1.96 * OR_STEP * s, 4),
                       hi=round(-OR_STEP * b + 1.96 * OR_STEP * s, 4),
                       p_value=float(f"{2 * _norm_sf(abs(b / s)):.3g}") if s > 0 else None,
                       covariates=";".join(names[2:]), status="ok")
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def build_per_person(args, internal_dir):
    depth = load_depth(args.depth_table, args.people_limit)
    log(f"{len(depth)} people with regional depth")
    pp = attach_cohort(depth, args.cohort)
    pp, meds = add_normalised(pp)
    n_depth = len(pp)
    pp_all = pp[["person_id", "platform"] + [depth_col(r, q, sc) for r in REGIONS for q in MAPQS
                                             for sc in SCALES]].copy()
    if args.include_related:
        log("--include-related: no kinship filter")
    else:
        keep = C.greedy_unrelated(pp["person_id"], C.load_relatedness_pairs(args.relatedness_table))
        pp = pp[pp["person_id"].isin(keep)].reset_index(drop=True)
        log(f"unrelated: {len(pp)} of {n_depth}")
    pids = list(pp["person_id"])

    t1 = C.load_table1(args.table1)
    t1 = t1[t1["person_id"].isin(set(pids))].copy()
    t1["gene_start"] = pd.to_numeric(t1["gene_start"], errors="coerce")
    log(f"Table 1: {len(t1)} rows for {t1['person_id'].nunique()} people")
    genes = [g.strip() for g in args.expected_genes.split(",") if g.strip()]
    parts = [hla_outcomes(t1, pids, genes), dq_outcomes(t1, pids, args.dq_anomalous_kb),
             contig_outcomes(t1, pids, args.frag_contigs)]
    del t1
    parts.append(kir_outcomes(args.kir_outroot, pids, internal_dir, args.workers))
    for df in parts:
        pp = pp.merge(df, on="person_id", how="left")
    for outcome, _ in OUTCOMES + CONTINUOUS:
        pp[outcome] = pd.to_numeric(pp[outcome], errors="coerce")
    return pp, genes, n_depth, meds, pp_all


EXCLUDE_COLS = ["person_id", "region", "mapq", "scale", "cutoff", "value", "below_cutoff", "exclude"]


def exclude_people_table(pp_all, prop, scale, mapq=0, exclude_regions=REGIONS):
    """INTERNAL coverage filter. pp_all: every person with depth (normalised, before the kinship
    filter). prop: proposed_table() output. One row per person x region in REGIONS; `exclude` is
    True when the person is below the region proposal in any of `exclude_regions`."""
    reg = prop[(prop["role"] == "region_proposal") & (prop["scale"] == scale)
               & (prop["mapq"].astype(int) == int(mapq))]
    cut = {}
    for r in reg.itertuples(index=False):
        v = r.proposed_cutoff if r.status == "proposed" else ""
        cut[r.region] = float(v) if str(v) != "" else None
    rows = []
    below_any = pd.Series(False, index=pp_all.index)
    per_region = {}
    for region in REGIONS:
        col = depth_col(region, mapq, scale)
        val = pp_all[col] if col in pp_all.columns else pd.Series(np.nan, index=pp_all.index)
        c = cut.get(region)
        below = (val < c).fillna(False) if c is not None else pd.Series(False, index=pp_all.index)
        per_region[region] = (val, c, below)
        if region in exclude_regions:
            below_any |= below
    for region in REGIONS:
        val, c, below = per_region[region]
        rows.append(pd.DataFrame({
            "person_id": pp_all["person_id"].values, "region": region, "mapq": mapq, "scale": scale,
            "cutoff": "" if c is None else c, "value": val.round(4).values,
            "below_cutoff": below.astype(bool).values, "exclude": below_any.astype(bool).values}))
    out = pd.concat(rows, ignore_index=True)[EXCLUDE_COLS]
    return out.sort_values(["person_id", "region"]).reset_index(drop=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth-table", default=DEFAULT_DEPTH, help="58's INTERNAL per_person_depth.tsv")
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--cohort", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--include-related", action="store_true",
                    help="skip the unrelated filter (default: greedy unrelated at kin >= 0.0442)")
    ap.add_argument("--kir-outroot", default=C.DEFAULT_KIR_OUTROOT)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    ap.add_argument("--people-limit", type=int,
                    help="smoke test: first N person_ids of the depth table (sorted)")
    ap.add_argument("--workers", type=int, default=8, help="KIR GTF parsing processes")
    ap.add_argument("--expected-genes", default=",".join(EXPECTED_GENES),
                    help="genes expected on every haplotype (DRB3/4/5 deliberately excluded)")
    ap.add_argument("--dq-anomalous-kb", type=float, default=100.0,
                    help="same-contig DQA1-DQB1 spacing that marks a structurally anomalous contig "
                         "(50e: 98%% of same-contig haplotypes lie in 16.6-24.2 kb)")
    ap.add_argument("--frag-contigs", type=int, default=7,
                    help="mhc_fragmented = a haplotype with >= this many MHC contigs (cohort mean "
                         "4.3 Revio / 4.9 Sequel IIe in 50b)")
    ap.add_argument("--bins", default=",".join(f"{x:g}" for x in DEFAULT_BIN_EDGES),
                    help="bin edges on the absolute rel_depth scale")
    ap.add_argument("--cutoffs", default=",".join(f"{x:g}" for x in DEFAULT_CUTOFFS),
                    help="cutoff grid on the absolute rel_depth scale")
    ap.add_argument("--norm-bins", default=",".join(f"{x:g}" for x in DEFAULT_NORM_BIN_EDGES),
                    help="bin edges on the cohort-normalised scale")
    ap.add_argument("--norm-cutoffs", default=",".join(f"{x:g}" for x in DEFAULT_NORM_CUTOFFS),
                    help="cutoff grid on the cohort-normalised scale")
    ap.add_argument("--proposal-scale", choices=list(SCALES), default="norm",
                    help="scale the cutoff rule is primarily applied on (both are reported)")
    ap.add_argument("--min-events-removed", type=float, default=0.50,
                    help="X: a proposed cutoff must remove at least this fraction of the events")
    ap.add_argument("--max-people-removed", type=float, default=0.10,
                    help="Y: a proposed cutoff may remove at most this fraction of people, overall "
                         "and in every ancestry")
    ap.add_argument("--min-ancestry-n", type=int, default=100,
                    help="M: ancestries with fewer people are not used for the Y constraint")
    ap.add_argument("--rel-clip", type=float, default=REL_CLIP)
    ap.add_argument("--exclude-region", default=",".join(REGIONS),
                    help="regions whose proposal feeds `exclude` in exclude_people.tsv (comma list; "
                         "default both: a person is excluded if below in any)")
    ap.add_argument("--exclude-mapq", type=int, default=0, choices=list(MAPQS),
                    help="MAPQ floor of the region proposal used for exclude_people.tsv")
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args(argv)
    exclude_regions = tuple(r.strip() for r in args.exclude_region.split(",") if r.strip())
    bad = set(exclude_regions) - set(REGIONS)
    if bad:
        ap.error(f"--exclude-region: unknown region(s) {sorted(bad)}; choose from {list(REGIONS)}")

    internal_dir = os.path.expanduser(args.internal_dir)
    internal_path = os.path.join(internal_dir, "per_person_depth_outcome.tsv")
    _disc.assert_internal_path_allowed(internal_path)  # fail before any work, not after
    for p in (args.depth_table, args.table1, args.cohort):
        if not os.path.exists(p):
            raise SystemExit(f"FATAL: input not found: {p}")
    if not os.path.isdir(args.kir_outroot):
        raise SystemExit(f"FATAL: KIR outroot not found: {args.kir_outroot}")

    t0 = time.time()
    edges = {"abs": tuple(float(x) for x in args.bins.split(",")),
             "norm": tuple(float(x) for x in args.norm_bins.split(","))}
    cutoffs = {"abs": tuple(round(float(x), 4) for x in args.cutoffs.split(",")),
               "norm": tuple(round(float(x), 4) for x in args.norm_cutoffs.split(","))}
    lock = acquire_lock(internal_dir)
    try:
        pp, genes, n_depth, meds, pp_all = build_per_person(args, internal_dir)
        _disc.write_tsv(pp, internal_path, _disc.INTERNAL)

        mode = args.disclosure
        out = os.path.expanduser(args.out_dir)
        os.makedirs(out, exist_ok=True)
        _disc.write_tsv(outcome_by_bin(pp, edges, mode), os.path.join(out, "outcome_by_depth_bin.tsv"), mode)
        _disc.write_tsv(gene_absent_by_bin(pp, genes, edges, mode),
                        os.path.join(out, "hla_gene_absent_by_bin.tsv"), mode)
        _disc.write_tsv(cutoff_scan(pp, cutoffs, mode), os.path.join(out, "cutoff_scan.tsv"), mode)
        _disc.write_tsv(logistic_table(pp, mode, args.rel_clip), os.path.join(out, "logistic_or.tsv"), mode)
        _disc.write_tsv(contigs_linear(pp, mode, args.rel_clip), os.path.join(out, "contigs_linear.tsv"), mode)
        prop = proposed_table(pp, cutoffs, args.min_events_removed, args.max_people_removed,
                              args.min_ancestry_n, mode, args.proposal_scale)
        _disc.write_tsv(prop, os.path.join(out, "proposed_cutoff.tsv"), mode)
        excl = exclude_people_table(pp_all, prop, args.proposal_scale, args.exclude_mapq,
                                    exclude_regions)
        _disc.write_tsv(excl, os.path.join(internal_dir, "exclude_people.tsv"), _disc.INTERNAL)
        n_excl = int(excl.drop_duplicates("person_id")["exclude"].sum())
        log(f"exclude_people.tsv: {n_excl} of {excl['person_id'].nunique()} people excluded "
            f"(regions {','.join(exclude_regions)}, MAPQ>={args.exclude_mapq}, {args.proposal_scale})")

        run = [dict(metric="people_with_depth", value=_cnt(n_depth, mode)),
               dict(metric="people_analysed", value=_cnt(len(pp), mode)),
               dict(metric="kinship_filter", value="none" if args.include_related else "greedy_kin_ge_0.0442")]
        for outcome, _region in OUTCOMES:
            a = pp[outcome].notna()
            run.append(dict(metric=f"assessable_{outcome}", value=_cnt(int(a.sum()), mode)))
            run.append(dict(metric=f"events_{outcome}", value=_cnt(int(pp.loc[a, outcome].sum()), mode)))
        run += [dict(metric="expected_genes", value=",".join(genes)),
                dict(metric="dq_anomalous_kb", value=args.dq_anomalous_kb),
                dict(metric="frag_contigs", value=args.frag_contigs),
                dict(metric="bins_abs", value=args.bins), dict(metric="cutoffs_abs", value=args.cutoffs),
                dict(metric="bins_norm", value=args.norm_bins),
                dict(metric="cutoffs_norm", value=args.norm_cutoffs),
                dict(metric="proposal_scale", value=args.proposal_scale),
                dict(metric="exclude_regions", value=",".join(exclude_regions)),
                dict(metric="exclude_mapq", value=args.exclude_mapq),
                dict(metric="people_excluded_by_coverage", value=_cnt(n_excl, mode))]
        run += [dict(metric=f"cohort_median_rel_{r}_q{q}_{pl}", value=m if ok else f"{m} (all-platform)")
                for r, q, pl, m, ok in meds]
        run += [
                dict(metric="min_events_removed_X", value=args.min_events_removed),
                dict(metric="max_people_removed_Y", value=args.max_people_removed),
                dict(metric="min_ancestry_n_M", value=args.min_ancestry_n),
                dict(metric="rule", value=rule_text(args.min_events_removed, args.max_people_removed,
                                                    args.min_ancestry_n, args.proposal_scale)),
                dict(metric="rel_clip_regression", value=args.rel_clip),
                dict(metric="people_limit", value=args.people_limit or ""),
                dict(metric="disclosure", value=mode),
                dict(metric="wall_seconds", value=int(time.time() - t0))]
        _disc.write_tsv(pd.DataFrame(run), os.path.join(out, "run_summary.tsv"), mode)
        log(f"done: {len(pp)} people; tables in {out}")
    finally:
        os.remove(lock)
    return 0


if __name__ == "__main__":
    sys.exit(main())
