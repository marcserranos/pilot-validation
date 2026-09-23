#!/usr/bin/env python3
"""Per-ancestry allele discovery / saturation curves (Pakistan Genome Resource, Nature 2026,
Fig 3e style; Marc's proposal, endorsed by Cole on the 2026-09-22 call -- see
sprints/CALL_SUMMARY_2026-09-22.md Sec 6).

QUESTION
--------
As we add more people to each ancestry group, how many DISTINCT HLA protein (2-field) alleles
have we found so far, and how far from saturated is each ancestry? Cole's expectation from the
call: the AFR-ancestry curve is the LEAST saturated (still climbing steeply at current N), because
AFR haplotype diversity in AoU outstrips what the reference catalogue and this cohort's other
ancestry groups show. This script tests that formally rather than eyeballing it.

METHOD
------
1. Reuse 24_novelty_by_field.py's own pipeline (imported as a sibling module, never re-derived) to
   build the fully-labeled per-(person, hap, gene) call table: `add_call_labels` (field_class:
   known / f2_protein novel-protein-altering / f3_synonymous / f4_noncoding), `load_refdata` +
   `match_sequences` + `classify_sequence` (seq_class: known vs novel_protein, needs the actual
   translated protein/CDS, not just the consensus string, because a "new" field in the consensus
   name does not by itself distinguish two DIFFERENT novel proteins seen in two different people),
   and `allele_ids` (assigns each clean call a `prot_id`: the 2-field IPD-IMGT name for known
   alleles, or `<gene>_prot_<sha8>` for a novel protein -- this project's one and only definition
   of "distinct novel protein allele", also used by 34_novel_protein_recurrence.py).
2. `build_people` (24's own cohort machinery, in turn reused by 29's LD script) gives, per person,
   the AoU-predicted ancestry (`anc_pred`, plain argmax over p_<ANC> -- the "AoU predicted
   ancestry" the task asked for as the primary scheme) and a probability-thresholded version
   (`anc_strict`, threshold configurable; run twice, at the project default 0.9 and at 0.95 for
   the stricter sensitivity check the call asked for) plus the greedy-unrelated flag (kinship <
   0.0442, third-degree-or-closer removed).
3. Sampling UNIT = person (both haplotypes), per the task spec -- not haplotype (04/29's unit).
   For each ancestry and each gene group (pooled classical genes A/B/C/DRB1/DQA1/DQB1/DPA1/DPB1,
   and each classical gene individually for the supplement grid), build one allele-identity set
   per person (`prot_id`s seen on either haplotype, restricted to `keep_clean` calls). 25 random
   orderings of the people (matching the Pakistan Fig 3e convention exactly -- 25 subsamples per
   point) are then walked once each, and at every step (= subsample of size N) we record, for the
   *cumulative population added so far*:
     - number of distinct alleles seen at all (>= 1 carrier)      -- "distinct protein alleles"
     - ... restricted to alleles with field_class == f2_protein AND seq_class == novel_protein
                                                                    -- "distinct NOVEL alleles"
     - number of alleles that have reached >= 1 / >= 2 / >= 10 carriers (all-alleles; Pakistan
       Fig 3e's own three carrier thresholds)
   Mean + 2.5/97.5 percentile bands across the 25 permutations are what gets plotted and exported.
4. Extrapolation: fit the Clench (1979) / Michaelis-Menten saturating curve S(N) = Smax*N/(b+N) to
   each ancestry's mean >=1-carrier curve, reusing `fit_clench_asymptote`/`clench_curve` from
   04_allele_saturation.py verbatim (same estimator already used and documented for haplotype-unit
   curves; the derivative dS/dN = Smax*b/(b+N)^2 is the closed-form marginal discovery rate this
   module cites for "expected new alleles per additional person"). Expected new alleles per next
   1,000 PEOPLE = clench_curve(N_now+1000) - clench_curve(N_now). Chao2 (incidence) is reported
   alongside as a second, independent richness estimator -- 04's module-level note on why two
   different estimator families are always reported side by side, never silently preferred.
5. Equal-N comparison: N_min = the largest N such that every "well-powered" ancestry (>=
   MIN_PEOPLE_PER_ANCESTRY people) has reached it. At N_min, the 25 permutation replicates ARE an
   empirical sampling distribution for "distinct alleles found after N_min people" -- report
   mean + 95% band per ancestry, plus a bootstrap difference test (AFR vs every other ancestry):
   with a higher-replicate-count (`--n-bootstrap-test`, default 500) independent re-run of the same
   permutation machinery, the two-sided p-value is 2*min(P(diff<=0), P(diff>=0)) over paired
   (AFR - other) differences at step N_min.
6. "% of allele space explored, using IPD-IMGT richness" (Marc's rework of the old
   imgt_coverage_today metric, which measured haplotype-level catalogue hit-rate, not allele-space
   coverage): per gene x ancestry, (# distinct 2-field IPD-IMGT-catalogued alleles observed in this
   ancestry) / (# 2-field alleles IPD-IMGT/HLA lists for that gene, from `load_refdata`'s own
   CDSseq index -- same refdata snapshot used everywhere else in this project, so the two counts
   are comparable). Caveat, stated in the README: IPD-IMGT/HLA's catalogue itself is EUR-biased in
   its submission history, so a LOW score for AFR/SAS/etc. partly reflects "the reference catalogue
   under-represents this ancestry" as much as "this ancestry is under-sampled here" -- the two are
   not separable from this number alone.

DISCLOSURE
----------
No participant-level row is ever written off the VM. Every exported curve point is a
mean/2.5/97.5-percentile summary across >= 25 random subsamples of >= MIN_PEOPLE_PER_ANCESTRY
people; any such statistic derived from fewer than 20 people is written as the literal string
"<20" and must never be treated as 0 downstream.

USAGE
-----
    pixi run -e spechla -- python3 scripts/hla_popgen/39_saturation_by_ancestry.py \\
        [--table1 PATH] [--cohort-membership PATH] [--refdata DIR] \\
        [--relatedness-table PATH] [--out-dir DIR] [--n-permutations 25] \\
        [--n-bootstrap-test 500] [--strict-threshold 0.95] [--limit N]
"""
import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _viz_common as vc

try:
    from scipy.optimize import curve_fit  # noqa: F401 -- imported by 04's module, checked here too
    HAVE_SCIPY_OPT = True
except ImportError:
    HAVE_SCIPY_OPT = False


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_m24 = None
_m04 = None


def m24():
    global _m24
    if _m24 is None:
        _m24 = _load_module("24_novelty_by_field.py", "novelty_by_field")
    return _m24


def m04():
    global _m04
    if _m04 is None:
        _m04 = _load_module("04_allele_saturation.py", "allele_saturation")
    return _m04


SUPPRESS_BELOW = 20
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DEFAULT_OUTROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
DEFAULT_PEOPLE_OUTROOT = os.path.join(DEFAULT_OUTROOT, "people")
# ENVIRONMENT quirk #35: never point this at ~/mnt (stale manual mount, hangs I/O).
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
DEFAULT_REFDATA = os.path.expanduser("~/tools/Immuannot_refdata")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/39_saturation_by_ancestry")

KIN_MIN = 0.0442
STRICT_MIN_DEFAULT = 0.95   # task: "stricter-probability sensitivity version (>=0.95)"
CLASSICAL_GENES_BARE = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
ANCESTRY_ORDER = vc.ANCESTRY_ORDER
ANCESTRY_COLORS = vc.ANCESTRY_COLORS

THRESHOLDS = [1, 2, 10]
N_PERMUTATIONS_DEFAULT = 25          # Pakistan Fig 3e: 25 random subsamples per point
N_BOOTSTRAP_TEST_DEFAULT = 200       # higher-replicate re-run, for the equal-N difference test
MIN_PEOPLE_PER_ANCESTRY = 100        # "well-powered" floor for the common-N comparison
SCHEMES = ["pred", "strict"]         # anc_pred (primary) + anc_strict (sensitivity, >=0.95)


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
def suppress(n, threshold=SUPPRESS_BELOW):
    if pd.isna(n):
        return n
    n = int(n)
    return "0" if n == 0 else (f"<{threshold}" if n < threshold else str(n))


def suppress_df(df, cols, threshold=SUPPRESS_BELOW):
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = out[c].apply(lambda v: suppress(v, threshold))
    return out


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


# ---------------------------------------------------------------------------
# Step 1-2: build the fully-labeled call table + people table (mirrors 24_novelty_by_field.main(),
# stopping short of its report/cluster generation -- we only need per-call identity + ancestry).
# ---------------------------------------------------------------------------
def build_labeled_calls(args):
    m = m24()
    print("[39] loading table1 ...", file=sys.stderr, flush=True)
    t1 = m.load_table1(args.table1, args.limit)
    t1 = m.add_call_labels(t1)
    t1["gene_b"] = t1["gene"].map(m.gene_bare)
    print(f"  {len(t1)} rows, {t1['person_id'].nunique()} people", file=sys.stderr, flush=True)

    # build_people(..., strict_min=args.strict_threshold) gives us BOTH anc_pred (plain argmax,
    # no threshold -- the task's primary "AoU predicted ancestry" scheme) and anc_strict
    # (probability >= args.strict_threshold, default 0.95 -- the task's sensitivity scheme) in one
    # call; no need to re-derive strict_ancestry a second time.
    people, n_removed = m.build_people(
        t1["person_id"].unique(), args.cohort_membership, args.relatedness_table,
        args.kin_min, args.strict_threshold, args.skip_relatedness)
    people = people.rename(columns={"anc_strict": "anc_strict95"})

    c1 = t1[t1["copy_index"] == 1].sort_values(["person_id", "hap", "gene", "contig"])
    calls = c1.drop_duplicates(["person_id", "hap", "gene"]).merge(
        people, on="person_id", how="left")
    calls["unrelated"] = calls["unrelated"].fillna(False).astype(bool)
    for c in ("anc_pred", "anc_strict95"):
        calls[c] = calls[c].fillna("UNASSIGNED")

    genes_needed = set(t1.loc[t1["depth"].isin([2, 3]), "gene_b"])
    print(f"[39] loading refdata ({len(genes_needed)} genes) ...", file=sys.stderr, flush=True)
    ref = m.load_refdata(os.path.expanduser(args.refdata), genes_needed)
    print(f"  {ref.n_records} known CDS records from {len(ref.files)} files",
          file=sys.stderr, flush=True)
    print("[39] matching depth-2/3 calls to cds.fa.gz ...", file=sys.stderr, flush=True)
    matched, mstats = m.match_sequences(t1, os.path.expanduser(args.outroot), args.threads)
    cache = {}
    seq_rows = []
    for r in matched:
        gb = m.gene_bare(r["gene"])
        key = (gb, r["seq"], r["depth"])
        if key not in cache:
            cache[key] = m.classify_sequence(r["seq"], gb, r["depth"], ref)
        info = cache[key]
        seq_rows.append({"person_id": r["person_id"], "hap": r["hap"], "gene": r["gene"],
                         "nearest_allele": r["nearest_allele"], **info})
    seq_df = pd.DataFrame(seq_rows, columns=["person_id", "hap", "gene", "nearest_allele",
                                             "cds_core", "protein", "frameshift",
                                             "premature_stop", "cds_known", "protein_known",
                                             "mapped_prot2", "mapped_cds3", "seq_class"])
    calls = calls.merge(seq_df, on=["person_id", "hap", "gene"], how="left")
    calls = m.allele_ids(calls)

    ref_catalogue_size = imgt_catalogue_size(ref, CLASSICAL_GENES_BARE)
    return calls, ref_catalogue_size, mstats, n_removed


def imgt_catalogue_size(ref, genes_bare):
    """Distinct 2-field IPD-IMGT protein names per gene, from the same CDSseq snapshot used for
    novelty matching everywhere else in this project (ref.name_prot: allele name -> protein)."""
    m = m24()
    counts = {}
    for g in genes_bare:
        names = [n for n, p in ref.name_prot.items() if m.gene_bare(n.split("*")[0]) == g]
        two_field = {m.n_field_name(n, 2) for n in names}
        two_field.discard(None)
        counts[g] = len(two_field)
    return counts


# ---------------------------------------------------------------------------
# Step 3: person -> allele-identity sets, per ancestry / gene-group / category
# ---------------------------------------------------------------------------
def build_person_unit_sets(calls, ancestry_by_person, ancestry, genes_bare, category):
    """Returns dict person_id -> set(prot_id) for `genes_bare` pooled, restricted to `category` in
    {"all", "novel"}, `keep_clean` calls only, and (if `ancestry` is not None) to that ancestry."""
    sub = calls[calls["gene_b"].isin(genes_bare) & calls["keep_clean"] & calls["prot_id"].notna()]
    if category == "novel":
        sub = sub[(sub["field_class"] == "f2_protein") & (sub["seq_class"] == "novel_protein")]
    if ancestry is not None:
        sub = sub[sub["person_id"].map(ancestry_by_person) == ancestry]
    by_person = {}
    for pid, prot_id in zip(sub["person_id"], sub["prot_id"]):
        by_person.setdefault(pid, set()).add(prot_id)
    return by_person


def rarefaction_threshold_curves(unit_sets, thresholds, n_permutations, seed):
    """unit_sets: list of sets (one per person, both haplotypes pooled). Returns dict:
    'distinct': (n_permutations, m) array of cumulative distinct-allele counts,
    thresholds[k]: (n_permutations, m) array of #alleles that have reached >=k carriers so far.
    """
    m = len(unit_sets)
    rng = np.random.default_rng(seed)
    distinct = np.zeros((n_permutations, m))
    thr = {k: np.zeros((n_permutations, m)) for k in thresholds}
    idx = np.arange(m)
    for p in range(n_permutations):
        order = rng.permutation(idx)
        carrier_count = {}
        n_distinct = 0
        n_reached = {k: 0 for k in thresholds}
        for step, i in enumerate(order):
            for allele in unit_sets[i]:
                c = carrier_count.get(allele, 0)
                if c == 0:
                    n_distinct += 1
                c += 1
                carrier_count[allele] = c
                if c in n_reached:
                    n_reached[c] += 1
            distinct[p, step] = n_distinct
            for k in thresholds:
                thr[k][p, step] = n_reached[k]
    out = {"distinct": distinct}
    out.update(thr)
    return out


def summarize_curve(arr):
    """arr: (n_permutations, m). Returns steps(1..m), mean, lo2.5, hi97.5."""
    m = arr.shape[1]
    steps = np.arange(1, m + 1)
    mean = arr.mean(axis=0)
    lo = np.percentile(arr, 2.5, axis=0)
    hi = np.percentile(arr, 97.5, axis=0)
    return steps, mean, lo, hi


# ---------------------------------------------------------------------------
# Step 4: extrapolation (Clench/Michaelis-Menten, reused verbatim from 04)
# ---------------------------------------------------------------------------
def fit_and_extrapolate(steps, mean, lo, hi, unit_sets_for_chao2):
    mod = m04()
    fit = mod.fit_clench_asymptote(steps, mean, lo, hi)
    n_now = int(steps[-1])
    per_1000 = np.nan
    if fit["fit_ok"]:
        per_1000 = mod.clench_curve(n_now + 1000, fit["s_max"], fit["b"]) - mean[-1]
    S_obs, m_units, Q = mod.incidence_freqs(unit_sets_for_chao2)
    chao2 = mod.chao2(S_obs, Q, m_units)
    return {
        "n_now": n_now, "s_obs_now": float(mean[-1]),
        "clench_s_max": fit["s_max"], "clench_b": fit["b"], "clench_fit_ok": fit["fit_ok"],
        "clench_s_max_lo": fit.get("s_max_lo", np.nan), "clench_s_max_hi": fit.get("s_max_hi", np.nan),
        "expected_new_per_1000": float(per_1000) if per_1000 == per_1000 else np.nan,
        "chao2_richness": float(chao2), "chao2_undetected": float(chao2 - S_obs),
    }


# ---------------------------------------------------------------------------
# Equal-N comparison + AFR-vs-rest bootstrap difference test
# ---------------------------------------------------------------------------
def equal_n_comparison(by_anc_units, n_min, n_bootstrap, seed):
    rows = []
    reps = {}
    for anc, units in sorted(by_anc_units.items()):
        if len(units) < n_min:
            continue
        curves = rarefaction_threshold_curves(units, THRESHOLDS, n_bootstrap, seed=seed + hash(anc) % 1000)
        vals = curves["distinct"][:, n_min - 1]
        reps[anc] = vals
        rows.append({"ancestry": anc, "n_min": n_min, "n_available": len(units),
                     "mean_distinct_at_n_min": float(vals.mean()),
                     "lo2_5": float(np.percentile(vals, 2.5)),
                     "hi97_5": float(np.percentile(vals, 97.5))})
    diffs = []
    if "AFR" in reps:
        for anc, vals in sorted(reps.items()):
            if anc == "AFR":
                continue
            d = reps["AFR"] - vals
            p_le0 = float((d <= 0).mean())
            p_ge0 = float((d >= 0).mean())
            p_two_sided = 2 * min(p_le0, p_ge0)
            diffs.append({"comparison": f"AFR - {anc}", "n_min": n_min,
                         "mean_diff": float(d.mean()), "lo2_5": float(np.percentile(d, 2.5)),
                         "hi97_5": float(np.percentile(d, 97.5)),
                         "p_two_sided_bootstrap": min(p_two_sided, 1.0)})
    return pd.DataFrame(rows), pd.DataFrame(diffs)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_main_panel(curve_data, out_stem, title_suffix=""):
    """Main-figure-ready panel (89 mm wide): per-ancestry >=1-carrier discovery curves, shaded
    2.5/97.5 percentile bands, direct end-of-line labels instead of a legend."""
    with vc.nature_style():
        fig, ax = plt.subplots(figsize=(vc.mm(89), vc.mm(70)))
        for anc in ANCESTRY_ORDER:
            if anc not in curve_data:
                continue
            steps, mean, lo, hi = curve_data[anc]
            color = ANCESTRY_COLORS.get(anc, "#333333")
            ax.plot(steps, mean, color=color, linewidth=0.9, zorder=3)
            ax.fill_between(steps, lo, hi, color=color, alpha=0.18, linewidth=0, zorder=2)
            ax.annotate(anc, (steps[-1], mean[-1]), xytext=(3, 0), textcoords="offset points",
                        fontsize=5.5, color=color, va="center", fontweight="bold")
        ax.set_xlabel("cohort size (people)", fontsize=7)
        ax.set_ylabel("cumulative distinct protein alleles", fontsize=7)
        if title_suffix:
            ax.set_title(title_suffix, fontsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        return vc.save_fig(fig, out_stem)


def fig_supplement_grid(curve_by_gene_anc, genes_bare, out_stem):
    with vc.nature_style():
        ncols = 4
        nrows = int(np.ceil(len(genes_bare) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(vc.mm(183), vc.mm(45 * nrows)))
        axes = np.atleast_1d(axes).ravel()
        for ax, gene in zip(axes, genes_bare):
            for anc in ANCESTRY_ORDER:
                key = (gene, anc)
                if key not in curve_by_gene_anc:
                    continue
                steps, mean, lo, hi = curve_by_gene_anc[key]
                color = ANCESTRY_COLORS.get(anc, "#333333")
                ax.plot(steps, mean, color=color, linewidth=0.6)
                ax.fill_between(steps, lo, hi, color=color, alpha=0.15, linewidth=0)
            ax.set_title(gene, fontsize=6)
            ax.tick_params(labelsize=5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        for ax in axes[len(genes_bare):]:
            ax.axis("off")
        handles = [plt.Line2D([0], [0], color=ANCESTRY_COLORS[a], lw=1.2, label=a)
                  for a in ANCESTRY_ORDER]
        fig.legend(handles=handles, loc="lower center", ncol=len(ANCESTRY_ORDER), fontsize=5.5,
                  frameon=False, bbox_to_anchor=(0.5, -0.02))
        fig.supxlabel("cohort size (people)", fontsize=6.5)
        fig.supylabel("cumulative distinct protein alleles", fontsize=6.5)
        return vc.save_fig(fig, out_stem)


def fig_thresholds_panel(curve_data_by_thr, ancestry, out_stem):
    """Pakistan Fig 3e style, single ancestry: 3 lines (>=1, >=2, >=10 carriers)."""
    thr_colors = {1: "#333333", 2: "#B4472E", 10: "#0072B2"}
    with vc.nature_style():
        fig, ax = plt.subplots(figsize=(vc.mm(89), vc.mm(65)))
        for k in THRESHOLDS:
            if k not in curve_data_by_thr:
                continue
            steps, mean, lo, hi = curve_data_by_thr[k]
            ax.plot(steps, mean, color=thr_colors[k], linewidth=0.9)
            ax.fill_between(steps, lo, hi, color=thr_colors[k], alpha=0.15, linewidth=0)
            ax.annotate(f">= {k} carrier{'s' if k > 1 else ''}", (steps[-1], mean[-1]),
                       xytext=(3, 0), textcoords="offset points", fontsize=5.5,
                       color=thr_colors[k], va="center")
        ax.set_xlabel("cohort size (people)", fontsize=7)
        ax.set_ylabel("cumulative distinct alleles", fontsize=7)
        ax.set_title(ancestry, fontsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        return vc.save_fig(fig, out_stem)


def fig_imgt_richness(df, out_stem):
    with vc.nature_style():
        genes = sorted(df["gene"].unique())
        fig, ax = plt.subplots(figsize=(vc.mm(183), vc.mm(60)))
        x = np.arange(len(genes))
        width = 0.8 / len(ANCESTRY_ORDER)
        for i, anc in enumerate(ANCESTRY_ORDER):
            d = df[df["ancestry"] == anc].set_index("gene")
            y = [d.loc[g, "pct_explored"] if g in d.index else np.nan for g in genes]
            ax.bar(x + i * width, y, width=width, color=ANCESTRY_COLORS.get(anc, "#333"),
                  label=anc)
        ax.set_xticks(x + width * (len(ANCESTRY_ORDER) - 1) / 2)
        ax.set_xticklabels(genes, fontsize=6)
        ax.set_ylabel("% of IPD-IMGT 2-field alleles observed", fontsize=7)
        ax.legend(fontsize=5.5, frameon=False, ncol=len(ANCESTRY_ORDER))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        return vc.save_fig(fig, out_stem)


# ---------------------------------------------------------------------------
# Render-only mode: rebuild figures + README from previously-exported aggregate TSVs (curves.tsv,
# extrapolation_fits.tsv, equal_n_comparison.tsv, afr_vs_rest_bootstrap_test.tsv,
# imgt_allele_space_explored.tsv, summary.json) -- no AoU data, no m24/m04 heavy pipeline needed.
# Meant to run LOCALLY after pulling back only those aggregate files from the VM (VM_CHANNEL.md:
# "preferred -- faster iteration").
# ---------------------------------------------------------------------------
def render_only(args):
    out_dir = args.out_dir
    curve_df = pd.read_csv(os.path.join(out_dir, "curves.tsv"), sep="\t")
    fit_path = os.path.join(out_dir, "extrapolation_fits.tsv")
    fit_df = pd.read_csv(fit_path, sep="\t") if os.path.exists(fit_path) else pd.DataFrame()
    eq_path = os.path.join(out_dir, "equal_n_comparison.tsv")
    eq_df = pd.read_csv(eq_path, sep="\t") if os.path.exists(eq_path) else pd.DataFrame()
    diff_path = os.path.join(out_dir, "afr_vs_rest_bootstrap_test.tsv")
    diff_df = pd.read_csv(diff_path, sep="\t") if os.path.exists(diff_path) else pd.DataFrame()
    imgt_path = os.path.join(out_dir, "imgt_allele_space_explored.tsv")
    imgt_df = pd.read_csv(imgt_path, sep="\t") if os.path.exists(imgt_path) else pd.DataFrame()
    summary_path = os.path.join(out_dir, "summary.json")
    with open(summary_path) as fh:
        summary = json.load(fh)

    pred = curve_df[curve_df["scheme"] == "pred"]

    def _curve_tuple(sub):
        sub = sub.sort_values("n")
        return (sub["n"].to_numpy(), sub["mean_distinct"].to_numpy(),
                sub["lo2_5"].to_numpy(), sub["hi97_5"].to_numpy())

    main_panel_data = {}
    pooled_all = pred[(pred["gene_group"] == "classical_pooled") & (pred["category"] == "all")]
    for anc, sub in pooled_all.groupby("ancestry"):
        main_panel_data[anc] = _curve_tuple(sub)

    thresholds_panel_data = {}
    for k in THRESHOLDS:
        cat = f"carriers_ge_{k}"
        sub_all = pred[(pred["gene_group"] == "classical_pooled") & (pred["category"] == cat)]
        for anc, sub in sub_all.groupby("ancestry"):
            thresholds_panel_data.setdefault(anc, {})[k] = _curve_tuple(sub)

    supplement_data = {}
    for gene in CLASSICAL_GENES_BARE:
        sub_gene = pred[(pred["gene_group"] == gene) & (pred["category"] == "all")]
        for anc, sub in sub_gene.groupby("ancestry"):
            supplement_data[(gene, anc)] = _curve_tuple(sub)

    fig_main_panel(main_panel_data,
                   os.path.join(out_dir, "fig1_main_saturation_panel"),
                   "distinct protein alleles, classical genes pooled")
    if supplement_data:
        fig_supplement_grid(supplement_data, CLASSICAL_GENES_BARE,
                            os.path.join(out_dir, "fig2_supplement_grid_by_gene"))
    for anc, thr_data in thresholds_panel_data.items():
        fig_thresholds_panel(thr_data, anc, os.path.join(out_dir, f"fig3_thresholds_{anc}"))
    if len(imgt_df):
        fig_imgt_richness(imgt_df[imgt_df["scheme"] == "pred"],
                          os.path.join(out_dir, "fig4_imgt_allele_space_explored"))

    write_readme(out_dir, summary, fit_df, eq_df, diff_df, imgt_df)
    print(f"[39] render-only done -> {out_dir}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run(args):
    ensure_dir(args.out_dir)
    t0 = time.time()
    calls, ref_catalogue_size, mstats, n_removed = build_labeled_calls(args)

    curve_rows, fit_rows, imgt_rows = [], [], []
    equal_n_rows, diff_rows = [], []
    main_panel_data = {}
    thresholds_panel_data = {}
    supplement_data = {}

    for scheme, anc_col in (("pred", "anc_pred"), ("strict95", "anc_strict95")):
        keep = calls[calls["unrelated"]]
        ancestry_by_person = dict(zip(keep["person_id"], keep[anc_col]))

        # --- pooled classical genes: main analysis ---
        by_anc_units_all, by_anc_units_novel = {}, {}
        for anc in ANCESTRY_ORDER:
            person_sets_all = build_person_unit_sets(
                calls, ancestry_by_person, anc, CLASSICAL_GENES_BARE, "all")
            person_sets_novel = build_person_unit_sets(
                calls, ancestry_by_person, anc, CLASSICAL_GENES_BARE, "novel")
            units_all = list(person_sets_all.values())
            units_novel = list(person_sets_novel.values())
            if len(units_all) < MIN_PEOPLE_PER_ANCESTRY:
                continue
            by_anc_units_all[anc] = units_all
            by_anc_units_novel[anc] = units_novel

            curves = rarefaction_threshold_curves(units_all, THRESHOLDS, args.n_permutations,
                                                  seed=args.seed)
            steps, mean, lo, hi = summarize_curve(curves["distinct"])
            if scheme == "pred":
                main_panel_data[anc] = (steps, mean, lo, hi)
                thresholds_panel_data[anc] = {
                    k: summarize_curve(curves[k]) for k in THRESHOLDS}
            for i, s in enumerate(steps):
                curve_rows.append({"scheme": scheme, "ancestry": anc, "gene_group": "classical_pooled",
                                  "category": "all", "n": int(s),
                                  "mean_distinct": round(float(mean[i]), 2),
                                  "lo2_5": round(float(lo[i]), 2), "hi97_5": round(float(hi[i]), 2)})
                for k in THRESHOLDS:
                    _, mk, lk, hk = summarize_curve(curves[k])
                    curve_rows.append({"scheme": scheme, "ancestry": anc, "gene_group": "classical_pooled",
                                      "category": f"carriers_ge_{k}", "n": int(s),
                                      "mean_distinct": round(float(mk[i]), 2),
                                      "lo2_5": round(float(lk[i]), 2), "hi97_5": round(float(hk[i]), 2)})

            fit = fit_and_extrapolate(steps, mean, lo, hi, units_all)
            fit.update({"scheme": scheme, "ancestry": anc, "gene_group": "classical_pooled",
                       "category": "all"})
            fit_rows.append(fit)

            if units_novel:
                curves_nov = rarefaction_threshold_curves(units_novel, [1], args.n_permutations,
                                                          seed=args.seed + 7)
                steps_n, mean_n, lo_n, hi_n = summarize_curve(curves_nov["distinct"])
                for i, s in enumerate(steps_n):
                    curve_rows.append({"scheme": scheme, "ancestry": anc, "gene_group": "classical_pooled",
                                      "category": "novel", "n": int(s),
                                      "mean_distinct": round(float(mean_n[i]), 2),
                                      "lo2_5": round(float(lo_n[i]), 2), "hi97_5": round(float(hi_n[i]), 2)})
                fit_n = fit_and_extrapolate(steps_n, mean_n, lo_n, hi_n, units_novel)
                fit_n.update({"scheme": scheme, "ancestry": anc, "gene_group": "classical_pooled",
                             "category": "novel"})
                fit_rows.append(fit_n)

        # --- per-gene supplement (pred scheme only, to bound VM time) ---
        if scheme == "pred":
            for gene in CLASSICAL_GENES_BARE:
                for anc in ANCESTRY_ORDER:
                    person_sets = build_person_unit_sets(calls, ancestry_by_person, anc, [gene], "all")
                    units = list(person_sets.values())
                    if len(units) < MIN_PEOPLE_PER_ANCESTRY:
                        continue
                    curves = rarefaction_threshold_curves(units, [1], args.n_permutations,
                                                          seed=args.seed + 13)
                    steps, mean, lo, hi = summarize_curve(curves["distinct"])
                    supplement_data[(gene, anc)] = (steps, mean, lo, hi)
                    for i, s in enumerate(steps):
                        curve_rows.append({"scheme": scheme, "ancestry": anc, "gene_group": gene,
                                          "category": "all", "n": int(s),
                                          "mean_distinct": round(float(mean[i]), 2),
                                          "lo2_5": round(float(lo[i]), 2), "hi97_5": round(float(hi[i]), 2)})

        # --- equal-N comparison + AFR-vs-rest test (pred scheme, the primary scheme) ---
        if scheme == "pred" and by_anc_units_all:
            n_min = min(len(v) for v in by_anc_units_all.values())
            eq_df, diff_df = equal_n_comparison(by_anc_units_all, n_min, args.n_bootstrap_test,
                                               args.seed + 99)
            eq_df["scheme"] = scheme
            diff_df["scheme"] = scheme
            equal_n_rows.append(eq_df)
            diff_rows.append(diff_df)

        # --- IMGT allele-space-explored richness, per gene x ancestry ---
        keep_c = calls[calls["unrelated"] & calls["keep_clean"] & calls["prot_id"].notna()
                       & (calls["seq_class"].fillna("known") != "novel_protein")]
        keep_c = keep_c.copy()
        keep_c["ancestry"] = keep_c["person_id"].map(ancestry_by_person)
        for (gene, anc), grp in keep_c.groupby(["gene_b", "ancestry"]):
            if gene not in CLASSICAL_GENES_BARE or anc not in ANCESTRY_ORDER:
                continue
            n_observed = grp["prot_id"].nunique()
            cat_size = ref_catalogue_size.get(gene, 0)
            imgt_rows.append({"scheme": scheme, "gene": gene, "ancestry": anc,
                             "n_distinct_known_alleles_observed_disp": suppress(n_observed),
                             "imgt_catalogue_size_2field": cat_size,
                             "pct_explored": round(100.0 * n_observed / cat_size, 2)
                             if cat_size else np.nan})

    # ---- write TSVs ----
    curve_df = pd.DataFrame(curve_rows)
    curve_df.to_csv(os.path.join(args.out_dir, "curves.tsv"), sep="\t", index=False)

    fit_df = pd.DataFrame(fit_rows)
    fit_df.to_csv(os.path.join(args.out_dir, "extrapolation_fits.tsv"), sep="\t", index=False)

    imgt_df = pd.DataFrame(imgt_rows)
    imgt_df.to_csv(os.path.join(args.out_dir, "imgt_allele_space_explored.tsv"), sep="\t",
                   index=False)

    if equal_n_rows:
        eq_all = pd.concat(equal_n_rows, ignore_index=True)
        eq_all.to_csv(os.path.join(args.out_dir, "equal_n_comparison.tsv"), sep="\t", index=False)
    else:
        eq_all = pd.DataFrame()
    if diff_rows:
        diff_all = pd.concat(diff_rows, ignore_index=True)
        diff_all.to_csv(os.path.join(args.out_dir, "afr_vs_rest_bootstrap_test.tsv"), sep="\t",
                        index=False)
    else:
        diff_all = pd.DataFrame()

    # ---- figures ----
    fig_main_panel(main_panel_data,
                  os.path.join(args.out_dir, "fig1_main_saturation_panel"),
                  "distinct protein alleles, classical genes pooled")
    if supplement_data:
        fig_supplement_grid(supplement_data, CLASSICAL_GENES_BARE,
                           os.path.join(args.out_dir, "fig2_supplement_grid_by_gene"))
    for anc, thr_data in thresholds_panel_data.items():
        fig_thresholds_panel(thr_data, anc,
                            os.path.join(args.out_dir, f"fig3_thresholds_{anc}"))
    if len(imgt_df):
        fig_imgt_richness(imgt_df[imgt_df["scheme"] == "pred"],
                         os.path.join(args.out_dir, "fig4_imgt_allele_space_explored"))

    summary = {
        "n_people_unrelated": suppress(int(calls[calls["unrelated"]]["person_id"].nunique())),
        "n_people_removed_relatedness": n_removed,
        "n_permutations": args.n_permutations,
        "n_bootstrap_test": args.n_bootstrap_test,
        "strict_threshold": args.strict_threshold,
        "classical_genes": CLASSICAL_GENES_BARE,
        "min_people_per_ancestry": MIN_PEOPLE_PER_ANCESTRY,
        "match_stats": {k: v for k, v in mstats.items()} if isinstance(mstats, dict) else None,
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    write_readme(args.out_dir, summary, fit_df, eq_all, diff_all, imgt_df)
    print(f"[39] done -> {args.out_dir} ({summary['elapsed_sec']}s)", file=sys.stderr, flush=True)


def write_readme(out_dir, summary, fit_df, eq_df, diff_df, imgt_df):
    L = []
    L.append("# Per-ancestry allele discovery / saturation curves\n\n")
    L.append("_Generated by `39_saturation_by_ancestry.py`. Spec: Marc's proposal, endorsed by "
            "Cole on the 2026-09-22 call (sprints/CALL_SUMMARY_2026-09-22.md Sec 6), styled after "
            "the Pakistan Genome Resource (Nature, 2026) Fig 3e curves._\n\n")
    L.append("## Question\n\n")
    L.append("As more people are added to each ancestry group, how many distinct HLA protein "
            "(2-field) alleles have been found so far, and how far from saturated is each "
            "ancestry? Call expectation: AFR is the least saturated.\n\n")
    L.append("## Method\n\n")
    L.append("Sampling unit = **person** (both haplotypes). 25 random orderings of each "
            "ancestry's unrelated people (Pakistan Fig 3e convention); mean + 2.5/97.5 percentile "
            "bands plotted/exported at every cohort size. Identity = `prot_id` from "
            "24_novelty_by_field.py's own `allele_ids()` (2-field IPD-IMGT name for known "
            "alleles, `<gene>_prot_<sha8>` for a novel protein), restricted to clean "
            "(non-artifact, non-frameshift) calls. Saturation model: Clench (1979) / "
            "Michaelis-Menten S(N)=Smax*N/(b+N), reused verbatim from "
            "`04_allele_saturation.py::fit_clench_asymptote`; Chao2 (incidence) reported "
            "alongside as an independent richness estimator. Equal-N comparison at N_min = the "
            "largest N every well-powered ancestry (>=100 people) reaches; AFR-vs-rest "
            f"difference tested via {summary['n_bootstrap_test']}-replicate bootstrap, two-sided "
            "percentile p-value.\n\n")
    L.append(f"- Unrelated people (all ancestries, primary scheme): "
            f"**{summary['n_people_unrelated']}**\n")
    L.append(f"- Removed for relatedness: {summary['n_people_removed_relatedness']}\n")
    L.append(f"- Permutations per curve point: {summary['n_permutations']}\n\n")

    L.append("## Extrapolation (classical genes pooled, all alleles, `pred` ancestry scheme)\n\n")
    if len(fit_df):
        sub = fit_df[(fit_df["scheme"] == "pred") & (fit_df["gene_group"] == "classical_pooled")
                     & (fit_df["category"] == "all")]
        L.append("| Ancestry | N now | Distinct now | Expected new / next 1,000 people | "
                "Clench S_max | Chao2 richness |\n|---|---|---|---|---|---|\n")
        for _, r in sub.sort_values("ancestry").iterrows():
            L.append(f"| {r['ancestry']} | {r['n_now']} | {r['s_obs_now']:.0f} | "
                    f"{r['expected_new_per_1000']:.1f} | {r['clench_s_max']:.0f} | "
                    f"{r['chao2_richness']:.0f} |\n")
    L.append("\n")

    L.append("## Equal-N comparison and AFR-vs-rest test\n\n")
    if len(eq_df):
        L.append("| Ancestry | N_min | Mean distinct | 95% band |\n|---|---|---|---|\n")
        for _, r in eq_df.sort_values("ancestry").iterrows():
            L.append(f"| {r['ancestry']} | {r['n_min']} | {r['mean_distinct_at_n_min']:.1f} | "
                    f"[{r['lo2_5']:.1f}, {r['hi97_5']:.1f}] |\n")
    L.append("\n")
    if len(diff_df):
        L.append("| Comparison | Mean diff | 95% band | p (two-sided bootstrap) |\n"
                "|---|---|---|---|\n")
        for _, r in diff_df.iterrows():
            L.append(f"| {r['comparison']} | {r['mean_diff']:.1f} | "
                    f"[{r['lo2_5']:.1f}, {r['hi97_5']:.1f}] | {r['p_two_sided_bootstrap']:.4f} |\n")
    L.append("\n")

    L.append("## % of allele space explored (IPD-IMGT richness), classical genes\n\n")
    L.append("Defined as: (distinct 2-field IPD-IMGT-catalogued alleles observed in this "
            "ancestry) / (2-field alleles IPD-IMGT/HLA lists for that gene, from the same "
            "refdata CDSseq snapshot used for novelty matching everywhere in this project). "
            "**Caveat: IPD-IMGT/HLA's own catalogue is itself EUR-biased in submission history — "
            "a low score for a given ancestry partly reflects reference under-representation, "
            "not only under-sampling here; the two are not separable from this number alone.**\n\n")
    if len(imgt_df):
        piv = imgt_df[imgt_df["scheme"] == "pred"].pivot_table(
            index="gene", columns="ancestry", values="pct_explored").round(1)
        # Avoid pandas.DataFrame.to_markdown() -- pulls in the optional 'tabulate' dependency,
        # not guaranteed present in the pixi env. Build the markdown table by hand instead.
        cols = list(piv.columns)
        L.append("| gene | " + " | ".join(str(c) for c in cols) + " |\n")
        L.append("|---|" + "---|" * len(cols) + "\n")
        for gene, row in piv.iterrows():
            cells = [f"{v:.1f}" if pd.notna(v) else "n/a" for v in row]
            L.append(f"| {gene} | " + " | ".join(cells) + " |\n")
        L.append("\n")

    L.append("## Figures\n\n")
    L.append("- `fig1_main_saturation_panel` — main-figure candidate (89 mm), classical genes "
            "pooled, per-ancestry curves with direct end labels.\n")
    L.append("- `fig2_supplement_grid_by_gene` — per-classical-gene grid, all ancestries "
            "overlaid.\n")
    L.append("- `fig3_thresholds_<ANC>` — Pakistan Fig 3e style per-ancestry panel: >=1/>=2/>=10 "
            "carrier thresholds.\n")
    L.append("- `fig4_imgt_allele_space_explored` — % of IPD-IMGT catalogue observed, by gene x "
            "ancestry.\n\n")

    L.append("## Caveats\n\n")
    L.append("- Novel-protein-allele curves (category=`novel` in `curves.tsv`) use the strict "
            "S01/S02 definition (field_class==f2_protein AND seq_class==novel_protein, clean "
            "calls only); recurrence within a novel cluster follows 34_novel_protein_recurrence's "
            "convention but is not itself re-derived here.\n")
    L.append("- Per-gene supplement and equal-N/bootstrap test use the `pred` (AoU predicted, "
            "argmax) ancestry scheme only, to bound VM runtime; `strict95` results for the pooled "
            "classical-gene curve are in `curves.tsv`/`extrapolation_fits.tsv` "
            "(scheme==`strict95`) as the sensitivity check.\n")
    L.append("- Clench S_max brackets (`clench_s_max_lo/hi` in `extrapolation_fits.tsv`) are "
            "curve-based sensitivity brackets (refit to the 2.5/97.5 percentile bands), not a "
            "calibrated bootstrap CI — see 04_allele_saturation.py's own caveat on this "
            "estimator.\n")

    L.append("\n## Distilled\n\n")
    L.append("- Sampling unit: person (both haplotypes); 25 permutations/point, Pakistan Fig 3e "
            "convention.\n")
    L.append("- Identity = 24's own `prot_id` (2-field IPD-IMGT name or `<gene>_prot_<sha8>` for "
            "novel).\n")
    L.append("- Saturation model: Clench/Michaelis-Menten (reused from 04), Chao2 reported "
            "alongside.\n")
    L.append("- AFR-vs-rest tested by bootstrap difference at equal N, not eyeballed.\n")
    L.append("- IPD-IMGT allele-space-explored caveat: catalogue itself is EUR-biased.\n")

    with open(os.path.join(out_dir, "README.md"), "w") as fh:
        fh.writelines(L)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--outroot", default=DEFAULT_PEOPLE_OUTROOT,
                    help="where <person_id>/immuannot_output/... live (for match_sequences)")
    ap.add_argument("--refdata", default=DEFAULT_REFDATA)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=STRICT_MIN_DEFAULT)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS_DEFAULT)
    ap.add_argument("--n-bootstrap-test", type=int, default=N_BOOTSTRAP_TEST_DEFAULT)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, default=None, help="limit table1 rows (smoke test)")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--render-only", action="store_true",
                    help="Skip the AoU pipeline entirely; rebuild figures + README from "
                         "previously-exported aggregate TSVs already in --out-dir (curves.tsv, "
                         "extrapolation_fits.tsv, equal_n_comparison.tsv, "
                         "afr_vs_rest_bootstrap_test.tsv, imgt_allele_space_explored.tsv, "
                         "summary.json). For local rendering after pulling back only the "
                         "aggregate files from the VM.")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.render_only:
        render_only(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
