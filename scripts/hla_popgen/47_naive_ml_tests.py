#!/usr/bin/env python3
"""S04 WS-D -- naive ML / simple tests on allele-carriage feature matrices (HLA two-field +
KIR gene presence), aggregate-only outputs.

Marc's ask (ORCHESTRATOR_HANDOFF.md sect2 WS-D): easy-to-understand ML that reveals the
*structure* of the data (batch effects, ancestry leakage, HLA<->KIR co-evolution) rather than
a disease-prediction claim. Every task here is a positive-control / sanity-check, not a novel
science result -- see the README for how to read each number.

## Tasks (each a separate binary/one-vs-rest classification)

1. **ancestry (one-vs-rest)**: HLA/KIR carriage -> AFR/AMR/EAS/EUR/MID/SAS, one binary model per
   ancestry (same convention as 42_repertoire_baseline.py's ancestry positive control). Expected:
   strongly predictable (HLA is *the* classical ancestry-informative locus) -- a confound warning
   for any future disease-carriage association, not a finding in itself.
2. **sequencing platform**: HLA/KIR carriage -> revio/sequel2e/sequel2 (`--platform-col`,
   default `platform`, from cohort_membership.tsv per SCHEMA.md Table 4; CLI-overridable because
   the column's exact name/values were not re-verified on every VM instance -- the S03 42-script
   README documents an instance where an expected column was entirely absent). Predictable would
   mean a technical/batch artifact riding on the calls, not biology.
3. **cA/cB from HLA**: person-level KIR A/B content (does either haplotype carry a B-content gene,
   43_kir_full_aggregate.py's B_CONTENT_GENES rule) predicted from HLA carriage alone. HLA and KIR
   sit on different chromosomes (6 vs 19) with no known direct genetic linkage, so this is a
   negative control: a strong positive result would be surprising and worth a second look (or a
   confound, e.g. ancestry correlating with both).

Two models per task: L1-logistic regression (`liblinear`) and a depth-3 decision tree (readable
by eye -- print its splits). 5-fold stratified CV, pooled out-of-fold AUROC with a permutation
baseline (>=20 label shuffles by default) so every AUROC ships with a null distribution, not a
bare number. Top |coefficient| / feature_importance_ features are reported ONLY for
carrier-count >=20 (the feature-construction step already drops rarer alleles before any model
sees them, so this is enforced structurally, not by a late filter that could leak a rare allele
name next to a small count).

## PCA (aggregate-only)

PCA of the same carriage matrix. Per-person coordinates NEVER leave the VM (disclosure hard
rule -- WS-D brief, ORCHESTRATOR_HANDOFF.md sect3). Only two aggregate views are written:
  - per-(ancestry x platform) group centroids (PC1/PC2 mean + n, groups with n<20 dropped);
  - a 2D binned density grid (`--pca-bins`, default 20x20) over PC1/PC2, cell counts <20
    suppressed to the literal string "<20" (never 0 -- a true empty bin stays 0).

## Data columns that need VM verification

`cohort_membership.tsv`'s `platform` column is documented in SCHEMA.md Table 4 as sourced from
`immuannot_cohort_full.tsv`, with values in `{revio, sequel2e, sequel2}` per 41_kir_pilot.py's
ASSEMBLY_PLATFORMS -- but 42_repertoire_baseline.py's README records a VM instance where
`hla_calls_rich.tsv` was "missing the expected person/allele columns entirely." This script
therefore never hardcodes a column name for anything beyond `person_id`/`ancestry_pred`: it
takes `--platform-col` as a CLI arg and, if the column is absent, prints every available column
in `cohort_membership.tsv` and exits (loud, not a silent skip) rather than guessing.

## Running (VM, real mode)

    cd ~/s04 && PYTHONPATH=~/s04:~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 47_naive_ml_tests.py \\
        --table1 ~/pipeline_outputs/hla_calls_rich.tsv \\
        --cohort-membership ~/pipeline_outputs/cohort_membership.tsv \\
        --relatedness-table ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv \\
        --kir-outroot ~/pipeline_outputs_kir \\
        --kir-pilot-script ~/s03/41_kir_pilot.py \\
        --platform-col platform \\
        --out-dir ~/s04/results/47

If `--platform-col` is wrong, rerun with the column name the script prints.

## Local synthetic dry run

    python3 scripts/hla_popgen/47_naive_ml_tests.py --synthetic --n-people 400 \\
        --out-dir /tmp/naive_ml_synthetic

Runs the identical pipeline against an in-memory synthetic cohort with a fixed seed and a
planted ancestry<->HLA-allele association (so the ancestry task should show real AUROC lift)
and NO planted HLA<->KIR association (so the cA/cB task should hover near chance).
"""
import argparse
import gzip
import importlib.util
import json
import os
import random
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
MIN_CARRIERS_FOR_FEATURE = 20
CLASSICAL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                    "HLA-DRB1"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def suppressed(n):
    """1-19 -> '<20' string; a true 0 stays '0' (never conflate the two -- see
    feedback_suppressed_counts_are_not_zero.md)."""
    n = int(n)
    if n == 0:
        return "0"
    return "<%d" % SUPPRESS_BELOW if n < SUPPRESS_BELOW else str(n)


# ---------------------------------------------------------------------------
# KIR module reuse (41_kir_pilot.py's parser/novelty/unrelatedness -- never re-derived).
# ---------------------------------------------------------------------------
def load_kir_module(script_path):
    spec = importlib.util.spec_from_file_location("kir_pilot_47", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Feature matrix construction
# ---------------------------------------------------------------------------
def hla_two_field(consensus):
    """First two colon-delimited fields after 'GENE*', stopping before any 'new' token. Mirrors
    _viz_common.to_nfield(consensus, 2) but reproduced here (stdlib+pandas only dependency,
    matching 43_kir_full_aggregate.py's own no-matplotlib-import discipline for an aggregate
    export script)."""
    if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
        return None
    s = str(consensus).strip()
    if not s or s.upper() in {"UNDETERMINED", "NA", ""}:
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = []
    for f in s.split(":"):
        if f == "":
            continue
        if f.strip().lower() == "new":
            break
        fields.append(f)
    if len(fields) < 2:
        return None
    return ":".join(fields[:2])


def build_hla_carriage(table1_df, min_carriers=MIN_CARRIERS_FOR_FEATURE):
    """Person x (gene, 2-field allele) boolean carriage matrix, classical genes only, restricted
    to alleles carried by >=min_carriers unrelated people (enforced by the caller passing an
    already-unrelated-filtered table1_df). Returns (DataFrame indexed by person_id, kept_cols,
    n_dropped_rare_alleles)."""
    df = table1_df[table1_df["gene"].isin(CLASSICAL_GENES)].copy()
    df["allele2"] = df["consensus"].map(hla_two_field)
    df = df.dropna(subset=["allele2"])
    df["feat"] = df["gene"] + "*" + df["allele2"]
    carriers = df.groupby("feat")["person_id"].nunique()
    keep_feats = sorted(carriers[carriers >= min_carriers].index)
    n_dropped = int((carriers < min_carriers).sum())
    people = sorted(table1_df["person_id"].unique())
    mat = pd.DataFrame(False, index=people, columns=keep_feats)
    sub = df[df["feat"].isin(keep_feats)][["person_id", "feat"]].drop_duplicates()
    for pid, feat in sub.itertuples(index=False):
        mat.at[pid, feat] = True
    mat.index.name = "person_id"
    return mat, n_dropped


def build_kir_features(pids, kir_outroot, kir_mod):
    """Person x KIR-gene presence (either haplotype) boolean matrix, plus a person-level cA/cB
    label ('cB' if either haplotype carries a B_CONTENT_GENES member, else 'cA'; 'missing' if
    both haplotypes failed to parse -- excluded from the cA/cB task by the caller)."""
    genes = kir_mod.KIR_GENES
    presence = pd.DataFrame(False, index=pids, columns=genes)
    content = pd.Series("missing", index=pids, dtype=object)
    for pid in pids:
        any_hap = False
        is_cb = False
        for hap in ("hap1", "hap2"):
            gtf_path = os.path.join(kir_outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
            if not os.path.exists(gtf_path):
                continue
            try:
                rows = kir_mod.parse_hap_gtf(gtf_path)
            except (OSError, EOFError):
                continue
            genes_called = {r["gene"] for r in rows if r["gene"] in genes}
            if not genes_called:
                continue
            any_hap = True
            for g in genes_called:
                presence.at[pid, g] = True
            if kir_mod.B_CONTENT_GENES.intersection(genes_called):
                is_cb = True
        if any_hap:
            content.loc[pid] = "cB" if is_cb else "cA"
    presence.index.name = "person_id"
    return presence, content


# ---------------------------------------------------------------------------
# Modeling: L1-LR + depth-3 tree, 5-fold CV AUROC, permutation baseline.
# ---------------------------------------------------------------------------
def run_binary_task(X, y, feature_names, n_splits=5, n_perms=20, seed=0):
    """X: 2D bool/float ndarray, y: 0/1 ndarray. Returns dict with cv AUROC (LR/tree), permuted
    null AUROC array, and top |coef|/importance features (already restricted upstream to
    >=20-carrier features, so nothing further to mask here)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    from sklearn.tree import DecisionTreeClassifier

    n_pos, n = int(y.sum()), len(y)
    if n_pos < MIN_CARRIERS_FOR_FEATURE or (n - n_pos) < MIN_CARRIERS_FOR_FEATURE:
        return None  # underpowered / disclosive class size -- caller skips this task entirely

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    lr = LogisticRegression(penalty="l1", solver="liblinear", C=1.0, max_iter=2000)
    tree = DecisionTreeClassifier(max_depth=3, random_state=seed)

    lr_oof = cross_val_predict(lr, X, y, cv=skf, method="predict_proba")[:, 1]
    tree_oof = cross_val_predict(tree, X, y, cv=skf, method="predict_proba")[:, 1]
    lr_auc = roc_auc_score(y, lr_oof)
    tree_auc = roc_auc_score(y, tree_oof)

    rng = np.random.default_rng(seed)
    perm_aucs = []
    for _ in range(n_perms):
        y_perm = rng.permutation(y)
        try:
            oof = cross_val_predict(lr, X, y_perm, cv=skf, method="predict_proba")[:, 1]
            perm_aucs.append(roc_auc_score(y_perm, oof))
        except ValueError:
            continue  # a degenerate fold (single class) under permutation -- skip that draw

    lr.fit(X, y)
    coefs = lr.coef_.ravel()
    top_idx = np.argsort(-np.abs(coefs))[:15]
    top_lr = [(feature_names[i], float(coefs[i])) for i in top_idx if coefs[i] != 0]

    tree.fit(X, y)
    imp = tree.feature_importances_
    top_idx_t = np.argsort(-imp)[:15]
    top_tree = [(feature_names[i], float(imp[i])) for i in top_idx_t if imp[i] > 0]

    perm_arr = np.array(perm_aucs) if perm_aucs else np.array([np.nan])
    p_value = float((np.sum(perm_arr >= lr_auc) + 1) / (len(perm_arr) + 1)) if perm_aucs else None

    return {
        "n_total": n, "n_pos": n_pos,
        "lr_auroc": float(lr_auc), "tree_auroc": float(tree_auc),
        "n_perms": len(perm_aucs),
        "perm_auroc_mean": float(np.nanmean(perm_arr)),
        "perm_auroc_p95": float(np.nanpercentile(perm_arr, 95)) if perm_aucs else None,
        "perm_pvalue_lr": p_value,
        "top_lr_features": top_lr,
        "top_tree_features": top_tree,
    }


# ---------------------------------------------------------------------------
# PCA -- aggregate-only exports.
# ---------------------------------------------------------------------------
def run_pca(X, people, ancestry_of, platform_of, out_dir, n_bins=20):
    from sklearn.decomposition import PCA
    Xc = X.astype(float)
    Xc = Xc - Xc.mean(axis=0, keepdims=True)
    n_comp = min(2, Xc.shape[1]) if Xc.shape[1] > 0 else 0
    if n_comp < 2 or Xc.shape[0] < 2:
        log("[47] PCA skipped: not enough features/people.")
        return
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(Xc)  # NEVER written per-person -- aggregated below only.

    anc = np.array([ancestry_of.get(p, "UNASSIGNED") for p in people])
    plat = np.array([platform_of.get(p, "UNKNOWN") for p in people])

    # Centroids per (ancestry, platform), n>=20 only.
    rows = []
    for a in sorted(set(anc)):
        for pl in sorted(set(plat)):
            mask = (anc == a) & (plat == pl)
            n = int(mask.sum())
            if n < SUPPRESS_BELOW:
                continue
            rows.append({"ancestry": a, "platform": pl, "n": n,
                        "pc1_mean": float(coords[mask, 0].mean()),
                        "pc2_mean": float(coords[mask, 1].mean())})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "pca_group_centroids.tsv"),
                              sep="\t", index=False)

    # 2D binned density grid over PC1/PC2, cells 1-19 suppressed to '<20', true 0 stays 0.
    x, y = coords[:, 0], coords[:, 1]
    if len(x) >= 2 and np.ptp(x) > 0 and np.ptp(y) > 0:
        hist, xedges, yedges = np.histogram2d(x, y, bins=n_bins)
        grid_rows = []
        for i in range(hist.shape[0]):
            for j in range(hist.shape[1]):
                c = int(hist[i, j])
                grid_rows.append({
                    "pc1_bin_lo": float(xedges[i]), "pc1_bin_hi": float(xedges[i + 1]),
                    "pc2_bin_lo": float(yedges[j]), "pc2_bin_hi": float(yedges[j + 1]),
                    "count": suppressed(c),
                })
        pd.DataFrame(grid_rows).to_csv(os.path.join(out_dir, "pca_density_grid.tsv"),
                                       sep="\t", index=False)
    log(f"[47] PCA: explained_variance_ratio_ = {pca.explained_variance_ratio_.tolist()}")


# ---------------------------------------------------------------------------
# Synthetic cohort generator (local dry run, no VM/no data/no network).
# ---------------------------------------------------------------------------
def make_synthetic(n_people=400, seed=20260926):
    rng = np.random.default_rng(seed)
    pids = [f"synthP{i:05d}" for i in range(n_people)]
    ancestries = rng.choice(ANCESTRY_ORDER, size=n_people, p=[0.2, 0.15, 0.15, 0.3, 0.1, 0.1])
    platforms = rng.choice(["revio", "sequel2e", "sequel2"], size=n_people, p=[0.6, 0.3, 0.1])

    genes = CLASSICAL_GENES
    n_alleles_per_gene = 6
    table1_rows = []
    kir_content = {}
    for pid, anc in zip(pids, ancestries):
        for gene in genes:
            for hap in ("hap1", "hap2"):
                # Planted ancestry <-> allele association: allele index skewed by ancestry.
                anc_idx = ANCESTRY_ORDER.index(anc)
                weights = np.ones(n_alleles_per_gene)
                weights[anc_idx % n_alleles_per_gene] += 3.0
                weights /= weights.sum()
                a_idx = rng.choice(n_alleles_per_gene, p=weights)
                consensus = f"{gene}*{a_idx + 1:02d}:01"
                table1_rows.append({"person_id": pid, "hap": hap, "gene": gene,
                                    "consensus": consensus})
        # No planted HLA<->KIR association: cA/cB assigned independently of HLA/ancestry.
        kir_content[pid] = "cB" if rng.random() < 0.4 else "cA"

    table1 = pd.DataFrame(table1_rows)
    cohort = pd.DataFrame({"person_id": pids, "ancestry_pred": ancestries, "platform": platforms})
    return table1, cohort, kir_content


class _FakeKirModule:
    """Minimal stand-in exposing only what build_kir_features needs, for the synthetic path
    (avoids depending on real hap*.gtf.gz files / importlib-loading 41_kir_pilot.py for a
    pure-synthetic smoke test)."""
    KIR_GENES = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR3DL1", "KIR3DS1"]
    B_CONTENT_GENES = {"KIR2DL2", "KIR2DS1"}


def build_kir_features_synthetic(pids, kir_content, seed=0):
    rng = np.random.default_rng(seed)
    genes = _FakeKirModule.KIR_GENES
    presence = pd.DataFrame(False, index=pids, columns=genes)
    content = pd.Series(index=pids, dtype=object)
    for pid in pids:
        c = kir_content[pid]
        content.loc[pid] = c
        for g in genes:
            base_p = 0.5
            if c == "cB" and g in _FakeKirModule.B_CONTENT_GENES:
                base_p = 0.9
            presence.at[pid, g] = rng.random() < base_p
    presence.index.name = "person_id"
    return presence, content


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def run_pipeline(table1_df, cohort_df, kir_presence, kir_content, out_dir, platform_col,
                 n_perms, pca_bins, status_path):
    os.makedirs(out_dir, exist_ok=True)
    ancestry_of = dict(zip(cohort_df["person_id"].astype(str),
                          cohort_df["ancestry_pred"].astype(str).str.upper()))
    if platform_col not in cohort_df.columns:
        log(f"[47] FATAL: --platform-col {platform_col!r} not in cohort_membership.tsv. "
            f"Available columns: {list(cohort_df.columns)}")
        sys.exit(2)
    platform_of = dict(zip(cohort_df["person_id"].astype(str),
                           cohort_df[platform_col].astype(str)))

    hla_mat, n_dropped_hla = build_hla_carriage(table1_df)
    people = sorted(set(hla_mat.index) & set(kir_presence.index))
    hla_mat = hla_mat.loc[people]
    kir_mat = kir_presence.loc[people]
    combined = pd.concat([hla_mat, kir_mat], axis=1)
    feature_names = list(combined.columns)
    X_all = combined.values.astype(float)

    results = []

    # Task 1: ancestry one-vs-rest.
    y_anc = np.array([ancestry_of.get(p, "UNASSIGNED") for p in people])
    for anc in ANCESTRY_ORDER:
        y = (y_anc == anc).astype(int)
        r = run_binary_task(X_all, y, feature_names, n_perms=n_perms)
        if r is None:
            log(f"[47] task=ancestry_{anc}: skipped (underpowered, n_pos or n_neg < 20)")
            continue
        r["task"] = f"ancestry_{anc}"
        results.append(r)

    # Task 2: sequencing platform, one-vs-rest per observed platform value.
    y_plat = np.array([platform_of.get(p, "UNKNOWN") for p in people])
    for plat in sorted(set(y_plat)):
        y = (y_plat == plat).astype(int)
        r = run_binary_task(X_all, y, feature_names, n_perms=n_perms)
        if r is None:
            log(f"[47] task=platform_{plat}: skipped (underpowered)")
            continue
        r["task"] = f"platform_{plat}"
        results.append(r)

    # Task 3: cA/cB from HLA only (KIR features excluded from X for this task by design).
    content_of = kir_content if isinstance(kir_content, dict) else kir_content.to_dict()
    y_content = np.array([content_of.get(p, "missing") for p in people])
    valid = y_content != "missing"
    if valid.sum() >= 2 * MIN_CARRIERS_FOR_FEATURE:
        y = (y_content[valid] == "cB").astype(int)
        Xh = hla_mat.loc[[p for p, v in zip(people, valid) if v]].values.astype(float)
        r = run_binary_task(Xh, y, list(hla_mat.columns), n_perms=n_perms)
        if r is not None:
            r["task"] = "cB_from_HLA"
            results.append(r)
        else:
            log("[47] task=cB_from_HLA: skipped (underpowered)")
    else:
        log("[47] task=cB_from_HLA: skipped (too few people with valid KIR content)")

    # Write results TSV + top-features TSV (features already >=20-carrier filtered).
    metric_rows, feat_rows = [], []
    for r in results:
        metric_rows.append({
            "task": r["task"], "n_total": r["n_total"], "n_pos": suppressed(r["n_pos"]),
            "lr_auroc": round(r["lr_auroc"], 4), "tree_auroc": round(r["tree_auroc"], 4),
            "n_perms": r["n_perms"],
            "perm_auroc_mean": round(r["perm_auroc_mean"], 4) if r["perm_auroc_mean"] == r["perm_auroc_mean"] else "",
            "perm_auroc_p95": round(r["perm_auroc_p95"], 4) if r["perm_auroc_p95"] is not None else "",
            "perm_pvalue_lr": round(r["perm_pvalue_lr"], 4) if r["perm_pvalue_lr"] is not None else "",
        })
        for feat, coef in r["top_lr_features"]:
            feat_rows.append({"task": r["task"], "model": "l1_logreg", "feature": feat,
                              "weight": round(coef, 4)})
        for feat, imp in r["top_tree_features"]:
            feat_rows.append({"task": r["task"], "model": "tree_depth3", "feature": feat,
                              "weight": round(imp, 4)})
    pd.DataFrame(metric_rows).to_csv(os.path.join(out_dir, "naive_ml_metrics.tsv"),
                                     sep="\t", index=False)
    pd.DataFrame(feat_rows).to_csv(os.path.join(out_dir, "naive_ml_top_features.tsv"),
                                   sep="\t", index=False)

    run_pca(X_all, people, ancestry_of, platform_of, out_dir, n_bins=pca_bins)

    with open(status_path, "w") as f:
        f.write(f"done n_people={len(people)} n_features={len(feature_names)} "
                f"n_hla_features={hla_mat.shape[1]} n_kir_features={kir_mat.shape[1]} "
                f"n_dropped_rare_hla_alleles={n_dropped_hla} n_tasks_run={len(results)}\n")
    log(f"[47] wrote {len(results)} task rows -> {out_dir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--cohort-membership",
                    default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table",
                    default=os.path.expanduser(
                        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/"
                        "aux/relatedness/samples_relatedness.tsv"))
    ap.add_argument("--kir-outroot", default=os.path.expanduser("~/pipeline_outputs_kir"))
    ap.add_argument("--kir-pilot-script",
                    default=os.path.expanduser("~/s03/41_kir_pilot.py"))
    ap.add_argument("--platform-col", default="platform",
                    help="Column in cohort_membership.tsv holding sequencing platform. If "
                         "absent, the script prints all available columns and exits.")
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s04/results/47"))
    ap.add_argument("--limit", type=int, default=None,
                    help="Restrict to the first N people (sorted, deterministic) for a smoke test.")
    ap.add_argument("--n-perms", type=int, default=20)
    ap.add_argument("--pca-bins", type=int, default=20)
    ap.add_argument("--synthetic", action="store_true",
                    help="Local dry run: generate an in-memory synthetic cohort instead of "
                         "reading VM paths.")
    ap.add_argument("--n-people", type=int, default=400, help="--synthetic only.")
    ap.add_argument("--seed", type=int, default=20260926)
    args = ap.parse_args()

    t0 = time.perf_counter()
    os.makedirs(args.out_dir, exist_ok=True)
    status_path = os.path.join(args.out_dir, "STATUS.txt")

    if args.synthetic:
        table1_df, cohort_df, kir_content_true = make_synthetic(args.n_people, args.seed)
        pids = sorted(cohort_df["person_id"])
        kir_presence, kir_content = build_kir_features_synthetic(pids, kir_content_true, args.seed)
        run_pipeline(table1_df, cohort_df, kir_presence, kir_content, args.out_dir,
                    args.platform_col, args.n_perms, args.pca_bins, status_path)
        log(f"[47] synthetic run done in {time.perf_counter()-t0:.0f}s")
        return

    from _viz_common import load_table1, load_cohort_membership  # local import: only needed live
    table1_df = load_table1(args.table1)
    if args.limit:
        keep = sorted(table1_df["person_id"].unique())[:args.limit]
        table1_df = table1_df[table1_df["person_id"].isin(keep)]
    cohort_df = load_cohort_membership(args.cohort_membership)

    kir_mod = load_kir_module(args.kir_pilot_script)
    pairs = kir_mod.load_relatedness_pairs(args.relatedness_table)
    all_pids = sorted(set(table1_df["person_id"]) | set(cohort_df["person_id"].astype(str)))
    kept, removed = kir_mod.greedy_unrelated(all_pids, pairs, kin_min=kir_mod.KIN_MIN)
    log(f"[47] unrelated subset: {len(kept)}/{len(all_pids)} ({len(removed)} relatives dropped)")
    table1_df = table1_df[table1_df["person_id"].isin(kept)]
    unrelated_pids = sorted(kept & set(table1_df["person_id"]))

    kir_presence, kir_content = build_kir_features(unrelated_pids, args.kir_outroot, kir_mod)

    run_pipeline(table1_df, cohort_df, kir_presence, kir_content, args.out_dir,
                args.platform_col, args.n_perms, args.pca_bins, status_path)
    log(f"[47] done in {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
