#!/usr/bin/env python3
"""07 -- Which embedding setting is actually best? (report: reports/07_embedding_settings/)

We adopted SCEPTR's beta-chain variant on TRBV + CDR3 after a two-model comparison whose
only criterion was a V-gene structure check -- and that check is circular once the V gene is
an input. This sweeps the settings properly and scores each one on what we care about, with
two independent kinds of evidence:

  CLONOTYPE LEVEL, external label. VDJdb TCRs of known epitope: do TCRs recognising the same
  epitope embed closer together? Pairwise AUROC and 5-nearest-neighbour balanced accuracy.
  This is the only axis with a real biological label, but note SCEPTR's training repertoires
  may overlap public TCR databases, so absolute values are optimistic for every SCEPTR
  variant equally; the ranking is the result.

  PERSON LEVEL, no label needed. Two poolings, both scored with z-scored split-half
  identifiability (does a person's own second half come back as their nearest neighbour,
  features standardised first; report 06) and out-of-fold age R^2, plus log-depth R^2 as
  the nuisance axis:
    mean pooling        the embedding alone -- the column that actually compares embeddings
    mean + TRBV usage   block-balanced (repfig.balanced_concat), report 06's best family;
                        TRBV usage is identical for every setting, so differences here are
                        diluted by construction
  (Until 2026-10-04 this axis used k=100 CLR cluster profiles, which scored 0.5% in report
  06 and left every setting near the floor; those columns are gone.)

  Epitope scores carry 95% CIs from an epitope-level bootstrap (repfig.epitope_bootstrap,
  500 replicates) and paired CIs for the difference from b_sceptr on the same replicates.

Settings swept (whatever the installed sceptr exposes; failures are reported, not fatal):
  input          b_sceptr on TRBV+CDR3  vs  cdr3_only  vs  the paired-chain default fed
                 beta only -- i.e. how much the V gene and the model's chain assumption matter
  capacity       tiny (16-d) / small (32) / default (64) / large (128)
  training       mlm_only (no autocontrastive stage), synthetic_data (trained on OLGA-
                 generated sequences), shuffled_data (chain pairing randomised)
  architecture   average_pooling instead of <cls>, blosum tokenisation
  baselines      CDR3 3-mer composition (SVD) and TRBV one-hot -- not learned at all

Reported against dimension and wall time, because a 16-d embedding that matches a 128-d one
is the better choice for everything downstream.

Usage (from aleix/RNA-seq/, after 02 and 06):
  pixi run python3 -u scripts/07_embedding_settings.py [--n-people 1000] [--per-person 250]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial.distance import pdist  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import repfig as R  # noqa: E402

VARIANTS = [
    ("b_sceptr (TRBV+CDR3)", "b_sceptr", True),
    ("default (beta only)", "default", True),
    ("cdr3_only", "cdr3_only", False),
    ("tiny (16-d)", "tiny", True),
    ("small (32-d)", "small", True),
    ("large (128-d)", "large", True),
    ("average_pooling", "average_pooling", True),
    ("blosum", "blosum", True),
    ("mlm_only", "mlm_only", True),
    ("synthetic_data", "synthetic_data", True),
    ("shuffled_data", "shuffled_data", True),
]
MIN_TCR_PER_EPITOPE = 30
MAX_TCR_PER_EPITOPE = 300
SVD_DIM = 100


def epitope_scores(X, epi, rng, max_pairs=2_000_000):
    """Pairwise AUROC (same- vs different-epitope) and 5-NN balanced accuracy, plus the
    distance matrix and per-epitope recall that the bootstrap reuses."""
    from scipy.spatial.distance import cdist
    D = cdist(X, X)
    same = epi[:, None] == epi[None, :]
    iu = np.triu_indices(len(epi), 1)
    take = rng.choice(len(iu[0]), size=min(max_pairs, len(iu[0])), replace=False)
    i, j = iu[0][take], iu[1][take]
    auc = roc_auc_score(same[i, j], -D[i, j] + rng.normal(0, 1e-9, len(i)))
    Dk = D + np.diag(np.full(len(epi), np.inf))
    nn = np.argsort(Dk + rng.uniform(0, 1e-9, Dk.shape), axis=1)[:, :5]
    pred = np.array([pd.Series(epi[r]).value_counts().idxmax() for r in nn])
    cats = sorted(set(epi))
    recall = [np.mean(pred[epi == e] == e) for e in cats]
    return auc, float(np.mean(recall)), D, dict(zip(cats, recall))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="~/pipeline_outputs/rnaseq/cohort_full.tsv")
    ap.add_argument("--pheno-dir", default="~/pipeline_outputs/rnaseq/pheno")
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas")
    ap.add_argument("--vdjdb", default=None)
    ap.add_argument("--n-people", type=int, default=1000)
    ap.add_argument("--per-person", type=int, default=250,
                    help="top clonotypes per person by read support")
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/07_embedding_settings")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    import importlib
    m06 = importlib.import_module("06_person_representation")
    from embed_cdr3s import accepted_trbv, embed_sceptr, sceptr_model
    from prep_vdjdb_pool import DEFAULT_URL
    load_vdjdb = importlib.import_module("04_antigen_specificity").load_vdjdb

    # ---------------- cohort slice ----------------
    clono = pd.read_pickle(os.path.join(os.path.expanduser(args.cache_dir),
                                        "trb_clonotypes_full.pkl"))
    clono["research_id"] = clono["research_id"].astype(str)
    clono["trbv"] = clono["trbv"].astype(str)
    clono = clono[clono["trbv"] != ""]
    people = R.load_people(args.cohort, args.pheno_dir).dropna(subset=["age", "sex"])
    people = people[people["age"] >= 18]
    n_cl = clono.groupby("research_id").size()
    people = people[people["research_id"].isin(n_cl[n_cl >= args.per_person].index)]
    if args.n_people and len(people) > args.n_people:
        people = people.sample(args.n_people, random_state=args.seed)
    people = people.sort_values("research_id").reset_index(drop=True)
    pid = pd.Series(np.arange(len(people)), index=people["research_id"])
    clono = clono[clono["research_id"].isin(pid.index)]
    # rank by reads, ties broken by the same seeded per-person key as embed_cdr3s.py
    clono["_tb"] = pd.util.hash_pandas_object(
        clono["research_id"] + "|" + clono["trbv"] + "|" + clono["cdr3aa"].astype(str),
        index=False).to_numpy()
    clono = (clono.sort_values(["reads", "_tb"], ascending=[False, True], kind="stable")
                  .groupby("research_id", observed=True).head(args.per_person)
                  .drop(columns="_tb").reset_index(drop=True))
    person_idx = pid.reindex(clono["research_id"]).to_numpy()
    n_people = len(people)
    age = people["age"].to_numpy(float)
    depth = np.log(clono.groupby("research_id")["reads"].sum()
                        .reindex(people["research_id"]).to_numpy())
    half = rng.random(len(clono)) < 0.5
    print(f"{n_people:,} people x up to {args.per_person} clonotypes = {len(clono):,} "
          f"sequences per variant", file=sys.stderr)

    # ---------------- VDJdb benchmark set ----------------
    vdj = load_vdjdb(os.path.expanduser(args.vdjdb) if args.vdjdb else None, DEFAULT_URL,
                     os.path.expanduser(args.cache_dir))
    uni = vdj[vdj.groupby(["trbv", "cdr3aa"])["epitope"].transform("nunique") == 1]
    cnt = uni["epitope"].value_counts()
    bench = (uni[uni["epitope"].isin(cnt[cnt >= MIN_TCR_PER_EPITOPE].index)]
             .sample(frac=1, random_state=args.seed)
             .groupby("epitope").head(MAX_TCR_PER_EPITOPE).reset_index(drop=True))
    # Same TCRs for every setting (TRBV accepted by b_sceptr), so bootstrap replicates are
    # shared and differences between settings are paired.
    bench = bench[bench["trbv"].isin(accepted_trbv(sceptr_model("v+cdr3"), bench["trbv"]))]
    cnt = bench["epitope"].value_counts()
    bench = bench[bench["epitope"].isin(cnt[cnt >= MIN_TCR_PER_EPITOPE].index)]
    bench = bench.reset_index(drop=True)
    boot = R.epitope_bootstrap(bench["epitope"].to_numpy(), seed=args.seed)
    ref_reps = {}
    vcode_all, _ = pd.factorize(clono["trbv"])
    print(f"Benchmark: {len(bench):,} TCRs, {bench['epitope'].nunique()} epitopes",
          file=sys.stderr)

    # ---------------- evaluate each setting ----------------
    def evaluate(name, Xc, pi, hm, vc, Xb, bidx, seconds, dim):
        """Xc: cohort clonotype vectors, with person index `pi`, split-half mask `hm` and
        TRBV codes `vc`. Xb: benchmark vectors for bench rows `bidx`."""
        row = dict(setting=name, dim=dim, seconds=seconds)
        if Xb is not None and len(Xb) > 50:
            epi = bench["epitope"].to_numpy()[bidx]
            auc, knn, D, rec = epitope_scores(Xb, epi, rng)
            paired = len(bidx) == len(bench)
            b_ = boot if paired else R.epitope_bootstrap(epi, seed=args.seed)
            ra, rk = R.boot_scores(D, rec, b_)
            row.update(epitope_auroc=auc, epitope_knn5=knn)
            row["epitope_auroc_lo"], row["epitope_auroc_hi"] = R.ci(ra)
            row["epitope_knn5_lo"], row["epitope_knn5_hi"] = R.ci(rk)
            if paired and not ref_reps:
                ref_reps.update(name=name, auc=ra, knn=rk)
            if paired and ref_reps:
                row["auroc_diff_vs_ref_lo"], row["auroc_diff_vs_ref_hi"] = R.ci(ra - ref_reps["auc"])
                row["knn5_diff_vs_ref_lo"], row["knn5_diff_vs_ref_hi"] = R.ci(rk - ref_reps["knn"])
                row["reference_setting"] = ref_reps["name"]
            del D

        def pooled(m):
            mean = m06.mean_pool(Xc[m], pi[m], n_people)
            trbv = m06.profile(vc[m], vc.max() + 1, pi[m], n_people)
            return mean, trbv

        (mA, tA), (mB, tB), (mF, tF) = pooled(hm), pooled(~hm), pooled(np.ones(len(Xc), bool))
        row["id_meanpool_z"], row["id_meanpool_z_pct_rank"] = m06.identifiability(
            R.zscore_like(mA, mF), R.zscore_like(mB, mF))
        cA, cB = R.balanced_concat([mA, tA], [mF, tF]), R.balanced_concat([mB, tB], [mF, tF])
        cF = R.balanced_concat([mF, tF], [mF, tF])
        row["id_combined_z"], _ = m06.identifiability(R.zscore_like(cA, cF),
                                                      R.zscore_like(cB, cF))
        row["age_r2_meanpool"], row["age_mae_meanpool"] = m06.cv_regress(mF, age, args.seed)
        row["depth_r2_meanpool"], _ = m06.cv_regress(mF, depth, args.seed)
        row["age_r2_combined"], _ = m06.cv_regress(cF, age, args.seed)
        print(f"  {name:24s} dim {dim:4d}  {seconds:6.0f}s  epitope AUROC "
              f"{row.get('epitope_auroc', np.nan):.3f}  ID(mean,z) "
              f"{100 * row['id_meanpool_z']:5.1f}%  ID(+TRBV,z) {100 * row['id_combined_z']:5.1f}%"
              f"  age R2 {row['age_r2_meanpool']:+.3f}", file=sys.stderr)
        return row

    rows = []
    for label, fn, use_v in VARIANTS:
        try:
            from sceptr import variant
            model = getattr(variant, fn)()
        except Exception as e:
            print(f"  !! {label}: unavailable ({type(e).__name__}: {e})", file=sys.stderr)
            continue
        try:
            if use_v:
                good = accepted_trbv(model, pd.concat([clono["trbv"], bench["trbv"]]))
                cm = clono["trbv"].isin(good).to_numpy()
                bm = bench["trbv"].isin(good).to_numpy()
                if cm.mean() < 0.5:
                    raise RuntimeError(f"only {cm.mean():.0%} of TRBV symbols accepted")
                dfc = pd.DataFrame({"TRBV": clono["trbv"][cm].to_numpy(),
                                    "CDR3B": clono["cdr3aa"][cm].to_numpy()})
                dfb = pd.DataFrame({"TRBV": bench["trbv"][bm].to_numpy(),
                                    "CDR3B": bench["cdr3aa"][bm].to_numpy()})
            else:
                cm = np.ones(len(clono), bool)
                bm = np.ones(len(bench), bool)
                dfc = pd.DataFrame({"CDR3B": clono["cdr3aa"].to_numpy()})
                dfb = pd.DataFrame({"CDR3B": bench["cdr3aa"].to_numpy()})
            t0 = time.time()
            Xc, _ = embed_sceptr(model, dfc)
            secs = time.time() - t0
            Xb, _ = embed_sceptr(model, dfb)
        except Exception as e:
            print(f"  !! {label}: failed ({type(e).__name__}: {e})", file=sys.stderr)
            continue
        rows.append(evaluate(label, Xc, person_idx[cm], half[cm], vcode_all[cm], Xb,
                             np.flatnonzero(bm), secs, Xc.shape[1]))
        del Xc, Xb

    # ---------------- non-learned baselines ----------------
    from sklearn.decomposition import TruncatedSVD
    from scipy.sparse import csr_matrix
    for label in ("CDR3 3-mer (SVD)", "TRBV one-hot"):
        t0 = time.time()
        if label.startswith("CDR3"):
            seqs = pd.concat([clono["cdr3aa"], bench["cdr3aa"]]).to_numpy()
            vocab, r_, c_ = {}, [], []
            for i, s in enumerate(seqs):
                for j in range(len(s) - 2):
                    c_.append(vocab.setdefault(s[j:j + 3], len(vocab)))
                    r_.append(i)
            M = csr_matrix((np.ones(len(r_)), (r_, c_)), shape=(len(seqs), len(vocab)))
            M = M.multiply(1 / np.maximum(np.sqrt(M.multiply(M).sum(1)).A1, 1e-9)[:, None])
            Z = TruncatedSVD(SVD_DIM, random_state=args.seed).fit_transform(M.tocsr())
        else:
            cats = sorted(set(clono["trbv"]) | set(bench["trbv"]))
            idx = {v: i for i, v in enumerate(cats)}
            allv = pd.concat([clono["trbv"], bench["trbv"]]).map(idx).to_numpy()
            Z = np.eye(len(cats), dtype=np.float32)[allv]
        secs = time.time() - t0
        rows.append(evaluate(label, Z[:len(clono)], person_idx, half, vcode_all,
                             Z[len(clono):], np.arange(len(bench)), secs, Z.shape[1]))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(outdir, "embedding_settings.csv"), index=False)

    # ---------------- figure ----------------
    import cnsplots as cns
    mp = cns.multipanel(max_width=540)
    pal = R.nature()
    W, H, GAP = 105, 105, 34

    def group(n):
        if n in ("CDR3 3-mer (SVD)", "TRBV one-hot"):
            return "not learned"
        if any(t in n for t in ("tiny", "small", "large")):
            return "capacity"
        if any(t in n for t in ("mlm_only", "synthetic", "shuffled")):
            return "training ablation"
        if any(t in n for t in ("average_pooling", "blosum")):
            return "architecture"
        return "input"
    gcol = {"input": pal[0], "capacity": pal[1], "training ablation": pal[2],
            "architecture": pal[4], "not learned": "0.65"}
    res["group"] = res["setting"].map(group)

    def barh(ax, col, title, xlabel, ref=None, sub=None):
        r = res.dropna(subset=[col]).sort_values(col)
        y = np.arange(len(r))
        left = ref if ref is not None else 0
        ax.barh(y, r[col] - left, left=left, height=0.68,
                color=[gcol[g] for g in r["group"]], linewidth=0)
        if ref is not None:
            ax.axvline(ref, color="0.5", lw=0.5, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels(r["setting"], fontsize=4.5)
        ax.set_title(title, pad=11 if sub else 3)
        if sub:
            ax.text(0.5, 1.02, sub, transform=ax.transAxes, ha="center", va="bottom",
                    fontsize=5.5)
        R.style(ax, None, xlabel, None)

    ax = mp.panel("a", width=145, height=H + 25, margin_right=12, margin_bottom=GAP)
    barh(ax, "epitope_auroc", "Known-epitope structure", "Pairwise AUROC", ref=0.5,
         sub=f"{len(bench):,} VDJdb TCRs, {bench['epitope'].nunique()} epitopes, 95% CI")
    r_ = res.dropna(subset=["epitope_auroc"]).sort_values("epitope_auroc")
    ax.errorbar(r_["epitope_auroc"], np.arange(len(r_)),
                xerr=[np.maximum(r_["epitope_auroc"] - r_["epitope_auroc_lo"], 0),
                      np.maximum(r_["epitope_auroc_hi"] - r_["epitope_auroc"], 0)],
                fmt="none", ecolor="0.2", elinewidth=0.5, capsize=1.2)

    ax = mp.panel("b", width=W, height=H + 25, margin_right=12, margin_bottom=GAP)
    res["id_meanpool_z_pct"] = 100 * res["id_meanpool_z"]
    barh(ax, "id_meanpool_z_pct", "Split-half identity", "Correct (%)",
         sub="mean pooling, z-scored")
    ax.set_yticklabels([])

    ax = mp.panel("c", width=W, height=H + 25, margin_right=55, margin_bottom=GAP)
    barh(ax, "age_r2_meanpool", "Age prediction", "R² (mean pooling)")
    ax.set_yticklabels([])
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=g) for g, c in gcol.items()], frameon=False)
    cns.take_legend_out(title="Setting type", ax=ax)

    ax = mp.panel("d", width=W, height=H, margin_right=12)
    for g_, s_ in res.groupby("group"):
        ax.scatter(s_["dim"], s_["epitope_auroc"], s=16, linewidths=0, color=gcol[g_])
    for _, r_ in res.iterrows():
        if pd.notna(r_["epitope_auroc"]):
            ax.annotate(r_["setting"].split(" (")[0], (r_["dim"], r_["epitope_auroc"]),
                        fontsize=4, xytext=(2, 2), textcoords="offset points")
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_title("Cost vs quality", pad=11)
    ax.text(0.5, 1.02, "smaller and higher is better", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.5)
    R.style(ax, None, "Embedding dimension", "Epitope AUROC")

    ax = mp.panel("e", width=W, height=H, margin_right=12)
    ok = res.dropna(subset=["epitope_auroc"])
    ax.scatter(ok["epitope_auroc"], ok["id_meanpool_z"] * 100, s=16, linewidths=0,
               color=[gcol[g] for g in ok["group"]])
    rho = spearmanr(ok["epitope_auroc"], ok["id_meanpool_z"])[0] if len(ok) > 2 else np.nan
    ax.set_title("Do the two axes agree?", pad=11)
    ax.text(0.5, 1.02, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=5.5)
    R.style(ax, None, "Epitope AUROC (clonotype level)", "Identifiability (%)")

    ax = mp.panel("f", width=W, height=H, margin_right=12)
    r = res.dropna(subset=["id_combined_z"]).sort_values("id_combined_z")
    y = np.arange(len(r))
    ax.barh(y - 0.2, 100 * r["id_meanpool_z"], height=0.38, color=pal[3],
            linewidth=0, label="mean pooling")
    ax.barh(y + 0.2, 100 * r["id_combined_z"], height=0.38, color=pal[0],
            linewidth=0, label="mean + TRBV usage")
    ax.set_yticks(y)
    ax.set_yticklabels(r["setting"], fontsize=4.5)
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    ax.set_title("Pooling vs embedding", pad=3)
    R.style(ax, None, "Identifiability (%)", None)
    R.save("fig_embedding_settings", outdir)

    summary = pd.DataFrame({
        "metric": ["people", "clonotypes_per_variant", "benchmark_tcrs", "settings_ok",
                   "best_epitope_auroc", "best_identifiability", "best_age_r2"],
        "value": [n_people, len(clono), len(bench), len(res),
                  res.loc[res["epitope_auroc"].idxmax(), "setting"]
                  if res["epitope_auroc"].notna().any() else "none",
                  res.loc[res["id_meanpool_z"].idxmax(), "setting"],
                  res.loc[res["age_r2_meanpool"].idxmax(), "setting"]]})
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + res.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
