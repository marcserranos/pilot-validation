#!/usr/bin/env python3
"""02 -- Population-scale TRB repertoire atlas (report: reports/02_repertoire_atlas/).

What 7,922 whole-blood TRB repertoires recovered from bulk RNA-seq look like, and one
textbook check that they carry real immunology:

  a  cohort by genetic ancestry (n, % female, median age)
  b  TRB clonotypes recovered per person (full repertoire, not the top-500 embedding sample)
  c  CDR3 length distribution by ancestry
  d  TRBV gene usage across people (median, IQR, 5-95%)
  e  sharing spectrum: how many clonotypes are carried by exactly k people
  f  depth-corrected diversity (rarefied richness) vs age
  g  clonal expansion (fraction of TRB reads in the top 10 clonotypes) vs age

Known biology tested in f/g: TCR repertoire diversity falls and clonal expansion rises with
age (Britanova et al. 2014 J Immunol, n=65; Qi et al. 2014 PNAS). Richness is rarefied
(Hurlbert: expected distinct clonotypes in a fixed-size read subsample) so sequencing depth
and blood T-cell content cannot masquerade as diversity. Age effects are estimated by OLS
adjusted for sex and ancestry; a sensitivity model also adjusts for log total TRB reads.

Same clonotype definition as embed_cdr3s.py (imported, not re-implemented): cdr3.out,
CDR3_score >= 0.02, canonical junction, <= 30 aa, TRB, unique (TRBV gene, CDR3aa) per person
with reads summed. Age = 2020 - year_of_birth (draws span ~2017-2023: +/- ~3 y, fine for
decade-scale trends).

Outputs: figures + de-identified CSVs to --outdir (all plotted groups/bins >= 20 people);
VM-local clonotype cache (real research_ids, never commit) reused by scripts 03 and 04.

Usage (from aleix/RNA-seq/):
  pixi run python3 -u scripts/02_repertoire_atlas.py [--workers 16]
"""
import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import repfig as R  # noqa: E402
from embed_cdr3s import load_person_cdr3s  # noqa: E402

TOP_V_GENES = 24


def load_one(args):
    pheno_dir, rid = args
    df, source = load_person_cdr3s(pheno_dir, rid, 0.02, False, False, 30)
    if df is None:
        return rid, None
    df = df[df["chain"] == "TRB"]
    trbv = df["v_gene"].astype(str).str.split("*").str[0].replace({"nan": ""})
    out = (pd.DataFrame({"trbv": trbv.to_numpy(), "cdr3aa": df["cdr3aa"].to_numpy(),
                         "reads": df["reads"].to_numpy()})
             .groupby(["trbv", "cdr3aa"], as_index=False)["reads"].sum())
    out["reads"] = np.maximum(np.rint(out["reads"]), 1).astype(np.int64)
    return rid, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="~/pipeline_outputs/rnaseq/cohort_full.tsv")
    ap.add_argument("--pheno-dir", default="~/pipeline_outputs/rnaseq/pheno",
                    help="demographics cache (person.tsv) from query_overlap_phenotypes.py")
    ap.add_argument("--trust4-dir", default="~/pipeline_outputs/rnaseq",
                    help="dir with <research_id>/ TRUST4 outputs")
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas",
                    help="VM-local clonotype cache for scripts 03/04 (real research_ids)")
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/02_repertoire_atlas")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    cache_dir = os.path.expanduser(args.cache_dir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

    people = R.load_people(args.cohort, args.pheno_dir)
    print(f"Cohort: {len(people):,} people; demographics for "
          f"{people['age'].notna().sum():,} (age) / {people['sex'].notna().sum():,} (sex)",
          file=sys.stderr)

    # ---------------- load every full TRB repertoire ----------------
    cache = os.path.join(cache_dir, "trb_clonotypes_full.pkl")
    if os.path.exists(cache):
        clono = pd.read_pickle(cache)
        print(f"Loaded clonotype cache {cache}", file=sys.stderr)
    else:
        t0 = time.time()
        jobs = [(os.path.expanduser(args.trust4_dir), r) for r in people["research_id"]]
        parts, missing = [], []
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for i, (rid, df) in enumerate(ex.map(load_one, jobs, chunksize=16), 1):
                if df is None:
                    missing.append(rid)
                    continue
                df.insert(0, "research_id", rid)
                parts.append(df)
                if i % 1000 == 0:
                    print(f"  loaded {i:,}/{len(jobs):,} ({time.time() - t0:.0f}s)",
                          file=sys.stderr)
        if missing:
            print(f"!! {len(missing):,} people without cdr3.out skipped", file=sys.stderr)
        clono = pd.concat(parts, ignore_index=True)
        clono["research_id"] = clono["research_id"].astype("category")
        clono["trbv"] = clono["trbv"].astype("category")
        clono.to_pickle(cache)
        print(f"  {len(clono):,} TRB clonotypes in {time.time() - t0:.0f}s -> {cache}",
              file=sys.stderr)

    # ---------------- per-person metrics ----------------
    g = clono.groupby("research_id", observed=True)
    per = pd.DataFrame({"n_clonotypes": g.size(), "total_reads": g["reads"].sum()})
    top10 = (clono.sort_values("reads", ascending=False)
                  .groupby("research_id", observed=True)["reads"].head(10)
                  .groupby(clono["research_id"], observed=True).sum())
    per["top10_frac"] = top10 / per["total_reads"]
    depth = int(max(50, np.floor(per["total_reads"].quantile(0.05))))
    t0 = time.time()
    per["rarefied_richness"] = pd.Series({rid: R.rarefied_richness(r.to_numpy(), depth)
                                          for rid, r in g["reads"]})
    print(f"Rarefaction at {depth} reads: {per['rarefied_richness'].notna().sum():,} people "
          f"eligible ({time.time() - t0:.0f}s)", file=sys.stderr)
    per = per.reset_index().merge(people, on="research_id", how="left")

    # ---------------- cohort-level aggregates ----------------
    anc_col = R.ancestry_colors()
    anc = [a for a in R.ANCESTRY_ORDER if a in set(per["ancestry"])]

    by_anc = per.groupby("ancestry").agg(
        n=("research_id", "size"), female_frac=("sex", lambda s: s.dropna().eq("Female").mean()),
        median_age=("age", "median"), median_clonotypes=("n_clonotypes", "median"),
        median_total_reads=("total_reads", "median")).reindex(anc)
    by_anc.to_csv(os.path.join(outdir, "cohort_by_ancestry.csv"))

    clono["len"] = clono["cdr3aa"].str.len()
    lens = (clono.merge(people[["research_id", "ancestry"]], on="research_id")
                 .groupby(["ancestry", "len"]).size().rename("n").reset_index())
    lens["frac"] = lens["n"] / lens.groupby("ancestry")["n"].transform("sum")
    lens.to_csv(os.path.join(outdir, "cdr3_length_by_ancestry.csv"), index=False)

    vcalled = clono[clono["trbv"] != ""]
    usage = (vcalled.groupby(["research_id", "trbv"], observed=True).size()
                    .unstack(fill_value=0))
    usage = usage.div(usage.sum(axis=1), axis=0)
    vsum = pd.DataFrame({"median": usage.median(), "q25": usage.quantile(.25),
                         "q75": usage.quantile(.75), "q05": usage.quantile(.05),
                         "q95": usage.quantile(.95),
                         "n_people_using": (usage > 0).sum()})
    vsum = vsum[vsum["n_people_using"] >= R.MIN_PEOPLE].sort_values("median", ascending=False)
    vsum.to_csv(os.path.join(outdir, "trbv_usage.csv"))

    share = (vcalled.groupby(["trbv", "cdr3aa"], observed=True)["research_id"]
                    .nunique().rename("n_people").reset_index())
    share.to_pickle(os.path.join(cache_dir, "trb_sharing.pkl"))
    spectrum = share["n_people"].value_counts().sort_index().rename("n_clonotypes")
    spectrum.index.name = "k_people"
    spectrum.to_csv(os.path.join(outdir, "sharing_spectrum.csv"))

    # ---------------- age models ----------------
    m = per.dropna(subset=["age", "sex", "ancestry", "rarefied_richness"]).copy()
    m = m[m["age"] >= 18]
    m["age_decade"] = m["age"] / 10
    m["log_reads"] = np.log(m["total_reads"])
    ref = {"sex": "Male", "ancestry": "EUR"}
    X = R.design(m, ["age_decade"], ["sex", "ancestry"], ref)
    Xs = R.design(m, ["age_decade", "log_reads"], ["sex", "ancestry"], ref)
    fits = {
        "log_rarefied_richness": R.ols(np.log(m["rarefied_richness"]), X),
        "log_rarefied_richness+log_reads": R.ols(np.log(m["rarefied_richness"]), Xs),
        "top10_frac_pct": R.ols(100 * m["top10_frac"], X),
        "top10_frac_pct+log_reads": R.ols(100 * m["top10_frac"], Xs),
    }
    pd.concat(fits, names=["model", "term"]).to_csv(os.path.join(outdir, "age_models.csv"))
    rich = fits["log_rarefied_richness"].loc["age_decade"]
    rich_s = fits["log_rarefied_richness+log_reads"].loc["age_decade"]
    top = fits["top10_frac_pct"].loc["age_decade"]
    pct = lambda b: 100 * (np.exp(b) - 1)  # noqa: E731

    rich_bins = R.quantile_summary(m, "age_bin", "rarefied_richness", order=R.AGE_LABELS)
    top_bins = R.quantile_summary(m.assign(top10_pct=100 * m["top10_frac"]), "age_bin",
                                  "top10_pct", order=R.AGE_LABELS)
    pd.concat({"rarefied_richness": rich_bins, "top10_frac_pct": top_bins}).to_csv(
        os.path.join(outdir, "age_bins.csv"))

    # ---------------- figure ----------------
    import cnsplots as cns
    mp = cns.multipanel(max_width=540)
    blue = R.nature()[3]
    W, H, GAP = 100, 100, 30

    ax = mp.panel("a", width=W, height=H, margin_right=30, margin_bottom=GAP)
    y = np.arange(len(anc))[::-1]
    ax.barh(y, by_anc["n"], color=[anc_col[a] for a in anc], height=0.7, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels(anc)
    for yi, a in zip(y, anc):
        ax.text(by_anc.loc[a, "n"] + by_anc["n"].max() * 0.03, yi,
                f"{int(by_anc.loc[a, 'n']):,}", va="center", fontsize=5.5)
    ax.set_xlim(0, by_anc["n"].max() * 1.25)
    R.style(ax, "Cohort", "People", None)

    ax = mp.panel("b", width=W, height=H, margin_right=15, margin_bottom=GAP)
    lx = np.log10(per["n_clonotypes"])
    hidden = R.safe_hist(ax, lx, bins=np.linspace(lx.quantile(.001), lx.quantile(.999), 40),
                         color=blue)
    ax.axvline(np.log10(per["n_clonotypes"].median()), color="0.2", lw=0.6, ls="--")
    ticks = [t for t in (30, 100, 300, 1000, 3000, 10000, 30000)
             if lx.min() - 0.2 <= np.log10(t) <= lx.max() + 0.2]
    ax.set_xticks(np.log10(ticks))
    ax.set_xticklabels([f"{t:,}" for t in ticks])
    R.style(ax, "TRB clonotypes per person", "Clonotypes", "People")
    ax.text(0.97, 0.97, f"median\n{per['n_clonotypes'].median():,.0f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5)

    ax = mp.panel("c", width=W, height=H, margin_right=40, margin_bottom=GAP)
    for a in anc:
        s_ = lens[(lens["ancestry"] == a) & lens["len"].between(8, 22)]
        ax.plot(s_["len"], 100 * s_["frac"], color=anc_col[a], lw=0.9, label=a)
    ax.legend(frameon=False)
    cns.take_legend_out(title="Ancestry", ax=ax)
    R.style(ax, "CDR3 length", "Length (aa)", "Clonotypes (%)")

    ax = mp.panel("d", width=440, height=85, margin_right=20, margin_bottom=GAP)
    vtop = vsum.head(TOP_V_GENES)
    R.pointrange(ax, 100 * vtop[["median", "q25", "q75", "q05", "q95"]],
                 xlabels=[v.replace("TRBV", "") for v in vtop.index], marker_size=10)
    R.style(ax, f"TRBV usage across {len(usage):,} people ({len(vtop)} most used genes)",
            "TRBV gene", "Clonotypes (%)")

    ax = mp.panel("e", width=W, height=H, margin_right=15)
    ax.loglog(spectrum.index, spectrum.to_numpy(), "o", ms=1.6, color=blue, mew=0)
    ax.minorticks_off()
    pub = spectrum[spectrum.index >= 2].sum() / spectrum.sum()
    ax.text(0.97, 0.97, f"{100 * pub:.1f}% of clonotypes\nshared by ≥2 people",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5)
    R.style(ax, "Sharing spectrum", "People sharing", "Clonotypes")

    for lab, summ, col, title, ylab, eff in (
            ("f", rich_bins, R.nature()[2], "Diversity vs age",
             f"Rarefied richness ({depth:,} reads)",
             f"{pct(rich['coef']):+.1f}% per decade "
             f"({pct(rich['lo']):+.1f} to {pct(rich['hi']):+.1f})"),
            ("g", top_bins, R.nature()[0], "Clonal expansion vs age",
             "Top-10 clonotype reads (%)",
             f"{top['coef']:+.2f} pp per decade ({top['lo']:+.2f} to {top['hi']:+.2f})")):
        ax = mp.panel(lab, width=W, height=H, margin_right=15)
        R.pointrange(ax, summ, colors=[col] * len(summ), xlabels=summ.index)
        ax.tick_params(axis="x", labelrotation=45)
        ax.set_title(title, pad=11)
        ax.text(0.5, 1.02, eff.replace("-", "−"), transform=ax.transAxes, ha="center", va="bottom",
                fontsize=5.5)
        R.style(ax, None, "Age (years)", ylab)
    R.save("fig_repertoire_atlas", outdir)

    # ---------------- summary ----------------
    rows = [("people", len(per)), ("trb_clonotypes_total", len(clono)),
            ("median_clonotypes_per_person", per["n_clonotypes"].median()),
            ("median_total_trb_reads", per["total_reads"].median()),
            ("rarefaction_depth_reads", depth),
            ("people_in_age_models", len(m)),
            ("frac_unique_clonotypes_shared_ge2", pub),
            ("richness_pct_per_decade", pct(rich["coef"])),
            ("richness_pct_per_decade_lo", pct(rich["lo"])),
            ("richness_pct_per_decade_hi", pct(rich["hi"])),
            ("richness_p", rich["p"]),
            ("richness_pct_per_decade_adj_log_reads", pct(rich_s["coef"])),
            ("top10_pp_per_decade", top["coef"]), ("top10_pp_per_decade_lo", top["lo"]),
            ("top10_pp_per_decade_hi", top["hi"]), ("top10_p", top["p"]),
            ("hist_people_in_suppressed_bins", hidden)]
    summary = pd.DataFrame(rows, columns=["metric", "value"])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
