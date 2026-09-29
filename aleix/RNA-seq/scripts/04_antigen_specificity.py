#!/usr/bin/env python3
"""04 -- Known antigen specificities: does the embedding see them, and can we find them in
the cohort? (report: reports/04_antigen_specificity/)

PART 1, BENCHMARK (public data only). VDJdb human TRB TCRs with a known epitope. Does SCEPTR
place TCRs that recognize the same epitope closer together than TCRs that don't?
  a  pairwise AUROC (same- vs different-epitope pairs, score = -distance)
  b  per-epitope AUROC, SCEPTR (TRBV+CDR3) vs a non-learned CDR3 3-mer baseline
  (also written: leave-one-out 5-nearest-neighbour epitope classification, balanced accuracy)
Compared: b_sceptr on TRBV+CDR3 (our production setting), SCEPTR cdr3_only, CDR3 3-mer
composition (cosine), TRBV identity alone. SCEPTR's pretraining is unsupervised (no epitope
labels, per its paper), but public sequences may overlap its training repertoires -- treat
absolute numbers as optimistic, the method ranking as the result.

PART 2, THE COHORT. Match every person's full TRB repertoire (02's cache) to VDJdb TCRs of
known specificity: exact (same TRBV gene + identical CDR3) and near (same TRBV gene, same
length, <= 1 amino-acid difference -- the usual "same specificity group" radius).
  c  % of people carrying >= 1 TCR matching each pathogen
  d  CMV-matched TCR carriage by age decade, EBV and influenza A alongside
  e  adjusted odds ratios for CMV carriage (age, sex, ancestry, log repertoire size)
  f  the CMV epitopes carried in the cohort, with their restricting HLA allele

Known biology tested in d/e: CMV seroprevalence rises steeply with age, and CMV drives large,
persistent CD8 expansions ("memory inflation"); EBV infects ~90-95% of adults early in life,
so its carriage should be far flatter with age (NHANES: Bate et al. 2010 Clin Infect Dis;
Emerson et al. 2017 Nat Genet for CMV-associated public TCRs). A matched TCR is evidence of
exposure only if it is not something everyone makes anyway, so a sensitivity model drops the
highest-generation-probability quartile of VDJdb TCRs (OLGA), which are public "by chance".
Carriage here is sensitivity-limited (bulk RNA-seq, whole blood) -- not a serostatus call;
the trends, not absolute prevalence, are the result.

Panel f is the bridge to the HLA phase: VDJdb's CMV TCRs are dominated by a few restricting
alleles (e.g. HLA-A*02:01 for NLVPMVATV). Whether a person can carry these TCRs depends on
their HLA type -- which is why any ancestry difference in (e) cannot be read as a difference
in CMV exposure until HLA is in the model.

Usage (from aleix/RNA-seq/, after 02; needs sceptr, olga):
  pixi run python3 -u scripts/04_antigen_specificity.py [--vdjdb path/to/vdjdb.slim.txt|.zip]
"""
import argparse
import io
import os
import sys
import time
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial.distance import cdist  # noqa: E402

import repfig as R  # noqa: E402

MIN_TCR_PER_EPITOPE = 30
MAX_TCR_PER_EPITOPE = 300
PATHOGEN_MIN_TCR = 50
CANON = r"C[ACDEFGHIKLMNPQRSTVWY]{3,}[FW]"


def load_vdjdb(path, url, cache_dir):
    if path is None:
        path = os.path.join(cache_dir, "vdjdb.zip")
        if not os.path.exists(path):
            import urllib.request
            print(f"Downloading VDJdb {url}", file=sys.stderr)
            urllib.request.urlretrieve(url, path)
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            name = next(n for n in zf.namelist() if n.endswith("vdjdb.slim.txt"))
            raw = pd.read_csv(io.BytesIO(zf.read(name)), sep="\t", dtype=str)
    else:
        raw = pd.read_csv(path, sep="\t", dtype=str)
    v = raw[(raw["species"] == "HomoSapiens") & (raw["gene"] == "TRB")].copy()
    v["score"] = pd.to_numeric(v["vdjdb.score"], errors="coerce").fillna(0)
    v = v[v["score"] >= 1]
    v["trbv"] = v["v.segm"].fillna("").str.split(",").str[0].str.split("*").str[0]
    v["cdr3aa"] = v["cdr3"].fillna("")
    v = v[(v["trbv"].str.startswith("TRBV")) & v["cdr3aa"].str.fullmatch(CANON)
          & (v["cdr3aa"].str.len() <= 30)]
    v["pathogen"] = v["antigen.species"]
    v["epitope"] = v["antigen.epitope"]
    v["hla"] = v["mhc.a"].fillna("").str.replace("HLA-", "", regex=False)
    v = v[["trbv", "cdr3aa", "epitope", "pathogen", "hla", "score"]].drop_duplicates(
        ["trbv", "cdr3aa", "epitope"])
    print(f"VDJdb human TRB, score >= 1: {len(v):,} TCR-epitope pairs, "
          f"{v['epitope'].nunique():,} epitopes", file=sys.stderr)
    return v.reset_index(drop=True)


# ------------------------------------------------------------------ Part 1: benchmark

def kmer_matrix(seqs, k=3):
    from scipy.sparse import csr_matrix
    vocab, rows, cols = {}, [], []
    for i, s in enumerate(seqs):
        for j in range(len(s) - k + 1):
            cols.append(vocab.setdefault(s[j:j + k], len(vocab)))
            rows.append(i)
    X = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(seqs), len(vocab)))
    norm = np.sqrt(X.multiply(X).sum(axis=1)).A1
    return X.multiply(1 / np.maximum(norm, 1e-9)[:, None]).tocsr()


def distances(method, bench):
    if method == "SCEPTR (TRBV + CDR3)":
        return cdist(bench["emb_v"], bench["emb_v"])
    if method == "SCEPTR (CDR3 only)":
        return cdist(bench["emb_c"], bench["emb_c"])
    if method == "CDR3 3-mer":
        X = kmer_matrix(bench["cdr3aa"])
        return 1 - (X @ X.T).toarray()
    if method == "TRBV gene only":
        v = bench["trbv"].to_numpy()
        return (v[:, None] != v[None, :]).astype(float)
    raise ValueError(method)


def benchmark(bench, methods, rng):
    from sklearn.metrics import roc_auc_score
    y = bench["epitope"].to_numpy()
    same = y[:, None] == y[None, :]
    iu = np.triu_indices(len(y), 1)
    take = rng.choice(len(iu[0]), size=min(2_000_000, len(iu[0])), replace=False)
    pi, pj = iu[0][take], iu[1][take]
    epis = sorted(set(y))
    res, per_epi = [], {}
    for m in methods:
        D = distances(m, bench)
        auc = roc_auc_score(same[pi, pj], -D[pi, pj] + rng.normal(0, 1e-9, len(pi)))
        Dk = D + np.diag(np.full(len(y), np.inf))
        nn = np.argsort(Dk + rng.uniform(0, 1e-9, Dk.shape), axis=1)[:, :5]
        pred = [pd.Series(y[r]).value_counts().idxmax() for r in nn]
        recall = [np.mean(np.asarray(pred)[y == e] == e) for e in epis]
        res.append((m, auc, float(np.mean(recall))))
        per = {}
        for e in epis:
            a = y == e
            s_in = D[np.ix_(a, a)][np.triu_indices(a.sum(), 1)]
            s_out = D[np.ix_(a, ~a)].ravel()
            s_out = rng.choice(s_out, size=min(len(s_out), 50_000), replace=False)
            lab = np.r_[np.ones(len(s_in)), np.zeros(len(s_out))]
            per[e] = roc_auc_score(lab, -np.r_[s_in, s_out] + rng.normal(0, 1e-9, len(lab)))
        per_epi[m] = per
        print(f"  {m}: pairwise AUROC {auc:.3f}, 5-NN balanced accuracy "
              f"{np.mean(recall):.3f}", file=sys.stderr)
    return (pd.DataFrame(res, columns=["method", "pairwise_auroc", "knn5_balanced_acc"]),
            pd.DataFrame(per_epi))


# ------------------------------------------------------------------ Part 2: cohort matching

def masked_keys(v, s, ids):
    """For each (V, CDR3) emit one key per position with that residue wildcarded, plus the
    exact key. Two same-V, same-length CDR3s within Hamming 1 share >= 1 key."""
    out = []
    df = pd.DataFrame({"id": ids, "v": v, "s": s})
    df["L"] = df["s"].str.len()
    for L, g in df.groupby("L"):
        for i in range(L):
            out.append(pd.DataFrame({"id": g["id"].to_numpy(),
                                     "key": (g["v"] + "|" + g["s"].str[:i] + "*"
                                             + g["s"].str[i + 1:]).to_numpy()}))
    return pd.concat(out, ignore_index=True)


def match_cohort(uniq, vdj, chunk=400_000):
    """uniq: unique cohort (trbv, cdr3aa) with 'uid'. vdj: VDJdb rows with 'vid'.
    Returns (uid, vid, exact) pairs."""
    exact = uniq.merge(vdj[["trbv", "cdr3aa", "vid"]], on=["trbv", "cdr3aa"])[["uid", "vid"]]
    exact["exact"] = True
    vk = masked_keys(vdj["trbv"], vdj["cdr3aa"], vdj["vid"]).rename(columns={"id": "vid"})
    vlen = set(zip(vdj["trbv"], vdj["cdr3aa"].str.len()))
    cand = uniq[[t in vlen for t in zip(uniq["trbv"], uniq["cdr3aa"].str.len())]]
    near = []
    for i in range(0, len(cand), chunk):
        c = cand.iloc[i:i + chunk]
        ck = masked_keys(c["trbv"], c["cdr3aa"], c["uid"]).rename(columns={"id": "uid"})
        near.append(ck.merge(vk, on="key")[["uid", "vid"]].drop_duplicates())
    near = pd.concat(near, ignore_index=True).drop_duplicates()
    near = near.merge(exact, on=["uid", "vid"], how="left")
    near["exact"] = near["exact"].fillna(False).astype(bool)
    return near


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vdjdb", default=None, help="vdjdb.slim.txt or release .zip; default: "
                    "download prep_vdjdb_pool.DEFAULT_URL into --cache-dir")
    ap.add_argument("--cohort", default="~/pipeline_outputs/rnaseq/cohort_full.tsv")
    ap.add_argument("--pheno-dir", default="~/pipeline_outputs/rnaseq/pheno")
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas")
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/04_antigen_specificity")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    cache_dir = os.path.expanduser(args.cache_dir)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    from prep_vdjdb_pool import DEFAULT_URL
    vdj = load_vdjdb(os.path.expanduser(args.vdjdb) if args.vdjdb else None, DEFAULT_URL,
                     cache_dir)
    vdj["vid"] = np.arange(len(vdj))

    # ---------------- Part 1: benchmark ----------------
    from embed_cdr3s import accepted_trbv, embed_sceptr, sceptr_model
    tcr_epis = vdj.groupby(["trbv", "cdr3aa"])["epitope"].transform("nunique")
    uni = vdj[tcr_epis == 1]
    counts = uni["epitope"].value_counts()
    keep_epi = counts[counts >= MIN_TCR_PER_EPITOPE].index
    bench = (uni[uni["epitope"].isin(keep_epi)].sample(frac=1, random_state=args.seed)
                .groupby("epitope").head(MAX_TCR_PER_EPITOPE).reset_index(drop=True))
    mv, mc = sceptr_model("v+cdr3"), sceptr_model("cdr3")
    bench = bench[bench["trbv"].isin(accepted_trbv(mv, bench["trbv"]))].reset_index(drop=True)
    counts = bench["epitope"].value_counts()
    bench = bench[bench["epitope"].isin(counts[counts >= MIN_TCR_PER_EPITOPE].index)]
    bench = bench.reset_index(drop=True)
    print(f"Benchmark: {len(bench):,} TCRs, {bench['epitope'].nunique()} epitopes",
          file=sys.stderr)
    bench_d = {"trbv": bench["trbv"], "cdr3aa": bench["cdr3aa"].tolist(),
               "emb_v": embed_sceptr(mv, pd.DataFrame({"TRBV": bench["trbv"],
                                                       "CDR3B": bench["cdr3aa"]}))[0],
               "emb_c": embed_sceptr(mc, pd.DataFrame({"CDR3B": bench["cdr3aa"]}))[0],
               "epitope": bench["epitope"]}
    methods = ["SCEPTR (TRBV + CDR3)", "SCEPTR (CDR3 only)", "CDR3 3-mer", "TRBV gene only"]
    bres, per_epi = benchmark(bench_d, methods, rng)
    epi_meta = bench.groupby("epitope").agg(pathogen=("pathogen", "first"),
                                            hla=("hla", lambda s: s.mode().iat[0]),
                                            n_tcr=("cdr3aa", "size"))
    per_epi = per_epi.join(epi_meta)
    bres.to_csv(os.path.join(outdir, "benchmark_methods.csv"), index=False)
    per_epi.to_csv(os.path.join(outdir, "benchmark_per_epitope.csv"))

    # ---------------- Part 2: cohort matching ----------------
    clono = pd.read_pickle(os.path.join(cache_dir, "trb_clonotypes_full.pkl"))
    clono = clono[clono["trbv"] != ""].copy()
    clono["trbv"] = clono["trbv"].astype(str)
    clono["research_id"] = clono["research_id"].astype(str)
    uniq = clono[["trbv", "cdr3aa"]].drop_duplicates().reset_index(drop=True)
    uniq["uid"] = np.arange(len(uniq))
    t0 = time.time()
    pairs = match_cohort(uniq, vdj)
    print(f"Matched {pairs['uid'].nunique():,} unique cohort clonotypes to VDJdb "
          f"({pairs['exact'].sum():,} exact pairs) in {time.time() - t0:.0f}s", file=sys.stderr)

    vdj["pgen"] = R.olga_pgen(vdj["cdr3aa"], args.workers)
    pgen_cut = vdj["pgen"].quantile(0.75)
    pairs = pairs.merge(vdj[["vid", "pathogen", "epitope", "hla", "pgen"]], on="vid")
    hits = (clono[["research_id", "trbv", "cdr3aa"]]
            .merge(uniq, on=["trbv", "cdr3aa"]).merge(pairs, on="uid"))

    people = R.load_people(args.cohort, args.pheno_dir)
    n_clono = clono.groupby("research_id").size().rename("n_clonotypes")
    people = people.merge(n_clono, left_on="research_id", right_index=True, how="inner")
    n_people = len(people)

    path_counts = vdj["pathogen"].value_counts()
    pathogens = [p for p in path_counts.index if path_counts[p] >= PATHOGEN_MIN_TCR]

    def carriers(h):
        return set(h["research_id"])

    rows = []
    for p in pathogens:
        h = hits[hits["pathogen"] == p]
        for level, hh in (("exact", h[h["exact"]]), ("near", h)):
            k = len(carriers(hh))
            rows.append((p, level, path_counts[p], hh["uid"].nunique(), k,
                         k if k >= R.MIN_PEOPLE else np.nan))
    carriage = pd.DataFrame(rows, columns=["pathogen", "match", "vdjdb_tcrs",
                                           "cohort_clonotypes_matched", "people_raw",
                                           "people"])
    carriage["pct_people"] = 100 * carriage["people"] / n_people
    carriage.drop(columns="people_raw").to_csv(os.path.join(outdir, "carriage_by_pathogen.csv"),
                                               index=False)

    epi_rows = []
    for (p, e), h in hits.groupby(["pathogen", "epitope"]):
        k = h["research_id"].nunique()
        if k >= R.MIN_PEOPLE:
            epi_rows.append((p, e, h["hla"].mode().iat[0], k, 100 * k / n_people))
    epi_car = (pd.DataFrame(epi_rows, columns=["pathogen", "epitope", "hla", "people",
                                               "pct_people"])
                 .sort_values("people", ascending=False))
    epi_car.to_csv(os.path.join(outdir, "carriage_by_epitope.csv"), index=False)

    # ---------------- epidemiology ----------------
    ref = {"sex": "Male", "ancestry": "EUR"}
    models, bins = {}, {}
    m = people.dropna(subset=["age", "sex", "ancestry"]).copy()
    m = m[m["age"] >= 18]
    m["age_decade"] = m["age"] / 10
    m["log_repertoire"] = np.log(m["n_clonotypes"])
    X = R.design(m, ["age_decade", "log_repertoire"], ["sex", "ancestry"], ref)
    for p in [q for q in ("CMV", "EBV", "InfluenzaA") if q in pathogens]:
        for variant, h in (("all", hits[hits["pathogen"] == p]),
                           ("low-Pgen", hits[(hits["pathogen"] == p)
                                             & (hits["pgen"] < pgen_cut)])):
            y = m["research_id"].isin(carriers(h)).astype(float)
            models[(p, variant)] = R.logit(y, X)
            if variant == "all":
                b = m.assign(y=y).groupby("age_bin", observed=True)["y"].agg(["sum", "size"])
                b = b[b["size"] >= R.MIN_PEOPLE].reindex(
                    [a for a in R.AGE_LABELS if a in b.index])
                pr, lo, hi = R.wilson(b["sum"], b["size"])
                bins[p] = pd.DataFrame({"n": b["size"], "pct": 100 * pr, "lo": 100 * lo,
                                        "hi": 100 * hi}, index=b.index)
    pd.concat(models, names=["pathogen", "variant", "term"]).to_csv(
        os.path.join(outdir, "carriage_models.csv"))
    pd.concat(bins, names=["pathogen", "age_bin"]).to_csv(
        os.path.join(outdir, "carriage_by_age.csv"))

    cmv_hla = (hits[hits["pathogen"] == "CMV"][["research_id", "uid", "hla"]]
               .drop_duplicates())
    hla_share = (cmv_hla.groupby("hla")["research_id"].nunique().sort_values(ascending=False))
    hla_share = hla_share[hla_share >= R.MIN_PEOPLE]
    hla_share.rename("people").to_csv(os.path.join(outdir, "cmv_matches_by_hla.csv"))

    # ---------------- figure ----------------
    import cnsplots as cns
    from matplotlib.ticker import FixedLocator, NullLocator
    mp = cns.multipanel(max_width=540)
    pal = R.nature()
    W, H, GAP = 100, 100, 34
    mcol = dict(zip(methods, [pal[3], pal[1], pal[5], "0.6"]))
    short = {"SCEPTR (TRBV + CDR3)": "SCEPTR V+CDR3", "SCEPTR (CDR3 only)": "SCEPTR CDR3",
             "CDR3 3-mer": "CDR3 3-mer", "TRBV gene only": "TRBV only"}

    def subtitle(ax, text):
        ax.text(0.5, 1.02, text, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=5.5)

    ax = mp.panel("a", width=W, height=H, margin_right=15, margin_bottom=GAP)
    y = np.arange(len(bres))[::-1]
    ax.barh(y, bres["pairwise_auroc"] - 0.5, left=0.5, height=0.65, linewidth=0,
            color=[mcol[k] for k in bres["method"]])
    ax.axvline(0.5, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([short[k] for k in bres["method"]])
    ax.set_xlim(0.4, 1.0)
    ax.set_title("Same-epitope TCR pairs", pad=11)
    subtitle(ax, f"{len(bench):,} VDJdb TCRs, {bench['epitope'].nunique()} epitopes")
    R.style(ax, None, "Pairwise AUROC", None)

    ax = mp.panel("b", width=W, height=H, margin_right=65, margin_bottom=GAP)
    top_path = per_epi["pathogen"].value_counts().head(5).index.tolist()
    pcol = {"CMV": pal[0], "EBV": pal[2], "InfluenzaA": pal[1], "SARS-CoV-2": pal[4],
            "HIV-1": pal[6], "HomoSapiens": pal[3]}
    for p in top_path + ["other"]:
        s_ = (per_epi[per_epi["pathogen"] == p] if p != "other"
              else per_epi[~per_epi["pathogen"].isin(top_path)])
        if len(s_):
            ax.scatter(s_["CDR3 3-mer"], s_["SCEPTR (TRBV + CDR3)"], s=12, linewidths=0,
                       color=pcol.get(p, pal[8]),
                       label=p.replace("HomoSapiens", "Self (human)"), zorder=3)
    lim = [min(0.45, per_epi[["CDR3 3-mer", "SCEPTR (TRBV + CDR3)"]].min().min()), 1.0]
    ax.plot(lim, lim, color="0.6", lw=0.5, ls="--", zorder=0)
    ax.set(xlim=lim, ylim=lim)
    ax.legend(frameon=False)
    cns.take_legend_out(title="Pathogen", ax=ax)
    better = (per_epi["SCEPTR (TRBV + CDR3)"] > per_epi["CDR3 3-mer"]).mean()
    ax.set_title("Per epitope", pad=11)
    subtitle(ax, f"SCEPTR ahead for {100 * better:.0f}% of epitopes")
    R.style(ax, None, "AUROC, CDR3 3-mer", "AUROC, SCEPTR V+CDR3")

    ax = mp.panel("c", width=W, height=H, margin_right=15, margin_bottom=GAP)
    cp = carriage.pivot(index="pathogen", columns="match", values="pct_people")
    cp = cp.reindex([p for p in pathogens if p in cp.index]).head(6)
    y = np.arange(len(cp))[::-1]
    ax.barh(y + 0.18, cp["near"].fillna(0), height=0.36, color=pal[3], linewidth=0,
            label="≤1 aa, same V")
    ax.barh(y - 0.18, cp["exact"].fillna(0), height=0.36, color=pal[1], linewidth=0,
            label="Exact")
    ax.set_yticks(y)
    ax.set_yticklabels([p.replace("HomoSapiens", "Self (human)") for p in cp.index])
    ax.legend(frameon=False, loc="lower right", fontsize=5.5, handlelength=0.8)
    ax.set_title("Cohort carriage", pad=11)
    subtitle(ax, "people with ≥1 matching TCR")
    R.style(ax, None, "People (%)", None)

    ax = mp.panel("d", width=W, height=H, margin_right=15)
    for p, colr in (("CMV", pal[0]), ("EBV", pal[2]), ("InfluenzaA", pal[1])):
        if p not in bins:
            continue
        b = bins[p]
        xx = np.array([R.AGE_LABELS.index(a) for a in b.index])
        ax.errorbar(xx, b["pct"], yerr=[b["pct"] - b["lo"], b["hi"] - b["pct"]], color=colr,
                    lw=0.9, ms=3, marker="o", capsize=0, label=p)
    ax.set_xticks(np.arange(len(R.AGE_LABELS)))
    ax.set_xticklabels(R.AGE_LABELS, rotation=45)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left")
    R.style(ax, "Carriage by age", "Age (years)", "People with ≥1 match (%)")

    ax = mp.panel("e", width=W, height=H, margin_right=15)
    t = models.get(("CMV", "all"))
    if t is not None:
        terms = {"age_decade": "Age (+10 y)", "sex[Female]": "Female vs male"}
        terms.update({f"ancestry[{a}]": f"{a} vs EUR" for a in R.ANCESTRY_ORDER
                      if f"ancestry[{a}]" in t.index})
        R.forest(ax, t, terms, color=pal[0], xlabel="Odds ratio (95% CI)")
        ts = models[("CMV", "low-Pgen")].loc[list(terms)]
        yy = np.arange(len(terms))[::-1] - 0.28
        ax.hlines(yy, np.exp(ts["lo"]), np.exp(ts["hi"]), color=pal[0], lw=0.6, alpha=0.6)
        ax.scatter(np.exp(ts["coef"]), yy, s=10, marker="D", facecolors="white",
                   edgecolors=pal[0], linewidths=0.8, zorder=3)
        ax.set_xscale("log")
        lo_, hi_ = ax.get_xlim()
        ticks = [v for v in (0.125, 0.25, 0.5, 1, 2, 4, 8) if lo_ <= v <= hi_]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticklabels([f"{v:g}" for v in ticks])
    ax.set_title("CMV carriage, adjusted", pad=11)
    subtitle(ax, "● all VDJdb TCRs   ◇ low-P$_{gen}$ TCRs only")

    ax = mp.panel("f", width=W, height=H, margin_right=15)
    cmv_epi = epi_car[epi_car["pathogen"] == "CMV"].head(6)
    y = np.arange(len(cmv_epi))[::-1]
    ax.barh(y, cmv_epi["pct_people"], color=pal[0], height=0.65, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{e}  {h}" for e, h in zip(cmv_epi["epitope"], cmv_epi["hla"])],
                       fontsize=5.5)
    ax.set_title("CMV epitopes carried", pad=11)
    subtitle(ax, "with restricting HLA allele")
    R.style(ax, None, "People (%)", None)
    R.save("fig_antigen_specificity", outdir)

    cmv = models.get(("CMV", "all"))
    ebv = models.get(("EBV", "all"))
    rows = [("vdjdb_tcr_epitope_pairs", len(vdj)), ("benchmark_tcrs", len(bench)),
            ("benchmark_epitopes", bench["epitope"].nunique()),
            ("people", n_people), ("pgen_q75_cutoff", pgen_cut)]
    rows += [(f"pairwise_auroc_{r.method}", r.pairwise_auroc) for r in bres.itertuples()]
    rows += [(f"knn5_bal_acc_{r.method}", r.knn5_balanced_acc) for r in bres.itertuples()]
    if cmv is not None:
        rows += [("cmv_or_per_decade", np.exp(cmv.loc["age_decade", "coef"])),
                 ("cmv_or_per_decade_lo", np.exp(cmv.loc["age_decade", "lo"])),
                 ("cmv_or_per_decade_hi", np.exp(cmv.loc["age_decade", "hi"])),
                 ("cmv_age_p", cmv.loc["age_decade", "p"])]
    if ebv is not None:
        rows += [("ebv_or_per_decade", np.exp(ebv.loc["age_decade", "coef"])),
                 ("ebv_or_per_decade_lo", np.exp(ebv.loc["age_decade", "lo"])),
                 ("ebv_or_per_decade_hi", np.exp(ebv.loc["age_decade", "hi"]))]
    summary = pd.DataFrame(rows, columns=["metric", "value"])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
