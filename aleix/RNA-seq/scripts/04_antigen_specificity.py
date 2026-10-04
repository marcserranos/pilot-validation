#!/usr/bin/env python3
"""04 -- Known antigen specificities: does the embedding see them, and can we detect them
in the cohort? (report: reports/04_antigen_specificity/)

PART 1, BENCHMARK (public data only). VDJdb human TRB TCRs with a known epitope. Does SCEPTR
place TCRs recognising the same epitope closer together than TCRs that do not?
  a  pairwise AUROC (same- vs different-epitope pairs, score = -distance)
  b  per-epitope AUROC, SCEPTR (TRBV+CDR3) vs a non-learned CDR3 3-mer baseline
  (also written: 5-nearest-neighbour epitope classification, balanced accuracy vs chance)
Compared: b_sceptr on TRBV+CDR3 (our production setting), SCEPTR cdr3_only, CDR3 3-mer
composition (cosine), TRBV identity alone. SCEPTR's pretraining is unsupervised (no epitope
labels), but public sequences may overlap its training repertoires -- treat absolute numbers
as optimistic and the ranking of methods as the result.

PART 2, THE COHORT -- and the null that the naive version needs. Matching each person's
repertoire (02's cache) to VDJdb TCRs looks easy and is badly misleading: with ~1,100
clonotypes per person and 8,000 reference TCRs, most "matches" are recombination
coincidences, not evidence of exposure.
  c  the artifact, made explicit: observed carriage vs published US adult seroprevalence.
     A working assay would track the diagonal. Chance matching does not.
  d  the null: for every VDJdb TCR, decoy TCRs drawn from OLGA's generative model and
     matched on TRBV gene, CDR3 length and generation probability (repfig.matched_decoys),
     run through the identical matching pipeline. Shown as enrichment
     (observed / decoy-expected, 95% CI) per pathogen, with high-prevalence
     pathogens (EBV, CMV, influenza A) separated from internal negative controls (HIV-1,
     HCV, HBV, dengue, HTLV-1, and human self-peptides -- no exposure at all).
  e  excess matches (observed - expected) vs age for CMV and EBV, the biological test:
     CMV seroprevalence climbs steeply with age while EBV is near-universal from childhood
     (Bate et al. 2010 Clin Infect Dis; Dowd et al. 2013 PLoS One), so CMV excess should
     rise and EBV excess should not.
  f  CMV epitopes carried, with the HLA allele that restricts each -- the bridge to the HLA
     phase: whether a person can carry these TCRs at all depends on their HLA type, so no
     ancestry difference in carriage is interpretable until HLA is in the model.

Exact matching only (same TRBV gene, identical CDR3 amino acid sequence). A <= 1 aa variant
is also computed and is reported as a failure mode, not a result: it matches the majority of
people for every pathogen including the negative controls.

Usage (from aleix/RNA-seq/, after 02; needs sceptr, olga):
  pixi run python3 -u scripts/04_antigen_specificity.py [--decoy-reps 3] [--decoy-pool 400000]
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

# Published US adult seroprevalence, for reference only -- literature values, nothing here
# measures serostatus. EBV: Dowd 2013 PLoS One (~95% adults). CMV: Bate 2010 Clin Infect Dis
# (NHANES, ~50% age-adjusted). Influenza A: near-universal lifetime exposure. MCPyV: ~65%.
# HBV ever-infected ~4.3%, HCV ~1.0%, HIV-1 ~0.4%, HTLV-1 <0.1%, dengue <1% (continental US).
SEROPREVALENCE = {"EBV": 95.0, "InfluenzaA": 98.0, "CMV": 50.0, "MCPyV": 65.0, "HBV": 4.3,
                  "HCV": 1.0, "HIV-1": 0.4, "DENV": 0.5, "HTLV-1": 0.05}
POSITIVES = ["EBV", "InfluenzaA", "CMV"]
NEGATIVE_CONTROLS = ["HIV-1", "HCV", "HBV", "DENV", "HTLV-1", "HomoSapiens"]

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
    boot = R.epitope_bootstrap(y, seed=int(rng.integers(1 << 31)))
    res, per_epi, reps = [], {}, {}
    for m in methods:
        D = distances(m, bench)
        auc = roc_auc_score(same[pi, pj], -D[pi, pj] + rng.normal(0, 1e-9, len(pi)))
        Dk = D + np.diag(np.full(len(y), np.inf))
        nn = np.argsort(Dk + rng.uniform(0, 1e-9, Dk.shape), axis=1)[:, :5]
        pred = [pd.Series(y[r]).value_counts().idxmax() for r in nn]
        recall = [np.mean(np.asarray(pred)[y == e] == e) for e in epis]
        reps[m] = R.boot_scores(D, dict(zip(epis, recall)), boot)
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
    out = pd.DataFrame(res, columns=["method", "pairwise_auroc", "knn5_balanced_acc"])
    # 95% CIs from an epitope-level bootstrap (epitopes are the independent units), and
    # paired differences against the first (production) method on the same replicates.
    ref = methods[0]
    for col, k in (("pairwise_auroc", 0), ("knn5_balanced_acc", 1)):
        out[f"{col}_lo"], out[f"{col}_hi"] = zip(*[R.ci(reps[m][k]) for m in out["method"]])
        out[f"{col}_diff_vs_ref_lo"], out[f"{col}_diff_vs_ref_hi"] = zip(
            *[R.ci(reps[m][k] - reps[ref][k]) for m in out["method"]])
    out["reference_method"] = ref
    out["bootstrap_replicates"] = len(boot["plan"])
    return out, pd.DataFrame(per_epi)


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
    ap.add_argument("--decoy-reps", type=int, default=3,
                    help="matched decoy TCR sets per real VDJdb TCR (the null)")
    ap.add_argument("--decoy-pool", type=int, default=250_000,
                    help="synthetic TCRs drawn from OLGA to sample decoys from")
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
    cohort_keys = clono[["research_id", "trbv", "cdr3aa"]]

    people = R.load_people(args.cohort, args.pheno_dir)
    n_clono = clono.groupby("research_id").size().rename("n_clonotypes")
    people = people.merge(n_clono, left_on="research_id", right_index=True, how="inner")
    n_people = len(people)

    t0 = time.time()
    vdj["pgen"] = R.olga_pgen(vdj["cdr3aa"], args.workers)
    print(f"Pgen for {len(vdj):,} VDJdb TCRs: {time.time() - t0:.0f}s", file=sys.stderr)

    # Decoy null: same TRBV gene, same CDR3 length, same generation probability, but drawn
    # from OLGA's recombination model instead of from an antigen-selected repertoire.
    t0 = time.time()
    pool_dec = R.olga_generate(args.decoy_pool, seed=args.seed)
    pool_dec = pool_dec.drop_duplicates(["trbv", "cdr3aa"])
    pool_dec = pool_dec[~pool_dec.set_index(["trbv", "cdr3aa"]).index.isin(
        vdj.set_index(["trbv", "cdr3aa"]).index)]
    pool_dec["pgen"] = R.olga_pgen(pool_dec["cdr3aa"], args.workers)
    print(f"Decoy pool {len(pool_dec):,} synthetic TCRs + Pgen: {time.time() - t0:.0f}s",
          file=sys.stderr)
    decoys = R.matched_decoys(vdj[["trbv", "cdr3aa", "pgen", "vid"]], pool_dec,
                              args.decoy_reps, seed=args.seed)
    # How often the Pgen match had to be loosened, per pathogen: a decoy matched only on
    # V and length (tol = inf) controls recombination less well than one within 0.25.
    tol_by_tcr = decoys.drop_duplicates("vid_real").set_index("vid_real")["tol"]
    tv = vdj[["vid", "pathogen"]].copy()
    tv["tol"] = tv["vid"].map(tol_by_tcr)
    tv.loc[tv["vid"].isin(decoys.attrs.get("unmatched_vids", [])), "tol"] = -1.0
    tv["tier"] = tv["tol"].map({0.25: "<=0.25", 0.5: "<=0.5", 1.0: "<=1.0",
                                np.inf: "V+length only", -1.0: "no decoy"}).fillna("Pgen = 0")
    tiers = ["<=0.25", "<=0.5", "<=1.0", "V+length only", "no decoy", "Pgen = 0"]
    tol_tab = (pd.crosstab(tv["pathogen"], tv["tier"]).reindex(columns=tiers, fill_value=0))
    tol_tab.loc["ALL"] = tol_tab.sum()
    tol_tab = tol_tab.div(tol_tab.sum(axis=1), axis=0).mul(100).round(2)
    tol_tab.insert(0, "vdjdb_tcrs", pd.crosstab(tv["pathogen"], tv["tier"]).sum(axis=1)
                   .reindex(tol_tab.index).fillna(len(tv)).astype(int))
    tol_tab.to_csv(os.path.join(outdir, "decoy_tolerance.csv"))
    print("Decoy Pgen tolerance actually used (% of VDJdb TCRs):\n"
          + tol_tab.loc[["ALL"]].to_string(), file=sys.stderr)
    decoys = decoys.merge(vdj[["vid", "pathogen"]].rename(columns={"vid": "vid_real"}),
                          on="vid_real")

    def exact_hits(ref):
        """People x pathogen exact matches for a reference TCR set."""
        return cohort_keys.merge(ref[["trbv", "cdr3aa", "pathogen"]].drop_duplicates(),
                                 on=["trbv", "cdr3aa"])

    obs = exact_hits(vdj)
    obs_people = obs.groupby("pathogen")["research_id"].nunique()
    obs_count = obs.groupby(["research_id", "pathogen"]).size().rename("n").reset_index()

    exp_people, exp_count = [], []
    for rep, d in decoys.groupby("rep"):
        h = exact_hits(d)
        exp_people.append(h.groupby("pathogen")["research_id"].nunique())
        exp_count.append(h.groupby(["research_id", "pathogen"]).size().rename("n").reset_index())
    exp_people = pd.concat(exp_people, axis=1).reindex(obs_people.index).fillna(0)
    exp_mean = exp_people.mean(axis=1)
    exp_sd = exp_people.std(axis=1)

    # near (<=1 aa) matching, kept only to document that it fails
    near_pairs = (match_cohort(uniq, vdj)
                  .merge(vdj[["vid", "pathogen"]], on="vid")[["uid", "pathogen"]]
                  .drop_duplicates())
    near = (cohort_keys.merge(uniq, on=["trbv", "cdr3aa"])
                       .merge(near_pairs, on="uid"))
    near_people = near.groupby("pathogen")["research_id"].nunique()

    path_counts = vdj["pathogen"].value_counts()
    pathogens = [p for p in path_counts.index if path_counts[p] >= PATHOGEN_MIN_TCR]
    car = pd.DataFrame({"vdjdb_tcrs": path_counts.reindex(pathogens),
                        "observed_people": obs_people.reindex(pathogens).fillna(0),
                        "expected_people": exp_mean.reindex(pathogens).fillna(0),
                        "expected_sd": exp_sd.reindex(pathogens).fillna(0),
                        "near_people": near_people.reindex(pathogens).fillna(0)})
    car.index.name = "pathogen"
    car = car[car["observed_people"] >= R.MIN_PEOPLE]
    for c in ("observed", "expected", "near"):
        car[f"pct_{c}"] = 100 * car[f"{c}_people"] / n_people
    k = car["observed_people"].to_numpy()
    e = np.maximum(car["expected_people"].to_numpy(), 0.5)
    car["enrichment"] = k / e
    # Poisson-style CI on the observed count, holding the expectation fixed
    car["enr_lo"] = np.maximum(k - 1.96 * np.sqrt(k), 0.5) / e
    car["enr_hi"] = (k + 1.96 * np.sqrt(k)) / e
    car["seroprevalence_pct"] = [SEROPREVALENCE.get(p, np.nan) for p in car.index]
    car["role"] = ["positive" if p in POSITIVES else
                   ("negative control" if p in NEGATIVE_CONTROLS else "other")
                   for p in car.index]
    car.to_csv(os.path.join(outdir, "carriage_by_pathogen.csv"))

    epi_rows = []
    obs_epi = cohort_keys.merge(vdj[["trbv", "cdr3aa", "pathogen", "epitope", "hla"]]
                                .drop_duplicates(["trbv", "cdr3aa"]), on=["trbv", "cdr3aa"])
    for (p, ep), h in obs_epi.groupby(["pathogen", "epitope"]):
        n = h["research_id"].nunique()
        if n >= R.MIN_PEOPLE:
            epi_rows.append((p, ep, h["hla"].mode().iat[0], n, 100 * n / n_people))
    epi_car = (pd.DataFrame(epi_rows, columns=["pathogen", "epitope", "hla", "people",
                                               "pct_people"])
                 .sort_values("people", ascending=False))
    epi_car.to_csv(os.path.join(outdir, "carriage_by_epitope.csv"), index=False)

    # ---------------- excess matches vs age ----------------
    m = people.dropna(subset=["age", "sex", "ancestry"]).copy()
    m = m[m["age"] >= 18]
    m["age_decade"] = m["age"] / 10
    m["log_repertoire"] = np.log(m["n_clonotypes"])
    ref = {"sex": "Male", "ancestry": "EUR"}
    X = R.design(m, ["age_decade", "log_repertoire"], ["sex", "ancestry"], ref)
    models, bins = {}, {}
    for p in [q for q in ("CMV", "EBV", "InfluenzaA") if q in car.index]:
        o = (obs_count[obs_count["pathogen"] == p].set_index("research_id")["n"]
             .reindex(m["research_id"]).fillna(0).to_numpy())
        ex = np.zeros(len(m))
        for ec in exp_count:
            ex += (ec[ec["pathogen"] == p].set_index("research_id")["n"]
                   .reindex(m["research_id"]).fillna(0).to_numpy())
        ex /= len(exp_count)
        excess = o - ex
        models[p] = R.ols(excess, X)
        d = m.assign(excess=excess).groupby("age_bin", observed=True)["excess"]
        agg = pd.DataFrame({"n": d.size(), "mean": d.mean(),
                            "se": d.std() / np.sqrt(d.size())})
        bins[p] = agg[agg["n"] >= R.MIN_PEOPLE].reindex(
            [a for a in R.AGE_LABELS if a in agg.index])
    pd.concat(models, names=["pathogen", "term"]).to_csv(
        os.path.join(outdir, "excess_vs_age_models.csv"))
    pd.concat(bins, names=["pathogen", "age_bin"]).to_csv(
        os.path.join(outdir, "excess_by_age.csv"))

    # ---------------- figure ----------------
    import cnsplots as cns
    from matplotlib.ticker import FixedLocator, NullLocator
    mp = cns.multipanel(max_width=540)
    pal = R.nature()
    W, H, GAP = 105, 100, 34
    mcol = dict(zip(methods, [pal[3], pal[1], pal[5], "0.6"]))
    short = {"SCEPTR (TRBV + CDR3)": "SCEPTR V+CDR3", "SCEPTR (CDR3 only)": "SCEPTR CDR3",
             "CDR3 3-mer": "CDR3 3-mer", "TRBV gene only": "TRBV only"}
    rcol = {"positive": pal[0], "negative control": pal[3], "other": "0.7"}

    def sub(ax, t):
        ax.text(0.5, 1.02, t, transform=ax.transAxes, ha="center", va="bottom", fontsize=5.5)

    def logticks(ax, axis="x"):
        lo, hi = (ax.get_xlim() if axis == "x" else ax.get_ylim())
        t = [v for v in (0.01, 0.1, 0.25, 0.5, 1, 2, 4, 10, 25, 100) if lo <= v <= hi]
        ax_ = ax.xaxis if axis == "x" else ax.yaxis
        ax_.set_major_locator(FixedLocator(t))
        ax_.set_minor_locator(NullLocator())
        (ax.set_xticklabels if axis == "x" else ax.set_yticklabels)([f"{v:g}" for v in t])

    ax = mp.panel("a", width=W, height=H, margin_right=12, margin_bottom=GAP)
    y = np.arange(len(bres))[::-1]
    ax.barh(y, bres["pairwise_auroc"] - 0.5, left=0.5, height=0.65, linewidth=0,
            color=[mcol[k] for k in bres["method"]])
    ax.errorbar(bres["pairwise_auroc"], y,
                xerr=[np.maximum(bres["pairwise_auroc"] - bres["pairwise_auroc_lo"], 0),
                      np.maximum(bres["pairwise_auroc_hi"] - bres["pairwise_auroc"], 0)],
                fmt="none", ecolor="0.2", elinewidth=0.6, capsize=1.5)
    ax.axvline(0.5, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([short[k] for k in bres["method"]])
    ax.set_xlim(0.4, 0.75)
    ax.set_title("Same-epitope TCR pairs", pad=11)
    sub(ax, f"{len(bench):,} VDJdb TCRs, {bench['epitope'].nunique()} epitopes")
    R.style(ax, None, "Pairwise AUROC", None)

    ax = mp.panel("b", width=W, height=H, margin_right=58, margin_bottom=GAP)
    top_path = per_epi["pathogen"].value_counts().head(5).index.tolist()
    pcol = {"CMV": pal[0], "EBV": pal[2], "InfluenzaA": pal[1], "SARS-CoV-2": pal[4],
            "HIV-1": pal[6], "HomoSapiens": pal[3], "HCV": pal[8]}
    for p in top_path + ["other"]:
        s_ = (per_epi[per_epi["pathogen"] == p] if p != "other"
              else per_epi[~per_epi["pathogen"].isin(top_path)])
        if len(s_):
            ax.scatter(s_["CDR3 3-mer"], s_["SCEPTR (TRBV + CDR3)"], s=11, linewidths=0,
                       color=pcol.get(p, "0.7"),
                       label=p.replace("HomoSapiens", "Self"), zorder=3)
    lim = [min(0.45, per_epi[["CDR3 3-mer", "SCEPTR (TRBV + CDR3)"]].min().min()), 1.0]
    ax.plot(lim, lim, color="0.6", lw=0.5, ls="--", zorder=0)
    ax.set(xlim=lim, ylim=lim)
    ax.legend(frameon=False)
    cns.take_legend_out(title="Pathogen", ax=ax)
    better = (per_epi["SCEPTR (TRBV + CDR3)"] > per_epi["CDR3 3-mer"]).mean()
    ax.set_title("Per epitope", pad=11)
    sub(ax, f"SCEPTR ahead for {100 * better:.0f}% of epitopes")
    R.style(ax, None, "AUROC, CDR3 3-mer", "AUROC, SCEPTR V+CDR3")

    ax = mp.panel("c", width=W, height=H, margin_right=12, margin_bottom=GAP)
    cc = car.dropna(subset=["seroprevalence_pct"])
    for role, g_ in cc.groupby("role"):
        ax.scatter(g_["seroprevalence_pct"], g_["pct_near"], s=13, marker="^",
                   linewidths=0, color=rcol[role], alpha=0.85)
        ax.scatter(g_["seroprevalence_pct"], g_["pct_observed"], s=13, linewidths=0,
                   color=rcol[role], label=role)
    ax.plot([0.03, 100], [0.03, 100], color="0.6", lw=0.5, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    logticks(ax, "x")
    logticks(ax, "y")
    ax.legend(frameon=False, fontsize=5, loc="upper left")
    ax.set_title("Carriage vs seroprevalence", pad=11)
    sub(ax, "● exact   ▲ ≤1 aa   (dashed = agreement)")
    R.style(ax, None, "Published seroprevalence (%)", "People with ≥1 match (%)")

    ax = mp.panel("d", width=W, height=H, margin_right=12)
    ce = car.sort_values("enrichment", ascending=False).iloc[::-1]
    y = np.arange(len(ce))
    for i, r in enumerate(ce.itertuples()):
        c = rcol[r.role]
        ax.plot([r.enr_lo, r.enr_hi], [i, i], color=c, lw=1.0)
        ax.scatter([r.enrichment], [i], s=15, color=c, zorder=3)
    ax.axvline(1, color="0.6", lw=0.5, ls="--")
    ax.set_xscale("log")
    logticks(ax, "x")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{p.replace('HomoSapiens', 'Self')} "
                        f"({r.pct_observed:.0f}% vs {r.pct_expected:.0f}%)"
                        for p, r in zip(ce.index, ce.itertuples())], fontsize=5)
    ax.set_title("Enrichment over null", pad=11)
    sub(ax, "observed vs decoy-expected carriage")
    R.style(ax, None, "Enrichment (95% CI)", None)

    ax = mp.panel("e", width=W, height=H, margin_right=12)
    for p, colr in (("CMV", pal[0]), ("EBV", pal[2])):
        if p not in bins:
            continue
        b_ = bins[p]
        xx = np.array([R.AGE_LABELS.index(a) for a in b_.index])
        ax.errorbar(xx, b_["mean"], yerr=1.96 * b_["se"], color=colr, lw=0.9, ms=3,
                    marker="o", capsize=0, label=p)
    ax.axhline(0, color="0.6", lw=0.5, ls="--")
    ax.set_xticks(np.arange(len(R.AGE_LABELS)))
    ax.set_xticklabels(R.AGE_LABELS, rotation=45)
    ax.legend(frameon=False, fontsize=5.5)
    if "CMV" in models:
        c_ = models["CMV"].loc["age_decade"]
        sub(ax, f"CMV {c_['coef']:+.2f} per decade (P = {c_['p']:.2g})".replace("-", "−"))
    ax.set_title("Excess matches vs age", pad=11)
    R.style(ax, None, "Age (years)", "Observed − expected matches")

    ax = mp.panel("f", width=W, height=H, margin_right=12)
    cmv_epi = epi_car[epi_car["pathogen"] == "CMV"].head(6)
    y = np.arange(len(cmv_epi))[::-1]
    ax.barh(y, cmv_epi["pct_people"], color=pal[0], height=0.62, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{e_}  {h}" for e_, h in zip(cmv_epi["epitope"], cmv_epi["hla"])],
                       fontsize=5)
    ax.set_title("CMV epitopes carried", pad=11)
    sub(ax, "with restricting HLA allele")
    R.style(ax, None, "People with ≥1 exact match (%)", None)
    R.save("fig_antigen_specificity", outdir)

    rows = [("vdjdb_tcr_epitope_pairs", len(vdj)), ("benchmark_tcrs", len(bench)),
            ("benchmark_epitopes", bench["epitope"].nunique()), ("people", n_people),
            ("decoy_reps", args.decoy_reps), ("decoy_pool", len(pool_dec))]
    rows += [(f"pairwise_auroc_{r.method}", r.pairwise_auroc) for r in bres.itertuples()]
    rows += [(f"pairwise_auroc_{r.method}_lo", r.pairwise_auroc_lo) for r in bres.itertuples()]
    rows += [(f"pairwise_auroc_{r.method}_hi", r.pairwise_auroc_hi) for r in bres.itertuples()]
    rows += [(f"knn5_bal_acc_{r.method}", r.knn5_balanced_acc) for r in bres.itertuples()]
    rows += [(f"decoy_tol_pct_{t}", tol_tab.loc["ALL", t]) for t in tiers]
    for p in car.index:
        rows += [(f"pct_observed_{p}", car.loc[p, "pct_observed"]),
                 (f"pct_expected_{p}", car.loc[p, "pct_expected"]),
                 (f"pct_near_{p}", car.loc[p, "pct_near"]),
                 (f"enrichment_{p}", car.loc[p, "enrichment"])]
    for p, t in models.items():
        rows += [(f"excess_per_decade_{p}", t.loc["age_decade", "coef"]),
                 (f"excess_per_decade_p_{p}", t.loc["age_decade", "p"])]
    summary = pd.DataFrame(rows, columns=["metric", "value"])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
