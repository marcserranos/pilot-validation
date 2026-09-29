"""Shared helpers for the numbered report scripts (02+): palette, ancestry colours, person
table, disclosure-safe aggregation, small statistics, figure saving.

Disclosure rule used everywhere: any plotted group, bin or cell must pool >= MIN_PEOPLE
distinct participants (AoU's n < 20 rule, applied to people, not to rows). Helpers here
either enforce it or return the counts needed to enforce it.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")

import cnsplots as cns  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

MIN_PEOPLE = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
REF_YEAR = 2020  # age proxy = REF_YEAR - year_of_birth; AoU draws span ~2017-2023, so +/- ~3 y
AGE_BINS = [18, 30, 40, 50, 60, 70, 120]
AGE_LABELS = ["18–29", "30–39", "40–49", "50–59", "60–69", "70+"]
PANEL = 120


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def nature():
    return list(cns.palettes("Nature"))


def ancestry_colors():
    return dict(zip(ANCESTRY_ORDER, nature()[:len(ANCESTRY_ORDER)]))


def load_people(cohort_path, pheno_dir):
    """research_id, ancestry, sex (Female/Male/NaN), age (REF_YEAR - year_of_birth), age_bin.
    Fails loudly if the demographics cache is missing (it is built by
    query_overlap_phenotypes.py in a Workbench notebook)."""
    cohort = pd.read_csv(os.path.expanduser(cohort_path), sep="\t", dtype=str)
    person_path = os.path.join(os.path.expanduser(pheno_dir), "person.tsv")
    if not os.path.exists(person_path):
        die(f"missing {person_path} -- run query_overlap_phenotypes.py (notebook) first")
    person = pd.read_csv(person_path, sep="\t", dtype={"research_id": str})
    df = cohort[["research_id", "ancestry"]].merge(person, on="research_id", how="left")
    sex = df["sex_at_birth"].astype(str).str.lower()
    df["sex"] = np.where(sex == "female", "Female", np.where(sex == "male", "Male", None))
    df["age"] = REF_YEAR - pd.to_numeric(df["year_of_birth"], errors="coerce")
    df["age_bin"] = pd.cut(df["age"], AGE_BINS, right=False, labels=AGE_LABELS)
    return df[["research_id", "ancestry", "sex", "age", "age_bin"]]


def save(stem, outdir):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        cns.savefig(os.path.join(outdir, f"{stem}.{ext}"))
    plt.close("all")
    print(f"  wrote {stem}.svg/.pdf/.png", file=sys.stderr)


def style(ax, title=None, xlabel=None, ylabel=None):
    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)


def colorbar(ax, mappable, label):
    cax = ax.inset_axes([1.06, 0.15, 0.05, 0.7])
    cb = plt.colorbar(mappable, cax=cax)
    cb.set_label(label)
    cb.outline.set_linewidth(0.5)
    return cb


def quantile_summary(df, group, value, min_n=MIN_PEOPLE, order=None):
    """Per-group n, median, IQR, 5-95% of a per-person value; groups under min_n dropped."""
    g = df.dropna(subset=[group, value]).groupby(group, observed=True)[value]
    out = pd.DataFrame({"n": g.size(), "median": g.median(), "q25": g.quantile(.25),
                        "q75": g.quantile(.75), "q05": g.quantile(.05), "q95": g.quantile(.95)})
    out = out[out["n"] >= min_n]
    if order is not None:
        out = out.reindex([o for o in order if o in out.index])
    return out


def pointrange(ax, summary, colors=None, xlabels=None, marker_size=18):
    """Median dot, IQR thick bar, 5-95% thin bar -- a box plot without outlier points."""
    x = np.arange(len(summary))
    cols = colors if colors is not None else [nature()[3]] * len(summary)
    for i, (_, r) in enumerate(summary.iterrows()):
        ax.plot([i, i], [r.q05, r.q95], color=cols[i], lw=0.6, solid_capstyle="butt")
        ax.plot([i, i], [r.q25, r.q75], color=cols[i], lw=2.4, solid_capstyle="butt")
        ax.scatter([i], [r["median"]], s=marker_size, color="white", edgecolors=cols[i],
                   linewidths=1.0, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels if xlabels is not None else summary.index, rotation=0)
    ax.set_xlim(-0.6, len(summary) - 0.4)


def safe_hist(ax, values, bins, color, min_n=MIN_PEOPLE, **kw):
    """Histogram of per-person values; bins with < min_n people are not drawn."""
    counts, edges = np.histogram(values, bins=bins)
    keep = counts >= min_n
    ax.bar(edges[:-1][keep], counts[keep], width=np.diff(edges)[keep], align="edge",
           color=color, linewidth=0, **kw)
    return int(counts[~keep].sum())


def wilson(k, n, z=1.96):
    """95% Wilson interval for a proportion (vectorised)."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    p = np.where(n > 0, k / np.maximum(n, 1), np.nan)
    den = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return p, centre - half, centre + half


def design(df, numeric, categorical, ref):
    """Dummy-coded design matrix with intercept. ref: {col: reference level}."""
    parts = [pd.Series(1.0, index=df.index, name="Intercept")]
    for c in numeric:
        parts.append(df[c].astype(float).rename(c))
    for c in categorical:
        levels = [lv for lv in pd.unique(df[c].dropna()) if lv != ref[c]]
        for lv in sorted(levels, key=str):
            parts.append((df[c] == lv).astype(float).rename(f"{c}[{lv}]"))
    return pd.concat(parts, axis=1)


def ols(y, X):
    """OLS with classical SEs. Returns DataFrame coef, se, lo, hi, p."""
    Xv, yv = X.to_numpy(float), np.asarray(y, float)
    beta, *_ = np.linalg.lstsq(Xv, yv, rcond=None)
    resid = yv - Xv @ beta
    dof = len(yv) - Xv.shape[1]
    sigma2 = resid @ resid / dof
    cov = sigma2 * np.linalg.inv(Xv.T @ Xv)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    q = stats.t.ppf(0.975, dof)
    return pd.DataFrame({"coef": beta, "se": se, "lo": beta - q * se, "hi": beta + q * se,
                         "p": p}, index=X.columns)


def logit(y, X, iters=50):
    """Logistic regression by IRLS with Wald SEs. Returns coef table on the log-odds scale."""
    Xv, yv = X.to_numpy(float), np.asarray(y, float)
    beta = np.zeros(Xv.shape[1])
    for _ in range(iters):
        eta = np.clip(Xv @ beta, -30, 30)
        mu = 1 / (1 + np.exp(-eta))
        w = np.maximum(mu * (1 - mu), 1e-10)
        H = Xv.T @ (Xv * w[:, None])
        step = np.linalg.solve(H, Xv.T @ (yv - mu))
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    eta = np.clip(Xv @ beta, -30, 30)
    mu = 1 / (1 + np.exp(-eta))
    H = Xv.T @ (Xv * (mu * (1 - mu))[:, None])
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    z = beta / se
    return pd.DataFrame({"coef": beta, "se": se, "lo": beta - 1.96 * se,
                         "hi": beta + 1.96 * se, "p": 2 * stats.norm.sf(np.abs(z))},
                        index=X.columns)


def forest(ax, table, labels, transform=np.exp, ref_line=1.0, color=None, xlabel=""):
    """Horizontal effect plot from a coef table (rows in `labels` order, top to bottom)."""
    color = color or nature()[3]
    rows = table.loc[list(labels)]
    y = np.arange(len(rows))[::-1]
    est, lo, hi = transform(rows["coef"]), transform(rows["lo"]), transform(rows["hi"])
    ax.hlines(y, lo, hi, color=color, lw=1.0)
    ax.scatter(est, y, s=16, color=color, zorder=3)
    ax.axvline(ref_line, color="0.6", lw=0.6, ls="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[k] for k in rows.index])
    ax.set_xlabel(xlabel)


def rarefied_richness(reads, depth):
    """Expected number of distinct clonotypes in a random subsample of `depth` reads, drawn
    without replacement (Hurlbert 1971). reads: integer read counts per clonotype. NaN if the
    repertoire has fewer than `depth` reads. Exact, via log-binomials."""
    from scipy.special import gammaln
    n = np.asarray(reads, float)
    N = n.sum()
    if N < depth:
        return np.nan

    def lnC(a, b):
        return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)

    rest = N - n
    miss = np.where(rest >= depth, np.exp(lnC(rest, depth) - lnC(N, depth)), 0.0)
    return float(np.sum(1 - miss))


# ------------------------------------------------------------------ OLGA generation probability
_PGEN = None


def _init_olga():
    global _PGEN
    import olga
    import olga.generation_probability as pg
    import olga.load_model as lm
    d = os.path.join(os.path.dirname(olga.__file__), "default_models", "human_T_beta")
    g = lm.GenomicDataVDJ()
    g.load_igor_genomic_data(os.path.join(d, "model_params.txt"),
                             os.path.join(d, "V_gene_CDR3_anchors.csv"),
                             os.path.join(d, "J_gene_CDR3_anchors.csv"))
    m = lm.GenerativeModelVDJ()
    m.load_and_process_igor_model(os.path.join(d, "model_marginals.txt"))
    _PGEN = pg.GenerationProbabilityVDJ(m, g)


def _pgen_chunk(seqs):
    return [_PGEN.compute_aa_CDR3_pgen(s) for s in seqs]


def olga_pgen(seqs, workers):
    """Generation probability of each CDR3 amino-acid sequence under OLGA's default human
    TRB model, marginalised over V and J (Sethna et al. 2019). Parallel over processes."""
    from concurrent.futures import ProcessPoolExecutor
    seqs = list(seqs)
    if not seqs:
        return np.array([])
    chunks = [seqs[i:i + 500] for i in range(0, len(seqs), 500)]
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_olga) as ex:
        return np.concatenate([np.asarray(c, float) for c in ex.map(_pgen_chunk, chunks)])


def olga_generate(n, seed=0):
    """Draw n productive TRB CDR3s from OLGA's human generative model.
    Returns DataFrame(trbv, cdr3aa) -- the recombination null: what the V(D)J machinery
    makes by chance, before any antigen selection."""
    import olga
    import olga.load_model as lm
    import olga.sequence_generation as sq
    d = os.path.join(os.path.dirname(olga.__file__), "default_models", "human_T_beta")
    g = lm.GenomicDataVDJ()
    g.load_igor_genomic_data(os.path.join(d, "model_params.txt"),
                             os.path.join(d, "V_gene_CDR3_anchors.csv"),
                             os.path.join(d, "J_gene_CDR3_anchors.csv"))
    m = lm.GenerativeModelVDJ()
    m.load_and_process_igor_model(os.path.join(d, "model_marginals.txt"))
    gen = sq.SequenceGenerationVDJ(m, g)
    np.random.seed(seed)
    rows = [gen.gen_rnd_prod_CDR3() for _ in range(n)]
    return pd.DataFrame({"trbv": [g.genV[r[2]][0].split("*")[0] for r in rows],
                         "cdr3aa": [r[1] for r in rows]})


def matched_decoys(real, pool, n_per, seed=0, tol=0.25):
    """For each row of `real` (trbv, cdr3aa, pgen), draw n_per decoys from `pool`
    (trbv, cdr3aa, pgen) with the SAME TRBV gene and CDR3 length and log10 Pgen within
    `tol`, widening the tolerance if needed. Controls the three properties that drive
    chance matching, so any excess of real over decoy matches is antigen-driven, not
    recombination-driven. Returns DataFrame(rep, trbv, cdr3aa, vid_real)."""
    rng = np.random.default_rng(seed)
    pool = pool[pool["pgen"] > 0].copy()
    pool["L"] = pool["cdr3aa"].str.len()
    pool["lp"] = np.log10(pool["pgen"])
    idx = {k: (g["lp"].to_numpy(), g["cdr3aa"].to_numpy(), g["trbv"].to_numpy())
           for k, g in pool.groupby(["trbv", "L"], observed=True)}
    real = real[real["pgen"] > 0]
    out = []
    unmatched = 0
    for r in real.itertuples():
        key = (r.trbv, len(r.cdr3aa))
        if key not in idx:
            unmatched += 1
            continue
        lp, seqs, vs = idx[key]
        target = np.log10(r.pgen)
        for t in (tol, 2 * tol, 4 * tol, np.inf):
            ok = np.abs(lp - target) <= t if np.isfinite(t) else np.ones(len(lp), bool)
            if ok.sum() >= 1:
                break
        pick = rng.choice(np.flatnonzero(ok), size=n_per, replace=ok.sum() < n_per)
        for rep, j in enumerate(pick):
            out.append((rep, vs[j], seqs[j], r.vid))
    if unmatched:
        print(f"  decoys: {unmatched:,} real TCRs had no same-V/length pool entry",
              file=sys.stderr)
    return pd.DataFrame(out, columns=["rep", "trbv", "cdr3aa", "vid_real"])
