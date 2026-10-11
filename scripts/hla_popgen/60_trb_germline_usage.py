#!/usr/bin/env python3
"""Is "TRBV gene never used" in a person's T-cell repertoire explained by germline copy number 0?

Sprint S05, WS-B (SPRINT.md board). A per-person TRBV usage table built from RNA-seq repertoires
(Aleix's W1 V-usage counts) shows some TRBV genes with zero usage in a sizeable share of people
(TRBV4-3 the leading case). The deletion/insertion polymorphism spanning TRBV4-3, TRBV3-2 and
TRBV6-2/6-3 is known prior art (READING.md section B), so the question here is not whether it exists
but how much of "never used" it explains, by ancestry, and what is left over.

Copy-number input (`--cn-table`, `--cn-source {auto,58,65}`):
  58  58_regional_depth.py INTERNAL trb_cn_per_person.tsv: person_id, gene, rel_depth, copies,
      assessable. Depth over genes on GRCh38 primary chr7 only. 58 marks paralog-confounded genes
      (median MAPQ20/MAPQ0 depth < 0.8) as not assessable; their copies are reported but flagged.
      TRBV4-3, TRBV3-2, TRBV6-3 and TRBV7-8 are NOT on primary chr7: they sit on the ~21 kb TRB
      insertion, which GRCh38 carries only on an ALT contig, so 58 cannot measure them directly.
      58 mode therefore adds a labelled PROXY: reads from the insertion's copies map (at MAPQ 0) onto
      their primary-assembly paralogs TRBV4-2 / TRBV6-2 / TRBV3-1 and lift their depth. Per person,
      proxy_rel = median rel_depth over those paralogs (`--proxy-genes`); insertion present when
      proxy_rel >= `--proxy-min-rel` (1.5: no insertion ~1.0, one copy ~1.5, two copies ~2.0);
      proxy copies of the insertion genes = 0 / 1 / 2 by `--proxy-cuts` (1.25, 1.75; halfway between
      those expected values). Rows built this way carry cn_method = "proxy_58_paralog_excess".
  65  65_trb_assembly_genotype.py INTERNAL per_person_gene.tsv: person_id, gene (or 'feature'),
      copies (0/1/2), optional method / status / hap1_present / hap2_present (presence per haplotype
      from each person's assemblies). Rows with status != 'ok' (a haplotype's TRB locus incomplete,
      so copies is only a lower bound) are dropped; cn_method = the method column, or "assembly_65".
      65 merges near-identical genes into one class (e.g. 'TRBV6-2/6-3'); usage for such a class is
      the sum over its member genes. Derived features ('<class>[copy2]') join no usage column and
      are used only in the haplotype concordance table.
  auto  58 if the table has rel_depth and assessable columns, else 65.

Usage input (`--usage`, Aleix's per-person TRBV usage table; path and format unknown in advance):
  - .pkl with a dict {"V": person x gene counts, "J": ...} (the W1 cache vj_usage.pkl), or a pickled
    DataFrame; .tsv/.csv/.parquet in long form (`--id-col`, `--gene-col`, `--value-col`) or wide form
    (id column + one column per gene). Gene names lose any allele suffix ("TRBV4-3*01" -> "TRBV4-3").
  - Without `--usage`, files matching `--usage-glob` under `--usage-root` (default
    ~/workspace/aleix_rnaseq/, depth <= 4) are searched; exactly one match is used, otherwise the
    script exits listing the candidates.
  - `--value-kind counts|fraction|auto`. With counts, a zero is judged only where a chance zero is
    unlikely: P(chance zero) = (1 - p_g)^N, p_g = the gene's pooled usage over all people, N = the
    person's total TRBV count; never_used is assessed only when P < `--max-chance-zero` (0.01). This
    is the same rule as the W1 zero-usage estimator (pooled p_g is conservative: people lacking the
    gene only lower it). With fractions the totals are unknown and the filter is skipped (logged).

Cohort: people in both tables, joined to cohort_membership.tsv for ancestry; unrelated only by
default (_s05_common.greedy_unrelated, kin >= 0.0442); `--include-related` skips it. The id in the
usage table is assumed to be the AoU person_id (research_id).

PUBLIC outputs (`--out-dir`, aggregate only; participant counts through _disclosure, 1-19 -> "<20";
rates blank when their numerator, denominator or the complementary count is 1-19):
    contingency.tsv     gene x ancestry x copies (0/1/2): people, never used, used, % never used
    diagnostic.tsv      gene x ancestry: CN0 as a test for never-used. TP = CN0 & never used, FP = CN0
                        & used, FN = CN>=1 & never used, TN = CN>=1 & used; sensitivity (share of
                        never-used explained by CN0), specificity, PPV, and the residual: never-used
                        rate among CN>=1, next to the chance-zero rate expected there
    dosage.tsv          gene x ancestry x copies: usage fraction quantiles; Kendall tau-b of usage
                        fraction on copies (all people, and CN>=1 only = dosage beyond presence)
    haplotype_genes.tsv pairs of `--haplotype-genes` (default TRBV4-3, TRBV3-2, TRBV6-3, TRBV7-8,
                        TRBV6-2, TRBV12-4): copy-number concordance (3x3 counts, % identical, Cohen's
                        kappa), never-used concordance, and whether the pair's calls are identical by
                        construction (both from the same proxy)
    run_summary.tsv     inputs, formats detected, people per stage, parameters
INTERNAL (`--internal-dir`, never in the repo, never pulled):
    joined_per_person.tsv one row per person x gene: copies, method, rel_depth, usage, never_used,
                        chance-zero probability, ancestry

Usage (VM, from the repo; pandas/numpy/scipy in the spechla pixi env):
    cd ~/repos/pilot-validation
    pixi run -e spechla -- python3 -u scripts/hla_popgen/60_trb_germline_usage.py \\
        --cn-table ~/s05/internal/58/trb_cn_per_person.tsv --usage <path from the search> \\
        --out-dir ~/s05/results/60 --internal-dir ~/s05/internal/60
    Smoke test first: --people-limit 300 --out-dir ~/s05/results/60_smoke
        --internal-dir ~/s05/internal/60_smoke (own paths, ENVIRONMENT.md quirk #22).

Tests: scripts/hla_popgen/tests/test_60_trb_germline_usage.py (synthetic; no AoU data).
"""
import argparse
import fnmatch
import os
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _disclosure as _disc  # noqa: E402
import _s05_common as C  # noqa: E402

DEFAULT_CN_TABLE = os.path.expanduser("~/s05/internal/58/trb_cn_per_person.tsv")
DEFAULT_USAGE_ROOT = os.path.expanduser("~/workspace/aleix_rnaseq")
DEFAULT_USAGE_GLOB = "*usage*"
USAGE_EXTS = (".pkl", ".pickle", ".tsv", ".csv", ".parquet", ".tsv.gz", ".csv.gz")
DEFAULT_OUT_DIR = os.path.expanduser("~/s05/results/60")
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s05/internal/60")
ANCESTRY_ORDER = C.ANCESTRY_ORDER
INSERTION_GENES = ("TRBV4-3", "TRBV3-2", "TRBV6-3", "TRBV7-8")
PROXY_GENES = ("TRBV4-2", "TRBV6-2", "TRBV3-1")
HAPLOTYPE_GENES = ("TRBV4-3", "TRBV3-2", "TRBV6-3", "TRBV7-8", "TRBV6-2", "TRBV12-4")
PROXY_METHOD = "proxy_58_paralog_excess"
ID_CANDIDATES = ("person_id", "research_id", "sample_id", "rid", "id")
GENE_CANDIDATES = ("gene", "v_gene", "trbv", "v")
VALUE_CANDIDATES = ("count", "n", "clonotypes", "n_clonotypes", "usage", "frac", "fraction", "value")


def log(msg):
    print(f"[60 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def norm_trb(g):
    return str(g).split("*")[0].strip()


# ---------------------------------------------------------------------------
# Lock (as 58/59)
# ---------------------------------------------------------------------------
def acquire_lock(internal_dir, name="60.lock"):
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
# Copy number
# ---------------------------------------------------------------------------
def detect_cn_source(cols):
    cols = set(cols)
    if {"rel_depth", "assessable"} <= cols:
        return "58"
    if {"person_id", "copies"} <= cols and ({"gene", "feature"} & cols):
        return "65"
    raise SystemExit(f"FATAL: --cn-table columns {sorted(cols)} match neither 58's trb_cn_per_person.tsv "
                     f"(person_id, gene, rel_depth, copies, assessable) nor 65's table (person_id, gene, "
                     f"copies[, method])")


def _as_bool(s):
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def load_cn(path, source, proxy_genes, proxy_min_rel, proxy_cuts, insertion_genes):
    """-> long table person_id, gene, copies, cn_assessable, cn_method, rel_depth (NaN for 65)."""
    d = pd.read_csv(path, sep="\t", comment="#", dtype={"person_id": str, "gene": str, "feature": str})
    if source == "auto":
        source = detect_cn_source(d.columns)
    if "gene" not in d.columns and "feature" in d.columns:
        d = d.rename(columns={"feature": "gene"})  # 65's per_person_gene.tsv names it 'feature'
    d = d.drop_duplicates(["person_id", "gene"], keep="last")
    d["gene"] = d["gene"].map(norm_trb)
    d["copies"] = pd.to_numeric(d["copies"], errors="coerce")
    if source == "58":
        d["rel_depth"] = pd.to_numeric(d["rel_depth"], errors="coerce")
        d["cn_assessable"] = _as_bool(d["assessable"])
        d["cn_method"] = "depth_58"
        out = d[["person_id", "gene", "copies", "cn_assessable", "cn_method", "rel_depth"]]
        proxy = insertion_proxy(d, proxy_genes, proxy_min_rel, proxy_cuts)
        if len(proxy):
            rows = []
            for g in insertion_genes:
                if g in set(out["gene"]):
                    continue  # measured directly: never overwritten by the proxy
                p = proxy.assign(gene=g)
                rows.append(p)
            if rows:
                out = pd.concat([out] + rows, ignore_index=True)
        return out.reset_index(drop=True), source, proxy
    if "method" in d.columns:
        d["cn_method"] = d["method"].fillna("assembly_65").astype(str)
    else:
        d["cn_method"] = "assembly_65"
    d["cn_assessable"] = d["copies"].isin([0, 1, 2])
    if "status" in d.columns:
        # 65: 'incomplete' = a haplotype's TRB locus is not complete, copies is a lower bound
        d["cn_assessable"] &= d["status"].astype(str) == "ok"
    n_drop = int((~d["cn_assessable"]).sum())
    if n_drop:
        log(f"65: {n_drop} person x gene rows not assessable (status != ok or copies not 0/1/2) dropped")
    d = d[d["cn_assessable"]]
    d["rel_depth"] = np.nan
    return (d[["person_id", "gene", "copies", "cn_assessable", "cn_method", "rel_depth"]]
            .reset_index(drop=True), source, pd.DataFrame())


def insertion_proxy(cn58, proxy_genes, min_rel, cuts):
    """Per person: median rel_depth over the paralogs present in 58's table -> insertion present /
    proxy copies. Returns person_id, copies, cn_assessable, cn_method, rel_depth, insertion_present."""
    sub = cn58[cn58["gene"].isin(proxy_genes)]
    if sub.empty:
        log(f"proxy: none of {proxy_genes} in the 58 table; no insertion proxy")
        return pd.DataFrame()
    m = sub.groupby("person_id")["rel_depth"].median()
    lo, hi = cuts
    copies = np.where(m < lo, 0, np.where(m < hi, 1, 2))
    return pd.DataFrame({"person_id": m.index, "copies": copies, "cn_assessable": True,
                         "cn_method": PROXY_METHOD, "rel_depth": m.to_numpy(),
                         "insertion_present": (m.to_numpy() >= min_rel).astype(int)})


# ---------------------------------------------------------------------------
# Usage table
# ---------------------------------------------------------------------------
def find_usage(root, pattern, max_depth=4):
    """Files under root (depth <= max_depth) matching pattern with a tabular extension."""
    root = os.path.expanduser(root)
    hits, tabular = [], []
    if not os.path.isdir(root):
        return hits, tabular
    base = root.rstrip(os.sep).count(os.sep)
    for d, dirs, files in os.walk(root):
        if d.count(os.sep) - base >= max_depth:
            dirs[:] = []
        for f in files:
            if f.lower().endswith(USAGE_EXTS):
                p = os.path.join(d, f)
                tabular.append(p)
                if fnmatch.fnmatch(f.lower(), pattern.lower()):
                    hits.append(p)
    return sorted(hits), sorted(tabular)


def _pick(cols, wanted, given=None, what="column"):
    if given:
        if given not in cols:
            raise SystemExit(f"FATAL: usage table has no {what} {given!r}; columns: {list(cols)[:30]}")
        return given
    low = {str(c).lower(): c for c in cols}
    for w in wanted:
        if w in low:
            return low[w]
    return None


def _read_any(path):
    p = path.lower()
    if p.endswith((".pkl", ".pickle")):
        return pd.read_pickle(path)
    if p.endswith(".parquet"):
        return pd.read_parquet(path)
    sep = "," if p.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(path, sep=sep, comment="#")


def load_usage(path, id_col=None, gene_col=None, value_col=None, layout="auto"):
    """-> wide DataFrame (index person_id as str, columns TRBV genes, values counts or fractions)
    and a description string."""
    obj = _read_any(path)
    desc = []
    if isinstance(obj, dict):
        key = next((k for k in ("V", "v", "TRBV") if k in obj), None)
        if key is None:
            raise SystemExit(f"FATAL: pickled dict in {path} has keys {list(obj)}; expected 'V'")
        obj = obj[key]
        desc.append(f"pickle dict['{key}']")
    df = obj
    if not isinstance(df, pd.DataFrame):
        raise SystemExit(f"FATAL: {path} is not a table ({type(df).__name__})")
    cols = list(df.columns)
    gcol = _pick(cols, GENE_CANDIDATES, gene_col, "gene column")
    trb_cols = [c for c in cols if str(c).upper().startswith("TRBV")]
    if layout == "long" or (layout == "auto" and gcol is not None and not trb_cols):
        icol = _pick(cols, ID_CANDIDATES, id_col, "id column")
        vcol = _pick(cols, VALUE_CANDIDATES, value_col, "value column")
        if icol is None or gcol is None or vcol is None:
            raise SystemExit(f"FATAL: long-format usage table needs id/gene/value columns "
                             f"(--id-col/--gene-col/--value-col); columns: {cols[:30]}")
        d = df[[icol, gcol, vcol]].copy()
        d[gcol] = d[gcol].map(norm_trb)
        d = d[d[gcol].str.upper().str.startswith("TRBV")]
        d[vcol] = pd.to_numeric(d[vcol], errors="coerce").fillna(0)
        w = d.pivot_table(index=icol, columns=gcol, values=vcol, aggfunc="sum", fill_value=0)
        desc.append(f"long: id={icol}, gene={gcol}, value={vcol}")
    else:
        icol = _pick(cols, ID_CANDIDATES, id_col, "id column")
        w = df.set_index(icol) if icol is not None else df
        w = w[[c for c in w.columns if str(c).upper().startswith("TRBV")]]
        w = w.apply(pd.to_numeric, errors="coerce").fillna(0)
        w.columns = [norm_trb(c) for c in w.columns]
        w = w.T.groupby(level=0).sum().T  # allele-level columns -> gene
        desc.append(f"wide: id={'index' if icol is None else icol}, {w.shape[1]} TRBV columns")
    w.index = w.index.astype(str)
    w.index.name = "person_id"
    w = w.groupby(level=0).sum()
    return w, "; ".join(desc)


def class_members(name):
    """65 merges near-identical genes into one class named like 'TRBV6-2/6-3' -> ['TRBV6-2',
    'TRBV6-3']. Derived features ('...[copy2]') have no usage counterpart -> []."""
    if "[" in name:
        return []
    parts = name.split("/")
    if len(parts) == 1:
        return [name]
    stem = parts[0].rsplit("V", 1)[0] + "V" if "V" in parts[0] else ""
    return [parts[0]] + [p if p.upper().startswith("TRB") else stem + p for p in parts[1:]]


def add_class_usage(w, cn_genes):
    """For each merged CN class whose members are all usage columns, add a summed usage column under
    the class name, so usage and copy number refer to the same set of genes."""
    added = []
    for g in sorted(cn_genes):
        mem = class_members(g)
        if len(mem) > 1 and g not in w.columns and all(m in w.columns for m in mem):
            w[g] = w[mem].sum(axis=1)
            added.append(g)
    return w, added


def detect_value_kind(w):
    v = w.to_numpy(dtype=float)
    if np.allclose(v, np.round(v)) and np.nanmax(v.sum(axis=1)) > 1.5:
        return "counts"
    return "fraction"


def usage_long(w, kind, max_chance_zero, min_total, total_cols=None):
    """Wide usage -> long person_id, gene, usage_value, usage_total, usage_frac, never_used,
    chance_zero, usage_assessable. Totals over `total_cols` (default all columns) so merged-class
    columns are not double counted."""
    tot = w[total_cols].sum(axis=1) if total_cols is not None else w.sum(axis=1)
    frac = w.div(tot.replace(0, np.nan), axis=0)
    rows = []
    if kind == "counts":
        pooled = w.sum(axis=0) / float(tot.sum())
    for g in w.columns:
        r = pd.DataFrame({"person_id": w.index, "gene": g, "usage_value": w[g].to_numpy(),
                          "usage_total": tot.to_numpy(), "usage_frac": frac[g].to_numpy()})
        r["never_used"] = (r["usage_value"] <= 0).astype(int)
        if kind == "counts":
            r["chance_zero"] = (1 - float(pooled[g])) ** r["usage_total"]
            r["usage_assessable"] = (r["chance_zero"] < max_chance_zero) & (r["usage_total"] >= min_total)
        else:
            r["chance_zero"] = np.nan
            r["usage_assessable"] = r["usage_total"] > 0
        rows.append(r)
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Disclosure helpers
# ---------------------------------------------------------------------------
def _small(*vals):
    return any(0 < int(v) < _disc.SUPPRESS_BELOW for v in vals)


def _cnt(n, mode):
    return _disc.mask(int(n), _disc.PARTICIPANT, mode).display


def rate(n, d, mode, digits=2):
    """mask_rate + complementary rule: in PUBLIC the rate is also blank when d - n is 1-19 (the
    complement is a subtraction away from a published n and d)."""
    r = _disc.mask_rate(int(n), int(d), mode, pct_decimals=digits)
    if mode == _disc.PUBLIC and _small(int(d) - int(n)):
        r = dict(r, pct="", lo="", hi="")
    return r


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def _strata(d):
    yield "ALL", d
    for a in ANCESTRY_ORDER:
        s = d[d["ancestry"] == a]
        if len(s):
            yield a, s


def gene_meta(j):
    m = j.groupby("gene").agg(cn_method=("cn_method", "first"),
                              cn_assessable=("cn_assessable", "mean")).reset_index()
    m["cn_assessable"] = m["cn_assessable"] >= 0.5
    return dict(zip(m["gene"], zip(m["cn_method"], m["cn_assessable"])))


def contingency(j, mode):
    meta = gene_meta(j)
    rows = []
    for g, dg in j.groupby("gene"):
        for anc, s in _strata(dg):
            for k in (0, 1, 2):
                b = s[s["copies"] == k]
                n, nv = len(b), int(b["never_used"].sum())
                r = rate(nv, n, mode)
                row = dict(gene=g, cn_method=meta[g][0], cn_assessable=meta[g][1], ancestry=anc,
                           copies=k, n_people=_cnt(n, mode), n_never_used=_cnt(nv, mode),
                           n_used=_cnt(n - nv, mode), pct_never_used=r["pct"], lo=r["lo"], hi=r["hi"])
                if mode == _disc.PUBLIC and _small(nv) != _small(n - nv):
                    # exactly one half of the split is 1-19: blank the other half, since n minus it
                    # would reveal the masked one
                    row["n_used" if _small(nv) else "n_never_used"] = ""
                rows.append(row)
    return pd.DataFrame(rows)


def diagnostic(j, mode):
    meta = gene_meta(j)
    rows = []
    for g, dg in j.groupby("gene"):
        for anc, s in _strata(dg):
            cn0 = s["copies"] == 0
            nv = s["never_used"] == 1
            tp, fp = int((cn0 & nv).sum()), int((cn0 & ~nv).sum())
            fn, tn = int((~cn0 & nv).sum()), int((~cn0 & ~nv).sum())
            sens, spec, ppv = rate(tp, tp + fn, mode), rate(tn, tn + fp, mode), rate(tp, tp + fp, mode)
            resid = rate(fn, fn + tn, mode)
            ge1 = s[~cn0]
            exp0 = ge1["chance_zero"].mean() if len(ge1) and ge1["chance_zero"].notna().any() else np.nan
            hide_exp = mode == _disc.PUBLIC and _small(len(ge1))
            rows.append(dict(gene=g, cn_method=meta[g][0], cn_assessable=meta[g][1], ancestry=anc,
                             n_people=_cnt(len(s), mode), TP_cn0_never=_cnt(tp, mode),
                             FP_cn0_used=_cnt(fp, mode), FN_cn_ge1_never=_cnt(fn, mode),
                             TN_cn_ge1_used=_cnt(tn, mode),
                             sensitivity_pct=sens["pct"], sensitivity_lo=sens["lo"], sensitivity_hi=sens["hi"],
                             specificity_pct=spec["pct"], ppv_pct=ppv["pct"], ppv_lo=ppv["lo"], ppv_hi=ppv["hi"],
                             residual_never_used_cn_ge1_pct=resid["pct"], residual_lo=resid["lo"],
                             residual_hi=resid["hi"],
                             expected_chance_zero_cn_ge1_pct=("" if (hide_exp or not np.isfinite(exp0))
                                                              else round(100 * exp0, 3))))
    return pd.DataFrame(rows)


def _kendall(x, y):
    from scipy.stats import kendalltau
    if len(x) < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan, np.nan
    t, p = kendalltau(x, y)
    return float(t), float(p)


def dosage(j, mode):
    meta = gene_meta(j)
    rows = []
    for g, dg in j.groupby("gene"):
        dg = dg[dg["usage_frac"].notna()]
        for anc, s in _strata(dg):
            t_all, p_all = _kendall(s["copies"].to_numpy(), s["usage_frac"].to_numpy())
            ge1 = s[s["copies"] >= 1]
            t_ge1, p_ge1 = _kendall(ge1["copies"].to_numpy(), ge1["usage_frac"].to_numpy())
            hide = mode == _disc.PUBLIC and _small(len(s))
            hide1 = mode == _disc.PUBLIC and (_small(len(ge1)) or _small(int((ge1["copies"] == 1).sum()))
                                               or _small(int((ge1["copies"] == 2).sum())))
            for k in (0, 1, 2):
                b = s[s["copies"] == k]["usage_frac"]
                hb = len(b) == 0 or (mode == _disc.PUBLIC and _small(len(b)))
                q = (np.nan,) * 4 if hb else (b.quantile(.25), b.median(), b.quantile(.75), b.mean())
                rows.append(dict(
                    gene=g, cn_method=meta[g][0], cn_assessable=meta[g][1], ancestry=anc, copies=k,
                    n_people=_cnt(len(b), mode),
                    usage_pct_q25="" if hb else round(100 * q[0], 4),
                    usage_pct_median="" if hb else round(100 * q[1], 4),
                    usage_pct_q75="" if hb else round(100 * q[2], 4),
                    usage_pct_mean="" if hb else round(100 * q[3], 4),
                    kendall_tau_b_all="" if (hide or not np.isfinite(t_all)) else round(t_all, 4),
                    kendall_p_all="" if (hide or not np.isfinite(p_all)) else float(f"{p_all:.3g}"),
                    kendall_tau_b_cn_ge1="" if (hide1 or not np.isfinite(t_ge1)) else round(t_ge1, 4),
                    kendall_p_cn_ge1="" if (hide1 or not np.isfinite(p_ge1)) else float(f"{p_ge1:.3g}")))
    return pd.DataFrame(rows)


def cohen_kappa(a, b, cats=(0, 1, 2)):
    a, b = np.asarray(a), np.asarray(b)
    n = len(a)
    if n == 0:
        return np.nan
    po = float(np.mean(a == b))
    pe = sum(float(np.mean(a == c)) * float(np.mean(b == c)) for c in cats)
    return np.nan if pe >= 1 else (po - pe) / (1 - pe)


def haplotype_genes(cn, j, genes, mode):
    """Pairwise CN concordance (from the CN table: does not need usage) and never-used concordance
    (from the joined table) among the haplotype genes."""
    present = [g for g in genes if g in set(cn["gene"])]
    meth = cn.groupby("gene")["cn_method"].first().to_dict()
    ass = cn.groupby("gene")["cn_assessable"].mean().to_dict()
    w = cn[cn["gene"].isin(present)].pivot_table(index="person_id", columns="gene", values="copies",
                                                 aggfunc="first")
    w = w.join(cn[["person_id", "ancestry"]].drop_duplicates("person_id").set_index("person_id"))
    u = j[j["gene"].isin(present)].pivot_table(index="person_id", columns="gene", values="never_used",
                                               aggfunc="first")
    rows = []
    for a, b in combinations(present, 2):
        by_constr = meth.get(a) == meth.get(b) == PROXY_METHOD
        for anc, s in _strata(w.reset_index()):
            s = s.dropna(subset=[a, b])
            x, y = s[a].astype(int).to_numpy(), s[b].astype(int).to_numpy()
            n = len(s)
            same = int((x == y).sum())
            r = rate(same, n, mode)
            row = dict(gene_a=a, gene_b=b, method_a=meth.get(a), method_b=meth.get(b),
                       cn_assessable_a=ass.get(a, 0) >= 0.5, cn_assessable_b=ass.get(b, 0) >= 0.5,
                       identical_by_construction=by_constr, ancestry=anc, n_people=_cnt(n, mode),
                       pct_identical_copies=r["pct"], lo=r["lo"], hi=r["hi"])
            kap = cohen_kappa(x, y)
            row["cohen_kappa"] = "" if (not np.isfinite(kap) or (mode == _disc.PUBLIC and _small(n))) \
                else round(kap, 4)
            for i in (0, 1, 2):
                for k in (0, 1, 2):
                    row[f"n_{i}{k}"] = _cnt(int(((x == i) & (y == k)).sum()), mode)
            # never-used concordance
            if a in u.columns and b in u.columns:
                uu = u.loc[u.index.isin(s["person_id"]), [a, b]].dropna()
                na, nb = uu[a].astype(int).to_numpy(), uu[b].astype(int).to_numpy()
                both, ea = int(((na == 1) & (nb == 1)).sum()), int((na == 1).sum())
                rc = rate(int((na == nb).sum()), len(uu), mode)
                rb = rate(both, ea, mode)
                row.update(n_people_usage=_cnt(len(uu), mode), pct_never_used_concordant=rc["pct"],
                           n_never_used_a=_cnt(ea, mode), n_never_used_both=_cnt(both, mode),
                           pct_b_never_given_a_never=rb["pct"])
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def resolve_usage_path(args):
    if args.usage:
        if not os.path.exists(args.usage):
            raise SystemExit(f"FATAL: --usage {args.usage} not found")
        return args.usage
    hits, tabular = find_usage(args.usage_root, args.usage_glob)
    if len(hits) == 1:
        log(f"usage table found by search: {hits[0]}")
        return hits[0]
    listing = hits if hits else tabular
    msg = [f"FATAL: {'several' if hits else 'no'} usage table(s) matching {args.usage_glob!r} under "
           f"{args.usage_root} (depth <= 4). Pass --usage <path> explicitly "
           f"(and --id-col/--gene-col/--value-col for a long table)."]
    if listing:
        msg.append(("Matching" if hits else "Tabular files found (none match)") + " candidates:")
        msg += [f"  {p}" for p in listing[:50]]
        if len(listing) > 50:
            msg.append(f"  ... and {len(listing) - 50} more")
    elif not os.path.isdir(os.path.expanduser(args.usage_root)):
        msg.append(f"{args.usage_root} does not exist or is not mounted (ENVIRONMENT.md quirk #35).")
    else:
        msg.append("No tabular files at all under that root.")
    raise SystemExit("\n".join(msg))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cn-table", default=DEFAULT_CN_TABLE,
                    help="58's INTERNAL trb_cn_per_person.tsv or 65's INTERNAL per-person table")
    ap.add_argument("--cn-source", choices=["auto", "58", "65"], default="auto")
    ap.add_argument("--proxy-genes", default=",".join(PROXY_GENES),
                    help="58 mode: primary-assembly paralogs whose excess depth proxies the insertion")
    ap.add_argument("--proxy-min-rel", type=float, default=1.5,
                    help="58 mode: insertion present when median paralog rel_depth >= this")
    ap.add_argument("--proxy-cuts", default="1.25,1.75",
                    help="58 mode: proxy copies 0/1/2 of the insertion genes at these rel_depth cuts")
    ap.add_argument("--insertion-genes", default=",".join(INSERTION_GENES),
                    help="58 mode: genes given proxy copies (only those absent from the 58 table)")
    ap.add_argument("--usage", help="Aleix's per-person TRBV usage table (pkl/tsv/csv/parquet)")
    ap.add_argument("--usage-root", default=DEFAULT_USAGE_ROOT)
    ap.add_argument("--usage-glob", default=DEFAULT_USAGE_GLOB)
    ap.add_argument("--usage-layout", choices=["auto", "long", "wide"], default="auto")
    ap.add_argument("--id-col")
    ap.add_argument("--gene-col")
    ap.add_argument("--value-col")
    ap.add_argument("--value-kind", choices=["auto", "counts", "fraction"], default="auto")
    ap.add_argument("--max-chance-zero", type=float, default=0.01,
                    help="counts only: judge a zero only where (1-p_gene)^N_total is below this")
    ap.add_argument("--min-total", type=float, default=0,
                    help="counts only: people with fewer total TRBV counts are not assessed")
    ap.add_argument("--haplotype-genes", default=",".join(HAPLOTYPE_GENES))
    ap.add_argument("--cohort", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--include-related", action="store_true")
    ap.add_argument("--people-limit", type=int, help="smoke test: first N joined person_ids (sorted)")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args(argv)

    internal_dir = os.path.expanduser(args.internal_dir)
    internal_path = os.path.join(internal_dir, "joined_per_person.tsv")
    _disc.assert_internal_path_allowed(internal_path)
    for p in (args.cn_table, args.cohort):
        if not os.path.exists(p):
            raise SystemExit(f"FATAL: input not found: {p}")
    usage_path = resolve_usage_path(args)
    t0 = time.time()
    lock = acquire_lock(internal_dir)
    try:
        cuts = tuple(float(x) for x in args.proxy_cuts.split(","))
        cn, source, proxy = load_cn(args.cn_table, args.cn_source,
                                    [g.strip() for g in args.proxy_genes.split(",") if g.strip()],
                                    args.proxy_min_rel, cuts,
                                    [g.strip() for g in args.insertion_genes.split(",") if g.strip()])
        log(f"CN table ({source}): {cn['person_id'].nunique()} people, {cn['gene'].nunique()} genes")
        w, udesc = load_usage(usage_path, args.id_col, args.gene_col, args.value_col, args.usage_layout)
        kind = detect_value_kind(w) if args.value_kind == "auto" else args.value_kind
        # merged-class columns are added AFTER the kind check and are excluded from totals below
        gene_cols = list(w.columns)
        w, merged = add_class_usage(w, set(cn["gene"]))
        if merged:
            log(f"usage summed over merged CN classes: {merged}")
        if kind == "fraction":
            log("usage values are fractions: totals unknown, chance-zero filter skipped")
        log(f"usage: {udesc}; {len(w)} people, {w.shape[1]} genes, values = {kind}")

        coh = C.load_cohort(args.cohort)[["person_id", "ancestry_pred"]].drop_duplicates("person_id")
        coh = coh.rename(columns={"ancestry_pred": "ancestry"})
        pids = sorted(set(cn["person_id"]) & set(w.index))
        n_both = len(pids)
        if args.people_limit:
            pids = pids[:args.people_limit]
        if not args.include_related:
            keep = C.greedy_unrelated(pids, C.load_relatedness_pairs(args.relatedness_table))
            pids = sorted(p for p in pids if p in keep)
        log(f"people in both tables: {n_both}; analysed: {len(pids)}")
        pset = set(pids)
        cn = cn[cn["person_id"].isin(pset)].merge(coh, on="person_id", how="left")
        cn["ancestry"] = cn["ancestry"].where(cn["ancestry"].isin(ANCESTRY_ORDER), "UNK")
        ul = usage_long(w.loc[[p for p in w.index if p in pset]], kind, args.max_chance_zero, args.min_total,
                        total_cols=gene_cols)
        j = cn.merge(ul, on=["person_id", "gene"], how="inner")
        j = j[j["copies"].isin([0, 1, 2])]
        _disc.write_tsv(j, internal_path, _disc.INTERNAL)
        ja = j[j["usage_assessable"]]

        mode = args.disclosure
        out = os.path.expanduser(args.out_dir)
        os.makedirs(out, exist_ok=True)
        _disc.write_tsv(contingency(ja, mode), os.path.join(out, "contingency.tsv"), mode)
        _disc.write_tsv(diagnostic(ja, mode), os.path.join(out, "diagnostic.tsv"), mode)
        _disc.write_tsv(dosage(j, mode), os.path.join(out, "dosage.tsv"), mode)
        hg = [g.strip() for g in args.haplotype_genes.split(",") if g.strip()]
        _disc.write_tsv(haplotype_genes(cn, ja, hg, mode), os.path.join(out, "haplotype_genes.tsv"), mode)

        cn_genes, u_genes = set(cn["gene"]), set(w.columns)
        run = [dict(metric="cn_source", value=source),
               dict(metric="cn_table", value=os.path.basename(args.cn_table)),
               dict(metric="usage_table", value=os.path.basename(usage_path)),
               dict(metric="usage_format", value=udesc), dict(metric="usage_value_kind", value=kind),
               dict(metric="people_in_both_tables", value=_cnt(n_both, mode)),
               dict(metric="people_analysed", value=_cnt(len(pids), mode)),
               dict(metric="kinship_filter", value="none" if args.include_related else "greedy_kin_ge_0.0442"),
               dict(metric="genes_joined", value=len(cn_genes & u_genes)),
               dict(metric="merged_classes_usage_summed", value=",".join(merged)),
               dict(metric="genes_cn_only", value=",".join(sorted(cn_genes - u_genes))[:2000]),
               dict(metric="genes_usage_only", value=",".join(sorted(u_genes - cn_genes))[:2000]),
               dict(metric="person_gene_pairs_joined", value=len(j)),
               dict(metric="person_gene_pairs_usage_assessable", value=int(j["usage_assessable"].sum())),
               dict(metric="max_chance_zero", value=args.max_chance_zero if kind == "counts" else "n/a"),
               dict(metric="min_total", value=args.min_total)]
        if source == "58":
            run += [dict(metric="proxy_genes", value=args.proxy_genes),
                    dict(metric="proxy_min_rel", value=args.proxy_min_rel),
                    dict(metric="proxy_cuts", value=args.proxy_cuts),
                    dict(metric="proxy_label", value=PROXY_METHOD + " (proxy: paralog excess depth, not a "
                                                                   "direct call)")]
            if len(proxy):
                pp = proxy[proxy["person_id"].isin(pset)]
                run.append(dict(metric="proxy_insertion_present_people",
                                value=_cnt(int(pp["insertion_present"].sum()), mode)))
        run += [dict(metric="people_limit", value=args.people_limit or ""),
                dict(metric="disclosure", value=mode),
                dict(metric="wall_seconds", value=int(time.time() - t0))]
        _disc.write_tsv(pd.DataFrame(run), os.path.join(out, "run_summary.tsv"), mode)
        log(f"done: {len(pids)} people; tables in {out}")
    finally:
        os.remove(lock)
    return 0


if __name__ == "__main__":
    sys.exit(main())
