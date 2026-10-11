#!/usr/bin/env python3
"""Shared loaders for the S05 scripts (50/51/52). Aggregate-only helpers; no plotting.

Conventions reused from earlier scripts, not reinvented:
  - unrelated = greedy removal of the highest-degree node at kin >= 0.0442 (41_kir_pilot.py,
    24_novelty_by_field.py), so "unrelated" means the same thing across hla_popgen.
  - ancestry_pred is lowercase on disk (quirk #30) -> upper-cased on load.
  - HLA GTF gene names are "HLA-"-prefixed, KIR bare (quirk #44) -> norm_gene() strips "HLA-".
  - Disclosure goes through _disclosure.py (PUBLIC default) when importable; a local fallback
    applies the same hard rule (participant counts 1-19 -> "<20") if it is not on the path.
"""
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

KIN_MIN = 0.0442
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
PROB_COLS = [f"p_{a.lower()}" for a in ANCESTRY_ORDER]
ASSEMBLY_PLATFORMS = ("revio", "sequel2e", "sequel2")
DEFAULT_TABLE1 = os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv")
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
    "samples_relatedness.tsv")
DEFAULT_KIR_OUTROOT = os.path.expanduser("~/pipeline_outputs_kir")

try:
    import _disclosure as _disc  # noqa: E402
except ImportError:  # flattened deploy without _disclosure on PYTHONPATH
    _disc = None


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def norm_gene(g):
    g = str(g)
    return g[4:] if g.startswith("HLA-") else g


def mask_n(n, mode="public"):
    """Participant/carrier count -> display string. 0 and >=20 exact; 1-19 -> '<20' (PUBLIC)."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return ""
    n = int(n)
    if _disc is not None:
        return _disc.suppressed(n, mode)
    if mode == "public" and 0 < n < 20:
        return "<20"
    return str(n)


def wilson(k, n, z=1.96):
    if n <= 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, max(0.0, c - h), min(1.0, c + h))


def rate_cells(k, n, mode="public", digits=2):
    """Rate with numerator/denominator masking: blank the rate if either side is 1-19 (PUBLIC)."""
    k, n = int(k), int(n)
    small = mode == "public" and ((0 < k < 20) or (0 < n < 20))
    p, lo, hi = wilson(k, n)
    if small or n == 0:
        return {"n": mask_n(k, mode), "d": mask_n(n, mode), "pct": "", "lo": "", "hi": ""}
    return {"n": str(k), "d": str(n), "pct": f"{100*p:.{digits}f}",
            "lo": f"{100*lo:.{digits}f}", "hi": f"{100*hi:.{digits}f}"}


def load_cohort(path):
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    df["ancestry_pred"] = df["ancestry_pred"].astype(str).str.upper()
    for c in PROB_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    present = [c for c in PROB_COLS if c in df.columns]
    df["p_max"] = df[present].max(axis=1)
    df["platform"] = df["platform"].astype(str).str.lower()
    return df


def load_relatedness_pairs(path):
    if not path or not os.path.exists(path):
        log(f"WARNING: relatedness table missing at {path!r}; NO kinship filter applied.")
        return []
    df = pd.read_csv(path, sep="\t", dtype={"i.s": str, "j.s": str})
    return list(zip(df["i.s"], df["j.s"], pd.to_numeric(df["kin"], errors="coerce")))


def greedy_unrelated(person_ids, pairs, kin_min=KIN_MIN):
    ids = set(map(str, person_ids))
    adj = defaultdict(set)
    for i, j, k in pairs:
        i, j = str(i), str(j)
        if i == j or not (k >= kin_min) or i not in ids or j not in ids:
            continue
        adj[i].add(j)
        adj[j].add(i)
    removed = []
    while True:
        live = [(len(v), n) for n, v in adj.items() if v]
        if not live:
            break
        maxdeg = max(d for d, _ in live)
        victim = min(n for d, n in live if d == maxdeg)
        for nb in adj[victim]:
            adj[nb].discard(victim)
        adj[victim] = set()
        removed.append(victim)
    return ids - set(removed)


TABLE1_USECOLS = ["person_id", "hap", "contig", "gene", "copy_index", "consensus",
                  "novelty_depth", "template_distance", "gene_start", "gene_end", "cds_mut",
                  "template_warning"]


def load_table1(path, limit_people=None):
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in TABLE1_USECOLS,
                     dtype={"person_id": str, "contig": str, "consensus": str, "cds_mut": str,
                            "template_warning": str})
    df["gene"] = df["gene"].map(norm_gene)
    df["copy_index"] = pd.to_numeric(df["copy_index"], errors="coerce").fillna(1).astype(int)
    df["novelty_depth"] = pd.to_numeric(df["novelty_depth"], errors="coerce")
    if limit_people:
        keep = sorted(df["person_id"].unique())[:limit_people]
        df = df[df["person_id"].isin(set(keep))]
    return df


def two_field(consensus):
    """'HLA-DQA1*05:01:01:02' -> 'DQA1*05:01'; novel/undetermined -> None."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    c = norm_gene(consensus)
    gene, rest = c.split("*", 1)
    f = rest.split(":")
    if len(f) < 2 or "new" in f[:2] or not f[0] or not f[1]:
        return None
    return f"{gene}*{f[0]}:{f[1]}"


def write_tsv(df, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    log(f"wrote {path} ({len(df)} rows)")
