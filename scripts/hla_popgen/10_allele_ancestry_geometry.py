#!/usr/bin/env python3
"""Per-allele ancestry-centroid geometry: positions ALLELES (not people) in ancestry-proportion
space by the average admixture composition of their carriers, operationalizing the QC heuristic
in research/NOVEL_LIT.md section 5 item 3 -- "a genuinely old, population-restricted allele should
cluster within one or a few ancestry groups... random scatter is more consistent with recurrent
independent assembly error than with real allelic identity by descent" -- and extending it to a
side-by-side comparison of several concrete approaches rather than one fixed choice:

  1. Per-allele ancestry centroid, known alleles (Table 1) + this cohort's novel alleles, plotted in
     a 3-component ternary simplex (VIZ_LIT.md 1.4/2.4 lineage, but one point per ALLELE, not per
     PERSON). NO MINIMUM CARRIER COUNT: an allele with exactly 1 carrier gets a centroid equal to
     that person's own admixture proportions -- deliberate, per explicit instruction. Marker alpha
     (fixed size) scales with carrier count instead, so a singleton renders faint and a common
     allele renders opaque -- overplotting at a simplex vertex (where many alleles cluster) becomes
     visually legible density, not a solid blob, without excluding any point.

  NOTE on novel-allele carrier identity: 03_novel_alleles.py's novel_alleles.tsv (Table 3) is
  aggregate-only by design (SCHEMA.md) and carries no person_id. Rather than modify that script to
  add a carrier-linkage side table -- 00-04 is documented (_viz_common.py's own docstring) as a
  different concurrent agent's ownership boundary, and this repo/VM has a second agent actively
  working in it -- this script independently recomputes the SAME (gene, sha1(observed CDS))
  clustering key directly from Table 1 + the raw per-person cds.fa.gz files (read-only, hash only,
  never exports a sequence), so novel_id values line up with Table 3's without importing or editing
  any 00-04 file. Table 3 itself IS still read (a data file, not code) purely for confidence_tier/
  nearest_allele metadata -- consuming another script's TSV output is the established, sanctioned
  pattern every 05-09 figure script already uses for Tables 1/2/4.
  2. The same centroid, computed two ways (--renorm): renormalize each carrier's proportions to the
     chosen component subset BEFORE averaging (matches the existing per-person plot_ternary in
     06_figures_structure.py) vs. average the full 6-way proportions first and renormalize once. Both
     are always computed side by side in the report table.
  3. Enrichment, computed two independent ways, both always reported (never used to filter what gets
     plotted): a continuous bootstrap (resample n_carriers people from the cohort at random, build a
     null distribution of centroids, report the observed centroid's distance from the cohort mean as
     a z-score/empirical p-value) and a discrete Fisher-exact/chi-square test on carriers' ancestry_
     pred labels vs. the rest of the cohort.
  4. Higher dimensions: a 4-component tetrahedron (3D barycentric projection, the "pyramid") and a
     PCA of the full 6-way per-allele centroid vectors (arbitrary-dimension generalization, reusing
     the same manual PCA already used for per-person dosage vectors in 06_figures_structure.py).
  5. Known alleles are circles; novel (this-cohort) alleles are triangles. The --label-top-known
     most common KNOWN alleles are text-labeled (novel_id values are opaque content hashes, not
     worth labeling, and labeling everything would overlap badly at this density).
  6. Multi-gene views (--mode multi): a combined ternary with EVERY allele across a gene list in one
     figure, colored by GENE identity instead of ancestry (to see where diversity concentrates
     across genes), plus a small-multiples grid of one ancestry-colored ternary per gene. Both reuse
     a single efficient pass over the raw data (one cds.fa.gz scan across all target genes, not one
     scan per gene).

Usage (fixtures): see 05_figures_frequency.py's docstring for the fixture pipeline, then run
03_novel_alleles.py (produces novel_alleles.tsv; --outroot must point at the person-directory root
so this script's own cds.fa.gz re-scan can find them), then:
    python3 scripts/hla_popgen/10_allele_ancestry_geometry.py --outroot /tmp/hla_fixtures/people \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --novel-alleles /tmp/hla_fixtures/novel_alleles.sample.tsv --cohort lr

Real run (VM): python3 scripts/hla_popgen/10_allele_ancestry_geometry.py --cohort lr --gene B
"""
import argparse
import gzip
import hashlib
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 -- registers the '3d' projection

import _viz_common as vc

_HEADER_KEY_RE = re.compile(r"^(.*)_(\d+)$")

DEFAULT_TERNARY = ["AFR", "EUR", "AMR"]
DEFAULT_TETRA = ["AFR", "EUR", "AMR", "EAS"]

# Regular tetrahedron vertices (equal pairwise distances) -- the 4-component barycentric target.
TETRA_VERTS = np.array([
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
])


# ---------------------------------------------------------------------------
# Carrier sets
# ---------------------------------------------------------------------------
def known_allele_carriers(table1, cohort_ids, gene_bare_name, field=2):
    """-> {allele_label: set(person_id)}, restricted to the given cohort. A person is counted once
    per allele even if it appears on both haplotypes or with copy_index>1 (Table 1's grain is finer
    than person -- groupby().apply(set) dedups naturally)."""
    sub = table1[(table1["gene_bare"] == gene_bare_name) &
                 (table1["person_id"].isin(cohort_ids))].copy()
    sub["allele_field"] = sub["consensus"].map(lambda a: vc.to_nfield(a, field))
    sub = sub.dropna(subset=["allele_field"])
    return {allele: set(grp["person_id"]) for allele, grp in sub.groupby("allele_field")}


def _parse_cds_fasta(path):
    """Independent reimplementation of 03_novel_alleles.py's parse_cds_fasta -- deliberately NOT
    imported from that file (00-04 is a different concurrent agent's ownership boundary; see module
    docstring). Parses one hap's cds.fa.gz: '>{contig}_{gene}_{i}' headers -> {rest: [(i, seq), ...]}
    (reference/IMMUANNOT_GTF_SPEC.md part D/E). Read-only; the sequence is used only to compute a
    hash and is never exported, kept, or written anywhere (matches the project's 'novel sequences
    stay on the VM' policy)."""
    index = defaultdict(list)
    if not os.path.exists(path):
        return index
    key, seq_chunks = None, []

    def flush():
        if key is None:
            return
        seq = "".join(seq_chunks).upper()
        m = _HEADER_KEY_RE.match(key)
        if m:
            index[m.group(1)].append((int(m.group(2)), seq))

    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                key = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    flush()
    return index


def novel_allele_carriers(table1, table3, cohort_ids, gene_bare_name, outroot):
    """-> ({novel_id: set(person_id)}, {novel_id: meta_dict}). Independently recomputes the
    (gene, sha1(observed CDS)) clustering key -- the SAME key novel_alleles.tsv (Table 3) itself
    uses -- directly from Table 1's novel rows + the raw per-person cds.fa.gz files, so novel_id
    values line up with Table 3's without importing or modifying 03_novel_alleles.py. Table 3 is
    still read (a data file, not code) purely for confidence_tier/nearest_allele metadata --
    consuming another script's TSV output is the established pattern every 05-09 figure script
    already uses for Tables 1/2/4. Restricted to the target gene only, for cost.

    Table 3's own `gene` column carries the prefixed form (e.g. 'HLA-G', confirmed against real
    fixture output -- novel_id is built as f'{gene}_nov_{sha1}' with that same prefixed value),
    unlike Table 1's derived `gene_bare` column -- exactly the bare-vs-prefixed mismatch
    _viz_common.py's own NONCLASSICAL_GENES comment warns silently drops rows with zero error."""
    gene_prefixed = vc.gene_display(gene_bare_name)
    cohort_set = set(cohort_ids)
    t = table1[(table1["gene_bare"] == gene_bare_name) &
               (table1["person_id"].isin(cohort_set))].copy()
    novel = t[t["is_novel"]]
    if novel.empty:
        return {}, {}

    # Ambiguity check identical to 03_novel_alleles.py's own rule: >1 Table-1 row sharing
    # (person_id, hap, contig, gene) means real copy_index>1, i.e. the cds.fa.gz join is genuinely
    # ambiguous (reference/IMMUANNOT_GTF_SPEC.md part D "Caveat") -- skip those, never guess.
    group_sizes = t.groupby(["person_id", "hap", "contig", "gene"]).size()

    clusters = defaultdict(set)  # sha1 -> set(person_id)
    for (person_id, hap), sub in novel.groupby(["person_id", "hap"]):
        cds_path = os.path.join(outroot, person_id, "immuannot_output", hap, "cds.fa.gz")
        fasta_index = None  # lazy -- don't touch disk for ambiguous-only haps
        for _, row in sub.iterrows():
            key = (row["person_id"], row["hap"], row["contig"], row["gene"])
            if group_sizes.loc[key] > 1:
                continue
            if fasta_index is None:
                fasta_index = _parse_cds_fasta(cds_path)
            rest = f"{row['contig']}_{row['gene']}"
            candidates = fasta_index.get(rest, [])
            if len(candidates) != 1:
                continue
            _i, seq = candidates[0]
            if not seq:
                continue
            sha1 = hashlib.sha1(seq.encode("ascii")).hexdigest()
            clusters[sha1].add(person_id)

    table3_by_id = table3.set_index("novel_id") if "novel_id" in table3.columns else None
    carriers, meta = {}, {}
    for sha1, pids in clusters.items():
        nid = f"{gene_prefixed}_nov_{sha1[:8]}"
        carriers[nid] = pids
        row = table3_by_id.loc[nid] if table3_by_id is not None and nid in table3_by_id.index else {}
        meta[nid] = {"confidence_tier": row.get("confidence_tier"),
                     "nearest_allele": row.get("nearest_allele")}
    return carriers, meta


# ---------------------------------------------------------------------------
# Centroid -- NO MINIMUM CARRIER COUNT. n=1 returns that person's own (renormalized) proportions.
# ---------------------------------------------------------------------------
def compute_centroid(carrier_ids, cohort_people, components, renorm="per_person"):
    """Returns (centroid_dict_or_None, n_used). None/0 only when literally no carrier has usable
    p_* data -- never gated on a minimum count otherwise."""
    cols = [f"p_{c.lower()}" for c in components]
    sub = cohort_people[cohort_people["person_id"].isin(carrier_ids)][["person_id"] + cols]
    sub = sub.dropna(subset=cols)
    n_used = len(sub)
    if n_used == 0:
        return None, 0
    if renorm == "per_person":
        row_sum = sub[cols].sum(axis=1)
        renormed = sub[cols].div(row_sum.where(row_sum != 0, np.nan), axis=0)
        centroid = renormed.mean(axis=0, skipna=True)
    else:  # per_centroid: average full proportions first, renormalize the subset average once
        mean_full = sub[cols].mean(axis=0)
        total = mean_full.sum()
        centroid = mean_full / total if total else mean_full
    return dict(zip(components, centroid.values)), n_used


def other_renorm(renorm):
    return "per_centroid" if renorm == "per_person" else "per_person"


# ---------------------------------------------------------------------------
# Enrichment -- two independent methods, both always computed, neither used to filter plotting.
# ---------------------------------------------------------------------------
def bootstrap_null(cohort_people, n_carriers, components, n_boot=1000, seed=0, renorm="per_person"):
    """Null distribution of centroids from resampling n_carriers people uniformly (with
    replacement) from the eligible cohort. At n_carriers=1 this is just 'one random person's own
    proportions' -- appropriately wide, so a singleton allele will rarely reach significance,
    without the centroid computation itself ever being blocked."""
    rng = np.random.default_rng(seed)
    cols = [f"p_{c.lower()}" for c in components]
    pool = cohort_people[cols].dropna()
    n_pool = len(pool)
    if n_pool == 0 or n_carriers <= 0:
        return np.empty((0, len(components)))
    vals = pool.values
    idx = rng.integers(0, n_pool, size=(n_boot, n_carriers))
    sampled = vals[idx]  # (n_boot, n_carriers, k)
    if renorm == "per_person":
        row_sum = sampled.sum(axis=2, keepdims=True)
        row_sum = np.where(row_sum == 0, np.nan, row_sum)
        renormed = sampled / row_sum
        centroids = np.nanmean(renormed, axis=1)
    else:
        mean_full = sampled.mean(axis=1)
        tot = mean_full.sum(axis=1, keepdims=True)
        centroids = mean_full / np.where(tot == 0, np.nan, tot)
    return centroids


def bootstrap_enrichment(centroid, null_centroids, cohort_mean, components):
    """Euclidean distance (raw barycentric coordinates -- a documented simplification; Aitchison/
    CLR geometry is the compositionally-correct alternative, not implemented here) of the observed
    centroid from the cohort-wide mean, vs. the same distance under the null."""
    obs = np.array([centroid[c] for c in components])
    base = np.array([cohort_mean[c] for c in components])
    obs_dist = float(np.linalg.norm(obs - base))
    if np.isnan(obs_dist):
        # A carrier's mass can legitimately fall entirely outside this component subset (e.g. a
        # 3-way AFR/EUR/AMR ternary for someone whose admixture is 100% EAS/MID/SAS) -- per_person
        # renorm correctly produces NaN (0/0) rather than a misleading value in that case. Without
        # this guard, `null_dists >= nan` is elementwise False everywhere, so .mean() silently
        # returns 0.0 -- a spuriously "significant" p-value on an undefined observation. Report NA
        # instead of guessing.
        return {"distance": np.nan, "z_score": np.nan, "p_value": np.nan}
    if len(null_centroids) == 0:
        return {"distance": obs_dist, "z_score": np.nan, "p_value": np.nan}
    null_dists = np.linalg.norm(null_centroids - base, axis=1)
    null_dists = null_dists[~np.isnan(null_dists)]
    if len(null_dists) == 0:
        return {"distance": obs_dist, "z_score": np.nan, "p_value": np.nan}
    mu = float(null_dists.mean())
    sd = float(null_dists.std(ddof=1)) if len(null_dists) > 1 else 0.0
    z = (obs_dist - mu) / sd if sd > 0 else np.nan
    p = float((null_dists >= obs_dist).mean())
    return {"distance": obs_dist, "z_score": z, "p_value": p}


def discrete_enrichment(carrier_ids, cohort_people):
    """Carriers' ancestry_pred label distribution vs. the rest of the cohort's -- a fully
    independent second read on "is this allele ancestry-restricted?" that never touches continuous
    proportions. Fisher exact (one ancestry vs. rest, safe at tiny n) + an overall chi-square."""
    labeled = cohort_people[cohort_people["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
    is_carrier = labeled["person_id"].isin(carrier_ids)
    carrier_counts = labeled.loc[is_carrier, "ancestry_pred"].value_counts()
    noncarrier_counts = labeled.loc[~is_carrier, "ancestry_pred"].value_counts()
    n_carriers, n_noncarriers = int(is_carrier.sum()), int((~is_carrier).sum())
    if n_carriers == 0 or n_noncarriers == 0:
        return {"chi2_p": np.nan, "fisher_p_by_ancestry": {a: np.nan for a in vc.ANCESTRY_ORDER}}
    fisher_p = {}
    for anc in vc.ANCESTRY_ORDER:
        a = int(carrier_counts.get(anc, 0))
        b = n_carriers - a
        c = int(noncarrier_counts.get(anc, 0))
        d = n_noncarriers - c
        try:
            _, p = scipy_stats.fisher_exact([[a, b], [c, d]])
        except ValueError:
            p = np.nan
        fisher_p[anc] = p
    table = np.array([[int(carrier_counts.get(a, 0)) for a in vc.ANCESTRY_ORDER],
                       [int(noncarrier_counts.get(a, 0)) for a in vc.ANCESTRY_ORDER]])
    try:
        _, chi2_p, _, _ = scipy_stats.chi2_contingency(table)
    except ValueError:
        chi2_p = np.nan
    return {"chi2_p": chi2_p, "fisher_p_by_ancestry": fisher_p}


# ---------------------------------------------------------------------------
# Assemble per-allele records for one component set (ternary, tetrahedron, or full 6-way)
# ---------------------------------------------------------------------------
def build_allele_records(table1, table3, cohort_people, gene_bare_name, components, outroot,
                          field=2, renorm="per_person", n_boot=1000, seed=0):
    cohort_ids = set(cohort_people["person_id"])
    known = known_allele_carriers(table1, cohort_ids, gene_bare_name, field=field)
    novel, novel_meta = novel_allele_carriers(table1, table3, cohort_ids, gene_bare_name, outroot)

    cohort_mean, _ = compute_centroid(cohort_ids, cohort_people, components, renorm=renorm)
    null_cache = {}

    def get_null(n):
        if n not in null_cache:
            null_cache[n] = bootstrap_null(cohort_people, n, components, n_boot=n_boot, seed=seed,
                                            renorm=renorm)
        return null_cache[n]

    def make_record(allele_id, kind, pids, tier=None):
        centroid, n_used = compute_centroid(pids, cohort_people, components, renorm=renorm)
        if centroid is None:
            return None
        centroid_alt, _ = compute_centroid(pids, cohort_people, components,
                                            renorm=other_renorm(renorm))
        enrich = bootstrap_enrichment(centroid, get_null(n_used), cohort_mean, components)
        disc = discrete_enrichment(pids, cohort_people)
        return {"allele_id": allele_id, "kind": kind, "n_carriers": n_used,
                "centroid": centroid, "centroid_alt": centroid_alt,
                "confidence_tier": tier, **enrich, **disc}

    records = []
    for label, pids in known.items():
        r = make_record(label, "known", pids)
        if r:
            records.append(r)
    for nid, pids in novel.items():
        r = make_record(nid, "novel", pids, tier=novel_meta[nid]["confidence_tier"])
        if r:
            records.append(r)
    return records, cohort_mean


# ---------------------------------------------------------------------------
# Multi-gene version: ONE pass over the raw data across every gene in `genes`, instead of calling
# the single-gene functions above once per gene (which would re-open every person's cds.fa.gz once
# per gene it has a novel call in, instead of once total). Each record carries an extra "gene" key.
# Every record's allele_id is gene-prefixed (e.g. "B*07:02") since alleles from different genes can
# otherwise collide on the bare 2-field string.
# ---------------------------------------------------------------------------
def build_records_for_genes(table1, table3, cohort_people, genes, components, outroot,
                             field=2, renorm="per_person", n_boot=1000, seed=0):
    cohort_ids = set(cohort_people["person_id"])
    cohort_mean, _ = compute_centroid(cohort_ids, cohort_people, components, renorm=renorm)
    null_cache = {}

    def get_null(n):
        if n not in null_cache:
            null_cache[n] = bootstrap_null(cohort_people, n, components, n_boot=n_boot, seed=seed,
                                            renorm=renorm)
        return null_cache[n]

    def make_record(allele_id, gene_bare_name, kind, pids, tier=None):
        centroid, n_used = compute_centroid(pids, cohort_people, components, renorm=renorm)
        if centroid is None:
            return None
        centroid_alt, _ = compute_centroid(pids, cohort_people, components,
                                            renorm=other_renorm(renorm))
        enrich = bootstrap_enrichment(centroid, get_null(n_used), cohort_mean, components)
        disc = discrete_enrichment(pids, cohort_people)
        return {"allele_id": allele_id, "gene": gene_bare_name, "kind": kind, "n_carriers": n_used,
                "centroid": centroid, "centroid_alt": centroid_alt,
                "confidence_tier": tier, **enrich, **disc}

    genes = set(genes)
    t = table1[table1["gene_bare"].isin(genes) & table1["person_id"].isin(cohort_ids)].copy()

    records = []

    # ---- known alleles, every gene at once ----
    t["allele_field"] = t["consensus"].map(lambda a: vc.to_nfield(a, field))
    known_sub = t.dropna(subset=["allele_field"])
    for (gene_bare_name, allele), grp in known_sub.groupby(["gene_bare", "allele_field"]):
        r = make_record(f"{gene_bare_name}*{allele}", gene_bare_name, "known", set(grp["person_id"]))
        if r:
            records.append(r)

    # ---- novel alleles, every gene, ONE cds.fa.gz open per (person, hap) across all their genes ----
    novel_t = t[t["is_novel"]]
    if not novel_t.empty:
        group_sizes = t.groupby(["person_id", "hap", "contig", "gene"]).size()
        clusters = defaultdict(set)  # (gene_prefixed, sha1) -> set(person_id)
        for (person_id, hap), sub in novel_t.groupby(["person_id", "hap"]):
            cds_path = os.path.join(outroot, person_id, "immuannot_output", hap, "cds.fa.gz")
            fasta_index = None
            for _, row in sub.iterrows():
                key = (row["person_id"], row["hap"], row["contig"], row["gene"])
                if group_sizes.loc[key] > 1:
                    continue
                if fasta_index is None:
                    fasta_index = _parse_cds_fasta(cds_path)
                rest = f"{row['contig']}_{row['gene']}"
                candidates = fasta_index.get(rest, [])
                if len(candidates) != 1:
                    continue
                _i, seq = candidates[0]
                if not seq:
                    continue
                sha1 = hashlib.sha1(seq.encode("ascii")).hexdigest()
                clusters[(row["gene"], sha1)].add(person_id)

        table3_by_id = table3.set_index("novel_id") if "novel_id" in table3.columns else None
        for (gene_prefixed, sha1), pids in clusters.items():
            gene_bare_name = vc.gene_bare(gene_prefixed)
            nid = f"{gene_prefixed}_nov_{sha1[:8]}"
            row = (table3_by_id.loc[nid]
                   if table3_by_id is not None and nid in table3_by_id.index else {})
            r = make_record(nid, gene_bare_name, "novel", pids,
                             tier=row.get("confidence_tier"))
            if r:
                records.append(r)

    return records, cohort_mean


# ---------------------------------------------------------------------------
# Plotting -- ONE encoding for confidence (alpha ~ carrier count, fixed marker size), not both --
# encoding the same variable in two visual channels at once made it unclear which was which and
# didn't add information. Never a hard exclusion of low-n points. Known alleles are circles, novel
# (this-cohort) alleles are triangles (an earlier bordered-circle design wasn't legible enough at
# real cohort density).
# ---------------------------------------------------------------------------
MARKER_SIZE = 36


def _alpha_scale(records):
    ns = np.array([r["n_carriers"] for r in records], dtype=float)
    scale = np.sqrt(ns / max(ns.max(), 1)) if len(ns) else ns
    return np.clip(0.15 + 0.7 * scale, 0.15, 0.85)


def _dominant_color(centroid, components):
    dominant = max(components, key=lambda c: centroid[c])
    return vc.ANCESTRY_COLORS.get(dominant, vc.MISSING_COLOR)


def _legend_handles():
    from matplotlib.lines import Line2D
    return [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="gray", markersize=8,
               label="known allele"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="gray", markersize=9,
               label="novel allele (this cohort)"),
    ]


def _label_top_known(ax, records, xy, top_n, fontsize=7):
    """Text-label the top_n known alleles by carrier count only -- skips novel alleles on purpose:
    their IDs are opaque content-hash strings (e.g. 'HLA-B_nov_27a0cb81'), not human-readable, and
    there can be hundreds of them, so labeling them would just add clutter with no payoff. Labeling
    every point (known or novel) would overlap badly at this density; capping to the top_n most
    common known alleles keeps labels legible since well-powered alleles tend to be more spread out
    than the dense cloud of rare/novel points near the vertices."""
    if top_n <= 0:
        return
    known_idx = [i for i, r in enumerate(records) if r["kind"] == "known"]
    known_idx.sort(key=lambda i: records[i]["n_carriers"], reverse=True)
    for i in known_idx[:top_n]:
        x, y = xy[i]
        ax.annotate(records[i]["allele_id"], (x, y), fontsize=fontsize, xytext=(3, 3),
                    textcoords="offset points")


def _scatter_one(ax, xy_or_xyz, color, is_novel, alpha, size=MARKER_SIZE):
    ax.scatter(*xy_or_xyz, s=size, alpha=alpha, color=color, marker="^" if is_novel else "o",
               linewidths=0)


def plot_ternary_alleles(records, components, out_path, label_top_known=12):
    verts = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])

    fig, ax = plt.subplots(figsize=(7, 6.5))
    ax.add_patch(plt.Polygon(verts, fill=False, edgecolor="black", linewidth=1))
    for (vx, vy), name in zip(verts, components):
        ax.text(vx, vy - 0.05 if vy == 0 else vy + 0.03, name, ha="center", fontsize=10,
                fontweight="bold")

    alphas = _alpha_scale(records)
    xy_all = []
    for r, alpha in zip(records, alphas):
        xy = np.array([r["centroid"][c] for c in components]) @ verts
        xy_all.append(xy)
        color = _dominant_color(r["centroid"], components)
        _scatter_one(ax, xy, color, r["kind"] == "novel", alpha)
    _label_top_known(ax, records, xy_all, label_top_known)

    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.15, 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(handles=_legend_handles(), fontsize=8, frameon=False, loc="upper right")
    ax.set_title(f"Per-allele ancestry centroid ({'/'.join(components)} simplex)\n"
                 "alpha ~ carrier count (no minimum-carrier-count floor)", fontsize=10)
    fig.tight_layout()
    vc.savefig(fig, out_path)


def plot_tetrahedron_alleles(records, components, out_path):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    for i, j in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
        xs, ys, zs = zip(TETRA_VERTS[i], TETRA_VERTS[j])
        ax.plot(xs, ys, zs, color="black", linewidth=0.8)
    for v, name in zip(TETRA_VERTS, components):
        ax.text(v[0], v[1], v[2], name, fontsize=10, fontweight="bold")

    alphas = _alpha_scale(records)
    for r, alpha in zip(records, alphas):
        xyz = np.array([r["centroid"][c] for c in components]) @ TETRA_VERTS
        color = _dominant_color(r["centroid"], components)
        _scatter_one(ax, xyz, color, r["kind"] == "novel", alpha)

    ax.set_axis_off()
    ax.legend(handles=_legend_handles(), fontsize=8, frameon=False, loc="upper right")
    ax.set_title(f"Per-allele ancestry centroid ({'/'.join(components)} tetrahedron)", fontsize=10)
    fig.tight_layout()
    vc.savefig(fig, out_path)


def plot_pca_centroids(records, out_path, label_top_known=12):
    X = np.array([[r["centroid"][c] for c in vc.ANCESTRY_ORDER] for r in records])
    if len(X) < 3:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "Too few alleles for a PCA of centroids", ha="center", va="center")
        ax.axis("off")
        vc.savefig(fig, out_path)
        return
    Xs = vc.standardize(X)
    scores, var_ratio = vc.pca_fit_transform(Xs, n_components=2)

    alphas = _alpha_scale(records)
    fig, ax = plt.subplots(figsize=(7, 6))
    for r, (x, y), alpha in zip(records, scores, alphas):
        color = _dominant_color(r["centroid"], vc.ANCESTRY_ORDER)
        _scatter_one(ax, (x, y), color, r["kind"] == "novel", alpha)
    _label_top_known(ax, records, scores, label_top_known)
    ax.set_xlabel(f"PC1 ({100 * var_ratio[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * var_ratio[1]:.1f}%)" if len(var_ratio) > 1 else "PC2")
    ax.legend(handles=_legend_handles(), fontsize=8, frameon=False, loc="upper right")
    ax.set_title("PCA of per-allele 6-way ancestry centroids\n(the arbitrary-dimension "
                 "generalization of the ternary/tetrahedron simplex figures)", fontsize=10)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Multi-gene plots: (1) every allele across a gene list in one ternary, colored by GENE identity
# instead of ancestry, to see where diversity concentrates across genes; (2) a small-multiples grid
# of one ancestry-colored ternary per gene (same style as the single-gene plot above), for
# comparing genes side by side. Both take the flat record list from build_records_for_genes().
# ---------------------------------------------------------------------------
def _gene_palette(genes):
    genes = sorted(genes)
    cmap = plt.get_cmap("tab10" if len(genes) <= 10 else "tab20")
    return {g: cmap(i % cmap.N) for i, g in enumerate(genes)}


def plot_gene_colored_ternary(records, components, out_path, point_size=10, alpha_max=0.5,
                               alpha_min=0.08):
    """The 'mega' plot: every allele from every gene in `records`, one point per allele, colored by
    gene rather than ancestry. Deliberately smaller points and a lower alpha ceiling than the
    single-gene ternary -- this is far denser (every gene's full allele set at once)."""
    verts = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    genes = sorted({r["gene"] for r in records})
    palette = _gene_palette(genes)

    fig, ax = plt.subplots(figsize=(9, 8.5))
    ax.add_patch(plt.Polygon(verts, fill=False, edgecolor="black", linewidth=1))
    for (vx, vy), name in zip(verts, components):
        ax.text(vx, vy - 0.05 if vy == 0 else vy + 0.03, name, ha="center", fontsize=11,
                fontweight="bold")

    raw_alphas = _alpha_scale(records)
    # Rescale the usual 0.15-0.85 alpha range down to alpha_min-alpha_max for this much denser
    # plot -- same relative confidence signal, just capped so overlaps read as density, not a
    # solid mass of color.
    alphas = alpha_min + (raw_alphas - 0.15) / (0.85 - 0.15) * (alpha_max - alpha_min)
    for r, alpha in zip(records, alphas):
        xy = np.array([r["centroid"][c] for c in components]) @ verts
        marker = "^" if r["kind"] == "novel" else "o"
        ax.scatter(*xy, s=point_size, alpha=alpha, color=palette[r["gene"]], marker=marker,
                   linewidths=0)

    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.15, 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[g],
                       markersize=8, label=vc.gene_display(g)) for g in genes]
    handles.append(Line2D([0], [0], marker="^", color="none", markerfacecolor="gray",
                           markersize=8, label="novel (any gene)"))
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper right", ncol=2)
    ax.set_title(f"All alleles across {len(genes)} genes, colored by gene "
                 f"({'/'.join(components)} simplex)\nalpha ~ carrier count (capped for density); "
                 "triangle = novel allele (this cohort)", fontsize=10)
    fig.tight_layout()
    vc.savefig(fig, out_path)


def plot_gene_grid_ternary(records, components, out_path, grid_genes, label_top_known=6,
                            ncols=3, point_size=None):
    """Small multiples: one ancestry-colored ternary per gene in `grid_genes`, same visual style
    (alpha ~ carrier count, circle=known/triangle=novel, top-N known alleles labeled) as the
    single-gene plot, for comparing genes side by side."""
    verts = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    n = len(grid_genes)
    nrows = -(-n // ncols)  # ceil division
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.6 * nrows), squeeze=False)
    axes = axes.flatten()
    size = point_size or (MARKER_SIZE * 0.6)

    for ax, gene in zip(axes, grid_genes):
        gene_records = [r for r in records if r["gene"] == gene]
        ax.add_patch(plt.Polygon(verts, fill=False, edgecolor="black", linewidth=1))
        for (vx, vy), name in zip(verts, components):
            ax.text(vx, vy - 0.05 if vy == 0 else vy + 0.03, name, ha="center", fontsize=8,
                    fontweight="bold")
        if gene_records:
            alphas = _alpha_scale(gene_records)
            xy_all = []
            for r, alpha in zip(gene_records, alphas):
                xy = np.array([r["centroid"][c] for c in components]) @ verts
                xy_all.append(xy)
                color = _dominant_color(r["centroid"], components)
                _scatter_one(ax, xy, color, r["kind"] == "novel", alpha, size=size)
            _label_top_known(ax, gene_records, xy_all, label_top_known, fontsize=6)
        ax.set_xlim(-0.15, 1.15)
        ax.set_ylim(-0.15, 1.0)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{vc.gene_display(gene)} (n={len(gene_records)} alleles)", fontsize=10)

    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle("Per-gene ancestry-centroid geometry (small multiples, ancestry-colored)",
                  fontsize=12)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _fmt_centroid(centroid, components):
    return ", ".join(f"{c}={centroid[c]:.3f}" for c in components)


def _fmt(x, fmt=".3f"):
    return "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else format(x, fmt)


def write_allele_report(path, gene, cohort_label, renorm, records_ternary, ternary_components,
                         records_tetra, tetra_components):
    lines = [f"# Per-allele ancestry-centroid geometry -- gene HLA-{gene}, cohort `{cohort_label}`\n",
              f"Ternary components: {'/'.join(ternary_components)}. Tetrahedron components: "
              f"{'/'.join(tetra_components)}. Primary centroid renorm: `{renorm}` (both variants "
              "always computed; the alt column shows the other one).\n",
              "**No minimum-carrier-count floor is applied to centroid computation** -- an allele "
              "with 1 carrier is shown with a centroid equal to that person's own admixture "
              "proportions. n_carriers plus marker alpha communicate confidence visually "
              "instead of a hard cutoff. Bootstrap enrichment (continuous) and Fisher/chi2 "
              "enrichment (discrete, on ancestry_pred labels) are two independent significance "
              "reads, reported side by side -- neither is used to decide what gets plotted.\n"]

    lines.append("\n## Ternary-space alleles (sorted by bootstrap enrichment p-value)\n")
    lines.append("| allele | kind | n_carriers | centroid | centroid (alt renorm) | bootstrap z "
                 "| bootstrap p | chi2 p | confidence_tier |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for r in sorted(records_ternary, key=lambda r: r.get("p_value") if r.get("p_value") == r.get("p_value") else 1.0):
        lines.append(
            f"| {r['allele_id']} | {r['kind']} | {r['n_carriers']} | "
            f"{_fmt_centroid(r['centroid'], ternary_components)} | "
            f"{_fmt_centroid(r['centroid_alt'], ternary_components)} | "
            f"{_fmt(r.get('z_score'))} | {_fmt(r.get('p_value'))} | {_fmt(r.get('chi2_p'))} | "
            f"{r.get('confidence_tier') or ''} |")

    lines.append("\n## Tetrahedron-space alleles\n")
    lines.append("| allele | kind | n_carriers | centroid | bootstrap p | chi2 p |")
    lines.append("|---|---|---|---|---|---|")
    for r in sorted(records_tetra, key=lambda r: r.get("p_value") if r.get("p_value") == r.get("p_value") else 1.0):
        lines.append(
            f"| {r['allele_id']} | {r['kind']} | {r['n_carriers']} | "
            f"{_fmt_centroid(r['centroid'], tetra_components)} | {_fmt(r.get('p_value'))} | "
            f"{_fmt(r.get('chi2_p'))} |")

    vc.write_report(path, lines)


# ---------------------------------------------------------------------------
# I/O for Table 3 (read as a data file only -- see novel_allele_carriers() docstring)
# ---------------------------------------------------------------------------
def load_novel_table(novel_alleles_path):
    if not os.path.exists(novel_alleles_path):
        sys.exit(f"FATAL: novel_alleles.tsv (Table 3) not found at {novel_alleles_path!r}. "
                 f"Run 03_novel_alleles.py first.")
    return pd.read_csv(novel_alleles_path, sep="\t")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run_single_gene(args, table1, table3, cohort_people, cohort_label, people_root, out_dir):
    gene = vc.gene_bare(args.gene)
    ternary_components = args.ternary_ancestries.split(",")
    tetra_components = args.tetrahedron_ancestries.split(",")

    print(f"Building ternary-space allele records ({'/'.join(ternary_components)}) ...",
          file=sys.stderr)
    records_ternary, _ = build_allele_records(
        table1, table3, cohort_people, gene, ternary_components, people_root, field=args.field,
        renorm=args.renorm, n_boot=args.n_boot, seed=args.seed)
    plot_ternary_alleles(records_ternary, ternary_components,
                          os.path.join(out_dir, f"ternary_{gene}.png"),
                          label_top_known=args.label_top_known)

    print(f"Building tetrahedron-space allele records ({'/'.join(tetra_components)}) ...",
          file=sys.stderr)
    records_tetra, _ = build_allele_records(
        table1, table3, cohort_people, gene, tetra_components, people_root, field=args.field,
        renorm=args.renorm, n_boot=args.n_boot, seed=args.seed)
    plot_tetrahedron_alleles(records_tetra, tetra_components,
                             os.path.join(out_dir, f"tetrahedron_{gene}.png"))

    print("Building full 6-way allele records for PCA ...", file=sys.stderr)
    records_full, _ = build_allele_records(
        table1, table3, cohort_people, gene, vc.ANCESTRY_ORDER, people_root, field=args.field,
        renorm=args.renorm, n_boot=args.n_boot, seed=args.seed)
    plot_pca_centroids(records_full, os.path.join(out_dir, f"pca_centroids_{gene}.png"),
                        label_top_known=args.label_top_known)

    write_allele_report(os.path.join(out_dir, "allele_ancestry_report.md"), gene, cohort_label,
                        args.renorm, records_ternary, ternary_components, records_tetra,
                        tetra_components)

    print(f"Wrote figures + report to {out_dir!r} "
          f"({len(records_ternary)} ternary alleles, {len(records_tetra)} tetrahedron alleles, "
          f"{len(records_full)} PCA alleles).", file=sys.stderr)


def run_multi_gene(args, table1, table3, cohort_people, people_root, out_dir):
    genes = [vc.gene_bare(g.strip()) for g in args.multi_genes.split(",") if g.strip()]
    grid_genes = ([vc.gene_bare(g.strip()) for g in args.grid_genes.split(",") if g.strip()]
                  if args.grid_genes else genes[:6])
    missing_grid = [g for g in grid_genes if g not in genes]
    if missing_grid:
        sys.exit(f"FATAL: --grid-genes {missing_grid} not present in --multi-genes {genes} -- "
                 f"the grid is drawn from the same single data pass as the combined plot, so grid "
                 f"genes must be a subset.")
    ternary_components = args.ternary_ancestries.split(",")

    print(f"Building multi-gene allele records for {len(genes)} genes "
          f"({'/'.join(ternary_components)}) -- single pass over the raw data ...", file=sys.stderr)
    records, _ = build_records_for_genes(
        table1, table3, cohort_people, genes, ternary_components, people_root, field=args.field,
        renorm=args.renorm, n_boot=args.n_boot, seed=args.seed)
    n_known = sum(1 for r in records if r["kind"] == "known")
    n_novel = sum(1 for r in records if r["kind"] == "novel")
    print(f"  {len(records)} total alleles ({n_known} known, {n_novel} novel) across {len(genes)} "
          f"genes.", file=sys.stderr)

    mega_path = os.path.join(out_dir, "ternary_all_genes.png")
    plot_gene_colored_ternary(records, ternary_components, mega_path,
                               point_size=args.mega_point_size)

    grid_path = os.path.join(out_dir, "ternary_gene_grid.png")
    plot_gene_grid_ternary(records, ternary_components, grid_path, grid_genes,
                            label_top_known=max(args.label_top_known // 2, 1))

    print(f"Wrote {mega_path!r} (all {len(genes)} genes, colored by gene) and {grid_path!r} "
          f"({len(grid_genes)}-gene small-multiples grid, ancestry-colored).", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    vc.add_common_args(ap)
    ap.add_argument("--mode", choices=["single", "multi"], default="single",
                     help="'single': the original one-gene ternary/tetrahedron/PCA figures "
                          "(--gene). 'multi': a combined all-genes ternary colored by gene "
                          "(--multi-genes) plus a small-multiples ancestry-colored grid "
                          "(--grid-genes), both from one shared data pass.")
    ap.add_argument("--gene", default="B", help="[--mode single] Bare gene name (default HLA-B).")
    ap.add_argument("--multi-genes", default=",".join(vc.CLASSICAL_GENES_BARE),
                     help="[--mode multi] Comma-separated gene list for the combined ternary "
                          "colored by gene identity. Default: all 8 classical genes.")
    ap.add_argument("--grid-genes", default=None,
                     help="[--mode multi] Comma-separated subset of --multi-genes to render as "
                          "the ancestry-colored small-multiples grid. Default: first 6 of "
                          "--multi-genes. Must be a subset (drawn from the same data pass).")
    ap.add_argument("--mega-point-size", type=int, default=10,
                     help="[--mode multi] Marker size (points^2) for the combined all-genes "
                          "ternary -- kept small since it is far denser than the single-gene plot.")
    ap.add_argument("--field", type=int, choices=[2, 4], default=2,
                     help="Allele resolution for known-allele grouping (default 2-field).")
    ap.add_argument("--ternary-ancestries", default=",".join(DEFAULT_TERNARY),
                     help="Comma-separated triple of ancestries for the ternary figure(s).")
    ap.add_argument("--tetrahedron-ancestries", default=",".join(DEFAULT_TETRA),
                     help="[--mode single] Comma-separated quadruple of ancestries for the "
                          "tetrahedron figure.")
    ap.add_argument("--renorm", choices=["per_person", "per_centroid"], default="per_person",
                     help="Primary centroid math: renormalize each carrier's subset proportions "
                          "before averaging (per_person, matches 06_figures_structure.py's "
                          "per-person ternary) vs. average full proportions first and renormalize "
                          "once (per_centroid). Both are always computed; this picks which is the "
                          "'main' column in the report and which drives the plotted position.")
    ap.add_argument("--label-top-known", type=int, default=12,
                     help="Text-label this many of the most common KNOWN alleles by carrier count "
                          "on the ternary/PCA figures (0 disables; halved per-gene in --mode multi's "
                          "grid). Novel alleles are never labeled -- their IDs are opaque "
                          "content-hash strings, not human-readable -- and labeling every point "
                          "would overlap badly at real cohort density.")
    ap.add_argument("--n-boot", type=int, default=1000,
                     help="Bootstrap resamples for the continuous enrichment null distribution.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--novel-alleles", default=None,
                     help="Override path to novel_alleles.tsv (Table 3, read as a data file only "
                          "-- see module docstring). Default: <outroot>/novel_alleles.tsv (this "
                          "script's --outroot is the AGGREGATE-file root, per _viz_common.py's "
                          "convention -- same as 05/06/07 -- not the person-directory root).")
    ap.add_argument("--people-root", default=None,
                     help="Root holding <person_id>/immuannot_output/hap{1,2}/cds.fa.gz, needed "
                          "only for the independent novel-allele carrier recomputation (see module "
                          "docstring). Default: <outroot>/people (RUNBOOK.md 'Step 1c' convention "
                          "-- 03_novel_alleles.py's OWN --outroot default, confusingly a different "
                          "root than this script's --outroot; not the same flag).")
    args = ap.parse_args()

    table1_path, _table2_path, cohort_path, _sr_path = vc.resolve_paths(args)
    novel_path = args.novel_alleles or os.path.join(args.outroot, "novel_alleles.tsv")
    people_root = args.people_root or os.path.join(args.outroot, "people")

    table1 = vc.load_table1(table1_path)
    cohort_df = vc.load_cohort_membership(cohort_path)
    cohort_people, cohort_label = vc.select_cohort_people(cohort_df, args.cohort, args.td_max)
    print(f"Cohort: {cohort_label} -- {len(cohort_people)} people.", file=sys.stderr)
    table3 = load_novel_table(novel_path)

    out_dir = args.out_dir or vc.default_out_dir("10_allele_ancestry_geometry", args.cohort,
                                                  args.td_max)
    vc.ensure_dir(out_dir)

    if args.mode == "multi":
        run_multi_gene(args, table1, table3, cohort_people, people_root, out_dir)
    else:
        run_single_gene(args, table1, table3, cohort_people, cohort_label, people_root, out_dir)


if __name__ == "__main__":
    main()
