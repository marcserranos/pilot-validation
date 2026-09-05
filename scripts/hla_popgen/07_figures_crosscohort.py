#!/usr/bin/env python3
"""Cross-cohort figures -- the ones VIZ_LIT.md Part 3 flags as this project's genuine paper
differentiators, because essentially no published HLA paper has three cohorts of overlapping
participants at different read-length/resolution technologies plus a template-distance accuracy
axis and true haplotype phase:

  1. Short-read vs long-read allele-frequency scatter, 2-field, log-log (VIZ_LIT.md 3.1) --
     identifies systematic SR miscalling/reference bias. Ranked #2 in VIZ_LIT.md Part 5.
  2. template_distance distribution by ancestry, plus the td-max SWEEP of "% exact (td==0) by
     ancestry" (VIZ_LIT.md 3.2) -- the headline reference-bias evidence figure. Ranked #5.
  3. Resolution cascade: 2-field -> 4-field fan-out by ancestry, summarized as conditional
     Shannon entropy H(4field | 2field) per gene per ancestry (VIZ_LIT.md 3.4). Ranked #6.
  4. Non-classical gene diversity map (MICA/MICB/TAP1/TAP2/DRB3-4-5/DQA2/DQB2/HLA-E/F/G), plus a
     gene x cohort coverage-availability matrix making the LR-only coverage advantage explicit
     (VIZ_LIT.md 3.5).
  5. Cis heterodimer (DQA1~DQB1, DPA1~DPB1) frequency by ancestry from Table 2, plus the
     pairable-fraction-by-ancestry number -- assembly fragmentation is itself a result
     (SCHEMA.md Table 2, VIZ_LIT.md 1.8/3.3). Pairable fraction is RE-DERIVED here directly from
     Table 1's (person_id, hap, contig, gene) grain (never modifies 01_extract_rich.py, which
     only prints the cohort-wide unpairable count to stderr and does not persist it per ancestry).

## On the three-cohort abstraction here (design note, also surfaced in the handoff report)

Every OTHER figure in this sub-project is reproducible against exactly one of {sr, lr, lr_td}.
These figures are inherently comparisons BETWEEN cohorts, so `--cohort` here selects only the
long-read SIDE of the comparison (`lr` or `lr_td<N>`); the short-read side is always the full `sr`
cohort (there is no "compare sr against itself" figure). `--cohort sr` is refused with a pointer
to this docstring rather than silently doing something arbitrary.

Usage (fixtures): see 05_figures_frequency.py's docstring for the fixture pipeline; then e.g.
    python3 scripts/hla_popgen/07_figures_crosscohort.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --table2 /tmp/hla_fixtures/hla_cis_pairs.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --sr-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --cohort lr

Real run (VM): python3 scripts/hla_popgen/07_figures_crosscohort.py --cohort lr
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

CLASSICAL = vc.CLASSICAL_GENES_BARE
NONCLASSICAL = vc.NONCLASSICAL_GENES
PAIRS = [("DQA1", "DQB1"), ("DPA1", "DPB1")]


def build_long_frame(table1, cohort_people, genes, resolution=2):
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)].copy()
    sub["allele_n"] = sub["consensus"].map(lambda a: vc.to_nfield(a, resolution))
    sub = sub.dropna(subset=["allele_n"])
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
    sub["allele_label"] = sub["gene_bare"] + "*" + sub["allele_n"]
    return sub


def build_sr_long_frame(sr_long, cohort_people, genes):
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    sub = sr_long[sr_long["person_id"].isin(person_ids) & sr_long["gene_bare"].isin(genes)].copy()
    sub["allele_n"] = sub["allele"].map(lambda a: vc.to_nfield(a, 2))
    sub = sub.dropna(subset=["allele_n"])
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
    sub["allele_label"] = sub["gene_bare"] + "*" + sub["allele_n"]
    return sub


def pooled_freq(long_df, genes):
    """{(gene, allele_label): (freq, n_copies)} pooled across ancestry."""
    out = {}
    for gene in genes:
        g = long_df[long_df["gene_bare"] == gene]
        total = len(g)
        if total == 0:
            continue
        counts = g["allele_label"].value_counts()
        for allele, n in counts.items():
            out[allele] = (n / total, n)
    return out


def freq_by_ancestry(long_df, genes):
    """{(ancestry, allele_label): (freq, n_copies)} -- same per-gene normalization as pooled_freq
    but stratified by ancestry. A pooled (cohort-wide) SR-vs-LR comparison blurs exactly the
    signal VIZ_LIT.md 3.1 is about: reference bias is fundamentally an ANCESTRY-dependent claim
    (the IPD reference itself is historically Euro-centric), not a whole-cohort average -- this is
    the stratification needed to actually test that claim rather than just eyeball a scatter."""
    out = {}
    for anc in vc.ANCESTRY_ORDER:
        sub_anc = long_df[long_df["ancestry_pred"] == anc]
        for gene in genes:
            g = sub_anc[sub_anc["gene_bare"] == gene]
            total = len(g)
            if total == 0:
                continue
            counts = g["allele_label"].value_counts()
            for allele, n in counts.items():
                out[(anc, allele)] = (n / total, n)
    return out


def compute_sr_lr_disagreement_by_ancestry(sr_freq_anc, lr_freq_anc, n_people_lr, n_people_sr,
                                            min_n=5, eur_common_thresh=0.15):
    """Per-ancestry SR-vs-LR agreement summary -- the headline scientific claim of this figure
    made quantitative and testable, not just visual: (a) mean absolute frequency difference per
    ancestry (expect EUR lowest -- SR calling is tuned/validated against a Euro-centric reference)
    and (b) a DIRECTIONAL bias check restricted to alleles common in the LR-derived EUR population
    (LR EUR frequency > eur_common_thresh): mean(sr_freq - lr_freq) for those alleles, computed
    within EACH ancestry. A positive value in a non-EUR ancestry means SR systematically
    OVER-calls EUR-common alleles there -- the specific reference-bias mechanism VIZ_LIT.md 3.1
    and the Mapping Bias literature (PMC4426377) describe, not just "SR and LR disagree more.\"

    n_people_lr/n_people_sr: dict {ancestry: n} -- carried through into the returned frame purely
    so the PLOT can mark thin ancestry groups (MID/SAS are genuinely the smallest AoU ancestry
    groups, in both the fixtures and real data -- a bug found in review: with n=6-8 people, MID/
    SAS bars dominated this figure visually while being the least powered estimates in it, exactly
    inverting how a reader should weight them). `n_alleles`/`n_eur_common_alleles` alone do NOT
    capture this -- they count alleles that individually cleared `min_n` COMBINED sr+lr copies,
    which says nothing about how many people/ancestry-copies underlie the AVERAGE across those
    alleles for a small ancestry group."""
    eur_common = {a for (anc, a), (f, n) in lr_freq_anc.items() if anc == "EUR" and f > eur_common_thresh}
    rows = []
    for anc in vc.ANCESTRY_ORDER:
        alleles = ({a for (a2, a) in sr_freq_anc if a2 == anc} |
                   {a for (a2, a) in lr_freq_anc if a2 == anc})
        diffs, eur_common_diffs = [], []
        for a in alleles:
            fx, nx = sr_freq_anc.get((anc, a), (0.0, 0))
            fy, ny = lr_freq_anc.get((anc, a), (0.0, 0))
            if nx + ny < min_n:
                continue
            diffs.append(abs(fx - fy))
            if a in eur_common:
                eur_common_diffs.append(fx - fy)
        rows.append(dict(
            ancestry=anc,
            mean_abs_diff=float(np.mean(diffs)) if diffs else np.nan,
            n_alleles=len(diffs),
            eur_common_allele_bias=float(np.mean(eur_common_diffs)) if eur_common_diffs else np.nan,
            n_eur_common_alleles=len(eur_common_diffs),
            n_people_lr=int(n_people_lr.get(anc, 0)),
            n_people_sr=int(n_people_sr.get(anc, 0)),
        ))
    return pd.DataFrame(rows)


def _thin_n_bar_panel(ax, disagree_df, value_col, min_cell_n, ylabel, title, zero_line=False):
    """Shared thin-N treatment for a single-bar-per-ancestry panel: bars below `min_cell_n`
    (using min(n_people_lr, n_people_sr) -- the more limited of the two cohorts' contributing
    ancestry counts, since the disagreement estimate needs BOTH to be well-powered) are greyed out
    and hatched, matching the heatmap convention in 05_figures_frequency.py's plot_heatmap_gene
    and 07's own non-classical diversity map, and every bar is annotated with its underlying n
    (lr/sr person counts) regardless of thin status -- so a reader never has to guess.

    Bug found in review: MID/SAS are the smallest ancestry groups in this project's cohorts (both
    the fixtures and real AoU data) and were previously rendered as ordinary, full-color,
    equal-weight bars -- with n=6-8 people they had the LARGEST-magnitude bars in both panels
    purely from sampling noise, visually dominating the exact claim (EUR agrees best, on n=95-110)
    the figure exists to make. This is the same class of bug the thin-N cell guard elsewhere in
    this project's heatmaps already exists to prevent, just never applied to a single-bar-per-
    ancestry chart before now."""
    n_min = np.minimum(disagree_df["n_people_lr"].values, disagree_df["n_people_sr"].values)
    thin = n_min < min_cell_n
    values = disagree_df[value_col].values
    colors = [("#DDDDDD" if t else vc.ANCESTRY_COLORS[a])
              for a, t in zip(disagree_df["ancestry"], thin)]
    edgecolors = [vc.ANCESTRY_COLORS[a] for a in disagree_df["ancestry"]]

    if zero_line:
        ax.axhline(0, color="grey", linewidth=1)
    bars = ax.bar(disagree_df["ancestry"], np.nan_to_num(values), color=colors,
                  edgecolor=edgecolors, linewidth=0.9)
    finite = values[np.isfinite(values)]
    span = (float(np.max(np.abs(finite))) if len(finite) else 1.0) or 1.0
    for bar, t, v, nl, ns in zip(bars, thin, values, disagree_df["n_people_lr"],
                                  disagree_df["n_people_sr"]):
        if t:
            bar.set_hatch("////")
        vv = 0.0 if not np.isfinite(v) else v
        va = "bottom" if vv >= 0 else "top"
        y = vv + 0.03 * span if vv >= 0 else vv - 0.03 * span
        ax.text(bar.get_x() + bar.get_width() / 2, y, f"lr n={int(nl)}\nsr n={int(ns)}",
                ha="center", va=va, fontsize=6.5)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    # Headroom so the n-count text never clips against the axes frame.
    ymin, ymax = ax.get_ylim()
    pad = 0.18 * (ymax - ymin if ymax > ymin else 1)
    ax.set_ylim(ymin - pad, ymax + pad)


def plot_sr_lr_disagreement_by_ancestry(disagree_df, min_cell_n, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    ax1, ax2 = axes
    _thin_n_bar_panel(ax1, disagree_df, "mean_abs_diff", min_cell_n,
                       "Mean |sr freq - lr freq|",
                       "SR-LR disagreement by ancestry\n(lower = better agreement)")
    _thin_n_bar_panel(ax2, disagree_df, "eur_common_allele_bias", min_cell_n,
                       "Mean (sr freq - lr freq), EUR-common alleles",
                       "Directional bias: does SR over-call EUR-common\nalleles in non-EUR ancestries?",
                       zero_line=True)
    fig.suptitle(f"Grey/hatched bars: fewer than {min_cell_n} people in the more limited of the "
                 f"sr/lr cohorts for that ancestry -- noise, not signal", fontsize=8.5, y=1.02)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 1. SR vs LR frequency scatter (VIZ_LIT.md 3.1)
# ---------------------------------------------------------------------------
def plot_sr_vs_lr_scatter(sr_freq, lr_freq, genes, out_path):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    cmap = plt.get_cmap("tab10")
    outliers = []
    for i, gene in enumerate(genes):
        alleles = {a for a in sr_freq if a.startswith(f"{gene}*")} | \
                  {a for a in lr_freq if a.startswith(f"{gene}*")}
        xs, ys, ns = [], [], []
        for a in alleles:
            fx, nx = sr_freq.get(a, (0.0, 0))
            fy, ny = lr_freq.get(a, (0.0, 0))
            if nx == 0 and ny == 0:
                continue
            xs.append(max(fx, 1e-4))
            ys.append(max(fy, 1e-4))
            ns.append(nx + ny)
            if fx > 0 and fy > 0 and (fx / fy > 2 or fy / fx > 2) and (nx + ny) >= 10:
                outliers.append((a, fx, fy))
        if not xs:
            continue
        ax.scatter(xs, ys, s=[8 + 0.4 * n for n in ns], alpha=0.6, color=cmap(i % 10),
                   label=gene, edgecolors="none")
    lims = [1e-4, 1.0]
    ax.plot(lims, lims, "--", color="grey", linewidth=1, label="y=x (no bias)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Allele frequency -- sr (AoU-native short-read)")
    ax.set_ylabel("Allele frequency -- lr (Immuannot long-read)")
    ax.set_title("Same-allele frequency, short-read vs long-read (2-field)\n"
                 "points off the diagonal = candidate SR miscalling / reference bias", fontsize=10)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)
    return sorted(outliers, key=lambda t: -max(t[1], t[2]))[:15]


# ---------------------------------------------------------------------------
# 2. template_distance by ancestry + td-max sweep (VIZ_LIT.md 3.2)
# ---------------------------------------------------------------------------
def plot_td_by_ancestry(table1, cohort_people, genes, min_cell_n, out_path):
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)].copy()
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER) & sub["template_distance"].notna()]
    # Thin-N is a PEOPLE count, not a row count: a violin built from few people but many
    # haplotype-call rows each (up to 16 rows/person across 8 classical genes x 2 haps) can look
    # just as "full" as one built from many people -- the row count alone would understate how
    # thin the underlying independent sample really is.
    n_people_by_anc = {anc: cohort_people.loc[cohort_people["ancestry_pred"] == anc,
                                               "person_id"].nunique()
                        for anc in vc.ANCESTRY_ORDER}

    fig, ax = plt.subplots(figsize=(9, 5.5))
    data, positions, colors, thin_flags = [], [], [], []
    for i, anc in enumerate(vc.ANCESTRY_ORDER):
        vals = sub.loc[sub["ancestry_pred"] == anc, "template_distance"].values
        if len(vals) == 0:
            continue
        data.append(vals)
        positions.append(i)
        thin = n_people_by_anc.get(anc, 0) < min_cell_n
        thin_flags.append(thin)
        colors.append("#BBBBBB" if thin else vc.ANCESTRY_COLORS[anc])
    if data:
        parts = ax.violinplot(data, positions=positions, showmedians=True, widths=0.8)
        for pc, c, thin in zip(parts["bodies"], colors, thin_flags):
            pc.set_facecolor(c)
            pc.set_alpha(0.35 if thin else 0.6)
            if thin:
                pc.set_hatch("////")
    ax.set_xticks(range(len(vc.ANCESTRY_ORDER)))
    ax.set_xticklabels([f"{anc}\n(n_people={n_people_by_anc.get(anc, 0)})"
                         for anc in vc.ANCESTRY_ORDER], fontsize=8)
    ax.set_ylabel("template_distance (edit distance to nearest IPD-IMGT/HLA reference allele)")
    ns = {anc: int((sub["ancestry_pred"] == anc).sum()) for anc in vc.ANCESTRY_ORDER}
    ax.set_title("template_distance by ancestry, classical genes pooled\n"
                 f"n (haplotype-call rows) = {ns}\n"
                 f"higher = more divergent from the (historically Euro-centric) IPD reference. "
                 f"Grey/hatched = fewer than {min_cell_n} people.",
                 fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)
    return sub


def plot_td_sweep(sub, td_max_values, n_people_by_anc, min_cell_n, out_path):
    """Companion sweep: % of haplotype calls at template_distance <= each cut, by ancestry --
    the 'movement as the filter tightens' figure SCHEMA.md insists on, and the raw-distribution
    complement the project's own commit history flagged as needed (a hard distance==0 cutoff was
    found unusable on real data).

    Thin-N treatment (audit fix, same principle as the sr/lr disagreement bars): a line built
    from few PEOPLE (n_people_by_anc, not the row count already shown in the legend) is drawn
    dashed and at reduced opacity so it doesn't visually compete with a well-powered line for the
    reader's eye, matching the same demotion used elsewhere for thin heatmap cells and bars."""
    fig, ax = plt.subplots(figsize=(8.5, 5))
    for anc in vc.ANCESTRY_ORDER:
        vals = sub.loc[sub["ancestry_pred"] == anc, "template_distance"].values
        if len(vals) == 0:
            continue
        pct = [100 * np.mean(vals <= t) for t in td_max_values]
        thin = n_people_by_anc.get(anc, 0) < min_cell_n
        ax.plot(td_max_values, pct, marker="o",
                color=("#BBBBBB" if thin else vc.ANCESTRY_COLORS[anc]),
                linestyle="--" if thin else "-", alpha=0.6 if thin else 1.0,
                label=f"{anc} (n_people={n_people_by_anc.get(anc, 0)}{', THIN' if thin else ''})")
    ax.set_xlabel("template_distance cut (<=)")
    ax.set_ylabel("% of haplotype calls retained")
    ax.set_title("Cumulative template_distance retention by ancestry\n"
                 f"(the movement here, not a single 0/1 threshold, is the result. Dashed/pale = "
                 f"fewer than {min_cell_n} people.)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 3. Resolution cascade (VIZ_LIT.md 3.4)
# ---------------------------------------------------------------------------
def resolution_cascade_entropy(table1, cohort_people, genes):
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)].copy()
    sub["a2"] = sub["consensus"].map(lambda a: vc.to_nfield(a, 2))
    sub["a4"] = sub["consensus"].map(lambda a: vc.to_nfield(a, 4))
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub.dropna(subset=["a2", "a4"])
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]

    n_people_by_anc = {anc: cohort_people.loc[cohort_people["ancestry_pred"] == anc,
                                               "person_id"].nunique()
                        for anc in vc.ANCESTRY_ORDER}

    rows = []
    for gene in genes:
        g = sub[sub["gene_bare"] == gene]
        for anc in vc.ANCESTRY_ORDER:
            ga = g[g["ancestry_pred"] == anc]
            if ga.empty:
                continue
            total = len(ga)
            cond_h = 0.0
            for a2, grp in ga.groupby("a2"):
                p_a2 = len(grp) / total
                sub_counts = grp["a4"].value_counts(normalize=True).values
                cond_h += p_a2 * vc.shannon_entropy(sub_counts)
            rows.append(dict(gene=gene, ancestry=anc, n=total,
                              n_people=n_people_by_anc.get(anc, 0),
                              n_2field_groups=ga["a2"].nunique(),
                              n_4field_subtypes=ga["a4"].nunique(),
                              conditional_entropy=cond_h))
    return pd.DataFrame(rows)


def plot_resolution_cascade(cascade_df, genes, min_cell_n, out_path):
    """Audit fix: this bar chart previously had NO thin-N treatment at all -- a gene/ancestry cell
    resting on a handful of PEOPLE (n_people, not the haplotype-call row count `n`, which can be
    up to ~16x larger per person) was plotted as an ordinary full-color, equal-weight bar, exactly
    the failure mode reported for the sr/lr disagreement figure. Fixed the same way: grey+hatch
    for thin cells, plus explicit Patch legend proxies (a single ax.bar() call per ancestry mixes
    thin/non-thin patches, so get_legend_handles_labels() could pick up a grey patch as that
    ancestry's representative color -- the same legend-fidelity bug fixed earlier in
    05_figures_frequency.py's plot_bar_gene)."""
    fig, ax = plt.subplots(figsize=(9, 5))
    width = 0.8 / len(vc.ANCESTRY_ORDER)
    x = np.arange(len(genes))
    for i, anc in enumerate(vc.ANCESTRY_ORDER):
        vals, thin = [], []
        for gene in genes:
            row = cascade_df[(cascade_df["gene"] == gene) & (cascade_df["ancestry"] == anc)]
            if row.empty:
                vals.append(np.nan)
                thin.append(True)
            else:
                vals.append(row["conditional_entropy"].iloc[0])
                thin.append(int(row["n_people"].iloc[0]) < min_cell_n)
        offs = x + (i - (len(vc.ANCESTRY_ORDER) - 1) / 2) * width
        colors = [vc.ANCESTRY_COLORS[anc] if not t else "#DDDDDD" for t in thin]
        bars = ax.bar(offs, np.nan_to_num(vals), width=width * 0.95, color=colors)
        for bar, t in zip(bars, thin):
            if t:
                bar.set_hatch("////")
    ax.set_xticks(x)
    ax.set_xticklabels(genes)
    ax.set_ylabel("H(4-field | 2-field) [nats]")
    ax.set_title("Resolution cascade: how much 4-field detail a 2-field call hides,\n"
                 f"by ancestry (higher = more ancestry-specific subtype diversity lost at "
                 f"2-field). Grey/hatched = fewer than {min_cell_n} people.", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=vc.ANCESTRY_COLORS[a]) for a in vc.ANCESTRY_ORDER]
    ax.legend(handles, vc.ANCESTRY_ORDER, fontsize=8, ncol=6, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, 1.15))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 4. Non-classical gene diversity map + coverage matrix (VIZ_LIT.md 3.5)
# ---------------------------------------------------------------------------
def plot_nonclassical_diversity(long_df_nc, genes, min_cell_n, out_path):
    rows = []
    for gene in genes:
        g = long_df_nc[long_df_nc["gene_bare"] == gene]
        for anc in vc.ANCESTRY_ORDER:
            ga = g[g["ancestry_pred"] == anc]
            n = len(ga)
            if n == 0:
                rows.append(dict(gene=gene, ancestry=anc, n=0, heterozygosity=np.nan))
                continue
            freqs = ga["allele_label"].value_counts(normalize=True).values
            rows.append(dict(gene=gene, ancestry=anc, n=n,
                              heterozygosity=vc.expected_heterozygosity(freqs)))
    df = pd.DataFrame(rows)
    mat = df.pivot(index="gene", columns="ancestry", values="heterozygosity").reindex(
        index=genes, columns=vc.ANCESTRY_ORDER)
    nmat = df.pivot(index="gene", columns="ancestry", values="n").reindex(
        index=genes, columns=vc.ANCESTRY_ORDER)

    # Bug found against the fixtures: imshow renders a NaN cell (a gene with 0 observations in
    # this ancestry) as fully transparent -- i.e. the plain white figure background shows through.
    # The thin-N hatch/text below was drawn in white (matching the convention used elsewhere in
    # this project's heatmaps, which works because those cells still have a real, non-white
    # viridis-colored backdrop). White-on-white is invisible: those rows rendered as blank space
    # with no hatch, no "n=0" text, no visual signal that anything was even attempted. Fix:
    # fillna(0) before imshow so every cell gets a real (low-end-of-colormap) color, same
    # convention 05_figures_frequency.py's heatmap already uses.
    fig, ax = plt.subplots(figsize=(7, 0.5 + 0.4 * len(genes)))
    im = ax.imshow(mat.fillna(0.0).values.astype(float), cmap="viridis", vmin=0)
    ax.set_xticks(range(len(vc.ANCESTRY_ORDER)))
    ax.set_xticklabels(vc.ANCESTRY_ORDER)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([vc.gene_display(g) for g in genes])
    for i in range(len(genes)):
        for j in range(len(vc.ANCESTRY_ORDER)):
            n = nmat.values[i, j]
            if pd.isna(n) or n < min_cell_n:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////",
                                            edgecolor="white", linewidth=0))
                ax.text(j, i, f"n={int(n) if not pd.isna(n) else 0}", ha="center", va="center",
                        fontsize=6, color="white")
    ax.set_title("Non-classical/framework gene diversity (expected heterozygosity) by ancestry\n"
                 "genes unavailable in AoU-native (sr) calls -- LR-only view", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.03)
    fig.tight_layout()
    vc.savefig(fig, out_path)
    return df


def plot_coverage_matrix(table1, sr_long, out_path):
    all_genes = CLASSICAL + NONCLASSICAL
    lr_genes = set(table1["gene_bare"].unique())
    sr_genes = set(sr_long["gene_bare"].unique())
    mat = np.zeros((len(all_genes), 2))
    for i, g in enumerate(all_genes):
        mat[i, 0] = 1.0 if g in sr_genes else 0.0
        mat[i, 1] = 1.0 if g in lr_genes else 0.0
    fig, ax = plt.subplots(figsize=(4, 0.35 * len(all_genes) + 1))
    ax.imshow(mat, cmap="Greens", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["sr", "lr"])
    ax.set_yticks(range(len(all_genes)))
    ax.set_yticklabels([vc.gene_display(g) for g in all_genes], fontsize=7)
    ax.set_title("Gene coverage by cohort", fontsize=10)
    for i in range(len(all_genes)):
        for j in range(2):
            ax.text(j, i, "yes" if mat[i, j] else "no", ha="center", va="center", fontsize=6)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 5. Cis heterodimer frequency + pairable fraction (VIZ_LIT.md 1.8/3.3, SCHEMA.md Table 2)
# ---------------------------------------------------------------------------
def compute_pairable_fraction(table1, cohort_people, pairs):
    """Re-derived directly from Table 1's (person_id, hap, contig, gene) grain -- independent of
    01_extract_rich.py's own stderr-only unpairable count, and broken down by ancestry (which
    that count is not). For each (person_id, hap) where BOTH genes of a pair are present anywhere
    in that hap file, records whether they share a contig (pairable) or not (fragmented)."""
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    sub = table1[table1["person_id"].isin(person_ids)]
    rows = []
    for (pid, hap), grp in sub.groupby(["person_id", "hap"]):
        anc = anc_map.get(pid)
        # anc can be pd.NA/NaN for the small number of people with a missing/malformed
        # ancestry_pred (ENVIRONMENT.md quirk #30) -- `NA not in list` raises TypeError rather
        # than returning True, because `in` evaluates `NA == x` per element and that returns NA,
        # not a bool (identical failure mode fixed in 06_figures_structure.py's dosage-matrix
        # filter). Guard with an explicit string check so NA/NaN short-circuits to "skip" instead
        # of reaching the `in` test.
        if not isinstance(anc, str) or anc not in vc.ANCESTRY_ORDER:
            continue
        by_gene = grp.groupby("gene_bare")["contig"].apply(set)
        for ga, gb in pairs:
            if ga not in by_gene.index or gb not in by_gene.index:
                continue
            shared = bool(by_gene[ga] & by_gene[gb])
            rows.append(dict(person_id=pid, hap=hap, pair=f"{ga}~{gb}", ancestry=anc,
                              pairable=shared))
    return pd.DataFrame(rows)


def plot_heterodimer_and_pairable(table2, pairable_df, cohort_people, pairs, min_cell_n, out_dir):
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    fig_paths = []
    report_lines = []
    for ga, gb in pairs:
        pair_label = f"{ga}~{gb}"
        sub = table2[table2["pair"] == pair_label].copy()
        sub["ancestry_pred"] = sub["person_id"].map(anc_map)
        sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
        if sub.empty:
            report_lines.append(f"- {pair_label}: 0 phased cis pairs found for this cohort.")
            continue
        counts = sub.groupby(["haplotype_label", "ancestry_pred"]).size().unstack(fill_value=0)
        counts = counts.reindex(columns=vc.ANCESTRY_ORDER, fill_value=0)
        totals = counts.sum(axis=0)
        top = counts.sum(axis=1).sort_values(ascending=False).head(8).index

        fig, ax = plt.subplots(figsize=(max(7, len(top) * 1.1), 5))
        width = 0.8 / len(vc.ANCESTRY_ORDER)
        x = np.arange(len(top))
        for i, anc in enumerate(vc.ANCESTRY_ORDER):
            n = totals.get(anc, 0)
            c = np.array([counts.loc[h, anc] if h in counts.index else 0 for h in top])
            p, lo, hi = vc.wilson_ci(c, np.full_like(c, n, dtype=float))
            thin = c < min_cell_n
            offs = x + (i - (len(vc.ANCESTRY_ORDER) - 1) / 2) * width
            colors = [vc.ANCESTRY_COLORS[anc] if not t else "#DDDDDD" for t in thin]
            ax.bar(offs, np.nan_to_num(p), width=width * 0.95, color=colors)
            # Same nan_to_num-vs-CI mismatch the coordinator flagged in 05_figures_frequency.py's
            # plot_bar_gene: lo/hi were computed from wilson_ci but never used here (dead CI, no
            # error bars at all) -- adding them the same safe way: mask NaN-p (no-data) cells out
            # of the errorbar call entirely instead of zero-filling them, which is exactly what
            # produced the negative-yerr crash in 05.
            has_data = ~np.isnan(p)
            err_lo = (p - lo)[has_data]
            err_hi = (hi - p)[has_data]
            assert np.all(err_lo >= -1e-9), f"{pair_label}/{anc}: negative err_lo {err_lo}"
            assert np.all(err_hi >= -1e-9), f"{pair_label}/{anc}: negative err_hi {err_hi}"
            if has_data.any():
                ax.errorbar(offs[has_data], p[has_data], yerr=[np.clip(err_lo, 0, None),
                                                                 np.clip(err_hi, 0, None)],
                             fmt="none", ecolor="black", elinewidth=0.6, capsize=1.5, alpha=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(top, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Cis-haplotype frequency")
        ax.set_title(f"{pair_label} phased cis-haplotype frequency by ancestry "
                     f"(true phase, not EM-imputed)", fontsize=10)
        # Same legend-fidelity bug as 05_figures_frequency.py's plot_bar_gene, fixed the same way:
        # explicit Patch proxies so a thin-N-grey bar can never masquerade as an ancestry's color
        # in the legend.
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=vc.ANCESTRY_COLORS[a]) for a in vc.ANCESTRY_ORDER]
        ax.legend(handles, vc.ANCESTRY_ORDER, fontsize=7, ncol=3, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        out_path = os.path.join(out_dir, f"heterodimer_{ga}_{gb}.png")
        vc.savefig(fig, out_path)
        fig_paths.append(out_path)

        pf = pairable_df[pairable_df["pair"] == pair_label]
        report_lines.append(f"\n### {pair_label}\n")
        report_lines.append(f"Phased cis-haplotype copies observed by ancestry: "
                             f"{dict(totals.astype(int))}")
        if not pf.empty:
            frac = pf.groupby("ancestry")["pairable"].mean()
            n_attempt = pf.groupby("ancestry").size()
            report_lines.append("\n| Ancestry | pairable fraction | n hap-instances with both genes present |")
            report_lines.append("|---|---|---|")
            for anc in vc.ANCESTRY_ORDER:
                if anc in frac.index:
                    report_lines.append(f"| {anc} | {frac[anc]:.1%} | {int(n_attempt[anc])} |")

    return fig_paths, report_lines


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    vc.add_common_args(ap)
    args = ap.parse_args()

    if args.cohort == "sr":
        sys.exit("FATAL: 07_figures_crosscohort.py's figures all compare sr against a long-read "
                 "side by construction -- pass --cohort lr or --cohort lr_td --td-max N to pick "
                 "the long-read side; sr is always the fixed short-read baseline (see this "
                 "script's module docstring, 'On the three-cohort abstraction here').")

    table1 = vc.load_table1(args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv"))
    table2 = vc.load_table2(args.table2 or os.path.join(args.outroot, "hla_cis_pairs.tsv"))
    sr_long = vc.load_sr_genotypes(args.sr_genotypes or os.path.join(args.outroot, "hla_genotypes.tsv"))
    cohort_df = vc.load_cohort_membership(
        args.cohort_membership or os.path.join(args.outroot, "cohort_membership.tsv"))

    lr_people, lr_label = vc.select_cohort_people(cohort_df, args.cohort, args.td_max)
    sr_people, sr_label = vc.select_cohort_people(cohort_df, "sr")
    print(f"LR side: {lr_label} -- {len(lr_people)} people. SR side: {sr_label} -- "
          f"{len(sr_people)} people.", file=sys.stderr)
    if len(lr_people) < 20 or len(sr_people) < 20:
        sys.exit("FATAL: too few people on one side of the comparison.")

    out_dir = args.out_dir or vc.default_out_dir("07_figures_crosscohort", args.cohort, args.td_max)
    vc.ensure_dir(out_dir)
    report = [f"# Cross-cohort figures -- sr vs `{lr_label}`\n",
              f"SR: {len(sr_people)} people. LR side: {len(lr_people)} people.\n"]
    fig_paths = []

    # ---- 1. SR vs LR frequency scatter ----------------------------------------------------
    print("Computing SR vs LR allele frequency scatter...", file=sys.stderr)
    lr_long2 = build_long_frame(table1, lr_people, CLASSICAL, resolution=2)
    sr_long2 = build_sr_long_frame(sr_long, sr_people, CLASSICAL)
    sr_freq = pooled_freq(sr_long2, CLASSICAL)
    lr_freq = pooled_freq(lr_long2, CLASSICAL)
    scatter_path = os.path.join(out_dir, "sr_vs_lr_frequency_scatter.png")
    outliers = plot_sr_vs_lr_scatter(sr_freq, lr_freq, CLASSICAL, scatter_path)
    fig_paths.append(scatter_path)
    report.append("\n## SR vs LR frequency outliers (>2x, combined n>=10)\n")
    report.append("| Allele | sr freq | lr freq | ratio |")
    report.append("|---|---|---|---|")
    for a, fx, fy in outliers:
        ratio = fy / fx if fx > 0 else float("inf")
        report.append(f"| {a} | {fx:.3%} | {fy:.3%} | {ratio:.2f}x |")

    # ---- 1b. SR-vs-LR disagreement BY ANCESTRY -- the pooled scatter above blurs exactly the
    # ancestry-dependent reference-bias signal VIZ_LIT.md 3.1 is about; this stratifies it.
    print("Computing SR-vs-LR disagreement by ancestry...", file=sys.stderr)
    sr_freq_anc = freq_by_ancestry(sr_long2, CLASSICAL)
    lr_freq_anc = freq_by_ancestry(lr_long2, CLASSICAL)
    n_people_lr = lr_people["ancestry_pred"].value_counts().to_dict()
    n_people_sr = sr_people["ancestry_pred"].value_counts().to_dict()
    disagree_df = compute_sr_lr_disagreement_by_ancestry(sr_freq_anc, lr_freq_anc,
                                                          n_people_lr, n_people_sr,
                                                          min_n=args.min_cell_n)
    disagree_path = os.path.join(out_dir, "sr_lr_disagreement_by_ancestry.png")
    plot_sr_lr_disagreement_by_ancestry(disagree_df, args.min_cell_n_people, disagree_path)
    fig_paths.append(disagree_path)
    report.append("\n## SR-vs-LR agreement by ancestry (headline reference-bias claim, made "
                   "quantitative)\n")
    report.append(f"Bars for an ancestry with fewer than {args.min_cell_n_people} people in the "
                   f"more limited of the sr/lr cohorts are greyed out and hatched in the figure "
                   f"-- treat those as noise, not signal. Note this uses the STRICTER "
                   f"--min-cell-n-people threshold ({args.min_cell_n_people}), not "
                   f"--min-cell-n ({args.min_cell_n}): each bar here averages ACROSS an "
                   f"ancestry's whole allele set, so it needs more independent people to be "
                   f"trustworthy than a single raw allele-copy count does (see "
                   f"MIN_CELL_N_PEOPLE's docstring in _viz_common.py).\n")
    report.append("| Ancestry | mean abs freq diff | n alleles | EUR-common-allele bias "
                   "(sr-lr) | n EUR-common alleles | n people (lr) | n people (sr) | thin? |")
    report.append("|---|---|---|---|---|---|---|---|")
    for _, r in disagree_df.iterrows():
        mad = "" if pd.isna(r["mean_abs_diff"]) else f"{r['mean_abs_diff']:.3%}"
        bias = "" if pd.isna(r["eur_common_allele_bias"]) else f"{r['eur_common_allele_bias']:+.3%}"
        thin = min(r["n_people_lr"], r["n_people_sr"]) < args.min_cell_n_people
        report.append(f"| {r['ancestry']} | {mad} | {int(r['n_alleles'])} | {bias} | "
                       f"{int(r['n_eur_common_alleles'])} | {int(r['n_people_lr'])} | "
                       f"{int(r['n_people_sr'])} | {'YES' if thin else ''} |")

    # ---- 2. template_distance by ancestry + sweep ------------------------------------------
    print("Computing template_distance by ancestry...", file=sys.stderr)
    td_path = os.path.join(out_dir, "template_distance_by_ancestry.png")
    td_sub = plot_td_by_ancestry(table1, lr_people, CLASSICAL, args.min_cell_n_people, td_path)
    fig_paths.append(td_path)
    td_sweep_vals = [int(v) for v in args.td_max_sweep.split(",")]
    sweep_path = os.path.join(out_dir, "template_distance_sweep.png")
    n_people_lr_map = lr_people["ancestry_pred"].value_counts().to_dict()
    plot_td_sweep(td_sub, td_sweep_vals, n_people_lr_map, args.min_cell_n_people, sweep_path)
    fig_paths.append(sweep_path)
    report.append("\n## template_distance: % of haplotype calls at or under each cut, by ancestry\n")
    report.append("| Ancestry | " + " | ".join(f"td<={t}" for t in td_sweep_vals) + " | n |")
    report.append("|---|" + "---|" * (len(td_sweep_vals) + 1))
    for anc in vc.ANCESTRY_ORDER:
        vals = td_sub.loc[td_sub["ancestry_pred"] == anc, "template_distance"].values
        if len(vals) == 0:
            continue
        cells = [f"{100 * np.mean(vals <= t):.1f}%" for t in td_sweep_vals]
        report.append(f"| {anc} | " + " | ".join(cells) + f" | {len(vals)} |")

    # ---- 3. Resolution cascade --------------------------------------------------------------
    print("Computing resolution cascade (2-field -> 4-field)...", file=sys.stderr)
    cascade_df = resolution_cascade_entropy(table1, lr_people, CLASSICAL)
    cascade_path = os.path.join(out_dir, "resolution_cascade.png")
    if not cascade_df.empty:
        plot_resolution_cascade(cascade_df, CLASSICAL, args.min_cell_n_people, cascade_path)
        fig_paths.append(cascade_path)
        report.append("\n## Resolution cascade: H(4-field | 2-field), nats\n")
        report.append(f"Grey/hatched bars in the figure: fewer than {args.min_cell_n_people} "
                       f"people for that ancestry (`n_people` below; `n` is haplotype-call ROWS, "
                       f"which can be much larger even for a thin ancestry -- do not read `n` "
                       f"alone as statistical power).\n")
        report.append("| Gene | Ancestry | n (rows) | n_people | n 2-field groups | "
                       "n 4-field subtypes | H(4|2) | thin? |")
        report.append("|---|---|---|---|---|---|---|---|")
        for _, r in cascade_df.iterrows():
            thin = int(r["n_people"]) < args.min_cell_n_people
            report.append(f"| {r['gene']} | {r['ancestry']} | {int(r['n'])} | "
                           f"{int(r['n_people'])} | {int(r['n_2field_groups'])} | "
                           f"{int(r['n_4field_subtypes'])} | {r['conditional_entropy']:.3f} | "
                           f"{'YES' if thin else ''} |")
    else:
        report.append("\nSkipped resolution cascade: no rows with both 2-field and 4-field "
                       "resolution in this cohort.\n")

    # ---- 4. Non-classical gene diversity + coverage ------------------------------------------
    print("Computing non-classical gene diversity map...", file=sys.stderr)
    lr_long_nc = build_long_frame(table1, lr_people, NONCLASSICAL, resolution=2)
    nc_path = os.path.join(out_dir, "nonclassical_diversity.png")
    nc_df = plot_nonclassical_diversity(lr_long_nc, NONCLASSICAL, args.min_cell_n, nc_path)
    fig_paths.append(nc_path)
    coverage_path = os.path.join(out_dir, "gene_coverage_by_cohort.png")
    plot_coverage_matrix(table1, sr_long, coverage_path)
    fig_paths.append(coverage_path)
    report.append("\n## Non-classical gene heterozygosity by ancestry (LR-only; sr has no calls "
                   "for these genes)\n")
    report.append("| Gene | Ancestry | n | Heterozygosity |")
    report.append("|---|---|---|---|")
    for _, r in nc_df.iterrows():
        het = "" if pd.isna(r["heterozygosity"]) else f"{r['heterozygosity']:.3f}"
        flag = " *thin*" if r["n"] < args.min_cell_n else ""
        report.append(f"| {r['gene']} | {r['ancestry']} | {int(r['n'])} | {het}{flag} |")

    # ---- 5. Cis heterodimer + pairable fraction -----------------------------------------------
    print("Computing cis-heterodimer frequency + pairable fraction by ancestry...", file=sys.stderr)
    pairable_df = compute_pairable_fraction(table1, lr_people, PAIRS)
    het_figs, het_report = plot_heterodimer_and_pairable(
        table2, pairable_df, lr_people, PAIRS, args.min_cell_n, out_dir)
    fig_paths += het_figs
    report.append("\n## Cis heterodimer frequency + pairable fraction by ancestry\n")
    report.extend(het_report)

    report.append("\n## Figures\n")
    for p in fig_paths:
        report.append(f"- `{p}`")

    vc.write_report(os.path.join(out_dir, "crosscohort_report.md"), report)


if __name__ == "__main__":
    main()
