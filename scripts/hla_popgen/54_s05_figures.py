#!/usr/bin/env python3
"""S05 figures from the committed PUBLIC aggregate tables of scripts 50/51/52/53 (no VM, no
participant data). One PNG per finding, written next to the tables they draw.

    python3 scripts/hla_popgen/54_s05_figures.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.path.join(REPO, "reports", "hla_popgen")
ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300", "#e34948"]
ANC_COL = dict(zip(ANC, SERIES))
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e0"

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.axisbelow": True, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
    "font.family": "DejaVu Sans"})


def num(s):
    return pd.to_numeric(s, errors="coerce")


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO))


def fig_waterfall():
    d = pd.read_csv(f"{R}/50_qc_foundation/waterfall.tsv", sep="\t")
    d = d[~d["step"].str.startswith("lost_")]
    keep = ["all_cohort_membership_rows", "in_lr", "unrelated_kin_lt_0.0442",
            "unrelated_and_ancestry_p_ge_0.90", "unrelated_and_ancestry_p_ge_0.95",
            "unrelated_and_ancestry_p_ge_0.98"]
    lab = {"all_cohort_membership_rows": "LR cohort rows (incl. sequel2)",
           "in_lr": "HLA called (Revio + Sequel IIe)", "unrelated_kin_lt_0.0442": "unrelated (kin < 0.0442)",
           "unrelated_and_ancestry_p_ge_0.90": "+ ancestry p ≥ 0.90",
           "unrelated_and_ancestry_p_ge_0.95": "+ ancestry p ≥ 0.95",
           "unrelated_and_ancestry_p_ge_0.98": "+ ancestry p ≥ 0.98"}
    d = d.set_index("step").loc[keep]
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    y = np.arange(len(d))[::-1]
    left = np.zeros(len(d))
    for a in ANC:
        v = num(d[a]).fillna(0).values
        ax.barh(y, v, left=left, color=ANC_COL[a], height=0.62, edgecolor="white", linewidth=1, label=a)
        left += v
    for yi, t in zip(y, num(d["total"]).values):
        ax.text(t + 120, yi, f"{int(t):,}", va="center", fontsize=8, color=INK)
    ax.set_yticks(y, [lab[k] for k in d.index])
    ax.set_xlabel("people")
    ax.set_xlim(0, 15000)
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=6, loc="lower center", bbox_to_anchor=(0.42, 1.0), fontsize=8, handlelength=1)
    save(fig, f"{R}/50_qc_foundation/fig_waterfall.png")


def fig_dq():
    d = pd.read_csv(f"{R}/50_qc_foundation/dq_ascertainment.tsv", sep="\t", dtype=str)
    a = d[(d["scope"] == "ALL") & (d["DQB1_group"] != "any")].copy()
    a["lab"] = "DQA1 " + a["DQA1_group"] + " / DQB1 " + a["DQB1_group"]
    a["pct"] = num(a["pct_split"])
    a["cross"] = a["DQA1_group"] != a["DQB1_group"]
    a = a.sort_values(["cross", "lab"])
    fig, ax = plt.subplots(figsize=(6.2, 2.3))
    y = np.arange(len(a))[::-1]
    ax.barh(y, a["pct"], color=[SERIES[1] if c else SERIES[0] for c in a["cross"]], height=0.6)
    for yi, p, n in zip(y, a["pct"], a["n_haps"]):
        if pd.isna(p):  # masked so a small one-contig count cannot be recovered
            ax.text(1.5, yi, f"masked (n={n})", va="center", fontsize=8, color=INK)
        else:
            ax.text(p + 1.5, yi, f"{p:.1f}%  (n={n})", va="center", fontsize=8, color=INK)
    ax.set_yticks(y, a["lab"])
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of haplotypes whose DQA1 and DQB1 sit on different contigs")
    ax.grid(axis="y", visible=False)
    ax.set_title("Cross-group DQ pairings occur almost only across contig breaks", loc="left", color=INK)
    save(fig, f"{R}/50_qc_foundation/fig_dq_ascertainment.png")


def fig_same_contig():
    d = pd.read_csv(f"{R}/50_qc_foundation/same_contig_matrix.tsv", sep="\t")
    d = d[d["stratum"] == "all"]
    genes = ["A", "C", "B", "MICA", "C4A", "DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "TAP1"]
    m = pd.DataFrame(np.nan, index=genes, columns=genes)
    for _, r in d.iterrows():
        if r["gene_a"] in genes and r["gene_b"] in genes:
            m.loc[r["gene_a"], r["gene_b"]] = m.loc[r["gene_b"], r["gene_a"]] = r["pct_same_contig"]
    for g in genes:
        m.loc[g, g] = 100
    fig, ax = plt.subplots(figsize=(4.9, 4.2))
    im = ax.imshow(m.values.astype(float), cmap="Blues", vmin=0, vmax=100)
    ax.set_xticks(range(len(genes)), genes, rotation=90)
    ax.set_yticks(range(len(genes)), genes)
    ax.grid(False)
    for i in range(len(genes)):
        for j in range(len(genes)):
            v = m.values[i, j]
            if i != j and not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                        color="white" if v > 60 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("% haplotypes with both genes on one contig")
    ax.set_title("Physical phasing level across the MHC", loc="left", color=INK)
    save(fig, f"{R}/50_qc_foundation/fig_same_contig_matrix.png")


def fig_novel_waterfall():
    d = pd.read_csv(f"{R}/50_qc_foundation/protein_novel_waterfall.tsv", sep="\t").set_index("gene")
    steps = ["distinct_protein_novel", "ge2_unrelated_carriers", "not_homopolymer_only",
             "no_template_warning", "ge2_platforms", "ge20_carriers_clean"]
    lab = ["distinct\nprotein-novel", "≥2 unrelated\ncarriers", "not homopolymer-\nindel only",
           "no template\nwarning", "seen on\n≥2 platforms", "≥20 carriers\n(clean)"]
    groups = {"classical (A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1)": ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"],
              "TAP1, TAP2, MICA, MICB": ["TAP1", "TAP2", "MICA", "MICB"]}
    fig, ax = plt.subplots(figsize=(8.6, 3.0))
    x = np.arange(len(steps))
    w = 0.38
    for k, (name, genes) in enumerate(groups.items()):
        v = d.loc[genes, steps].sum().values
        b = ax.bar(x + (k - .5) * w, v, w * 0.92, color=SERIES[k], label=name)
        for xi, vi in zip(x + (k - .5) * w, v):
            ax.text(xi, vi * 1.15 + 0.3, f"{int(vi)}", ha="center", fontsize=7.5, color=INK)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, 20000)
    ax.set_xticks(x, lab, fontsize=7.5)
    ax.set_ylabel("distinct sequences (log)")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Protein-novel allele robustness funnel (unrelated people)", loc="left", color=INK)
    save(fig, f"{R}/50_qc_foundation/fig_protein_novel_funnel.png")


def fig_platform():
    d = pd.read_csv(f"{R}/50_qc_foundation/novelty_platform_mh_or.tsv", sep="\t")
    genes = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
    d = d[d["gene"].isin(genes)]
    fig, ax = plt.subplots(figsize=(5.6, 2.6))
    y = np.arange(len(genes))[::-1]
    for k, (oc, name) in enumerate([("nov_protein", "protein-novel"), ("nov_noncoding", "non-coding novel")]):
        s = d[d["outcome"] == oc].set_index("gene").loc[genes, "mh_or_revio_vs_sequel2e"]
        ax.plot(s.values, y + (k - .5) * 0.25, "o", ms=6, color=SERIES[k], label=name, mec="white", mew=1)
    ax.axvline(1, color=MUTED, lw=1, ls="--")
    ax.set_yticks(y, genes)
    ax.set_xlim(0.3, 1.2)
    ax.set_xlabel("Mantel–Haenszel OR, Revio vs Sequel IIe (strata = ancestry)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Older platform → more 'novel' calls, within ancestry", loc="left", color=INK)
    save(fig, f"{R}/50_qc_foundation/fig_platform_novelty_or.png")


def fig_kir_framework():
    d = pd.read_csv(f"{R}/51_kir_structure_qc/framework_by_completeness.tsv", sep="\t").set_index("hap_class")
    cls = ["ALL", "complete", "multi_contig", "one_contig_open_end"]
    names = {"ALL": "all haplotypes", "complete": "complete (1 contig,\n3DL3…3DL2)",
             "multi_contig": "multi-contig", "one_contig_open_end": "1 contig,\nopen end"}
    genes = ["KIR3DL3", "KIR2DL4", "KIR3DP1", "KIR3DL2"]
    fig, ax = plt.subplots(figsize=(6.6, 2.6))
    x = np.arange(len(cls))
    w = 0.2
    for k, g in enumerate(genes):
        v = num(d.loc[cls, f"{g}_pct"]).values
        ax.bar(x + (k - 1.5) * w, v, w * 0.9, color=SERIES[k], label=g.replace("KIR", ""))
    ax.set_xticks(x, [f"{names[c]}\n(n={d.loc[c, 'n_haps']})" for c in cls], fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% haplotypes carrying gene")
    ax.grid(axis="x", visible=False)
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=8)
    save(fig, f"{R}/51_kir_structure_qc/fig_framework_by_completeness.png")


def fig_kir_pairs():
    d = pd.read_csv(f"{R}/51_kir_structure_qc/allelic_pair_events.tsv", sep="\t", dtype=str)
    d = d[d["scope"] == "all_haps"]
    fig, ax = plt.subplots(figsize=(6.0, 1.9))
    y = np.arange(len(d))[::-1]
    left = np.zeros(len(d))
    kinds = [("overlap", "same locus annotated twice"), ("same_contig", "distinct loci, same contig"),
             ("diff_contig", "different contigs")]
    for k, (c, name) in enumerate(kinds):
        v = num(d[f"{c}_share_pct"]).fillna(0).values
        ax.barh(y, v, left=left, color=[SERIES[7], SERIES[0], MUTED][k], height=0.55, edgecolor="white", label=name)
        left += v
    ax.set_yticks(y, [f"{p.replace('KIR', '').replace('+KIR', ' + ')}  (co-called {c}%)"
                      for p, c in zip(d["pair"], d["co_called_pct"])])
    for yi, l in zip(y, left):
        if l < 99:
            ax.text(l + 1, yi, "rest: masked", va="center", fontsize=7, color=INK2)
    ax.set_xlim(0, 115)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("% of co-call events")
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.4, 1.0), fontsize=7.5)
    save(fig, f"{R}/51_kir_structure_qc/fig_allelic_pairs.png")


def fig_atlas():
    d = pd.read_csv(f"{R}/52_mhc_haplotype_atlas/absent_cells_field1.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.45})
    for ax, pair, rows, cols, title in [
        (axes[0], "DQA1~DQB1", ["DQA1*01", "DQA1*02", "DQA1*03", "DQA1*04", "DQA1*05", "DQA1*06"],
         ["DQB1*02", "DQB1*03", "DQB1*04", "DQB1*05", "DQB1*06"], "DQA1 × DQB1 lineage"),
        (axes[1], "DRB1~DRBX", ["DRB1*01", "DRB1*15", "DRB1*16", "DRB1*03", "DRB1*11", "DRB1*12", "DRB1*13",
                                "DRB1*14", "DRB1*04", "DRB1*07", "DRB1*09", "DRB1*08", "DRB1*10"],
         ["DRBX*3", "DRBX*4", "DRBX*5"], "DRB1 lineage × paralog on haplotype")]:
        s = d[d["pair"] == pair].set_index(["allele_a", "allele_b"])
        m = np.full((len(rows), len(cols)), np.nan)
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                if (r, c) in s.index:
                    m[i, j] = s.loc[(r, c), "E_pooled"]
        ax.imshow(np.zeros_like(m), cmap="Greys", vmin=0, vmax=10, aspect="auto")
        im = ax.imshow(np.log10(m), cmap="Oranges", vmin=1, vmax=3.3, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(cols)):
                if not np.isnan(m[i, j]):
                    ax.text(j, i, f"{m[i, j]:.0f}", ha="center", va="center", fontsize=6.5,
                            color="white" if m[i, j] > 300 else INK)
        ax.set_xticks(range(len(cols)), [c.replace("DRBX*", "DRB").replace("DQB1*", "DQB1*\n") for c in cols], fontsize=7.5)
        ax.set_yticks(range(len(rows)), rows, fontsize=7.5)
        ax.grid(False)
        ax.set_title(title, loc="left", fontsize=9, color=INK)
    cb = fig.colorbar(im, ax=axes, shrink=0.8)
    cb.set_label("expected haplotypes if independent\n(cell observed = 0; log scale)")
    save(fig, f"{R}/52_mhc_haplotype_atlas/fig_absent_cells.png")


def fig_richness():
    d = pd.read_csv(f"{R}/52_mhc_haplotype_atlas/pooled_richness.tsv", sep="\t")
    genes = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
    pools = [("AFR", ANC_COL["AFR"]), ("EUR", ANC_COL["EUR"]), ("EAS", ANC_COL["EAS"]), ("AMR", ANC_COL["AMR"]),
             ("MIX_EUR45_AFR10_EAS45", INK2)]
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    x = np.arange(len(genes))
    for k, (p, col) in enumerate(pools):
        v = d[d["pool"] == p].set_index("gene").loc[genes, "richness_mean"].values
        ax.plot(x + (k - 2) * 0.12, v, "D" if p.startswith("MIX") else "o", ms=6, color=col, mec="white", mew=1,
                label="synthetic mix 45% EUR / 10% AFR / 45% EAS" if p.startswith("MIX") else p)
    ax.set_xticks(x, genes)
    ax.set_ylabel("2-field alleles observed")
    ax.set_ylim(0, 135)
    n = int(d["n_haps_equal"].iloc[0])
    ax.set_title(f"Equal-N richness ({n:,} haplotypes per pool, 50 resamples)", loc="left", color=INK)
    ax.grid(axis="x", visible=False)
    ax.legend(ncol=3, fontsize=7.5, loc="upper right")
    save(fig, f"{R}/52_mhc_haplotype_atlas/fig_equal_n_richness.png")


def fig_amr_sources():
    d = pd.read_csv(f"{R}/52_mhc_haplotype_atlas/amr_source_decomposition.tsv", sep="\t")
    genes = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
    cls = [("private_to_AMR", "not in AFR/EUR/EAS pools", SERIES[1]), ("only_AFR", "AFR only", ANC_COL["AFR"]),
           ("only_EUR", "EUR only", ANC_COL["EUR"]), ("only_EAS", "EAS only", ANC_COL["EAS"]),
           ("shared", "shared by ≥2 pools", "#c9c8c2")]
    d["c2"] = d["class"].where(~d["class"].str.startswith("shared"), "shared")
    t = d.groupby(["gene", "c2"])["pct_of_AMR_alleles"].sum().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    x = np.arange(len(genes))
    bottom = np.zeros(len(genes))
    for c, name, col in cls:
        v = t.reindex(genes)[c].fillna(0).values if c in t else np.zeros(len(genes))
        ax.bar(x, v, 0.6, bottom=bottom, color=col, edgecolor="white", linewidth=1, label=name)
        bottom += v
    ax.set_xticks(x, genes)
    ax.set_ylabel("% of AMR 2-field alleles")
    ax.set_ylim(0, 100)
    ax.grid(axis="x", visible=False)
    ax.legend(ncol=5, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    save(fig, f"{R}/52_mhc_haplotype_atlas/fig_amr_sources.png")


def fig_ase():
    p = f"{R}/53_hla_ase_pilot/ase_by_gene_class.tsv"
    if not os.path.exists(p):
        return
    d = pd.read_csv(p, sep="\t", dtype=str)
    d = d[d["pair_class"] != "mrna_identical"]
    d = d[num(d["median_abs_log2_obs"]).notna()]
    fig, ax = plt.subplots(figsize=(6.4, 2.4))
    lab = [f"HLA-{g}  {c.replace('protein_identical_mrna_diff', 'protein-identical').replace('protein_diff', 'protein-different')}"
           for g, c in zip(d["gene"], d["pair_class"])]
    y = np.arange(len(d))[::-1]
    ax.barh(y + 0.17, num(d["median_abs_log2_obs"]), 0.32, color=SERIES[0], label="observed |log2 hap1/hap2|")
    ax.barh(y - 0.17, num(d["median_abs_log2_sim_bias"]), 0.32, color="#c9c8c2", label="simulated (mapping bias)")
    ax.set_yticks(y, lab, fontsize=8)
    ax.set_xlabel("median |log2 allelic ratio|")
    ax.grid(axis="y", visible=False)
    ax.legend(fontsize=7.5, loc="lower right")
    save(fig, f"{R}/53_hla_ase_pilot/fig_ase.png")


def fig_structure_expression():
    p = f"{R}/55_structure_expression/structure_expression.tsv"
    if not os.path.exists(p):
        return
    d = pd.read_csv(p, sep="\t", dtype=str)
    panels = [("HLA-DRB5", "drb1_dqa1_same_contig_both_haps"), ("KIR2DL3", "kir_both_haps_complete"),
              ("KIR2DL1", "kir_both_haps_complete"), ("KIR2DL4", "kir_both_haps_complete"),
              ("KIR3DL1", "kir_both_haps_complete"), ("KIR2DS4", "kir_both_haps_complete")]
    fig, axes = plt.subplots(1, 6, figsize=(10.2, 2.5), gridspec_kw={"wspace": 0.5})
    for ax, (g, st) in zip(axes, panels):
        s = d[(d["gene"] == g) & (d["set"] == st) & (d["copies"].isin(["0", "1", "2"]))]
        rho = d[(d["gene"] == g) & (d["set"] == st) & (d["copies"] == "rho")]["tpm_median"]
        x = num(s["copies"]).values
        med, lo, hi = num(s["tpm_median"]).values, num(s["tpm_q25"]).values, num(s["tpm_q75"]).values
        ok = ~np.isnan(med)
        col = SERIES[0] if g == "HLA-DRB5" else (SERIES[2] if g in ("KIR2DL3", "KIR2DL1", "KIR2DL4") else SERIES[1])
        ax.vlines(x[ok], lo[ok], hi[ok], color=col, lw=2)
        ax.plot(x[ok], med[ok], "o", color=col, ms=6, mec="white", mew=1)
        ax.set_xticks([0, 1, 2])
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(0, None)
        ax.grid(axis="x", visible=False)
        ax.set_title(f"{g}\nρ = {float(rho.iloc[0]):.2f}" if len(rho) else g, fontsize=8.5, color=INK)
        ax.set_xlabel("copies called", fontsize=8)
    axes[0].set_ylabel("TPM (median, IQR)")
    save(fig, f"{R}/55_structure_expression/fig_structure_expression.png")


if __name__ == "__main__":
    fig_waterfall()
    fig_dq()
    fig_same_contig()
    fig_novel_waterfall()
    fig_platform()
    fig_kir_framework()
    fig_kir_pairs()
    fig_atlas()
    fig_richness()
    fig_amr_sources()
    fig_ase()
    fig_structure_expression()
