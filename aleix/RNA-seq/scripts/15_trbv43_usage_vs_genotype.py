#!/usr/bin/env python3
"""15 -- TRBV4-3 expression vs germline genotype (report: reports/15_trbv43_usage_vs_genotype/).

Report 10 (W1) found that 27.3% of people never use TRBV4-3 in their expressed TRB repertoire,
far above the ~0.04% expected by chance at its usage level, and read that as a germline
absence. Marc's script 65 genotypes the same thing from DNA: TRBV4-3 sits on a ~21.7 kb
insertion ("insertion A", with TRBV3-2 and a second TRBV6-2/6-3 copy) that GRCh38 carries only
on an ALT contig; 65 counts long reads spanning the breakpoint and calls 0/1/2 copies.

This report puts the two side by side. It uses committed, de-identified aggregate tables
only, so it runs off the VM:

  RNA (this workstream)  reports/10_vj_usage_ancestry/{zero_usage_by_ancestry,dosage_histogram,
                         codeletion_vs_absence}.csv; denominators from
                         reports/09_analysis_table/waterfall.csv (people with >= 823 TRB reads)
  DNA (script 65, main)  inputs/65_*.tsv, copied unchanged (see inputs/PROVENANCE.txt)

  a  % of people without TRBV4-3, RNA (never used) vs DNA (0 copies), by ancestry
  b  the same, RNA against DNA, with the identity line and the pooled ratio
  c  absence-allele frequency: RNA (implied from never-use under HWE) vs DNA (counted)
  d  cumulative share of people by TRBV4-3 usage, against the DNA 0-copy rate
  e  usage of other TRBV genes in non-users vs users, along the locus (position from TRUST4's
     hg38 reference, inputs/trust4_trb_gene_coords_hg38.tsv)

What it does not do: match people. The cohorts overlap (RNA 7,527 people at the read floor,
relatives included; DNA 11,856 unrelated long-read people; 8,327 have both assays) but the
comparison here is group-level. The per-person join is Marc's script 60.

Usage (from aleix/RNA-seq/):
  python3 scripts/15_trbv43_usage_vs_genotype.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import repfig as R  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPORTS = os.path.join(HERE, "..", "reports")
GENE = "TRBV4-3"
ANC = R.ANCESTRY_ORDER
BREAKPOINT = 142_346_200  # 65 probe build: insertion A sits here on GRCh38 primary chr7


def wilson_pct(k, n):
    p, lo, hi = R.wilson(k, n)
    return 100 * p, 100 * lo, 100 * hi


def num(x):
    """Masked cells in the 65 tables read as '<20*'; treat them as missing."""
    return pd.to_numeric(x, errors="coerce")


def load(args):
    z = pd.read_csv(os.path.join(args.w1, "zero_usage_by_ancestry.csv")).drop_duplicates("gene")
    z = z.set_index("gene").loc[GENE]
    wf = pd.read_csv(os.path.join(args.w0, "waterfall.csv")).set_index("step")
    n_rna = wf.loc["TRB reads >= 823"]
    g = pd.read_csv(os.path.join(args.inputs, "65_genotype_by_stratum_unrelated.tsv"), sep="\t")
    g = g[(g["feature"] == GENE)].copy()
    for c in ["n_called", "n_0copies", "n_1copy", "n_2copies", "insertion_allele_freq", "pct_0copies"]:
        g[c] = num(g[c])
    g = g.set_index("stratum")

    rows = []
    for a in ANC + ["ALL"]:
        col = "all" if a == "ALL" else a
        n = int(n_rna[col])
        frac = float(z["obs_zero_frac"] if a == "ALL" else z[f"obs_zero_{a}"])
        k = int(round(frac * n))
        r_p, r_lo, r_hi = wilson_pct(k, n)
        q = float(z["implied_absence_allele_freq"] if a == "ALL" else z[f"implied_absence_allele_freq_{a}"])
        # implied absence-allele freq = sqrt(excess zero fraction); CI by the same transform
        exp0 = float(z["exp_zero_frac"] if a == "ALL" else z[f"exp_zero_{a}"])
        q_lo = np.sqrt(max(r_lo / 100 - exp0, 0) / (1 - exp0))
        q_hi = np.sqrt(max(r_hi / 100 - exp0, 0) / (1 - exp0))
        d = g.loc[a] if a in g.index else None
        nd = d["n_called"] if d is not None else np.nan
        k0 = d["n_0copies"] if d is not None else np.nan
        if np.isfinite(nd) and np.isfinite(k0):
            d_p, d_lo, d_hi = wilson_pct(k0, nd)
            alt = 2 * k0 + d["n_1copy"]                     # absence alleles counted directly
            dq, dq_lo, dq_hi = [v / 100 for v in wilson_pct(alt, 2 * nd)]
        else:
            d_p = d_lo = d_hi = dq = dq_lo = dq_hi = np.nan
        rows.append(dict(ancestry=a, rna_people=n, rna_never_used_pct=r_p, rna_lo=r_lo, rna_hi=r_hi,
                         dna_people_called=nd, dna_0copy_pct=d_p, dna_lo=d_lo, dna_hi=d_hi,
                         rna_ratio_to_dna=r_p / d_p if np.isfinite(d_p) else np.nan,
                         rna_implied_absence_allele_freq=q, rna_q_lo=q_lo, rna_q_hi=q_hi,
                         dna_absence_allele_freq=dq, dna_q_lo=dq_lo, dna_q_hi=dq_hi,
                         dna_hwe_p=num(pd.Series([d["hwe_p"]]))[0] if d is not None else np.nan))
    comp = pd.DataFrame(rows)

    dos = pd.read_csv(os.path.join(args.w1, "dosage_histogram.csv"))
    dos = dos[dos["gene"] == GENE].sort_values(["usage_pct_hi", "usage_pct_lo"])
    mix = pd.read_csv(os.path.join(args.inputs, "65_af_mixture_fit_unrelated.tsv"), sep="\t")
    mix = mix[mix["feature"] == GENE].set_index("platform")
    cod = pd.read_csv(os.path.join(args.w1, "codeletion_vs_absence.csv"))
    probes = pd.read_csv(os.path.join(args.inputs, "65_trb_probe_mapping.tsv"), sep="\t", comment="#")
    coords = pd.read_csv(os.path.join(args.inputs, "trust4_trb_gene_coords_hg38.tsv"), sep="\t")
    coords = coords[coords["chrom"] == "chr7"]
    return comp, dos, mix, cod, probes, coords, int(n_rna["all"])


def cumulative(dos, n_total):
    """Share of people at or below each usage level. Bins the W1 script masked (< 20 people)
    or beyond its 4.15% axis are missing, so the curve ends below 100%."""
    zero = float(dos.loc[dos["usage_pct_hi"] == 0, "people"].sum())
    pos = dos[dos["usage_pct_hi"] > 0]
    x = np.concatenate([[0.0], pos["usage_pct_hi"].to_numpy()])
    y = 100 * np.cumsum(np.concatenate([[zero], pos["people"].to_numpy()])) / n_total
    return pd.DataFrame({"usage_pct_at_most": x, "cum_people_pct": y})


def annotate_codeletion(cod, probes, coords):
    """Usage of every other TRBV gene in TRBV4-3 non-users vs users, placed on the locus.
    Position from TRUST4's own hg38 reference (the coordinates every run used); insertion A's
    gene content and breakpoint from 65's probe mapping."""
    cod = cod[(cod["reference_gene"] == GENE) & (cod["median_usage_pct_gene_present"] >= 0.3)].copy()
    ins = probes[probes["note"].fillna("").str.contains("insertion A")]["gene"].tolist()
    on_ins = {g for name in ins for g in name.split("/")}      # TRBV6-2/6-3 -> both names
    pos = coords.set_index("gene")["start"]
    cod["start_hg38"] = cod["gene"].map(pos)

    def where(r):
        if r.gene in on_ins:
            return "second copy on insertion A"
        return "5′ of insertion" if r.start_hg38 < BREAKPOINT else "3′ of insertion"
    cod["position"] = cod.apply(where, axis=1)
    return cod.sort_values("start_hg38").reset_index(drop=True)


def figure(comp, cum, mix, cod, outdir):
    import cnsplots as cns
    col = R.ancestry_colors()
    pal = R.nature()
    RNA_C, DNA_C = pal[3], pal[0]
    mp = cns.multipanel(max_width=540)
    W, H, GAP = 120, 120, 36
    sub = comp[comp["ancestry"] != "ALL"].reset_index(drop=True)

    # a: paired bars by ancestry
    ax = mp.panel("a", width=2 * W + 50, height=H, margin_right=30, margin_bottom=GAP)
    order = list(sub.sort_values("dna_0copy_pct", ascending=False, na_position="last")["ancestry"]) + ["ALL"]
    c2 = comp.set_index("ancestry").loc[order]
    x = np.arange(len(order))
    w = 0.38
    ax.bar(x - w / 2, c2["rna_never_used_pct"], width=w, color=RNA_C, lw=0, label="RNA: never used")
    ax.errorbar(x - w / 2, c2["rna_never_used_pct"],
                yerr=[c2["rna_never_used_pct"] - c2["rna_lo"], c2["rna_hi"] - c2["rna_never_used_pct"]],
                fmt="none", ecolor="black", elinewidth=0.6, capsize=0)
    ax.bar(x + w / 2, c2["dna_0copy_pct"], width=w, color=DNA_C, lw=0, label="DNA: 0 copies")
    ok = c2["dna_0copy_pct"].notna().to_numpy()
    ax.errorbar((x + w / 2)[ok], c2["dna_0copy_pct"][ok],
                yerr=[(c2["dna_0copy_pct"] - c2["dna_lo"])[ok], (c2["dna_hi"] - c2["dna_0copy_pct"])[ok]],
                fmt="none", ecolor="black", elinewidth=0.6, capsize=0)
    for i in np.where(~ok)[0]:
        ax.text(x[i] + w / 2, 1.5, "masked", rotation=90, ha="center", va="bottom", fontsize=5, color="0.4")
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.legend(frameon=False, fontsize=5.5, loc="upper right")
    R.style(ax, f"People without {GENE}", None, "% of people (95% CI)")

    # b: RNA vs DNA, per ancestry
    ax = mp.panel("b", width=W, height=H, margin_right=15, margin_bottom=GAP)
    s = sub.dropna(subset=["dna_0copy_pct"])
    lim = (0, 62)
    ax.plot(lim, lim, color="0.6", lw=0.6, ls="--", zorder=1)
    allrow = comp.set_index("ancestry").loc["ALL"]
    ratio = allrow["rna_never_used_pct"] / allrow["dna_0copy_pct"]
    ax.plot(lim, [v * ratio for v in lim], color="0.25", lw=0.6, zorder=1)
    for _, r in s.iterrows():
        ax.errorbar(r.dna_0copy_pct, r.rna_never_used_pct,
                    xerr=[[r.dna_0copy_pct - r.dna_lo], [r.dna_hi - r.dna_0copy_pct]],
                    yerr=[[r.rna_never_used_pct - r.rna_lo], [r.rna_hi - r.rna_never_used_pct]],
                    fmt="o", ms=3.5, color=col[r.ancestry], elinewidth=0.6, zorder=3)
        ax.annotate(r.ancestry, (r.dna_0copy_pct, r.rna_never_used_pct), xytext=(4, -7),
                    textcoords="offset points", fontsize=5.5, color=col[r.ancestry])
    ax.text(0.04, 0.96, f"dashed: equal\nsolid: pooled ratio {ratio:.2f}", transform=ax.transAxes,
            va="top", fontsize=5.5)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    R.style(ax, "RNA against DNA", "DNA: 0 copies (%)", "RNA: never used (%)")

    # c: absence-allele frequency
    ax = mp.panel("c", width=W, height=H, margin_right=15, margin_bottom=GAP)
    s = sub.dropna(subset=["dna_absence_allele_freq"])
    lim = (0.3, 0.8)
    ax.plot(lim, lim, color="0.6", lw=0.6, ls="--", zorder=1)
    for _, r in s.iterrows():
        ax.errorbar(r.dna_absence_allele_freq, r.rna_implied_absence_allele_freq,
                    xerr=[[r.dna_absence_allele_freq - r.dna_q_lo], [r.dna_q_hi - r.dna_absence_allele_freq]],
                    yerr=[[r.rna_implied_absence_allele_freq - r.rna_q_lo],
                          [r.rna_q_hi - r.rna_implied_absence_allele_freq]],
                    fmt="o", ms=3.5, color=col[r.ancestry], elinewidth=0.6, zorder=3)
        ax.annotate(r.ancestry, (r.dna_absence_allele_freq, r.rna_implied_absence_allele_freq),
                    xytext=(4, -7), textcoords="offset points", fontsize=5.5, color=col[r.ancestry])
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    R.style(ax, "Absence-allele frequency", "DNA: counted", "RNA: implied (HWE)")

    # d: cumulative usage vs DNA rate
    ax = mp.panel("d", width=W, height=H, margin_right=15, margin_bottom=GAP)
    d_all = allrow["dna_0copy_pct"]
    m0 = 100 * float(mix.loc["revio", "w_0copy"]) if "revio" in mix.index else np.nan
    ax.axhspan(allrow["dna_lo"], allrow["dna_hi"], color=DNA_C, alpha=0.35, lw=0)
    ax.step(cum["usage_pct_at_most"], cum["cum_people_pct"], where="post", color=RNA_C, lw=1.0, zorder=3)
    ax.text(0.97, 0.20, f"DNA: 0 copies, {d_all:.1f}% (band: 95% CI)", transform=ax.transAxes,
            ha="right", color=DNA_C, fontsize=5.5)
    if np.isfinite(m0):
        ax.text(0.97, 0.13, f"mixture genotyper, Revio: {m0:.1f}%", transform=ax.transAxes, ha="right",
                color=DNA_C, fontsize=5.5)
    ax.text(0.97, 0.06, "RNA: people with usage ≤ x", transform=ax.transAxes, ha="right",
            color=RNA_C, fontsize=5.5)
    ax.set_xlim(-0.03, 1.0)
    ax.set_ylim(0, 70)
    R.style(ax, f"{GENE} usage, cumulative", "Usage, x (% of clonotypes)", "% of people")

    # e: other genes along the locus
    ax = mp.panel("e", width=W + 30, height=H, margin_right=15, margin_bottom=GAP)
    pc = {"5′ of insertion": pal[4], "second copy on insertion A": DNA_C, "3′ of insertion": "0.6"}
    kb = ((cod["start_hg38"] - BREAKPOINT) / 1e3).to_numpy()
    ax.axhline(1, color="0.6", lw=0.5, ls="--")
    ax.axvline(0, color=DNA_C, lw=0.6)
    ax.text(4, 0.33, "insertion A", color=DNA_C, fontsize=5.5, va="bottom")
    for k, r in zip(kb, cod.itertuples()):
        ax.scatter(k, r.ratio_absent_over_present, s=9, color=pc[r.position], zorder=3, lw=0)
        if r.position == "3′ of insertion" and (r.ratio_absent_over_present < 0.8 or r.ratio_absent_over_present > 1.15) \
                or r.position == "second copy on insertion A":
            ax.annotate(r.gene.replace("TRBV", "V"), (k, r.ratio_absent_over_present), xytext=(4, 0),
                        textcoords="offset points", fontsize=5, va="center")
    five = cod["position"] == "5′ of insertion"
    if five.any():
        ax.annotate(f"{cod.loc[five, 'gene'].iloc[0].replace('TRBV', 'V')}–{cod.loc[five, 'gene'].iloc[-1].replace('TRBV', 'V')}",
                    (kb[five.to_numpy()].min(), cod.loc[five, "ratio_absent_over_present"].max()),
                    xytext=(-2, 6), textcoords="offset points", fontsize=5, ha="left", color=pal[4])
    for lab, c in pc.items():
        ax.scatter([], [], s=9, color=c, label=lab, lw=0)
    ax.set_ylim(0.3, 1.5)
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    R.style(ax, "Other TRBV genes", "Position from insertion (kb, GRCh38)", "Median usage, non-users / users")
    R.save("fig_trbv43_usage_vs_genotype", outdir)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--w1", default=os.path.join(REPORTS, "10_vj_usage_ancestry"))
    ap.add_argument("--w0", default=os.path.join(REPORTS, "09_analysis_table"))
    ap.add_argument("--inputs", default=os.path.join(REPORTS, "15_trbv43_usage_vs_genotype", "inputs"))
    ap.add_argument("--outdir", default=os.path.join(REPORTS, "15_trbv43_usage_vs_genotype"))
    args = ap.parse_args()

    comp, dos, mix, cod, probes, coords, n_total = load(args)
    cum = cumulative(dos, n_total)
    cod = annotate_codeletion(cod, probes, coords)
    os.makedirs(args.outdir, exist_ok=True)
    comp.round(4).to_csv(os.path.join(args.outdir, "comparison_by_ancestry.csv"), index=False)
    cum.round(4).to_csv(os.path.join(args.outdir, "cumulative_usage.csv"), index=False)
    cod.round(4).to_csv(os.path.join(args.outdir, "codeletion_by_position.csv"), index=False)

    s = comp.dropna(subset=["dna_0copy_pct"])
    s = s[s["ancestry"] != "ALL"]
    allrow = comp.set_index("ancestry").loc["ALL"]
    first_bin = cum.iloc[1]
    summary = [
        ("rna_people", n_total),
        ("rna_never_used_pct_all", allrow["rna_never_used_pct"]),
        ("dna_people_called_all", allrow["dna_people_called"]),
        ("dna_0copy_pct_all", allrow["dna_0copy_pct"]),
        ("ratio_rna_over_dna_all", allrow["rna_ratio_to_dna"]),
        ("ratio_rna_over_dna_min", s["rna_ratio_to_dna"].min()),
        ("ratio_rna_over_dna_max", s["rna_ratio_to_dna"].max()),
        ("spearman_rank_rna_vs_dna_5_groups", s["rna_never_used_pct"].rank().corr(s["dna_0copy_pct"].rank())),
        ("pearson_rna_vs_dna_5_groups", np.corrcoef(s["rna_never_used_pct"], s["dna_0copy_pct"])[0, 1]),
        ("absence_allele_freq_rna_all", allrow["rna_implied_absence_allele_freq"]),
        ("absence_allele_freq_dna_all", allrow["dna_absence_allele_freq"]),
        ("max_abs_diff_absence_allele_freq_5_groups",
         (s["rna_implied_absence_allele_freq"] - s["dna_absence_allele_freq"]).abs().max()),
        ("cum_pct_usage_at_most_first_bin", first_bin["cum_people_pct"]),
        ("first_bin_upper_usage_pct", first_bin["usage_pct_at_most"]),
        ("dna_mixture_w0copy_revio", float(mix.loc["revio", "w_0copy"]) if "revio" in mix.index else np.nan),
        ("ratio_median_5prime_genes", cod.loc[cod["position"] == "5′ of insertion", "ratio_absent_over_present"].median()),
        ("ratio_median_3prime_genes", cod.loc[cod["position"] == "3′ of insertion", "ratio_absent_over_present"].median()),
        ("dna_mixture_w0copy_sequel2e", float(mix.loc["sequel2e", "w_0copy"]) if "sequel2e" in mix.index else np.nan),
    ]
    pd.DataFrame(summary, columns=["metric", "value"]).to_csv(os.path.join(args.outdir, "summary.csv"), index=False)
    print(pd.DataFrame(summary, columns=["metric", "value"]).to_string(index=False), file=sys.stderr)
    print(comp[["ancestry", "rna_never_used_pct", "dna_0copy_pct", "rna_ratio_to_dna",
                "rna_implied_absence_allele_freq", "dna_absence_allele_freq"]].round(3).to_string(index=False),
          file=sys.stderr)
    print(cod.groupby("position")["ratio_absent_over_present"].describe().round(3).to_string(), file=sys.stderr)
    figure(comp, cum, mix, cod, args.outdir)


if __name__ == "__main__":
    main()
