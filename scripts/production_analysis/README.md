# Production-run post-processing & supervisor report

Three scripts, built 2026-08-10 after the full-cohort Immuannot production run finished, to answer
Marc's ask: completeness/coverage/demographics, a confidence-threshold discrepancy sweep vs
AoU-native, a confidence distribution plot, and HLA-vs-ancestry clustering. All are cheap I/O over
small TSVs + per-person `.gtf.gz` files — **run these on the resized, cheap Workbench environment**
(`scripts/monitoring/README.md`'s "resize the compute" note), not the 96-core production VM.

## Run order

```bash
cd ~/repos/pilot-validation && pixi shell -e spechla   # pandas/matplotlib/sklearn/umap all live here
pixi install -e spechla   # first time only, to pick up umap-learn (added 2026-08-10)

python3 scripts/production_analysis/analyze_completeness_and_demographics.py
python3 scripts/production_analysis/analyze_confidence_vs_aou_native.py
python3 scripts/production_analysis/cluster_hla_by_ancestry.py
```

Each is independent (different inputs read, different output dir) — run them in any order, or
just the ones you need. Defaults assume the standard `~/pipeline_outputs/` layout from
`scripts/production_orchestrator/RESULTS_LOCATION.md`; override with `--calls`/`--cohort`/
`--timing`/`--aou-tsv`/`--outroot` if your paths differ.

## What each produces

| Script | Output dir | Figures | Answers |
|---|---|---|---|
| `analyze_completeness_and_demographics.py` | `production_analysis/completeness/` | `completeness_overview.png` (6-panel), `timing_stats.png` | Did the pipeline capture data correctly? Who got covered, by ancestry and platform? Basic operational stats. |
| `analyze_confidence_vs_aou_native.py` | `production_analysis/confidence/` | `baseline_concordance_by_gene.png`, `confidence_threshold_sweep.png`, `confidence_distribution.png` | How does AoU-native compare to Immuannot at full scale? How does that discrepancy change as confidence tightens, per gene? What does the confidence distribution look like per gene (mean/median/clustering)? |
| `cluster_hla_by_ancestry.py` | `production_analysis/clustering/` | `pca_pc1_pc2_by_ancestry.png`, `umap_by_ancestry.png` | Does HLA genotype structure correlate with ancestry under PCA/UMAP? |

All write a markdown report alongside the PNGs (`*_report.md` / `completeness_report.md`) with the
numeric tables behind each figure — paste these into the supervisor report directly, they're
already aggregate-only prose+tables, same discipline as every other analysis script in this repo.

## Design notes worth knowing before reading the figures

- **The confidence-threshold sweep chart's x-axis is reversed on purpose** — 0/strict sits on the
  right, so "raising the certainty threshold" reads left-to-right as discrepancy falls, matching
  how Marc described wanting to read it. Read the axis label, not just the numbers.
- **The sweep's discrepancy floor (right edge) is an estimate of AoU-native's real error rate
  against high-confidence ground truth** — not the discrepancy at threshold=0 alone; look at where
  the curve visibly levels off, same interpretation as the existing n=60 confidence-matched-truth
  work (`scripts/analyze_confidence_matched_truth.py`), now at full production scale.
- **The clustering script's PCA/UMAP scatter plots are one dot per real person** — no identifiers on
  them, but structurally person-level in a way most of this project's outputs deliberately aren't.
  This is standard AoU Workbench research-figure practice (AoU's own ancestry PCA plots are built
  the same way) — fine to build and view inside the Workbench. If either figure is ever shared
  outside the sponsor/supervisor review loop, treat it like any other research output under the
  still-open compliance question in `context/DECISIONS.md`.
- **The clustering script only uses complete-case people** (real calls at all 8 classical genes,
  both haplotypes) — no imputation. The report states exactly how many people that excludes; if the
  excluded fraction is large, that's itself a data point about the run's completeness (cross-check
  against `analyze_completeness_and_demographics.py`'s numbers), not a clustering-method problem.
- **`trim_tier == self_align_needed` (991 sequel2 people) will show ~0% completion everywhere** —
  Phase 2 was disabled for this production launch after Tier 3 failed testing
  (`context/EXPERIMENTS.md`, 2026-08-05). Expected, not a bug — the completeness script's markdown
  says this explicitly so it isn't mistaken for a problem when reviewing the report.

## Ideas not built yet — flagged, not implemented, in case they're wanted

Raised while designing the above, in the spirit of "what else would be worth showing a supervisor":

1. **Allele-frequency spectrum by ancestry** (top-K alleles per gene, faceted/stacked by ancestry
   group) — directly serves this project's Aim 1 (`context/TASK_CONTEXT.md`, ancestry-stratified
   HLA allele frequency) and hasn't been done yet at full cohort scale, only smaller pilots.
2. **A capstone "DRB1 is the hardest locus" summary figure** pulling together the six independent,
   already-documented lines of evidence for this (`context/DECISIONS.md`) into one compact visual —
   this project's single most-repeated finding, never presented as one unified figure.

Both are readily buildable from data these three scripts already load — say the word and either
can be added as a fourth script rather than folded into the existing ones (keeps each script's
output/runtime scoped to what it's already named for).

## Compliance reminder

Nothing in this directory downloads or exports raw per-person calls — every script here only ever
writes aggregate counts, rates, distributions, and (for the clustering script specifically)
unlabeled per-person points in a derived embedding space. Real allele calls and bare person_ids
never leave `~/pipeline_outputs/` or get written into any file under this directory. Same standing
rule as `scripts/production_orchestrator/RESULTS_LOCATION.md`.
