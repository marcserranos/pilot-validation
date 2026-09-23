# 38 -- HLA-A "deletion" short-read cross-validation (WS2, S03 call #8)

*Cole, call #8: "that is a big result... how is that plausible?" -- script 30
(`reports/hla_popgen/30_hla_sv/`) reports HLA-A bridged-absence "deletion" at 1.83% of bridged
haplotypes (n_bridged ~= 23k) and designates HLA-A a negative control (should be ~0). No
documented whole-gene HLA-A deletion haplotype exists in the literature (only null alleles; one
85-kb HLA-B deletion family is known). Prior: artifact.*

**Verdict: HLA-A "deletions" are 82.8% contradicted by AoU short-read genotypes -> the long-read
call is overwhelmingly an assembly artifact. The upper bound on a true HLA-A deletion rate in
this cohort is ~0.3% (1.83% x (1 - 0.828)), not the raw 1.83%.**

Script: `scripts/hla_popgen/38_hla_a_deletion_validation.py`. Run on the VM (full cohort, 11,856
unrelated people) via the S03 VM channel; see that script's docstring for a note on why the
dependency functions from scripts 16/24/30 are reproduced verbatim in this file rather than
dynamically imported (the VM's checked-out repo was on an unrelated branch missing those files).
Figure re-rendered locally from the committed aggregate tables here (`_viz_common.nature_style`).

## Q1 -- rerun on Table 1, biallelic test, within-person concordance, quality correlates

Bridged-deletion rerun reproduces script 30 exactly: 935,929 bridged absences, 3,032,899
unbridged, 23,699 haplotypes seen, 11,856 unrelated people -- confirms this script's independent
recomputation of `bridged_absences` matches the original.

**Per-ancestry HLA-A rate** (Wilson 95% CI; `hla_a_rate_by_ancestry.tsv`): overall 1.83%
(1.67-2.01%), ranging 1.46% (EUR) to 2.03% (AFR) -- overlapping CIs, no ancestry stands out.

**Biallelic (both-haplotype) loss** (`summary.json` `q1_hwe`): among the 11,420 people with both
HLA-A haplotypes bridged, 3 have both called deleted vs. **3.84 expected** under N*p^2
independence (p=0.018345 per haplotype). Poisson/binomial exact P = **1.0** -- observed is not an
excess over the independent-artifact null.

**Within-person hap1<->hap2 concordance** (`q1_within_person_concordance`): 2x2 table
[[11004, 269], [144, 3]], odds ratio **0.85**, Fisher P = **1.0**. No evidence of a shared
per-person cause (e.g. one bad assembly run) -- consistent with independent per-haplotype
artifacts, which is exactly what would ALSO produce an N*p^2-matching biallelic count. This is
why the HWE test alone is not sufficient evidence either way; see Q2 for the test that actually
distinguishes artifact from biology.

**Assembly-quality correlates** (`q1_quality_spearman`, Spearman rho vs. per-person
haplotype-deletion count 0/1/2): bridging-contig span (rho=-0.055, P=0.27, n=416) and gene count
(rho=-0.055, P=0.26, n=416) are not significant. `template_warning_rate` (rho=0.052,
P=2.5e-8), `n_genes_called` (rho=-0.058, P=7e-10) and `mean_template_distance` (rho=0.092,
P=1e-22) are nominally significant only because n=11,420 is huge -- all |rho|<0.1, i.e. no
practically meaningful relationship. Platform (`hla_a_rate_by_platform.tsv`): revio 3.53% vs.
sequel2e 4.77% of people with any HLA-A deletion call -- a mild difference, not a smoking gun (no
platform concentrates the artifact). `trim_tier` has a single value (`paf_region`) for this
cohort, so it is not an informative axis here.

## Q2/Q3 -- short-read contradiction test, calibration panel across all script-30 control genes

Core logic: for a person with exactly one LR HLA-A-deleted haplotype ("hemizygous carrier"), a
true deletion means short reads should look hemizygous/homozygous at that locus. If AoU short-read
genotyping instead reports **two distinct alleles**, the gene demonstrably exists on both
haplotypes -> the long-read call is wrong. `sr_calibration_panel.tsv` (n_used = ancestry-pooled
carriers/non-carriers with an SR genotype; Fisher exact carrier vs. non-carrier heterozygosity):

| gene | positive control | LR deletion % | SR-contradiction % (carriers) | non-carrier SR het % | Fisher P |
|---|---|---:|---:|---:|---:|
| **A** | no | 1.83 | **82.8** (n=412) | 91.5 | 4.2e-08 |
| B | no | 0.60 | 77.5 (n=138) | 95.9 | 1.8e-14 |
| C | no | 0.03 | 37.5 (n<20) | 91.5 | 2.0e-04 |
| DRA | no | 1.63 | 0.0 (n=367) | 0.0 | 1.0 (uninformative -- see caveats) |
| DQA1 | no | 1.18 | 72.9 (n=269) | 82.2 | 2.0e-04 |
| DQB1 | no | 2.00 | 80.4 (n=449) | 85.9 | 1.9e-03 |
| DPA1 | no | 1.09 | 45.1 (n=246) | 50.5 | 0.107 |
| DPB1 | no | 1.00 | 58.2 (n=225) | 81.1 | 8.1e-15 |
| DRB1 | no | 1.65 | 83.3 (n=360) | 93.2 | 4.1e-10 |
| DRB3 | **yes** | 49.3 | **0.82** (n=3777) | 54.4 | <1e-300 |
| DRB4 | **yes** | 69.6 | 0.0 (n<20) | 0.0 (n<20) | 1.0 (thin) |
| DRB5 | **yes** | 83.1 | **0.33** (n=2425) | 37.7 | 1.96e-104 |

**Multiple-testing correction:** 12 Fisher exact tests are reported in this panel (one per
script-30 control gene). Bonferroni-correcting across all 12 (alpha=0.05/12=0.00417): every
reported test still clears the corrected threshold **except DPA1** (raw p=0.107, already
non-significant uncorrected) and, trivially, the two uninformative rows (DRA and DRB4, both
Fisher P=1.0 -- see caveats below). This does not change the calibration story: HLA-A (p=4.2e-08)
and every other negative-control gene tested remain significant after correction, as do all three
positive controls.

**This is the calibration Cole asked for.** The three positive controls (DRB3/4/5 -- textbook,
DRB1-haplotype-group-determined absences) are SR-contradicted **<1%** of the time: short reads
overwhelmingly agree these genes are really gone. Every negative-control gene tested, including
HLA-A, is SR-contradicted **37-83%** of the time: short reads overwhelmingly disagree with the
long-read "deletion" call. HLA-A sits squarely with the other known-artifact genes (B, DQA1,
DQB1, DRB1 all >70%), not with the true deletions. Panel (a) of `fig_sr_validation` plots this
calibration directly.

**Fragmentation/misplacement check** (`summary.json` `fragmentation_evidence_hla_a`, Q2's second
ask): of 413 HLA-A hemizygous carriers checked, HLA-A is annotated on **more than one contig**
(any hap) in **123 (29.8%)** -- i.e. nearly a third of "deleted" calls have direct evidence in
their own assembly that the gene sequence exists elsewhere and was simply not placed on the
bridging contig. This is independent, mechanistic evidence for assembly
fragmentation/misplacement, not a true deletion.

## Q4 -- mechanism

No single mechanism fully explains the artifact (see Q1's quality correlates -- all weak) but
several pieces point the same direction: (1) panel (c) of the figure shows the bridging-contig
gene count for confirmed HLA-A "deletion" calls spans a wide range (2-30 genes, median 12) with
no concentration at the shortest/sparsest contigs -- so it is not simply "HLA-A gets deleted only
on badly fragmented contigs"; (2) `template_warning_rate` and `mean_template_distance` both
correlate weakly-but-significantly *positively* with deletion status (more warnings/divergence ->
slightly more likely to see a "deletion"), consistent with assembly/annotation noise as a
contributing factor without being the whole story; (3) the 29.8% direct fragmentation-evidence
rate (above) is the strongest single mechanistic signal recovered.

## Figure: `fig_sr_validation.{png,pdf}`

Nature-grade 3-panel supplement, rendered from the committed aggregate tables in this folder
(`sr_calibration_panel.tsv`, `summary.json`, `hla_a_deleted_contig_n_genes_hist.tsv`) via
`_viz_common.nature_style()`.

- **(a)** Calibration scatter -- LR bridged-deletion rate (x) vs. SR-contradiction rate in
  carriers (y), one point per script-30 control gene. HLA-A (orange, black edge) sits with the
  negative-control cluster (grey, top-left, high SR-contradiction), far from the positive
  controls (green, bottom-right, near-zero SR-contradiction).
- **(b)** Biallelic HLA-A loss: observed (3) vs. expected under independence (N*p^2 = 3.84).
- **(c)** Bridging-contig gene count for confirmed HLA-A "deletion" calls (n=426); hatched grey
  bars are small-cell-suppressed (`<20`, drawn at a nominal placeholder height, never as 0 or a
  real value) per the project's AoU disclosure rule.

## Caveats

- "Deleted" = script 30's strict *bridged absence* (a single contig spans the gene's position and
  the gene is not there), not raw absence -- see `30_hla_structural_variation.py`'s docstring.
- DRA's 0% SR-contradiction is **not evidence of a real deletion**: DRA's non-carrier SR
  heterozygosity rate is also 0%, meaning AoU's SR genotyping essentially never calls DRA
  heterozygous at all (very low diversity / poor SR typing) -- the contradiction test has no
  power for this gene. C and DRB4 have thin carrier counts (n<20) and are reported but should be
  read with wide uncertainty.
- The HWE/biallelic test (Q1) cannot by itself distinguish "no biology" from "independent
  per-haplotype artifacts," which is exactly why the SR-contradiction test (Q2) is the load-bearing
  result here, not the HWE test.
- All counts 1-19 are written as `<20` per the AoU small-cell rule; parsed as censored, never 0.
- Every number here is an aggregate computed on the VM from the full 11,856-person unrelated
  cohort; no participant-level row left the VM.

## Distilled

HLA-A "deletions" (1.83% of bridged haplotypes) are contradicted by AoU short-read genotypes in
**82.8%** of carriers -- comparable to every other genuinely-implausible negative-control gene
(B, C, DQA1, DQB1, DPB1, DRB1, all 37-83%) and starkly different from the true positive-control
deletions DRB3/4/5 (<1% contradicted). Roughly 30% of HLA-A "deletion" calls additionally show
direct evidence the gene sequence exists elsewhere in the same person's own assembly. No single
assembly-quality proxy fully explains it (all Spearman correlations are weak), but the pattern is
consistent with assembly/placement noise, not biology. **Recommended framing for the paper: true
HLA-A deletion rate < ~0.3-0.4% (upper bound); the 1.83% raw bridged-absence rate is not usable
as a biological estimate.**

---

## Supplement: restyled deletion-rate figure (script 38b)

The section below is unchanged from the companion agent's deliverable (`38b_deletion_supplement_fig.py`,
re-rendered from the already-committed `reports/hla_popgen/30_hla_sv/` tables, no VM access) --
kept here rather than duplicated, per the WS2 brief ("your README should link them").

### 38 -- deletion supplement figure (HLA-A / MHC gene deletions, style pass)

Source: script `38b_deletion_supplement_fig.py`, re-rendered from the already-committed tables
in `reports/hla_popgen/30_hla_sv/` (`deletion_rates_by_gene.tsv`,
`deletion_rates_by_gene_ancestry.tsv`). No VM access, no re-computation from raw calls -- every
number here was already public.

Cole (call #8): "simplify the colour scheme -- drop the orange/green distinction, keep
everything else." Script 30's original `fig_deletion_rates()` colour-coded genes green
(positive control) / grey (negative control) / orange (everything else). This version uses one
neutral colour for every gene; the two control groups are called out with a bracket + label
above their bars instead.

**Panel (a)** -- per-gene bridged-deletion rate, Wilson 95% CI, grouped by SCHEMA.md's
`gene_class`, controls marked by bracket + label rather than colour.

**Panel (b)** -- CNV-gene (DRB3/4/5/C4A/C4B) deletion frequency by ancestry, dot-plot with Wilson
95% CI, project ancestry palette.

Figures: `supp_deletions.png`, `supp_deletions.pdf`. Full detail in that script's own inline
documentation (see `38b_deletion_supplement_fig.py`'s module docstring).
