# 08 — Methods audit: filter steps, top-500 ties, BAM cap, environment

Script: `scripts/08_methods_audit.py` (plus `freeze_env.sh`). Run 2026-10-04 in the polish run
(5b83486); BAM-cap part rerun 2026-10-06 with RNA-SeQC2 metrics (b38b48f). No figure.

## 1. What the shared clonotype filter removes (`filter_steps.csv`)

All 7,922 people, 51.3 M `cdr3.out` rows:

| Step | Rows kept | % of previous |
|---|---|---|
| CDR3_score ≥ 0.02 | 42.4 M | 82.7 |
| canonical junction C…F/W | 39.1 M | 92.2 |
| length ≤ 30 aa | 39.1 M | 99.96 |
| unique clonotypes, all chains | 34.8 M | 89.1 |
| TRB | 9.44 M | 27.1 |

The median person ends with 1,137 TRB clonotypes.

## 2. Ties at the top-500 cut-off (`tie_summary.csv`, `trbv_pool_shift.csv`)

- 7,340 people have more than 500 TRB clonotypes. For 98.8% of them the cut-off falls inside a
  tie: 78% at 1 read, 20% at 2.
- Median tied block: 625 clonotypes, of which 199 are kept, so ~200 of 500 slots were decided
  by the tie-break.
- The old tie-break was alphabetical. It favoured early TRBV names: TRBV12-3 1.62×,
  TRBV10-3 1.77×, TRBV11-2 1.72× over-represented; TRBV7-9 0.65× and TRBV7-2 0.68× under.
  Total pool composition shift: 13.6 points.
- The seeded random tie-break changes 19% of each person's top 500. The old pool was
  reproduced exactly before the change (0 differences), so the comparison is clean.

## 3. The 12 GB BAM cap (`bam_cap_qc.csv`, `bam_cap_by_ancestry.csv`)

- 404 people in the long-read × RNA-seq overlap were excluded for BAM size (4.9%), including
  48.5% of the top depth decile. They are the deep tail:
  - mapped reads median 277 M vs 109 M kept (2.5×);
  - library complexity 2.1×;
  - genes detected 35,569 vs 30,626;
  - mapping rate 86.3% vs 90.9%;
  - duplicates 44.6% vs 39.8%;
  - all p < 1e-3.
- Excluded by ancestry: EAS 3.2%, SAS 4.2%, EUR 5.0%, AMR 5.1%, MID 5.4%, AFR 5.9%.
- Results describe typical-depth samples. TRUST4 was never run on the 404; including them
  would take ~4–6 h on one VM.

## 4. Software (`environment.txt`)

- TRUST4 v1.1.5-r573 (bioconda), reference FASTAs from repo commit 4032f90.
- STAR 2.7.11b in AoU's BAM headers.
- sceptr 1.2.0, olga 1.3.0, scikit-learn 1.9.1, umap-learn 0.5.12, cnsplots 0.7.0, torch
  2.14.0, pandas 2.3.3, scipy 1.18.0.
- numpy is present twice (conda 2.5.3, pip 2.4.6). trust4 is pinned in `pixi.toml`.
