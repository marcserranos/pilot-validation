# HLA-disease-allele association: real diagnosis vs real carrier status

Methodology: see this script's module docstring for the pairing rationale and the Mantel-Haenszel ancestry-adjustment formula. All numbers below are counts/statistics only -- no person_id, safe to commit (unlike the two per-person files this script reads).

**8 tests, Bonferroni threshold = 0.0063.**


| pair | n_cohort | n_carriers | n_diagnosed | n_both | OR (raw, 95% CI) | p (Fisher) | OR (ancestry-adj, CMH) | p (CMH) | sig? |
|---|---|---|---|---|---|---|---|---|---|
| B*27 / ank. spondylitis | 9355 | 423 | 28 | 9 | 10.2 (4.59-22.68) | 2.33e-06 | 9.5 | 4.35e-10 | Bonferroni |
| DQ2 / celiac | 9355 | 1227 | 38 | 5 | 1.0 (0.39-2.58) | 1.00e+00 | 0.92 | 9.47e-01 | no |
| DQ8 / T1D | 9355 | 1600 | 204 | 49 | 1.55 (1.12-2.15) | 1.09e-02 | 1.56 | 1.27e-02 | p<0.05 |
| Cw6 / psoriasis | 9355 | 1339 | 223 | 43 | 1.44 (1.03-2.02) | 4.14e-02 | 1.35 | 1.13e-01 | p<0.05 |
| DRB1*04 / RA | 9355 | 100 | 253 | 6 | 2.33 (1.01-5.36) | 5.36e-02 | 2.06 | 1.50e-01 | no |
| DRB1*15:01 / MS | 9355 | 738 | 68 | 10 | 2.03 (1.03-3.98) | 6.47e-02 | 1.91 | 9.52e-02 | no |
| DQB1*06:02 / narcolepsy | 9355 | 1551 | 25 | 2 | 0.44 (0.1-1.85) | 4.15e-01 | 0.41 | 4.06e-01 | no |
| B*51 / Behcet | 9355 | 766 | 6 | 2 | 5.62 (1.03-30.72) | 8.05e-02 | 4.4 | 2.03e-01 | no |

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/13_disease_allele_association/disease_allele_association_forest.png`
