# S06 log (append-only)

- [10-20] PLAN — sprint opened; scope, truth sets, metric list and queue in SPRINT.md. Reused: 56 (call score inputs, pairs),
  57b (per-call parent of origin), 58 (control depth, bedcov runner), 67 (pair/trio records, pruning curve idea).
- [10-20] FINDING (inspection) — no exon- or CDS-restricted depth exists in 58/59/67 (whole gene or window only). The exon
  coordinates are available on the VM (GENCODE v46 GTF; Immuannot allele reference with exon= / CDS=).
- [10-20] QUEUE — E7 added (Marc): trio-truth version of the gene-pair matrix; contig, not "physical", in titles.

- [10-23] RESTRUCTURE (Marc: finish sprint, nuanced QC, KIR too) — QC split into three outcomes. 76 (A missing call; B wrong allele among fully
  called child + fully called transmitting parent with 3-way class match / other-parent-only / neither; C novel). 76b pruning evidence
  (0.25% steps, exact bands, enrichment, NNP, family-grouped CV, transfer, power). KIR: 77 call table, 77b yak origin (174 trios,
  3,274 windows, 21% ambiguous), 77c trio/pair records (protein level), 78 chr19 exon/CDS depth (10 of 18 genes on GRCh38), 78b KIR pruning.
  Final report: scripts/hla_popgen/79_s06_final_report.py (local, git-excluded output reports/hla_popgen/79_s06_final/).
- [10-23] FINDING HLA — among fully called genotypes the true wrong-allele rate is ~0.3-0.4% (classical-8: 5 neither + 2 other-parent-only of 1,873;
  71 records: 8/2,025 with both people fully called vs 30/74 with a one-call parent). The earlier 1.7-1.9% 'trio error' was missing calls.
  Stream B has ~5 events: minimum detectable AUC 0.82, so 'nothing predicts' means 'not detectable'. Stream A (missing call, 158,730
  person-gene units, 12,531 incomplete): AUC 0.72 for MHC contigs and CDS bases<10x; 5% worst by cds_share_lt10 removes ~23% of incomplete (4.7x), 20% ~50%.
  Pairs: cds_share_lt10 AUC 0.81, x_opt ~6%.
- [10-23] FINDING KIR — protein-level trio 'neither' 128 of 4,886 named records; 85 of them sit on transmitting parents with an unnamed call (27.1% vs 0.9%
  with named parent): label instability, so the B_worst AUC 0.90 (tier_rank, cds_distance) is a label artefact; true KIR rate ~0.9% (43 events). Sequence-level
  CDS comparison is the definitive test (open).
- [10-23] OPEN: yak origin is per gene call (+-10 kb), no haplotype-level consistency; KIR 21% ambiguous windows; GRCh38 depth missing for 8 KIR genes.
