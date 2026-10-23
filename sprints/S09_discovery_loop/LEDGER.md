# S09 hypothesis ledger (pre-register BEFORE running; append results, never rewrite)

| id | mode | hypothesis | data | MVP / readout | kill criterion | status | result |
|---|---|---|---|---|---|---|---|
| S08-A | secure | HLA-F T237I/P205T are real inherited alleles (single flanking haplotype, no platform/centre skew) | LR 12.2k, 95 | flanking-background effective number vs matched controls | carriers >= control median | PASSED (S08) | 1.1-1.4 vs 7-8 effective backgrounds; controls F*01:25, G*01:01:33 identical pattern |
| S08-B | moonshot | HLA-F T237I/P205T expressed ~50% in heterozygous carriers (personal-reference RNA-seq ASE) | 27 het carriers with RNA-seq | variant-haplotype read fraction, bias-corrected | zero variant reads or imbalance beyond bias | RUNNING (96) | pending |
| S09-01 | secure | native short-read HLA calls vs LR template-distance-0 truth: concordance by gene x ancestry; falls with haplotype divergence; error structure (confusions, hidden homozygosity) | 535k native + 12k LR overlap | script 100: concordance by gene x ancestry x td tier, confusion map | flat within ancestry | RUNNING (100) | pending |
| S09-02 | secure | misclassification-corrected HLA-disease ORs differ from native (DQA1/DQB1 vs T1D/celiac) | native 535k, EHR, confusion matrix from overlap | OR native vs corrected | CIs agree for all alleles | proposed | |
| S09-03 | moonshot | HLA-F T237I/P205T x KIR3DS1 epistasis vs chronic-viral EHR codes (KIR2DS4 negative control) | LR KIR + HLA + EHR | interaction OR | CI includes 1 / cells <50 | proposed (power first) | |
| S09-04 | moonshot | HLA-A deletion shifts KIR gene content on the same people (NK ligand-dosage compensation) | LR SV calls x KIR | CMH per KIR gene, conditioned on C1/C2/Bw4 | no gene at q<0.05 | proposed | |
| S09-05 | heuristic | allele x phecode PheWAS with permutation null + negative-control phenotypes | native 535k, LR 12k | top HLA-prior phecodes first | no replication across ancestries | proposed | |
| S09-06 | naive | the ten 'just look' passes (template distance by gene/platform, KIR gene-count outliers, DRB3/4/5 co-occurrence gaps, ...) | LR | tables/plots | n/a (leads only) | proposed | |
| S09-07 | moonshot | parent-of-origin effect on HLA allele-specific expression | ~1000 trios x RNA | ASE ~ parent-of-origin | CI includes 0 after eQTL matching | proposed (check overlap) | |
| S09-08 | moonshot | balancing-selection halo (EHH) around HLA-F signature haplotypes, modulated by KIR3DS1 | LR phased haplotypes | EHH percentile vs matched random haplotypes | <=50th percentile in AFR and EAS | proposed | |
