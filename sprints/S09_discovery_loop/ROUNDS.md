# S09 think-tank rounds (Haiku proposals — hypotheses, not evidence)

## Round 1 (S08 scope, 5 personas) — consensus order
Transmission oracle (across platforms: same-pipeline errors transmit), simulated reference-bias null, RNA-seq ASE on personal references, primate outgroup (can only confirm antiquity),
per-run error barcode, mark-recapture error estimate. Self-corrections: trans-species "no match = artefact" was wrong; ESM scores plausibility, not selection; NKG2A "tone" from bulk RNA is a proxy.

## Round 2 (2026-10-24, 9 agents, modes: moonshot x6 themes, secure, naive, heuristic) — headline ideas
- **LR truth vs 535k native short-read (calib):** concordance by gene x ancestry; error vs divergence; 'hallucinated alleles' over deletions/nulls (HLA-A del, DRB4/5, C4); hidden homozygosity; misclassification-corrected HLA-disease ORs; per-call error-bar model; C1/C2 flips distorting KIR-ligand tests.
- **KIR x HLA x disease:** KIR3DL1 expression class x Bw4 dose vs infection; HLA-F T237I/P205T x KIR3DS1 (KIR2DS4 negative control); KIR content vs HLA-E/F across ancestry; HLA-A deletion vs KIR content; 'frustration index' (expression-weighted missing self); maternal KIR x fetal C2.
- **Strange correlations:** EHH/balancing halo at HLA-F haplotypes; gene-conversion breakpoint atlas; allele x lab-value pleiotropy network; C4 copy-number residuals vs DR and SLE; survivor haplotypes by age; parent-of-origin ASE; pathogen-epitope overlap; hypersensitivity (B*57:01/B*15:02) background dependence.
- **Functional:** HLA-F ASE; non-coding burden vs ASE; HLA-G isoform by haplotype; divergence vs NK/T composition; AF/ESM of alpha3 variants; NMD of stop-gained alleles.
- **Repertoire:** exact-allele public TCR clones; amino-acid divergence of the two alleles vs infection codes; groove electrostatics vs TRBV; TRBV x DR-DQ block.
- **Secure:** native benchmark; reference-bias dose-response; phased MHC haplotype atlas; trio Mendelian error rates; IPD submissions (F T237I/P205T); KIR gene-content frequencies with short-read false-negative rate.
- **Heuristic engines:** allele x phecode PheWAS (~540k tests, max-stat permutation, negative controls); residue-level scan conditioned on allele; haplotype-window scan; anomaly detection on divergence features; KIR x HLA epistasis mining; clustering + bandit with coalescent null.
- **Naive looks (10):** template-distance histograms by gene/platform; homozygosity per person; KIR gene-count outliers; novel burden vs divergence per ancestry; C4 CN vs DR/DQ; concordance heatmap; DRB1 vs DRB3/4/5 co-occurrence; warning rate vs sequencing date; related-pair sharing vs IBD; haplotype block-length distribution.

Next: critic + literature review (stronger model, web search) -> portfolio -> MVPs (ledger rows first).

## Critic + literature (Sonnet, web search; items marked [mem] unverified) — portfolio
Secure: S1 LR-truth benchmark of native calls (no AoU LR-vs-SR concordance paper found); S2 trio Mendelian error per gene; S3 IPD submission of F T237I/P205T (>=3 unrelated carriers, raw-read check, independent trio).
Exploratory: E1 personal-reference ASE at HLA-F carriers (Aguiar 2019 did GEUVADIS); E2 T237I x KIR3DS1 with KIR2DS4 control (KIR3DS1 LD with 3DL1 structure, few carriers); E3 divergence-scaled heterozygote advantage (Chowell 2018 — low novelty).
Moonshots: M1 non-coding/UTR burden -> ASE ('invisible eQTL'); M2 expression-weighted KIR-ligand 'frustration index'. Automated: A1 PheWAS (>=50 negative-control phenotypes; any negative-control hit = miscalibrated); A2 haplotype-window scan recovering the two known controls.
Cut: hallucinated-deletion claims (subset of S1), C4 CN, naive looks as standalone claims (keep as QC plots).
Traps: (1) circularity/ancestry — F alleles on one haplotype = haplotype tag until conditioned; (2) multiple testing — permutation not Bonferroni, suppressed != zero; (3) short-read error is allele/ancestry-differential, LR truth is not perfect (report trio LR-vs-LR as ceiling).
Papers: Garcia-Beltran 2016 + Burian 2016 (HLA-F/KIR3DS1), Aguiar 2019 (HLA ASE), Immuannot 2024, HLA*LA 2019, UKB HLA PheWAS 2020 [mem], Chowell 2018 [mem], Emerson 2017 [mem], Martin 2002 (KIR3DL1/Bw4/HIV) [mem].
