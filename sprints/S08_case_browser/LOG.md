# S08 case browser log (append-only; isolated from the S06 branch)

- [10-20] PLAN — Marc: stop working at the level of aggregate error rates, look at the bases. Built `80_case_browser*` (core,
  extract, viewer) + tests. Alignments are Immuannot's own (PAF cs replay); no yak/trioeval run; VM use = one read-only
  extraction into `~/case_browser`, deployed from a throwaway clone of branch `s08-case-browser`.
- [10-20] NOTE — a broad `git add` on the shared S06 branch swept early WIP of `80_*` into S06 commits c97667a / 1b936ef.
  Work moved to branch `s08-case-browser` (own worktree) to stay isolated; those two S06 commits hold stale copies.
- [10-20] RESULT (stratified sample of 316 cases; NOT population rates; counts only) — all 249 ladder alignments replay Immuannot's
  PAF cs exactly. HLA calls at <=10 CDS edits called "novel-protein" (n=145): 138 carry a frameshift, 94% of indels are 1 bp,
  ~half in homopolymer runs of 3-5; at ~50 edits the pattern flips to substitutions. Indels do not recur at fixed loci (462 distinct
  sites of 486 events, 23 recurrent), are uniform along the CDS, mostly G/C, deletions > insertions (263 vs 195), and 41/140
  indel-carrying cases restore the frame with a second indel. 42/145 cases are homopolymer/tandem indels only. Mendelian (30
  sampled mismatches): 17 within 6 edits of a parent haplotype (child's own call 0 edits in all 17), 12 at >=20 edits.
  Reading: low-distance "novel" calls are mostly consensus/assembly indel errors, not new alleles; consistent with the S05 demotion.
- [10-20] OPS — local copy built from a user-downloaded cases.json (classifier blocked agent download: PII). VM: throwaway clone
  ~/cb_code may remain (rm -rf); ~/case_browser holds cases.json, index.html, private id map.
- [10-21] RENUMBER — moved from 80 to the 90 series (`90_case_browser*`, `reports/hla_popgen/90_case_browser/`) so S06's 70-75+ never collides;
  76-89 left free as grace. Stale `80_*` copies still sit in S06 commits c97667a/1b936ef (not mine to rewrite).
- [10-21] ADDED (local, tested on synthetic data) — `stats` (cohort-wide composition of novel calls from hla_calls_rich + public IPD CDS; no BAM),
  `extract --exhaustive-snv CAP` (every novel call with no indel), edit-type / protein-altering filters, show-substitutions-only,
  Statistics tab, depth panel, short-read panel (UI + verdict rule `sr_verdict`; extraction not yet written, needs VM go-ahead).
  Exhaustiveness so far: the 316 cases are stratified EXAMPLES (<=3 per gene x distance bucket; KIR from 300 people), not all novel calls.
- [10-21] RESULT (VM run, whole cohort + 762 cases) — stats pass: 8,015 novel HLA calls with a CDS alignment (137 not rebuildable from cs). Novel-protein
  (7,635): 86% indel-only, 8% indels+substitutions, 6% substitution-only (449); 8% carry a protein-altering substitution; 93% have a frameshift;
  at 1 edit 91% are indel-only (72% of those in homopolymers/tandem runs). Novel-synonymous (227): 100% substitution-only. Novel-beyond-CDS (153): CDS identical.
  Substitution-only novel calls are 682 (F 145, G 211 capped at 60 each in the viewer, everything else complete).
- [10-21] RESULT (short reads, 567 cases, 1,097 sites; srWGS realigned to the called haplotype, ~1.5 s/case) — indels: 79% "not seen in short reads" and 8% contradicted
  (only 4% supported); protein-altering substitutions: 60% supported, 13% half-the-reads (other haplotype call may differ), 18% not seen; synonymous 68% supported.
  Of 439 substitution-only novel cases with short reads, 281 have every checked site supported (189 with a protein-altering substitution).
  Caveat: DRB3/DRB4 have no primary-assembly locus (not checked); verdicts are recomputed at build time from stored counts.
- [10-21] OPS — VM cleaned (code clone, logs, HTML removed); kept ~/case_browser/{cases.json,sr.jsonl,cohort_stats.json,case_map.private.tsv}. Local copy in the
  gitignored case_browser/ of the worktree.
- [10-23] THREAD — recurrent "novel" DQB1*02 A57D (codon 89 incl. signal peptide; <20 shortlist calls vs DQB1*02:02:01:01): literature subagent reports IPD-IMGT lists
  DQB1*02:180:01 (protein = 02:02:01:01 except D57; named "DQB1*02:180" since release 3.43, so present in 3.55). Likely a KNOWN allele, not novel; open question is why Immuannot
  (3.55 reference set) did not match it (partial/CDS-only entry for 02:180?). Needs: align carrier CDS+genomic to HLA29352 before any "novel" wording. Subagent claims are
  unverified secondary-source (agent did not open all pages); also lists Asp57 biology (T1D-protective, P9 pocket) as speculation. Other recurrent F/E signatures still to be checked the same way.
- [10-23] QC (local, shortlist+sr.jsonl; counts <20 withheld) — single-edit recurrent signatures in HLA-F (T237I, P205T, Y106C) and HLA-E (K195R): each always vs ONE ref label, balanced strands,
  MAPQ ~60, no CpG, mostly SR-supported; F/E alleles are absent from the only local IPD FASTA (3.38), so cross-allele check NOT possible -> must be checked against current IPD (hla_nuc) for
  named alleles before any novelty claim. DQB1 A57D: 8/9-like mostly homozygous, clean; but see IPD 02:180 above. DPB1 T16M: CpG transition, base already in IPD 3.38 (DPB1*296:01) -> known/hypermutable.
  DRB1 R262T: 3 different refs, all SR-skewed -> artefact. Many-edit recurrent B/C substitutions all sit on S/Q/N-suffixed or truncated ref templates (B*56:01:01:05S, B*44:02:01:02S, C*07:402, DRB5*01:131N,
  B*27:197) -> wrong-template artefact, NOT novelty; shortlist should drop suffix-ref and >~15-edit calls. Stop-gains: all upstream of last exon (NMD expected); only a handful (F, tier A,
  SR-supported) look clean; E/G/B/C/DRB5 nonsense are mostly SR 'contradicted/not seen'. Caveat: shortlist has no person-level zygosity; agent used SR het flag.
- [10-23] RESULT (public IPD-IMGT 3.66.0 protein alignments F_prot/E_prot, parsed directly; our codon numbering = file position +30, ref residues verified at every site) —
  recurrent shortlist signatures vs CURRENT IPD: HLA-F Y106C = F*01:25 (single substitution, KNOWN); HLA-E K195R on the E*01:03 background = E*01:152 (R128G+K195R, KNOWN);
  HLA-F T237I and P205T are in NO allele of the 3.66 alignment (140 F alleles) -> the only recurrent signatures still unnamed (exon 4, alpha3, outside the groove); caveat: '*' = unsequenced
  positions in partial entries not resolved. DQB1 A57D = DQB1*02:180 (per literature agent, not yet re-derived here). Net: of the recurrent single-edit leads, most "novelty" is the 3.55 snapshot being
  11 releases old; rerun Immuannot-equivalent matching against 3.66 (or mask by 3.66 CDS) before any novelty count. Literature agent: Immuannot paper reports ~41% of assemblies' HLA/KIR seqs novel vs v3.52 (concentrated in non-classical genes/pseudogenes) -> precedent. Submission standard: IPD requires phased full exons 2-3 (class I) / exon 2 (II), isolated from the other allele; names assigned by WHO committee.
- [10-23] VM (91 recurrence, 92 burden; one core, nice 10, <1 min each; outputs in ~/s08 on the VM, clone removed, raw/carrier files mode 600, nothing copied out) —
  RESULT 91: of 8,015 novel calls, 675 are substitution-only and rebuildable -> 441 distinct (gene, ref allele, substitution set) signatures; only 49 are seen in >= 2 people, the rest are singletons.
  The recurrent ones are strongly ANCESTRY-STRUCTURED among unrelated carriers (independent of any read-level evidence): HLA-G synonymous on G*01:01:03:01 is overwhelmingly EAS (the only
  signature with >= 20 persons besides HLA-F T237I; F T237I is mostly AFR with some AMR); F P205T mostly AFR; F Y106C (= F*01:25, known) entirely AMR; DQB1 A57D (=02:180?) AFR+AMR; E K195R (=E*01:152) AFR+AMR;
  G D130H AFR+AMR. Every recurrent key sits on ONE reference allele. Trio transmission test is nearly unusable: only 827 of 12,233 LR persons have a relative in the LR cohort, so carrier children with two
  LR parents are <20 per signature (transmission evidence must come from S06's 57b trio tables or a wider kin source). Interpretation: these are real germline variants carried by ancestry groups that IPD under-represents,
  so 'novel vs 3.55' mostly measures reference bias, not discovery.
  RESULT 92 (novel-call burden by platform x ancestry, 12,233 persons): indel-bearing novel-protein calls per person in classical genes are ~2x higher on Sequel II/IIe than on Revio at every ancestry
  (platform artefact -> consistent with indel 'novelty' being consensus error), whereas substitution-only novel-protein calls per person are roughly platform-flat and differ by ancestry
  (AFR/EAS/MID about 1.5x EUR/AMR/SAS). Novel calls per person overall (incl. synonymous/beyond-CDS) AFR > EAS ~ SAS ~ MID > AMR > EUR, as expected for a EUR-biased reference.
  Caveats: person-level relatedness available for LR persons is sparse; ancestry labels are AoU predictions; the burden table was read on the VM (not committed); pub-style tables with <20 suppression sit on the VM (~/s08/*_pub.tsv) pending Marc's go-ahead to commit.
  OPEN: (1) rebuild shortlist vs IPD 3.66 so 'novel' means novel today; (2) F T237I/P205T genomic coords + gnomAD/ALFA frequency; (3) transmission via 57b trios; (4) EHR path still unprobed.
- [10-23] RESULT (IPD-IMGT 3.66 hla_nuc.fasta, 47,419 CDS; local 93 on the 337-case shortlist, VM 91 --ipd-nuc on the cohort) — shortlist: 58 now exact matches to a named allele, 9 contain a partial entry, 270 still novel;
  after also dropping suffix-flagged reference alleles, >15-edit calls and non-tier-A short-read support, 98 remain (HLA-F 36, E 15, DQA1 12; 28 in the groove). Cohort: of 441 substitution-only signatures 54 are exact and
  8 partial matches to 3.66 (of the 49 seen in >= 2 people: 15 exact, 1 partial). Names: the large EAS HLA-G synonymous cluster = G*01:01:33; F Y106C = F*01:25; E K195R = E*01:152; DQB1 A57D = DQB1*02:180:02; G D130H = G*01:58;
  G G107S = G*01:64; G K170E = G*01:66. STILL UNNAMED + recurrent: HLA-F T237I and P205T (exact CDS in no 3.66 allele), plus small clusters (G A235V, G E190A [EUR; ref is a Q-suffix allele], F T279I, F synonymous).
  So: 'novel vs 3.55' overstated discovery; the defensible novelty set is the 3.66-unmatched, tier-A, non-suffix-reference calls. VM: clone and the IPD file removed again, outputs in ~/s08/r366 (raw files mode 600).
- [10-23] REPORT + THINK TANK — `94_s08_report.py` builds reports/hla_popgen/91_novel_recurrence/S08_report.html (uncommitted by design: Marc decides on committing result tables). Five Haiku personas (evolution, immunology, ML, skeptic, wildcard)
  proposed ideas, then critiqued each other. Consensus order: (1) transmission oracle (but across platforms — child long-read vs parental short-read — since same-pipeline errors transmit), (2) simulated reference-bias null, (3) RNA-seq ASE with
  personalised reference, (4) primate outgroup (can only confirm antiquity), (5) per-run error barcode for hidden batch effects, (6) mark-recapture per-platform error. G*01:01:33 is a known allele; only its soluble/membrane isoform question remains.
  Check trio x carrier x RNA-seq overlap before committing to any of them (power may be near zero for HLA-F T237I/P205T).
- [10-24] RESULT S08-B (script 96, VM, personal two-haplotype reference + simulated bias control; carriers <20 per signature, 24 analysable runs): HLA-F T237I and P205T are TRANSCRIBED in blood — variant-haplotype read fraction ~0.54 (T237I) and ~0.65 (P205T) vs ~0.58 for the known F*01:25 control; no carrier without variant reads;
  bias-corrected median log2 variant/other +0.26/+0.56/+0.47 (control behaves the same, so no allele-specific expression effect is claimed). With single flanking haplotype (95), no platform/centre skew, concordant relatives, and expression, the two HLA-F alleles pass every test we could run without new sequencing;
  still missing for an IPD submission: raw-read check at >=3 unrelated carriers and an independent trio. S08 report regenerated with the ASE section (not committed: results tables are Marc's call).
