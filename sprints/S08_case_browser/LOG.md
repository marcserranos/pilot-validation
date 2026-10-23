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
