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
