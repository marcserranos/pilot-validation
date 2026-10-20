# S08 case browser log (append-only; isolated from the S06 branch)

- [10-20] PLAN — Marc: stop working at the level of aggregate error rates, look at the bases. Built `80_case_browser*` (core,
  extract, viewer) + tests. Alignments are Immuannot's own (PAF cs replay); no yak/trioeval run; VM use = one read-only
  extraction into `~/case_browser`, deployed from a throwaway clone of branch `s08-case-browser`.
- [10-20] NOTE — a broad `git add` on the shared S06 branch swept early WIP of `80_*` into S06 commits c97667a / 1b936ef.
  Work moved to branch `s08-case-browser` (own worktree) to stay isolated; those two S06 commits hold stale copies.
