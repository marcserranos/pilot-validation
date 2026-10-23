# S09 — discovery loop (think tank → falsifiable MVP → critics/literature → think tank …)

Branch `s09-discovery-loop` (worktree `../pilot-validation-s09`), forked from S08 (`s08-case-browser`) on 2026-10-24. S08 answered "are the substitution-only
novel HLA calls real?" (report: `reports/hla_popgen/91_novel_recurrence/S08_report.html`, scripts 91-96). S09 is open-ended exploration of what the full dataset
can discover. New scripts take numbers **100+**.

## The loop
1. **Think tank** (Haiku personas, several MODES — see below) proposes ideas, each with mechanism, data, <1-day VM MVP, kill criterion.
2. **Critic + literature** (stronger model, web search) dedupes, checks novelty, builds a portfolio.
3. **MVP experiments** run on the VM (aggregates only, `<20` rule), results appended to the LEDGER.
4. Results feed the **next think-tank round** (orthogonal directions, critics, more literature) until something is publishable.

## Modes (marc, 2026-10-24: "some for the cover headline, some for secure claims, some naive exploration, some heuristic search")
- **moonshot** — cover-level headline, mechanism-first.  **secure** — near-certain, citable resources/benchmarks.
- **naive** — "just look" tables and plots, no hypothesis.  **heuristic** — automated scans with permutation nulls and negative controls.

## Rules
- Participant-level data never leaves the Workbench or enters the repo; aggregates only; counts 1-19 are `<20`. Local result tables are not committed without Marc's OK.
- Every experiment gets a LEDGER row BEFORE it runs (hypothesis, mode, data, kill criterion, multiple-testing family); results are appended, never rewritten.
- Negative-control phenotypes and known-allele positive controls in every test (S08 lesson: G*01:01:33, F*01:25 behave like real alleles; Sequel II/IIe indels do not).
- VM is shared: one core, `nice -n 10`, detached jobs with checkpoints, throwaway clone, clean up (see S08 `AGENT_BRIEF.md`).
- Think-tank outputs are Haiku proposals, not evidence; literature claims need the critic's web check.

Files: `LEDGER.md` (hypothesis ledger), `ROUNDS.md` (think-tank digests), `LOG.md` (append-only run log).

## Red-team checklist (adversarial reviewer, 2026-10-24) — every experiment must pass before we believe it
1 positive control recovered at the same power; 2 negative control / permuted labels null; 3 chance-level baseline reported beside every concordance (e.g. modal-G-group predictor);
4 one-per-kinship-group recount; 5 within-ancestry tests with interaction test before pooled claims; 6 platform/centre/depth strata; 7 read-level support for novel calls (MAPQ, strand, homopolymer, phase);
8 an independent non-IPD truth on a subset (HPRC / second annotator); 9 pre-registered family + LD/ancestry-preserving permutation null; 10 kill criterion + power statement, underpowered cells '<20' never zero.
Open follow-ups from the review: modal-baseline for the G-group metric (S09-01); second annotator on ~200 people; ordered-haplotype scoring; compare native allele frequencies in the 12k vs 535k; unique-MAPQ read support at HLA-F T237I/P205T sites.
