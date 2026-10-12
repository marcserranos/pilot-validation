# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-12 — S05 in progress: all scripts written (58–65), VM runs partly started

**Board:** [`sprints/S05_qc_phasing_kir_trb/SPRINT.md`](../sprints/S05_qc_phasing_kir_trb/SPRINT.md). Log: `LOG.md`; report
material: `METHODS_LEDGER.md` (one section per script, written as the work happens); reading: `READING.md`.
`main` holds everything (PRs #4–#12). The 10-06 scripts 50–57c own those numbers; new S05 scripts are 58–65.
Private meeting notes stay out of the repo (`sprints/_calls/` gitignored).

**VM state (big_run, n2-highmem-4, 4 vCPU / 31 GB).**
- Checkouts: `~/repos/pilot-validation-main` (main; 58 runs from it, don't touch until it ends), `~/repos/pv-65` (worktree for
  65), `~/repos/pilot-validation` (old branch; its `.pixi/envs/spechla` is the Python 3.8 env used for everything).
- Extra envs: `~/envs/yak`, `~/envs/mappy38` (use via `PYTHONPATH=~/envs/mappy38/lib/python3.8/site-packages`).
- Refs: `~/ref/gencode.v46.basic.annotation.gtf.gz`; `~/ref/trb/` (chr7_KI270803v1_alt, UCSC chr7 slice); 65 probes in `~/s05/results/65`.
- **58 full run started 2026-10-11 21:22 UTC** (12,233 people, 8 workers, ~0.22 people/s, ETA ~12:00 UTC 10-12). Progress:
  `tail -2 ~/s05/logs/58.log`. A keep-alive loop exits with it.
- Deletable: `~/s05/{results,internal}/65_smoke` (stopped attempt).

**Open, needs Marc:**
1. Four S05 decisions in `DECISIONS.md` (ancestry purity rule; demote classical protein-novel counts; KIR claims on complete
   haplotypes only; physical-haplotype evidence rule from trios).
2. Sequel II re-call (63): ~27 days on this app, ~3 days on a resized one — needs his go-ahead to resize after the smoke test.
3. `VM_CLEANUP_PLAN_v2.md` (S04) still pending.

**Literal next commands (VM, in this order; each script's `sprints/S05_qc_phasing_kir_trb/vm_parts/NN.sh` has the exact lines):**
1. When 58 ends: pull its PUBLIC tables (`~/s05/results/58/*.tsv`, ids check) to `reports/hla_popgen/58_regional_depth/`.
2. `git pull` in `~/repos/pilot-validation-main`; run 59 (cutoff + `exclude_people.tsv`), then 61 and 62 with the exclusion list,
   then 60 with 58's proxy.
3. 65 (`vm_parts/65.sh`): reads-mode smoke on 20 → set `--min-reads` (12 if coverage allows) → stratified 2,000 in assembly
   mode + reads mode on all → `--compare`; then rerun 60 with 65's table.
4. 64 (`vm_parts/64.sh`): IPD format checks, re-match-only run, then read support.
5. 63 (`vm_parts/63.sh`): read-only checks U1–U7, 2-person smoke; then ask Marc about the resize.

## As of 2026-10-06 — earlier S05 work (scripts 50–57c; committed 2026-10-11 on `s05-oct06-work`) (weekly update; 3 axes: QC/understanding, small experiments, strategy)

Access restored (AoU deactivation of 09-28/29 resolved). Branch still `s04-kir-recurrence-style-share` (uncommitted S05 work
on top of an uncommitted call-note scrub from 10-05 — commit both deliberately). Private report for Marc:
`reports/walkthrough_S05/index.html` (gitignored). Scripts 50–55 added; results in `reports/hla_popgen/5x_*/README.md`;
run log in EXPERIMENTS.md 2026-10-06.

**Headlines:** DQ G1/G2 claim survives an ascertainment test (exceptions = structurally anomalous contigs only); classical-HLA
protein-novel counts do not survive recurrence + homopolymer filters (demote pending read-level validation); KIR framework
deficit = fragmentation + contracted haplotypes; DRB5 presence/absence validated by RNA dosage; AMR > AFR richness explained by
admixture pooling; ASE with personal references is technically viable (53 pilot).
Evening: 56/56b relatives benchmark (confidence score works; A~B physical haplotypes unreliable) and 57/57b yak trio
benchmark (50 trios: class II switch/hamming ~0.2%, 0 parent mismatches for DQ/DR/DP same-contig pairs; class I contigs mix
parents; Mendelian 2-field typing mismatch 3.2%). yak lives in `~/envs/yak` (VM); network is open (bioconda/pypi/github).

**Literal next commands:**
1. Marc: decide exclusion list / purity rule, KIR complete-haplotype rule, novel-protein demotion (DECISIONS.md open items, S05).
2. VM: 53 results → pull `~/s05/results/53/ase_by_gene_class.tsv` (PUBLIC) if not yet pulled; INTERNAL per-person table stays at
   `~/s05/internal/53/`. Then scale ASE to all LR∩RNA on a resized app (16–32 vCPU; resize needs Marc's authorisation).
3. VM: fix Tier-3 self-align trim (991 sequel2, mostly AFR) and rerun HLA+KIR for them.
4. Novel alleles: the VM can download directly (network open) → fetch a current IPD-IMGT/HLA release and re-match; add HiFi
   read / short-read support per novel allele.
5. 57/57b done and integrated (report 57, S05 HTML section `#qb-trio`). Optional next steps:
   - genome-wide (not window) yak on a resized app;
   - trio-based switch rates per gene interval;
   - use the 179 trios as a truth set to tune the confidence score.
6. Marc: commit S05 + the 10-05 scrub; push/merge the branch (agent commits are classifier-blocked).
