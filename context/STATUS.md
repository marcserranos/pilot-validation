# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-11 — All of Us access restored; repo harmonized; planning before Sprint S05

**Access.** Marc lost Workbench access ~2026-09-28 (workspace showed Standard Tier only); UCSC
reinstated it, confirmed restored 2026-10-11. Assume the VM restarted and its disk state is
whatever it was on 2026-09-28 (quirk #14): remount, `ls`, verify before trusting anything.

**Branches.** PR #1 merged 2026-10-11: `main` now holds S03 + S04 and the harmonization pass and
is the canonical branch. `needle-view-cds-diversity-density`
are fully contained in `main` and can be deleted; keep `s04-kir-recurrence-style-share` (commit 5ebfcd4, gitignored walkthrough) until Marc decides.
Aleix's `aleix/hla-resolve-phase1` shares no history with `main` (separate root); leave it alone.

**Where things are.**
- Index of every result folder, current vs superseded: [`reports/hla_popgen/README.md`](../reports/hla_popgen/README.md).
- Last sprint (S04) plain-language findings: [`sprints/S04_kir_recurrence_style_share/FINDINGS_FOR_MARC.md`](../sprints/S04_kir_recurrence_style_share/FINDINGS_FOR_MARC.md); figure catalogue: `FIGURES_INDEX.md` in the same folder.
- Last supervisor meeting: 2026-10-06. Agreed focus: QC and phasing (incl. regional MHC/KIR coverage filter) and KIR; later, disease curation.
- Aleix's live state: `aleix/RNA-seq/context/STATUS.md` on his branch (W0 analysis table and W1 V/J usage done 10-08; TRBV4-3 germline-absence lead).

**Done since 2026-09-28 (S04 close-out, before access was lost):** 44 v4b full run (KIR/HLA
recurrence + catalogue completeness, valid numbers in 44/46 READMEs); 4 of 9 PUBLIC tables of 44 pulled.

**Open, needs Marc:**
1. Run `sprints/S04_kir_recurrence_style_share/VM_CLEANUP_PLAN_v2.md` (agent deletes are classifier-blocked).

**Literal next commands (VM, once Marc steers):**
1. Reattach the big_run app, remount gcsfuse, `ls` the bucket and `~/s04/results/44_public/`.
2. Pull the remaining 5 of 9 `44_public` tables, commit them under `reports/hla_popgen/44_kir_recurrence_saturation/`.
3. Then start the S05 board from the plan Marc approves.

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
