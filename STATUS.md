# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-11 — Experiments A, B, C done; D pipeline built, ready to launch

**Experiment A (AoU-native distribution stats) complete, including the ancestry extension.** Full-cohort (n=535,658) missingness/resolution/homozygosity computed. Two real findings: AoU-native is capped at 3-field resolution (never 4-field, answers the standing "is 4-field needed" question for this arm specifically), and DPA1's homozygosity varies ~2.5x by ancestry (afr 27% → eur 68%), tracking the expected Out-of-Africa diversity gradient — looks like real population structure, not a calling artifact, and reframes the earlier "confident SpecImmune-DPA1 divergence" finding. DQA1 (the standing AoU-outlier locus) shows no such gradient — consistent with its issue being method-specific, not population-genetic. Full detail: EXPERIMENTS.md 2026-07-10 entries.

**Experiment B (SpecHLA short-read padding sweep) complete.** Real fragment-size floor measured (median 419bp, p99 649bp) — nowhere near the LR sweep's 34kb floor, confirming they needed separate measurement. **pad2000–pad10000 is a strong new SpecHLA-specific default candidate** (n=1): ~4.3-4.6min vs ~18.5min+ at pad500000, zero call degradation down to 325bp. A genuine "should break" sanity check (truncating into the gene body, not just the margin) confirmed the pipeline *can* break under real starvation — found a real wrong-allele call and a total-failure point — validating rather than undermining the padding result. **Found and fixed a real bug in our own QC**: the mate-dropout metric used `samtools flagstat` (insensitive to the manipulation being tested); switched to the actual singleton-discard rate. This is a SpecHLA-specific recommendation — do not conflate with SpecImmune's pad100k.

**Experiment C (SpecImmune gene-panel restriction) investigated, tested, and rejected.** Read the actual source (local clone, not the VM) and found the restriction mechanism (`--HLA_fa`/`--HLA_exon_fa`) genuinely works — every typing stage discovers genes dynamically from the `--db` folder, confirmed by inspecting `individual_ref/`. Found and fixed two real upstream bugs (an `UnboundLocalError` that had made the feature untestable, and a gzip read-miscounting bug) — packaged as `patch_specimmune_for_gene_restriction.py`, kept in place on the VM as harmless independent fixes. **But the actual controlled comparison (both aligners, person 2522883) showed restriction is worse, not better**: 18-23% slower in both bwa and minimap2, plus reproducible miscalls at Gene A and DRB1 — both aligners converge on the *same* wrong alleles at A, and the same allele-family shift at DRB1. Likely cause: removing pseudogene/paralog decoy references (A and DRB1 have the closest relatives among the 8 classical genes) redirects ambiguous reads onto the wrong allele rather than eliminating wasted work. **Decision: not adopting gene-panel restriction.** Full writeup + presentation visual: EXPERIMENTS.md 2026-07-10 entries. VM cleaned up (`db_classical8`, `experiment_c_restricted_db*` removed); local Mac environment (Homebrew tools, pip packages, cloned SpecImmune repo) also cleaned up.

## Pick up here

**Experiment D is BUILT and ready to launch (2026-07-11).** Three new scripts, designed as a single unattended, resumable overnight pipeline (survives the VM's ~1h idle sleep — quirk #14 — by being idempotent; a morning re-run of the same command finishes whatever a sleep interrupted, skipping completed people via `expd.done` markers + per-gene completeness checks, never exit-code-trusting):

1. `build_experiment_d_cohort.py` — derives the cohort straight from the mounted v9 manifests (CRAM + LR + ancestry TSV), applies the `/revio/` filter (quirk #13), picks N per ancestry group deterministically. No notebook/BigQuery step. Chosen scope: **~8-10 × all 6 ancestry groups (~56 people)**, run across multiple nights via resumability. Refuses to overwrite an existing cohort.tsv without `--force` (protects an in-flight run).
2. `run_experiment_d.sh cohort.tsv` — the orchestrator. Per person: remount-if-needed → SR slice (pad10000, Exp B band) → LR slice (pad100k) → SpecHLA → SpecImmune (bwa, `--align_method_1 bwa`, gene-panel restriction NOT used per C) → 3-way compare. No `set -e`; per-step resume; MAX_ATTEMPTS=4 backstop; optional `ntfy.sh` phone pings (set `NTFY_TOPIC`).
3. `analyze_experiment_d.py cohort.tsv` — aggregates the accumulated per-person `comparison_log.csv` rows (run_label=`experiment_d`), joins ancestry, computes per-locus × per-ancestry concordance. Purpose-built to answer the two carried-over questions: DQA1 AoU-lone-outlier rate by ancestry, and DPA1 SpecImmune-divergence rate (confidence-weighted). Writes markdown + 2 PNG heatmaps under `~/pipeline_outputs/experiment_d/analysis/` — aggregate-only, but keep on the VM (egress caveat) and paste tables back.

**Literal launch sequence (from `~/repos/pilot-validation`, after `git pull`):**
```bash
cd ~/repos/pilot-validation && pixi shell -e spechla          # samtools env (all tool calls are also env-pinned internally)
python3 build_experiment_d_cohort.py --per-group 10           # review the printed per-ancestry counts + cohort.tsv
nohup bash run_experiment_d.sh ~/pipeline_outputs/experiment_d/cohort.tsv > ~/pipeline_outputs/experiment_d/nohup.out 2>&1 &
disown                                                          # survives SSH/tab close; re-run the same line each morning until done
# when enough people are done:
python3 analyze_experiment_d.py ~/pipeline_outputs/experiment_d/cohort.tsv
```
Per-person budget ~25 min (SpecHLA ~4.5 + SpecImmune bwa pad100k ~21 + slicing); ~56 people ≈ 23-33h wall, so multi-night. Each person doubles as the "≥1 more person" confirmation B's and C's n=1 padding optimizations still need — watch the concordance for any config-induced drift. Sanity-check person 1 (~25 min after launch, via the ntfy ping or `tail -f nohup.out`) before trusting the batch; per-person isolation means one bad person can't sink the run.

The supervisor conversation (DECISIONS → "AoU-native trust") remains open and worth having whenever convenient — still not a blocking prerequisite.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **Any automated SpecImmune call must use `pixi run -e specimmune`** (quirk #15) and verify its output file exists — never trust exit code 0 alone for this tool. Same discipline now applies to SpecHLA runs too (quirk #18) and to any QC metric — verify it actually responds to the thing being tested before trusting it (the flagstat-vs-singleton-rate lesson).
- **Check `reference/` before any external documentation research** (quirk #16) — the AoU org PDF and prior research notes are already saved there.
- **pad100k (SpecImmune-LR) and pad2000-10000 (SpecHLA-SR) are two separate, tool-specific recommendations** — don't conflate them when setting up Experiment D's pipeline.
- **Gene-panel restriction is a closed question, not a queued optimization** — don't re-attempt it for Experiment D without a specific reason to revisit (e.g., the "smarter panel that keeps A/DRB1's paralogs" idea floated in EXPERIMENTS.md, which would need its own validation, not a default assumption it'll work).
