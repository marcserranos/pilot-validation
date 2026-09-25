# VM handoff — for Marc, run in a JupyterLab terminal

> **Status 2026-09-25:** (a) done (0ec28ff), (b) done earlier (59fef06), (c) done (170-person summary in `41_kir_scoping/`), (d) done — full-cohort KIR run, see `KIR_FULL_RUN.md` and `reports/hla_popgen/43_kir_full_cohort/`. (e) partially: 1,500-person run pulled (b9ff562); the full 7,640-person run remains.

Auto-mode VM automation was refused by the permission classifier repeatedly this session
("Auto-Mode Bypass" at the websocket-terminal step). Everything computable locally is done;
these five steps need a human hand on the VM. Start each terminal session with:

```bash
exec > >(tee -a ~/.claude_session.log) 2>&1
```

General notes: VM repo `~/repos/pilot-validation` is on a DIFFERENT branch
(`needle-view-cds-diversity-density`) — do not `git checkout`/commit there for S03 work.
S03 scripts live in `~/s03/` (deploy any missing/updated local file there first via Jupyter's
file browser or `scp`/upload). Run with
`PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen`. Use `setsid nohup ... & disown`
for anything >40s. Assume the VM restarted since the last session — `ls ~/pipeline_outputs/`
to confirm the bucket mount is live before trusting any read.

---

## (a) Figure 1 v5 panels a/c, then local compose

Two VM export scripts, then one local compose. Confirm `~/s03/40a_admixture_bins.py`,
`~/s03/40b_novelty_rate_export.py`, `~/s03/40_figure1_v5.py`, `~/s03/24_novelty_by_field.py`
are present (upload from local `scripts/hla_popgen/` if not — 40a/40b both `_load_module`-style
reuse `24_novelty_by_field.build_people`, so 24 must sit in the same directory as 40a/40b).

```bash
cd ~/s03
mkdir -p results/40
PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen \
  python3 40a_admixture_bins.py --out results/40/panel_a_admixture_bins.tsv
PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen \
  python3 40b_novelty_rate_export.py --out results/40/panel_c_novelty_totals.tsv
```

Both use standard default paths (`~/pipeline_outputs/hla_calls_rich.tsv`,
`~/pipeline_outputs/cohort_membership.tsv`, default relatedness table) — no flags strictly
required beyond `--out`, but pass `--out` explicitly as above so the filenames match what
`compose` expects below. Each takes well under a minute (aggregation only, no long-read work).

**Copy back** `results/40/panel_a_admixture_bins.tsv` and `results/40/panel_c_novelty_totals.tsv`
to local `reports/hla_popgen/40_figure1_v5/` (these are bin/gene-level aggregates only, already
≥20-person gated — disclosure-safe to leave the VM per the scripts' own docstrings).

Then **locally** (not on the VM):

```bash
cd /Users/marcserrano/WORK/STANFORD/pilot-validation
python3 scripts/hla_popgen/40_figure1_v5.py compose \
  --panel-a-bins reports/hla_popgen/40_figure1_v5/panel_a_admixture_bins.tsv \
  --panel-c-totals reports/hla_popgen/40_figure1_v5/panel_c_novelty_totals.tsv \
  --out-dir reports/hla_popgen/40_figure1_v5
```

Check the regenerated `figure1_v5.png` (Read the PNG) for panel a/c now populated and no
overlap/hatching regressions before calling this done.

---

## (b) 37d/37e — held-out-fold EM rerun (critic #1 circularity fix)

`CRITIC_1.md` item 3 flagged the original 37d mask-and-rephase EM as circular (fit on a pool
that includes the truth set it's tested against). **Not yet run** — `scripts/hla_popgen/
37e_unrelated_fix_kfold_em.py` (untracked locally) is the fix (k-fold held-out EM + a naive
linkage-equilibrium baseline) but has no VM output yet; 37's README still only describes the
circular 37d result.

Deploy `37e_unrelated_fix_kfold_em.py` to `~/s03/` (needs `37_vm_run.py` and
`24_novelty_by_field.py` also present in `~/s03/` — same `_load_module`-from-own-directory
pattern as 40a/40b).

```bash
cd ~/s03
mkdir -p results/37
setsid nohup python3 -u 37e_unrelated_fix_kfold_em.py \
  --out-dir results/37 --resolution 2 --k-folds 5 \
  < /dev/null > results/37/37e_run.log 2>&1 &
disown
tail -f results/37/37e_run.log   # ctrl-C to stop watching; job keeps running
```

Also fixes a second critic finding (blocker #2): reruns on the CORRECT unrelated set (~11,856,
matching 29/38/39) instead of 37_vm_run.py's own buggy 13,252. Should finish in minutes (no
long-read work, just EM on genotype tables already in `~/pipeline_outputs/`).

**Copy back** whatever `results/37/*.tsv`/`.png`/`.pdf`/`.json` the script writes (aggregate-only
by construction — see the script's own `_sup()` disclosure gating) to
`reports/hla_popgen/37_dq_g1g2_signed_ld/`, and note in that README's Deliverable 3 section that
the held-out-fold result supersedes the original circular 37d numbers.

---

## (c) KIR extended pilot (170 people) — aggregation only

Already running unattended since ~2026-09-23 (`reports/hla_popgen/41_kir_scoping/README.md`
§4a). Check status, and once done, re-run the same command — it's idempotent/resumable and
just recomputes the aggregate summary from whatever's finished:

```bash
tail -20 ~/s03/results/41/run_extended.log
```

If it shows all 170 done (or you're ready to aggregate whatever's finished):

```bash
cd ~/repos/pilot-validation && pixi run -e specimmune -- python3 ~/s03/41_kir_pilot.py \
    --n-people 170 --seed 41 --threads 2 --outroot ~/s03/results/41/pipeline_outputs \
    --out-suffix .kir41 > ~/s03/results/41/run_extended_aggregate.log 2>&1
tail -5 ~/s03/results/41/run_extended_aggregate.log   # final aggregate-only Quality dict
```

**Copy back** `~/s03/results/41/pipeline_outputs/41_kir_pilot_summary.kir41.tsv` (aggregate-only,
per the README's disclosure note) to `reports/hla_popgen/41_kir_scoping/`, and fold the
170-person numbers into the README's §3/Distilled section (larger ancestry-stratified panel than
the original 20-person pilot).

---

## (d) KIR full-cohort run — needs your resize decision

**Not started — your call, not automated.** `41_kir_scoping/README.md` §5 (Go/no-go) recommends
GO, conditional on: (1) `--region`/`--pad`/`--out-suffix` flags added to
`run_production_orchestrator.py` (currently hardcoded to the HLA window — mechanical change,
mirrors `run_immuannot_person.py`'s existing `--region` flag) and (2) a spot-check of the
framework-gene misses (§3, 85–97.5% present across pilot haplotypes).

Cost/time estimate (§4, extrapolated from the 2-core pilot, flagged as a bigger extrapolation
leap than HLA's own estimate ever validated): **~14.5–16.2 h, ~$45** on a 96-core machine — vs
HLA's own 52–58h/~$160–200 on the same class of machine. Requires resizing the Workbench app
to 96 cores, which restarts the VM and kills any concurrently running jobs (i.e., don't do this
while (a)/(b)/(c)/(e) above are still running).

Once you decide to proceed: implement the orchestrator flags (small, mechanical — see
`run_immuannot_person.py`'s own `--region` for the pattern to mirror), do the framework-gene
spot-check, resize the app, then launch the full run with `--region chr19:54,600,000-54,920,000
--out-suffix .kir_full` (or similar) against a **separate `--outroot`**
(e.g. `~/pipeline_outputs_kir`, per README §4's explicit warning — reusing the HLA outroot risks
the resumability check silently skipping/mixing runs).

---

## (e) WS6 (42_repertoire_baseline) — pull results if the VM run finished

`sprints/S03_call8_figures_kir_prediction/LOG.md` (WS6 entries, ~[t+75m] through end) — the real
run hit and fixed 3 bugs (local gcsfuse path handling, BigQuery int64 person_id vs str, real
TRUST4 column names) and was resumed with this exact command (fresh log to avoid clobbering
`run5.log`):

```bash
cd ~/s03 && PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen setsid nohup python3 -u 42_repertoire_baseline.py \
  --results-bucket ~/workspace/aleix_rnaseq/repertoire_results \
  --manifest ~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq/manifest.tsv \
  --cohort-membership ~/pipeline_outputs/cohort_membership.tsv \
  --relatedness ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv \
  --hla-calls ~/pipeline_outputs/hla_calls_rich.tsv \
  --scratch ~/pipeline_outputs/rnaseq/ws6_scratch \
  --out-dir ~/s03/results/42 --n-jobs 2 --n-repeats 3 --max-kmers 1500 --with-xgboost \
  > ~/s03/results/42/run6.log 2>&1 & disown
```

Check whether it's already finished:

```bash
tail -30 ~/s03/results/42/run6.log
ls ~/s03/results/42/
```

If `42_repertoire_baseline_{metrics,top_features,qc}.tsv` and a figure exist, **copy them back**
to `reports/hla_popgen/42_repertoire_baseline/` (this report directory doesn't exist locally
yet — create it) and write a README following the same convention as 37–41 (question, method,
results, distilled section, caveats) using LOG.md's WS6 entries for the method narrative
(BenchRep-T protocol: VJ usage + gapped 4-mers, L1-LR + XGBoost, 3-fold stratified CV, pooled
OOF AUROC/AUPRC; EHR-depth covariates and positive controls per LOG.md `[t+75m] DECISION`). If
`run6.log` shows it's still running or died, `tail` the log for the error and re-launch with a
fresh log path (increment to `run7.log`) rather than reusing `run6.log`.
