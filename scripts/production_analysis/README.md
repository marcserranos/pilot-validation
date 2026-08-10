# Production-run post-processing & supervisor report

Five analysis scripts, built 2026-08-10 after the full-cohort Immuannot production run finished,
covering everything Marc asked for: completeness/coverage/demographics, a confidence-threshold
discrepancy sweep vs AoU-native, a confidence distribution plot, HLA-vs-ancestry clustering, an
allele-frequency-by-ancestry spectrum, and a capstone figure synthesizing this project's
"DRB1 is the hardest locus" finding. All are cheap I/O over small TSVs + per-person `.gtf.gz`
files — **run these on the resized, cheap Workbench environment**, not the 96-core production VM.

## ⚠ Known bug, already fixed — read before running anything

**`merge_fragments()` in `run_production_orchestrator.py` deduplicated `immuannot_calls.tsv` on
`person_id` alone, silently collapsing every person's 8+ gene rows down to just one (whichever
sorted alphabetically last — always one of `MICA`/`MICB`/`TAP1`/`TAP2`, never a classical gene).
Found 2026-08-10 against the real production output. Fixed in `merge_fragments()` (now keys on
`[person_id, gene]`), and a one-time repair script,
`scripts/production_orchestrator/rebuild_immuannot_calls.py`, rebuilds the canonical file from the
raw per-person `hap{1,2}.gtf.gz` files (deliberately kept, not pruned) — the per-worker fragments
this could otherwise re-merge from were already deleted by `merge_fragments()` itself.
**You don't need to run the repair manually — `run_all.sh` below detects this exact symptom and
repairs it automatically before running anything else, and every individual script also refuses to
run (loud `FATAL`, not a silent empty result) if it ever sees this again.** The repair reads local
disk files (not the gcsfuse bucket), so it's fast — but `run_all.sh` doesn't just assume that: it
first times a 200-person sample, extrapolates the real full-cohort time from a measured rate, and
only proceeds automatically if that estimate is under 15 minutes (`MAX_AUTO_REPAIR_MIN` env var to
change the threshold). If the real cohort turns out slower than expected, it stops and tells you
instead of silently running long. To check the timing yourself before running anything else:
```bash
pixi run -e spechla -- python3 scripts/production_orchestrator/rebuild_immuannot_calls.py --limit 200
```

## Run everything in one command

```bash
cd ~/repos/pilot-validation && git pull
pixi run -e spechla -- bash scripts/production_analysis/run_all.sh
```

That's it — `pixi run` activates the environment for just this command (no separate `pixi shell` +
`pixi install` steps needed). `run_all.sh`:
1. Checks `immuannot_calls.tsv` for the bug above; auto-repairs if needed (safe — the old file is
   renamed aside with a timestamp, never deleted).
2. Mounts gcsfuse if it isn't already (billing project `wb-cordial-leechee-9743`, confirmed
   2026-08-10 — override with `AOU_BILLING_PROJECT=...` if this ever runs in a different
   workspace). Non-fatal if it fails — just skips the one script that needs it.
3. Runs all 5 scripts, continuing past any individual failure rather than stopping the batch.
4. Prints one PASS/FAIL/SKIPPED summary table and exactly where every output landed.

Paste the terminal output back, or open the PNGs directly in the Jupyter file browser — no need to
run anything script-by-script or babysit each step.

## Running scripts individually (advanced / debugging)

```bash
pixi shell -e spechla   # wait for the (omni-hla-pilot:spechla) prompt
pixi install -e spechla   # first time only, to pick up umap-learn (added 2026-08-10)

python3 scripts/production_analysis/analyze_completeness_and_demographics.py
python3 scripts/production_analysis/analyze_confidence_vs_aou_native.py       # needs the gcsfuse mount -- see below
python3 scripts/production_analysis/cluster_hla_by_ancestry.py
python3 scripts/production_analysis/analyze_allele_frequency_by_ancestry.py
python3 scripts/production_analysis/summarize_drb1_evidence_capstone.py       # no pipeline data read at all
```

Each is independent (different inputs, different output dir) — run them in any order, or just the
ones you need. Defaults assume the standard `~/pipeline_outputs/` layout from
`scripts/production_orchestrator/RESULTS_LOCATION.md`; override with `--calls`/`--cohort`/
`--timing`/`--aou-tsv`/`--outroot` if your paths differ.

**Manual gcsfuse mount** (only needed for `analyze_confidence_vs_aou_native.py`; `run_all.sh`
does this automatically). Billing project is `wb-cordial-leechee-9743` (confirmed 2026-08-10 via
`gcloud config get-value project`, for the Stanford-pod workspace this production run actually
used — don't reuse the earlier pilot workspace's `wb-glacial-potato-8710`, they're different GCP
projects). Paste as two separate commands (quirk #2 — never chain the mount and a consumer command
together):
```bash
mkdir -p ~/mnt/aou-controlled
gcsfuse --billing-project wb-cordial-leechee-9743 --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled
```
Then verify it actually resolved before trusting it:
```bash
ls ~/mnt/aou-controlled/v9/wgs
```

## What each produces

| Script | Output dir | Figures | Answers | Needs the gcsfuse mount? |
|---|---|---|---|---|
| `analyze_completeness_and_demographics.py` | `production_analysis/completeness/` | `completeness_overview.png` (6-panel), `timing_stats.png` | Did the pipeline capture data correctly? Who got covered, by ancestry and platform? Basic operational stats. | No |
| `analyze_confidence_vs_aou_native.py` | `production_analysis/confidence/` | `baseline_concordance_by_gene.png`, `confidence_threshold_sweep.png`, `confidence_distribution.png` | How does AoU-native compare to Immuannot at full scale? How does that discrepancy change as confidence tightens, per gene? What does the confidence distribution look like per gene (mean/median/clustering)? | **Yes** — reads AoU-native's `hla_genotypes.tsv` from the mount |
| `cluster_hla_by_ancestry.py` | `production_analysis/clustering/` | `pca_pc1_pc2_by_ancestry.png`, `umap_by_ancestry.png` | Does HLA genotype structure correlate with ancestry under PCA/UMAP? | No |
| `analyze_allele_frequency_by_ancestry.py` | `production_analysis/allele_frequency/` | `allele_frequency_spectrum.png` (8-panel) | What does the actual allele-frequency spectrum look like per gene, per ancestry? (Aim 1) | No |
| `summarize_drb1_evidence_capstone.py` | `production_analysis/drb1_capstone/` | `drb1_evidence_capstone.png` (6-panel) | Why does this project keep concluding DRB1 is unreliable? (synthesizes 6 prior findings — see caveat below) | No (reads no pipeline data at all) |

All write a markdown report alongside the PNGs with the numeric tables behind each figure — paste
these into the supervisor report directly, they're already aggregate-only prose+tables, same
discipline as every other analysis script in this repo.

## Design notes worth knowing before reading the figures

- **The confidence-threshold sweep chart's x-axis is reversed on purpose** — 0/strict sits on the
  right, so "raising the certainty threshold" reads left-to-right as discrepancy falls, matching
  how Marc described wanting to read it. Read the axis label, not just the numbers.
- **The sweep's discrepancy floor (right edge) is an estimate of AoU-native's real error rate
  against high-confidence ground truth** — not the discrepancy at threshold=0 alone; look at where
  the curve visibly levels off, same interpretation as the existing n=60 confidence-matched-truth
  work (`scripts/analyze_confidence_matched_truth.py`), now at full production scale. The right
  edge often gets noisy from shrinking N — the shaded area behind the lines shows this; read a
  thin-N tail skeptically.
- **The clustering script's PCA/UMAP scatter plots are one dot per real person** — no identifiers on
  them, but structurally person-level in a way most of this project's outputs deliberately aren't.
  This is standard AoU Workbench research-figure practice (AoU's own ancestry PCA plots are built
  the same way) — fine to build and view inside the Workbench. If either figure is ever shared
  outside the sponsor/supervisor review loop, treat it like any other research output under the
  still-open compliance question in `context/DECISIONS.md`.
- **The clustering script only uses complete-case people** (real calls at all 8 classical genes,
  both haplotypes) — no imputation. The report states exactly how many people that excludes; if the
  excluded fraction is large, that's itself a data point about the run's completeness (cross-check
  against `analyze_completeness_and_demographics.py`'s numbers), not a clustering-method problem.
- **`trim_tier == self_align_needed` (991 sequel2 people) will show ~0% completion everywhere** —
  Phase 2 was disabled for this production launch after Tier 3 failed testing
  (`context/EXPERIMENTS.md`, 2026-08-05). Expected, not a bug — the completeness script's markdown
  says this explicitly so it isn't mistaken for a problem when reviewing the report.
- **The DRB1 capstone script is different in kind from the other four.** It reads no live pipeline
  data at all — it transcribes six already-published findings from earlier, smaller pilots
  (Experiment D n=60, the AoU callset validation report, Experiment F n=60), each hardcoded with an
  exact citation to its source file. It will not automatically reflect anything new the other four
  scripts find in the full production cohort. If the full-scale results change the DRB1 story
  (better or worse), that's a genuinely new, real finding worth adding as a 7th line — not
  something this script would surface on its own.
- **The allele-frequency script's ancestry variable is AoU's genetic-ancestry prediction**
  (`ancestry_pred`, from `immuannot_cohort_full.tsv`), same convention as every other ancestry
  figure in this project (`context/DECISIONS.md` — preferred over self-reported race).

## Step by step: from a freshly-resized VM to results in hand

Assumes the Workbench Cloud Environment has already been resized down per the earlier session's
confirmed procedure (stop → three-dot menu → Edit → change machine type → Update) — this section
picks up from there.

1. **Start the resized environment** from the Workbench UI and open its Jupyter terminal tab. It
   should boot noticeably faster than the 96-core VM did.
2. **Sanity-check you're on the right disk before doing anything else:**
   ```bash
   ls ~/pipeline_outputs/immuannot_calls.tsv ~/pipeline_outputs/immuannot_cohort_full.tsv
   ```
   Both should exist and be non-empty (`ls -la` to check size) — this confirms the resize
   preserved the persistent disk and you're not looking at a fresh/empty one.
3. **Pull the latest code:**
   ```bash
   cd ~/repos/pilot-validation && git pull
   ```
4. **Run everything with the one command from "Run everything in one command," above:**
   ```bash
   pixi run -e spechla -- bash scripts/production_analysis/run_all.sh
   ```
   This alone handles environment activation, the `umap-learn` install, the merge-bug auto-repair,
   the gcsfuse mount, and all 5 scripts — nothing else to do manually. If you'd rather run scripts
   one at a time (e.g. to debug a single failure), see "Running scripts individually" above instead
   — that section also covers the manual mount command and the billing-project caveat in full.
5. **Where results land:** everything goes under `~/pipeline_outputs/production_analysis/`, one
   subdirectory per script (`completeness/`, `confidence/`, `clustering/`, `allele_frequency/`,
   `drb1_capstone/`), each holding its PNG figure(s) + a `*_report.md`. Nothing writes outside that
   tree, and nothing here downloads or exports anything — to actually view the PNGs and put
   together the supervisor deck, either:
   - Open them directly in the Workbench's Jupyter file browser (image preview works there), or
   - Paste the markdown tables into your report as-is (already aggregate-only prose+tables), or
   - If you want the images themselves outside the Workbench: that's a real file-egress action
     (unlike the aggregate-only telemetry this project already treats as a separate, lesser
     question) — go through AoU's official reviewed download workflow for these PNGs specifically,
     don't improvise a new path for it.

## Compliance reminder

Nothing in this directory downloads or exports raw per-person calls — every script here only ever
writes aggregate counts, rates, distributions, and (for the clustering script specifically)
unlabeled per-person points in a derived embedding space. Real allele calls and bare person_ids
never leave `~/pipeline_outputs/` or get written into any file under this directory. Same standing
rule as `scripts/production_orchestrator/RESULTS_LOCATION.md`.
