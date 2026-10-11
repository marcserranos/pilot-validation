# Runbook — what to run on the VM, in order

> **Role:** the literal command sequence, what each step costs, and what to paste back.
> **Read:** when you're at the Workbench terminal and want to make progress.

Everything here runs on the **cheap resized VM**, not the 96-core production machine. No step
re-runs Immuannot. Total compute is minutes-to-an-hour, not hours-to-days.

## Step 0 — sanity + pull

```bash
ls ~/pipeline_outputs/immuannot_calls.tsv ~/pipeline_outputs/immuannot_cohort_full.tsv
```
Both must exist and be non-empty — this confirms you're on the persistent disk and not a fresh one.

```bash
cd ~/repos/pilot-validation && git pull
```

## Step 1 — recon (RUN THIS FIRST, ~2 min)

**This is the gate. Nothing downstream is trustworthy until it passes.**

```bash
pixi run -e spechla -- python3 scripts/hla_popgen/00_recon_vm.py --limit 50
```

It answers, empirically, the questions the source-code spec could not:

- **Do `hap{1,2}/cds.fa.gz` still exist on disk?** The entire novel-allele workstream depends on
  this. If it reports ~0%, stop and tell me — the workaround is a cheap partial re-run from the
  kept `hap{1,2}.trimmed.fa`, not a full re-run, but the plan changes.
- Do `mm2.ipd.cds.paf.gz` / `mm2.ipd.gen.paf.gz` / `gene.filtered.paf` survive? These hold the full
  ranked candidate-allele list — a bonus capability if present.
- Every GTF attribute key actually present, and any key the spec doesn't document.
- Whether any KIR gene appears (it must not — if it does, the trim window is wrong).
- What fraction of haplotypes span more than one contig (i.e. how much cis-pairing we lose).
- Whether `cds_distance` is quoted or unquoted — an unresolved contradiction in the spec.

**Paste back the JSON block it prints.** That single block settles the open questions.

## Step 1b — warning-token census (~1 min)

Recon measured `template_warning` on **95.3%** of transcript rows. That number makes any blanket
"no warning" filter useless (it would reject ~95% of calls), so the policy has to be token-aware.
This tells you which tokens to disqualify:

```bash
pixi run -e spechla -- python3 scripts/hla_popgen/00b_warning_census.py --limit 200
```

Paste back its two tables. Feed the answer to `03_novel_alleles.py --disqualifying-warnings`
(default: `inframe_stop` only).

## Step 1c — fix the Jupyter UI lag (one-time)

**Real incident, 2026-09-04:** `~/pipeline_outputs` has ~12,000 top-level `person_id`-named
directories. Opening it in the Workbench Jupyter file browser lags/crashes, because the UI tries to
list+stat every entry. The fix: move every `person_id` directory one level deeper, into
`~/pipeline_outputs/people/`, so the top level only has the aggregate `.tsv` files (`hla_calls_rich.tsv`,
`cohort_membership.tsv`, `immuannot_cohort_full.tsv`, `immuannot_calls.tsv`, `ancestry_preds.tsv`,
`hla_genotypes.tsv`, `novel_alleles_seqs.fa`) plus the one `people/` folder — a handful of entries the
UI can render instantly.

**Run this once**, before any script below — `00_recon_vm.py`/`00b_warning_census.py` above already
default to `--outroot ~/pipeline_outputs/people`, so if you're on a fresh VM that hasn't had this move
done yet, run it now, before re-running Step 1/1b too:

```bash
mkdir -p ~/pipeline_outputs/people
cd ~/pipeline_outputs
for d in */; do
  d="${d%/}"
  [[ "$d" =~ ^[0-9]+$ ]] && mv "$d" people/
done
cd -
```

Why this is safe at ~12,000-entry scale:
- No shell glob of all 12,000 names in one command (no `mv */ people/`, no brace expansion) — the
  `for d in */` loop expands the top-level directory listing only, one name at a time, so there's no
  `ARG_MAX`/command-line-length risk the way a single `mv <12000 names> people/` invocation would have.
- The `^[0-9]+$` regex test only ever matches pure-digit `person_id` basenames — it will never touch
  the aggregate `.tsv`/`.fa` files (non-directory, and don't match `*/`  anyway) or any other named
  directory (`people/` itself, or anything else that isn't a bare integer).
- **Idempotent/resumable**: if this is interrupted (dropped session, VM restart) and re-run, every
  already-moved `person_id` directory is simply no longer present at the top level, so `for d in */`
  won't see it again — the loop just picks up wherever it left off. Safe to re-run any number of times.

Verify afterward:
```bash
ls ~/pipeline_outputs | wc -l        # should now be small (a handful of .tsv/.fa files + people/)
ls ~/pipeline_outputs/people | wc -l # should be ~12,000-ish
```

## Step 2 — rich extraction (Tables 1 & 2)

Sanity-check on 200 people first — never launch the full pass blind, and **measure the rate before
committing to the full run**. Recon took ~6.5 s/person single-threaded; `01` is threaded, but
confirm the real rate rather than assuming it:

```bash
{ time pixi run -e spechla -- python3 scripts/hla_popgen/01_extract_rich.py --sample --limit 200 ; }
```

Then the full cohort. It is checkpointed and resumable, so a re-run continues rather than restarts.
**Launch it detached from the terminal and log to a file** — standing practice for any unattended
run past a few minutes (quirk #14/#22: a dropped browser session has silently killed a multi-hour
job here before, and printed-only output was lost when the VM idled out). `tmux` is not installed
on this VM and there is no sudo to add it (quirk #6) — use `nohup ... & disown` instead, exactly
the fallback quirk #14 itself names:

```bash
nohup pixi run -e spechla -- python3 scripts/hla_popgen/01_extract_rich.py \
    > ~/pipeline_outputs/01_extract.log 2>&1 &
disown
```
```bash
tail -f ~/pipeline_outputs/01_extract.log
```

`Ctrl-C` on the `tail` only stops watching — `disown` detaches the job from this shell, so it
survives a dropped session or a closed tab. Check it's still alive with
`ps -ef | grep 01_extract_rich`. Once the log shows it finished:

```bash
pixi run -e spechla -- python3 scripts/hla_popgen/01_extract_rich.py --validate
```

`--validate` asserts row-grain uniqueness on `(person_id, hap, contig, gene, copy_index)` — the
exact bug class that silently destroyed the last run's classical-gene calls (quirk #29). Do not
skip it.

Produces `~/pipeline_outputs/hla_calls_rich.tsv` (all 65 genes, contig-keyed, with distances,
novelty depth and codon-level diffs) and `hla_cis_pairs.tsv` (physically-phased heterodimers).

## Step 3 — cohorts (Table 4)

Needs the gcsfuse mount. Paste these as **two separate commands** (quirk #2 — never chain the mount
with a consumer):

```bash
mkdir -p ~/mnt/aou-controlled
```
```bash
gcsfuse --billing-project wb-cordial-leechee-9743 --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled
```
Verify it actually resolved before trusting it:
```bash
ls ~/mnt/aou-controlled/v9/wgs
```
```bash
pixi run -e spechla -- python3 scripts/hla_popgen/02_build_cohorts.py
```

Joins long-read calls, platform/trim tier, AoU genetic ancestry (**including the continuous 6-way
admixture proportions**, which the previous analysis round never used) and short-read membership.

## Step 4 — novel alleles

```bash
pixi run -e spechla -- python3 scripts/hla_popgen/03_novel_alleles.py
pixi run -e spechla -- python3 scripts/hla_popgen/04_allele_saturation.py
```

`03` clusters novel alleles by observed CDS identity, requires ≥2 unrelated carriers to pass QC,
and flags homopolymer-indel artifacts. `04` produces the discovery curves and Chao2/ACE richness
estimates — the "% of HLA allele space discovered, per gene per ancestry" headline.

## Step 5 — figures

```bash
pixi run -e spechla -- python3 scripts/hla_popgen/05_figures_frequency.py --cohort lr
pixi run -e spechla -- python3 scripts/hla_popgen/06_figures_structure.py --cohort lr
pixi run -e spechla -- python3 scripts/hla_popgen/07_figures_crosscohort.py --cohort lr
```

`07` always compares short-read against long-read, so `--cohort` selects only the **long-read side**
of the comparison — it accepts `lr` or `lr_td`, and deliberately refuses `sr` (comparing short-read
against itself is meaningless).

For the template-distance sweep, re-run 05/06 at each cut and compare:

```bash
for TD in 0 1 5; do
  pixi run -e spechla -- python3 scripts/hla_popgen/05_figures_frequency.py --cohort lr_td --td-max $TD
done
```

And the short-read cohort for the maximum-N frequency baseline:

```bash
pixi run -e spechla -- python3 scripts/hla_popgen/05_figures_frequency.py --cohort sr
```

## Where results land

Everything under `~/pipeline_outputs/` (per-person tables) and `reports/hla_popgen/` (aggregate
figures + markdown reports with the numeric tables behind each figure). Nothing here downloads or
exports anything.

**Novel-allele nucleotide sequences (`novel_alleles_seqs.fa`) stay on the VM** — they are
per-participant genomic sequence and never belong in a report or a repo. Getting the PNGs out for a
supervisor deck is a real file-egress action: use AoU's official reviewed download workflow, don't
improvise a path for it.

## Publishing results to the share bucket (do this at the END)

A share bucket already exists and collaborators (Cole) read from it. Full metadata:
`context/ENVIRONMENT.md` → "Our own buckets". Short version:

- Bucket: `gs://hla-calls-share-wb-cordial-leechee-9743` — workspace-owned, **no `--billing-project`
  flag** (unlike the AoU source bucket).
- `aggregate/` already holds the production run's `*.tsv`. **Wave 1 is done.** Append there:

```bash
gcloud storage cp ~/pipeline_outputs/*.tsv gs://hla-calls-share-wb-cordial-leechee-9743/aggregate/
```

- The per-person raw tree was deliberately **not** copied. ~81 MB/person × ~12k ≈ 970 GB, of which
  ~95% is alignment scratch (`mm2.*.paf.gz`, uncompressed `gene.filtered.paf`, logs). The valuable
  part is only `hap{1,2}.gtf.gz` + `hap{1,2}/cds.fa.gz`, ~3 MB/person. **`00_recon_vm.py` reports
  exact per-file-type sizes** — use its numbers, not an estimate, before copying anything bulky.
- Everything stays inside the VPC-SC perimeter, so this is **not egress**. It does not touch the
  open download/publication question in `DECISIONS.md`. Anyone with access must be individually
  Controlled-Tier approved — these are participant-level genotypes.
- Never publish `novel_alleles_seqs.fa` or any bare person_id to a report; per-person data stays in
  `~/pipeline_outputs/` and the bucket, not in the repo.

## If something fails

Every script fails loud with a remediation hint rather than continuing silently. The two most
common causes, both documented in `context/ENVIRONMENT.md`:

- **`FileNotFoundError` on `~/mnt/aou-controlled/...`** → the VM restarted and dropped the mount
  (quirk #11/#14). Remount, verify with `ls`, retry.
- **Anything ancestry-related returning zero matches** → `ancestry_pred` is lowercase in the real
  files (quirk #30). The scripts normalize on load; if you've written new code, do the same.
