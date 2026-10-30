# Handout — run the novel-allele ladder on the VM (paste this whole file to the agent)

You are a Claude Code agent driving Marc's All of Us Workbench VM through his own logged-in Chrome tab (computer-use). Marc hands you the tab; you never see or need credentials. Goal: run scripts 118/119 and produce the **exact-count** report plus a **cloud-safe** data file. Read `CLAUDE.md`, `context/ENVIRONMENT.md` (quirks #14, #31–34) and `sprints/S08_case_browser/RUNBOOK_novel_ladder.md` first; this handout is the short version plus the rules that bite.

## What is being computed
All "minimally novel" HLA alleles: Immuannot-novel calls that are rebuildable from the public reference allele, have no indel, and are not an exact CDS of a named allele in the newest IPD-IMGT/HLA release (3.66.0, 2026-10-14). Per allele: carrier counts (all / unrelated), platforms, transmission, effect, template flags, and a short-read confidence (LOD from srWGS reads realigned to the person's own haplotype), then 4 confidence buckets and saturation curves. Script docstrings in `scripts/hla_popgen/118_novel_ladder.py` / `119_novel_ladder_report.py` define everything.

**Status: 118/119 were written and tested only on synthetic data (no VM in the cloud session).** Expect to fix small things on first contact with the real files (column names in `hla_calls_rich.tsv`, `case_map.private.tsv`, the lrWGS manifest, `sr.jsonl`). Fix them in a **throwaway clone**, keep a list of every change, and hand the diff to Marc; do not push.

## Rules of engagement (from CLAUDE.md, restated)
1. **Assume the VM restarted.** Remount the bucket, verify with `ls` before trusting any AoU path (quirk #14).
2. **Chrome terminal:** start every terminal with `exec > >(tee -a ~/.claude_session.log) 2>&1`; read output from a Notebook/Console cell via `get_page_text` (strip ANSI per quirk #31); never screenshot the terminal to read text. **`type` and `Return` are always two separate tool calls.** The Workbench embeds JupyterLab in a cross-origin iframe: navigate to the inner JupyterLab URL first.
3. **`git commit` / `git push` typed via computer-use is often blocked.** Try at most 2–3 times, then give Marc the exact command. Prefer not to commit from the VM at all.
4. **Ask Marc before** killing a running process, restarting, or writing to a path another run may write to (quirk #23). Use fresh out-dirs (`~/s10`, `~/s10_cb`); **never write into `~/case_browser`** (a fresh extract renumbers case ids and would orphan its `sr.jsonl`).
5. **No job longer than one day**; use `--max-minutes` and resume (`sr.jsonl` is checkpointed).
6. **Exploratory relaxation (Marc, 2026-10-30):** for this exploratory run the `<20` blur is relaxed. You may use and report **exact counts** (including 1–19) to Marc, and the exact-count selector report is the main deliverable. What stays strict: the repo is public, so **no exact or small-count result is ever committed or pushed** (no commit messages, no copied tables, no figures with exact small counts). Everything that holds exact counts is gitignored (`*.private.*`, `ladder_data.*.json`, `S08_*LADDER*.html`, `~/s10`, `~/s10_cb`): keep it on the VM, never `git add -f`, never copy it into a clone. Do not move exact files off the VM on your own; if Marc wants one elsewhere he will say so (All of Us's own small-cell dissemination rule still applies to anything that leaves the Workbench, and that is Marc's call, not yours).
7. Do not quote or cite anything from `sprints/_calls/`.

## Steps (all outputs under `~/s10` and `~/s10_cb`)
1. Setup: fresh clone of branch `main-si057v` into `~/cb_code` (`git clone` works; auto-clone does not, quirk #3). Verify `ls ~/pipeline_outputs/hla_calls_rich.tsv ~/pipeline_outputs/cohort_membership.tsv ~/tools/Immuannot_refdata ~/case_browser/case_map.private.tsv`.
2. IPD: `curl -O https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_nuc.fasta` (43 MB). If the Workbench cannot reach GitHub, ask Marc to upload it. Check the release: `curl -s https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt | head -3` should say 3.66.0 (or newer: then pass that number as `--ipd-release` and tell Marc).
3. Short-read evidence for every substitution-only novel call (cap raised from the earlier 60 per gene):
   - `cd ~/cb_code && python3 scripts/hla_popgen/90_case_browser.py extract --out-dir ~/s10_cb --exhaustive-snv 5000` (see `--help`; the 10-21 run used the same entry point). Check the log line `exhaustive_snv_total:*` against ~675 substitution-only novel calls.
   - `python3 scripts/hla_popgen/90_case_browser.py sr --out-dir ~/s10_cb --max-cases 5000 --max-minutes 240` (~1.5 s/case, resumable; rerun the same command until `sr.jsonl` has a line for each case).
4. Ladder: `python3 scripts/hla_popgen/118_novel_ladder.py --ipd-nuc ~/hla_nuc.fasta --ipd-release 3.66.0 --out-dir ~/s10 --sr-dir ~/s10_cb`
   Sanity checks before trusting it: the funnel in the log should show ~8,015 novel rows in the LR cohort and ~675 `pool_calls`; `drop:indel_*` should be the large majority of novel-protein rows (86% indel-only per S08); the number of alleles flagged `ipd_exact` should be in the tens (S08 found 54 exact signatures); if any of these is wildly off, stop and tell Marc.
5. Exact report (stays on the VM): `python3 scripts/hla_popgen/119_novel_ladder_report.py --data ~/s10/ladder_data.private.json --out ~/s10/S08_NOVEL_LADDER.private.html` and have Marc open it from JupyterLab (File Browser → right-click → Open With → Notebook/HTML viewer, or "Trust HTML"). Check in the selector that the "Everything minimal" preset total matches the 118 log.
6. Optional masked build (only if Marc asks): `119 ... --data ~/s10/ladder_data.pub.json --out ~/s10/S08_REPORT.html`. Do **not** copy it anywhere until Marc agrees.
7. Cleanup (ask first): remove `~/cb_code` and `~/hla_nuc.fasta`. Keep `~/s10`, `~/s10_cb` (mode-600 files stay).

## What to report back to Marc (exact counts are fine in this exploratory run)
- Funnel totals from the 118 log; number of alleles in the pool; number per bucket (1–4); share protein-altering vs synonymous; exact histogram of carrier counts (alleles per number of unrelated carriers, 1, 2, 3 …), all vs unrelated, and per bucket; saturation numbers at a few sample sizes.
- Every code change you had to make (diff), every anomaly, runtime per step.
- Anything you could not verify. Do not guess numbers; do not tell Marc something ran if it did not.

## Hand-offs to Marc
`git commit`/`git push` of any fixes; deciding whether anything (including `ladder_data.pub.json`) may leave the VM or be published; deleting merged remote branches (see `context/STATUS.md`).
