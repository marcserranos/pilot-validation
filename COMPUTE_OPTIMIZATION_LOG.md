# Compute Optimization Log

*Tracks wall-clock/CPU-time per configuration as tweaks are tried, so cost/runtime estimates for the eventual n=25/100 scale-up (CURRENT_SPRINT.md objectives #3/#5) are based on real measurements, not guesses. Append new entries below rather than overwriting — the point is to compare configs over time.*

**Test VM**: `HLAcalling_pilot_v0_m`, 4 vCPU / 25GB RAM. All timings below are on this single VM — **not** representative of real parallel scaling (4 vCPUs means concurrent jobs on this VM cause core contention, not real parallelism; n=25/100 scaling will need separate/larger VMs, not more concurrent jobs here).

## SpecImmune — person 1017156, PacBio HiFi (Revio), chr6-sliced FASTQ (2,904 reads after QC)

### Config: defaults (`--align_method_1 bwa`, `--align_method_2 minimap2`) — 2026-07-08

**Read-binning phase only** (measured mid-run, before typing/consensus/annotation):
| Batch | Reads | CPU-sec | Real-sec |
|---|---|---|---|
| 1 | 2,362 | 5,824.012 | 1,498.632 |
| 2 | 542 | 1,763.509 | 470.864 |
| **Total** | **2,904** | **7,587.521 (~2.11 CPU-hr)** | **1,969.496 (~32.8 min wall)** |

**Full run (read-binning + typing/consensus/annotation + visualization attempt)**: **52 minutes wall-clock** total, on 4 vCPUs. Visualization step itself fails (missing `rsvg-convert`/`inkscape`, no sudo) but this is a small tail cost after the real work is done, not counted as wasted typing time.

### Config: `--align_method_1 minimap2` (everything else identical) — 2026-07-08

**Controlled test**: same person, same FASTQ, same `-j 4`, only `--align_method_1` changed from `bwa` to `minimap2`.

**Result: 16m0.794s wall-clock (real), 38m28.560s user, 1m6.415s sys** — vs. 52 min wall-clock for the full `bwa` run. **~3.25x faster.**

**But not a free win — real accuracy tradeoff found on this person.** Diffed the full result files (`diff` on `*.HLA.final.type.result.formatted.txt`):
- **A, B, DQA1, DQB1, DPA1**: identical between configs.
- **C**: same allele calls, only a minor read-support-count difference — not a real change.
- **DRB1**: haplotype 1 identical; haplotype 2 went from an unresolved `NA` (bwa) to a specific resolved call within the same allele family under minimap2 — arguably an improvement.
- **DPB1 — real regression under minimap2**: `bwa` config called both haplotypes homozygous at **100% read identity**. `minimap2` config dropped haplotype 1 to **96.1% identity** and collapsed haplotype 2 to `NA` with a completely different, unrelated ambiguous allele list. Genuine loss of calling confidence at this locus, not noise.

**Conclusion: do not default to `minimap2` yet.** 3.25x speedup is real and significant, but it isn't free — at least one classical gene (DPB1) got measurably worse for this person, while another (DRB1) got marginally better. Unknown whether the DPB1 regression is a one-off quirk of this person's specific DPB1 sequence/coverage or a systematic pattern across people. **Next step: run `--align_method_1 minimap2` on at least one more of the smoke-test people and check specifically whether DPB1 degrades again** — if it's consistent, `minimap2` is unsafe as a blanket default and should only be used with a DPB1-specific fallback/warning; if it was a one-off, the speedup may be safe to adopt broadly.

**Current default for the reusable pipeline script: keep `bwa` (the proven, fully-validated config)** until this is resolved — don't trade a 3.25x speedup for an unquantified per-locus accuracy risk without more evidence.

**Ruled out**: no way to restrict `make_db.py`'s database to just the 8 classical genes (HLA-A/B/C/DRB1/DQA1/DQB1/DPA1/DPB1) via CLI — only whole-family selection (`HLA`/`KIR`/`CYP`/`IG_TR`/`extend`). A custom filtered `--HLA_fa` database is theoretically possible but the input format isn't documented — not worth pursuing speculatively; the full ~30-locus panel it types is more compute than we strictly need, but the CLI doesn't offer a narrower option.

## SpecHLA — person 1017156, short-read, chr6-sliced FASTQ

Timing not yet captured with this level of detail — add on a future run.
