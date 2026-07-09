# Experiments Log

> **Role:** append-only record of every pipeline run and optimization test — what was run, the result, and the runtime. **Progress = this file read top to bottom.**
> **Edit:** **append** a dated entry per run; never rewrite past entries. Raw genotype/allele values stay in `SMOKE_TEST_PICKS.local.md` — reference summary stats here (e.g. "6/8 concordant"), not raw calls.
> **Read:** for what's been tried and how long it took.

## Test VM

`HLAcalling_pilot_v0_m`, 4 vCPU / 25 GB RAM. All timings are single-VM — **not** representative of parallel scaling (4 vCPUs means concurrent jobs contend for cores, not real parallelism; n=25/100 needs separate/larger VMs, not more concurrency here).

---

## 2026-07-08 — pipeline proven end-to-end; first bake-off + optimization test

**Milestone:** full access chain works — gcsfuse mount → chr6 slice (`chr6:29.5–33.5 Mb`; CRAM needs `-T` vs. local ref, BAM doesn't) → FASTQ → typing → comparison. Wrapped into `slice_and_fastq.sh` + `compare_hla_results.py` so it never needs re-deriving.

**Person 1017156 — full 3-way bake-off complete.** 6/8 classical genes concordant across AoU-native / SpecHLA-SR / SpecImmune-LR; 2 discordances (one where AoU-native is the lone outlier — 2 methods agree against it; one where SpecImmune-LR alone disagrees). n=1 — not conclusive about method superiority. Allele detail in `.local.md`.

**Candidate swap.** Originally-picked 2037879 & 1413564 had no real aligned lrWGS BAM despite passing the eligibility flag (→ ENVIRONMENT quirk #13). Replaced with 2522883 & 1253627 (files verified present: 21.2 GB / 20.2 GB) before committing compute.

### SpecImmune timings — default `bwa` config (`--align_method_1 bwa`, `--align_method_2 minimap2`)

| Person | Real | User | Sys | Notes |
|---|---|---|---|---|
| 1017156 | ~52 min | — | — | full run; read-binning alone ~32.8 min wall / ~2.11 CPU-hr (2,904 reads) |
| 2522883 | 49m50.867s | 171m20.521s | 1m13.252s | |
| 1253627 | 48m21.136s | 165m15.980s | 1m10.806s | |

All three cluster tightly → **~50 min is a consistent per-person baseline** on this VM, not a fluke. All hit the harmless cosmetic viz-PDF failure at the end (`No converter found for conversion to pdf`) — typing results already written by then.

### SpecImmune optimization test — `--align_method_1 minimap2` (person 1017156, everything else identical)

**16m0.794s real / 38m28.560s user / 1m6.415s sys — ~3.25× faster than bwa's 52 min.** But **not free** — full-result diff (`diff` on `*.HLA.final.type.result.formatted.txt`):
- A, B, DQA1, DQB1, DPA1: identical. C: identical calls, minor read-support-count diff only.
- DRB1: haplotype 2 *improved* (bwa `NA` → resolved call, same allele family).
- **DPB1 regression:** bwa called both haplotypes homozygous at 100% identity; minimap2 dropped hap1 to 96.1% and collapsed hap2 to `NA` with an unrelated ambiguous list. Real loss of calling confidence at this locus.

**Decision: keep `bwa` default** (see DECISIONS). Retest minimap2 on ≥1 more person, checking DPB1 specifically, before any adoption.

**Ruled out:** no CLI way to restrict `make_db.py` to just the 8 classical genes — only whole-family (`HLA`/`KIR`/`CYP`/`IG_TR`/`extend`). A custom `--HLA_fa` database is theoretically possible but its input format is undocumented — not worth speculative effort.

---

## 2026-07-09 — SpecHLA-SR timing baseline (first capture)

Default SpecHLA short-read config, chr6-sliced FASTQ, 4 vCPU VM.

| Person | Real | User | Sys |
|---|---|---|---|
| 2522883 | 21m8.218s | 60m46.152s | 0m22.426s |
| 1253627 | 19m40.039s | 56m54.781s | 0m24.992s |

**SpecHLA-SR (~20 min) is ~2.5× faster than SpecImmune-LR (~50 min)** on the same people/VM — first head-to-head runtime signal between the two tools. Two benign warnings reproduced and confirmed harmless: `sh: 11: less: not found` (the unlocated 2nd `less` call, see ENVIRONMENT known-issues) and `blastn … reduced to 4 to match available CPUs` (SpecHLA hardcodes `bwa mem -t 5` on this 4-core VM). Neither affected output — full allele tables produced for all 8 classical genes for both people.

---

## Queued experiments (not yet run)

- **minimap2 on a 2nd person** — resolve whether the DPB1 regression is one-off or systematic.
- **Region-padding narrowing.** Current window `chr6:29.5–33.5 Mb` (4 Mb) carries ~440 kb pad each side beyond the 8 classical genes (real span ~29.94–33.06 Mb, HLA-A→HLA-DPB1). Hypothesis: a narrower window speeds slicing/FASTQ/typing further, on top of any binning-aligner speedup — but risks read-pair dropout and boundary alignment artifacts, especially given HLA's structural variability person-to-person and the DRB1/DQB1 phasing ambiguity already seen even with the current generous padding. Method: same controlled approach (same person, only the window changed; diff full result files, not just runtime). Don't adopt a narrower default without this check.
