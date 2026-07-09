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

## 2026-07-09 — 3-way bake-off complete for all 3 smoke-test people (n=3)

First cross-person read across AoU-native-SR / SpecHLA-SR / SpecImmune-LR, 8 classical genes × 3 people = 24 locus-calls. (Raw genotypes in `SMOKE_TEST_PICKS.local.md`; only patterns here.)

**Concordance:** full 3-way agreement 13/24 (54%) at 2-field; per person 6/8, 3/8, 4/8.

**Read the discordances by DATA-TYPE independence, not by majority vote.** AoU-native and SpecHLA both run on the *same* short-read CRAM — independent algorithms, not independent data. So when they agree against SpecImmune, that is *not* two independent witnesses converging on truth; it can equally be a shared short-read blind spot (reads too short to span the polymorphic exons, mismapping between paralogous HLA genes, reference-bias dropout). Concordance alone cannot adjudicate a 2-vs-1 split where the 2 share a data type. Given long reads' stronger prior for this region, a SpecImmune divergence from the short-read pair deserves *more* benefit of the doubt, not less. Split:
- **Short-read-pair vs. SpecImmune (A, B, C, DRB1, DPA1, DQB1):** unadjudicated — shared short-read artifact OR SpecImmune error, can't tell from concordance. Sub-distinction matters: some divergences are *deeper resolution* (e.g. C 4th-field) = expected long-read behavior, trust it; at least one is a different allele *group* (A\*02:xxx vs A\*24:xx) = warrants scrutiny, since it's not explained by depth.
- **DQA1 — the one cross-technology finding (robust):** AoU is the lone outlier 3/3, with SpecHLA (SR) *and* SpecImmune (LR) — two *different* data types — agreeing against it at 2nd field every time. This survives the independence critique → AoU-native's DQA1 calls look genuinely off. Hypothesis: AoU's ensemble leans on HLA-HD alone for DQA1 (Polysolver/OptiType are class-I-oriented).
- DPB1 is the messiest locus for everyone; here SpecHLA itself emits suspicious high alleles (63:01, 105:01) — so odd calls aren't a SpecImmune-only phenomenon.

**Coverage note (corrected):** in AoU, HiFi long-read is ~8× vs. short-read ~30× — short-read is the *deeper* data here, not shallower. But lower long-read depth ≠ worse: long reads span whole HLA genes and phase across variants short reads can't link (the actual source of their advantage). Earlier "3,200 vs 1.2M reads" framing is **retracted** — read count is a bad proxy across ~100× read-length difference. Per-person/per-locus depth NOT yet measured; do not assert coverage effects until measured from the sliced BAMs.

**Adjudication pending (both cheap, data already on VM):** (a) per-locus `samtools depth`/`coverage` on `chr6.bam` (SR) + `chr6_LR.bam` (LR) per person — do the false-homozygous / haplotype-drop cases coincide with low LR depth? (b) SpecImmune per-call confidence (`One_guess` vs ambiguous top-pick) on each divergent locus. Weight divergences by evidence, not by vote.

**Caveat: n=3, no ancestry stratification.** Directional only — no method-superiority or AoU-bias claim yet.

---

## 2026-07-09 (cont.) — confidence adjudication + measured coverage change the read

Pulled per-call SpecImmune confidence (`Match_info` = allele|score|identity, `Reads_num`, `Step1_type`, ambiguous `Genotype` list) for the divergent loci, and measured real depth. **Most "SpecImmune outliers" are not real disagreements.**

- **Coverage (measured, `samtools coverage` over chr6:29.5–33.5M):** SR ~36–40× vs LR ~12× (both ~99.7% breadth). LR is ~3× shallower; a het locus averages ~6×/allele. Corrects the earlier ~8× guess; the retracted read-count comparison is superseded by this.
- **≥2 of the 6 "SpecImmune outliers" are `One_guess` tie-break artifacts:** the short-read-consensus allele was tied on score+identity AND was SpecImmune's own `Step1_type`; `One_guess` merely reported a different tied allele. **The compare script's use of `One_guess` manufactured these discordances** → fix it to score concordance against the ambiguous list / `Step1_type`, and surface identity + read-count + tie flags.
- **DB-version confound:** SpecImmune uses IPD-IMGT/HLA **3.64.0**, SpecHLA uses **3.38.0**. Some rare-allele divergences are likely newer-DB naming, not real calling differences. Standardize/track before treating rare calls as discordances.
- **Genuine confident SpecImmune divergence:** one DPA1 locus (0.999 identity, 21 reads, `Step1` agrees, homozygous) truly differs from AoU+SpecHLA — investigate (DB/Q-allele vs real LR signal).
- **DQA1 (AoU-outlier) reconfirmed + strengthened:** SpecImmune's DQA1 calls are high-identity confident (100% in one person) and agree with SpecHLA against AoU.
- **DPB1 hard for all** (identities ~0.96–0.99); SpecHLA emits the oddest calls here, not SpecImmune.

**Bottom line:** the naive top-pick scoring made SpecImmune look unreliable; it isn't. True concordance is materially higher than the raw 54% once ties + DB-version are handled. Two real residual signals: **DQA1 (AoU off)** and **one confident DPA1 divergence**; the rest are method/DB artifacts or genuinely-hard loci. (Per-locus detail with alleles in `SMOKE_TEST_PICKS.local.md`.)

**Action items:** (1) rework `compare_hla_results.py` SpecImmune handling (ambiguous-list concordance + confidence surfacing), re-run; (2) standardize/record the IMGT DB version per tool; (3) still measure per-gene depth for the DRB1/DPB1 hard loci.

---

## 2026-07-09 (cont.) — combined aligner x padding sweep, IN PROGRESS (unattended, `run_aligner_pad_sweep.sh`)

Supersedes the two queued items below (minimap2-on-2nd-person, region-padding-narrowing) with one combined, controlled experiment: person 2522883, 5 padding levels x 2 aligners (`--align_method_1 bwa`/`minimap2`, `--align_method_2` held at default) = 9 new SpecImmune runs + 1 reused baseline, run in series via `run_aligner_pad_sweep.sh`, launched unattended (`nohup`) while Marc is away for ~5h.

**Padding design (not a naive single-window shrink):** the 8 classical genes sit in 4 clusters (A | C+B | DRB1+DQA1+DQB1 | DPA1+DPB1) separated by up to ~1.3 Mb of non-classical sequence with zero classical-gene reads — carved multi-region `samtools view` windows per cluster, not one shrinking contiguous block (a naive "outermost-gene ± n" window barely shrinks, since it still has to span both ~1.2 Mb empty gaps). Gene coordinates pulled live from Ensembl REST (GRCh38, 2026-07-09), not hand-typed. Padding floor grounded in this person's real measured LR read-length distribution (`samtools view` + awk on `chr6_LR.bam`): median 16.5 kb, p99 34 kb, max 42.6 kb, n=3332 reads — "theoretical minimum" pad = p99 rounded to 40 kb (buffers the anchor-displacement risk: a divergent read's aligner-reported anchor can land in the flank even though the read's true reference footprint overlaps the gene).

| Config | Windows | Total bp | vs. current 4.00 Mb |
|---|---|---|---|
| fullpad (current) | 1 | 4.00 Mb | baseline |
| pad200k | 3 (clusters 3+4 merge) | 1.82 Mb | 2.2x smaller |
| pad100k | 4 | 1.02 Mb | 3.9x smaller |
| pad40k ("theoretical min") | 4 | 0.54 Mb | **7.3x smaller** |
| pad5k ("risky", below the read-length floor) | 4 | 0.26 Mb | 15.1x smaller — expect real degradation |

Actual wasted-read-count (not just bp) at each level is logged automatically by the sweep script (`samtools view -c` per carved window) — answers "how much unuseful data are we introducing" empirically once it runs, rather than by estimate.

**Mechanism note (why bwa vs minimap2 differ, for the record):** minimap2 seeds on exact-match minimizers; a read from a highly-divergent HLA allele shares few/no minimizers with GRCh38 → weak/missing seed → mis-chained or dropped. bwa-mem seeds shorter exact matches (SMEMs) then does full per-base Smith-Waterman extension, which tolerates that divergence much better — at large compute cost. This is a plausible mechanistic account of the earlier DPB1 regression, not yet independently confirmed against SpecImmune's exact invocation.

**Outputs:** `~/pipeline_outputs/<id>/sweep/progress.log` (live status, one line per event), `~/pipeline_outputs/<id>/sweep/comparisons.md` (all 10 human-readable matrices), `~/pipeline_outputs/<id>/comparison_log.csv` (master long-format log, all configs, machine-readable). **Comparison script rewritten 2026-07-09** to stop declaring a "verdict"/winner (see DECISIONS.md — vote-based outlier framing was wrong when 2 of 3 methods share the same underlying short-read data) — it now reports raw calls plus SpecImmune's own identity/read-support/tie/ambiguous-list-size per haplotype, for evidence-weighted reading downstream.

**Known risk, not yet resolved:** VM auto-sleeps after ~1h idle (quirk #14); unclear whether a backgrounded compute job counts as "activity" preventing that. A mid-sweep sleep would kill the in-flight run but not lose completed ones (incremental logging). Check next session whether the full sweep completed or stalled partway — `progress.log` timestamps will show exactly where.

## Queued experiments (not yet run)

(none — both prior queued items folded into the sweep above)
