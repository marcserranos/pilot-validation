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

**Ruled out (2026-07-08), RECONSIDERED (2026-07-09) — see Experiment C below:** no CLI way to restrict `make_db.py` to just the 8 classical genes — only whole-family (`HLA`/`KIR`/`CYP`/`IG_TR`/`extend`). A custom `--HLA_fa` database is theoretically possible but its input format is undocumented — deemed not worth speculative effort at the time. Revisited given the real padding payoff evidence below: since a chunk of SpecImmune's runtime is fixed per-gene overhead (index-building across ~30 loci) independent of input size, restricting to 8 genes could be a bigger lever than padding. Queued as Experiment C.

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

## 2026-07-09 (cont.) — combined aligner x padding sweep — COMPLETE, real results (second attempt, after env fix)

**First attempt's results were garbage (do not use)** — launched from a `spechla` shell, SpecImmune's pipeline needs `sniffles` (only in `specimmune` env), every gene failed silently while `main.py` still exited 0. Fixed by always invoking via `pixi run -e specimmune` + verifying the output file exists before declaring success (ENVIRONMENT.md quirk #15). **Rerun completed cleanly, all 10 configs produced real output** (confirmed: zero "result not found" warnings, verified file exists at `fullpad_minimap2` before trusting the rest, real per-config runtimes distinct from the fake ones).

### Runtime + read-count scaling (person 2522883)

| Config | Sliced reads | % of full-pad | Runtime bwa | Runtime minimap2 |
|---|---|---|---|---|
| fullpad | 3,332 | 100% | ~49m51s (historical, 2026-07-08) | 16m44s |
| pad200k | 1,867 | 56% | 28m18s | 15m51s |
| pad100k | 1,212 | 36% | 21m11s | 14m29s |
| pad40k ("theoretical min") | 616 | 18% | 15m34s | 12m06s |
| pad5k ("risky") | 329 | 10% | 12m34s | 10m52s |

**Answers the original "how much unuseful data are we introducing" question directly: at pad100k, ~64% of full-pad's reads fall outside the carved windows (pure waste); at pad40k, ~82%.** Runtime shrinks more slowly than read count (~3.2x at pad40k-bwa vs the 7.3x bp reduction) — a real chunk of SpecImmune's runtime is fixed overhead (indexing ~30 reference FASTAs per gene, tool startup) independent of input size. minimap2's speed edge over bwa **shrinks as padding narrows** (33min gap at fullpad → 1.7min gap at pad5k) since there's less binning work left for bwa's slowness to matter.

### Per-gene accuracy across the padding sweep

**Rock-solid at every pad level, both aligners — B, C, DQA1, DPA1, DPB1, DQB1:** identical calls, identical (or near-identical) identity scores, flat read counts, all the way down to pad5k. Padding doesn't touch 6 of 8 genes in this person.

**Two genes show real signal:**
- **DRB1** — reads decline cleanly and monotonically with padding (101→94→85→80→71→66→35→32→24→21 across the 10 configs, bwa/minimap2 interleaved) — the cleanest dose-response in the dataset. Notably, the call *improved* from fullpad to pad200k/pad100k: fullpad-bwa's messy `DRB1*04:90` (identity 0.952) became a clean `DRB1*04:05:01:01` (identity 0.99) once padding narrowed — the extra ~440kb of fullpad padding appears to have been adding binning noise, not signal, for this locus. But at **pad40k and pad5k, hap2's identity score becomes unparseable/NA** — the first real degradation signal, coinciding with reads thinning to 35 and 24.
- **Gene A** — stable through pad100k (9-10 reads, ambig=15), but at **pad40k/pad5k the ambiguous-candidate list explodes to 186** (reads only drop 9→8, but specificity collapses). Same pad level where DRB1 also starts degrading.

**Real, data-grounded floor: pad100k looks safe (zero degradation seen, 3.9x smaller than current default, ~21min bwa / ~14.5min minimap2). pad40k and below show measurable degradation on the two hardest genes (A, DRB1) even though no allele flipped outright.** n=1 person — directional, not a final answer; worth confirming on a 2nd/3rd person before adopting pad100k as default.

### bwa vs. minimap2 — no accuracy difference in this person, at any pad level

Every gene gives essentially identical calls between aligners at every padding level tested — **including DPB1, which does NOT reproduce the regression seen in person 1017156** (both aligners: identity ~0.999-1.0, same calls, stable across all 5 pad levels). **Update to the open minimap2-DPB1-regression question: at n=2 people, this looks more like a one-off specific to 1017156's DPB1 alleles/coverage than a systematic minimap2 problem** — though n=2 still isn't conclusive.

**Recommendation pending wider testing:** pad100k as a new default window (safe, 3.9x smaller, real runtime savings) is a stronger, better-evidenced lever than the aligner swap alone — and once padding is already narrow, the aligner choice matters far less (only ~6-7min difference by pad100k, vs ~33min at full pad).

Full per-gene, per-config raw values (allele calls, identity, read counts, ambiguous-list size, tie flags) in `~/pipeline_outputs/2522883/comparison_log.csv` and `~/pipeline_outputs/2522883/sweep/comparisons.md` on the VM (not duplicated here — no participant genotypes in this tracked file).

**Design recap:** person 2522883, 5 padding levels (fullpad/pad200k/pad100k/pad40k/pad5k) x 2 aligners (`--align_method_1 bwa`/`minimap2`, `--align_method_2` held at default), 9 new SpecImmune runs + 1 reused baseline. Padding uses carved multi-region `samtools view` windows per gene-cluster (A | C+B | DRB1+DQA1+DQB1 | DPA1+DPB1), not one shrinking contiguous block — a naive "outermost-gene ± n" window barely shrinks, since it still has to span two ~1.2 Mb empty non-classical gaps. Gene coordinates from Ensembl REST (GRCh38, 2026-07-09); padding floor grounded in this person's measured LR read-length distribution (median 16.5kb, p99 34kb, max 42.6kb, n=3332 reads) — pad40k = p99 rounded up ("theoretical minimum"), pad5k deliberately below that floor ("risky" control). bp per config: fullpad 4.00Mb, pad200k 1.82Mb (clusters 3+4 merge at this pad), pad100k 1.02Mb, pad40k 0.54Mb, pad5k 0.26Mb.

**Mechanism note (bwa vs minimap2, for the record):** minimap2 seeds on exact-match minimizers; a read from a highly-divergent HLA allele shares few/no minimizers with GRCh38 → weak/missing seed → mis-chained or dropped. bwa-mem seeds shorter exact matches (SMEMs) then does full per-base Smith-Waterman extension, tolerating divergence much better at large compute cost. Consistent with the earlier 1017156 DPB1 regression being real for that person, but the n=2 result above suggests it isn't universal.

**Outputs:** `~/pipeline_outputs/<id>/sweep/progress.log`, `~/pipeline_outputs/<id>/sweep/comparisons.md` (10 human-readable matrices), `~/pipeline_outputs/<id>/comparison_log.csv` (master long-format log). Comparison script (rewritten 2026-07-09) reports raw calls + confidence, no verdict — see DECISIONS.md.

## Roadmap locked 2026-07-09 — four experiments, in order, for the next session

Agreed with Marc after the aligner x padding sweep. Read this whole section before starting any of them — each has design notes and quirks that took real debugging to discover once already (don't rediscover). Order below is the intended sequence: A is cheap/parallel-anytime; B and C are compute-optimization investigations that should ideally land before D's big spend, but D is **not blocked waiting indefinitely** on them (see D's own notes). The two standalone confirmation items from the previous entry (pad100k floor, minimap2/DPB1 non-regression on a 3rd person) are folded into Experiment B/D's scope, not separate tasks.

### Experiment A — AoU-native distribution stats (cheap, no new compute, do first)

Goal: characterize the AoU-native HLA callset itself, before treating it as a fixed comparison baseline at ancestry scale (Experiment D).

Compute (all just pandas over the already-known TSV):
- **Per-locus missingness/NA rate** across the 8 classical genes. Cross-check against our own hardest loci (DRB1, DPB1) — if AoU's ensemble *also* struggles most there, that's independent corroboration those genes are just hard; if AoU is clean there while we aren't, that's a methodology-specific weakness worth flagging.
- **Resolution-depth distribution per locus** (2-field vs 3-field vs 4-field) — quantify, don't assume from the earlier "variable, 2-3 field" note in ENVIRONMENT.md.
- **Real non-null coverage of the TSV** — how many `research_id`s actually have a usable call. Never actually measured; only assumed "near-full" from file size.
- **Homozygosity rate per locus** as a sanity/reference-bias check.
- **Bonus (needs one new lookup first):** join against AoU's genetic-ancestry TSV for a first-cut allele-frequency-by-ancestry preview — a free early look at the project's actual Aim 1 deliverable. **The ancestry TSV's exact bucket path has never been confirmed** (unlike `hla_genotypes.tsv`, whose path was corrected this session via `gsutil ls` after being wrong for a whole session) — locate it the same way: `gsutil ls` exploration under `v9/wgs/` or wherever AoU's ancestry outputs live, don't guess the path from memory.

Quirks: the TSV is 331MiB and static — consider `gsutil cp`-ing it once to local persistent disk (same pattern as the reference genome) rather than depending on the gcsfuse mount every session; removes a whole class of session-start friction for repeated analysis. Corrected mount path: `v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv` (ENVIRONMENT.md).

### Experiment B — SpecHLA (short-read) padding sweep

Goal: test whether narrowing the short-read slicing window reduces binning noise the way it did for long-read DRB1 (see the 2026-07-09 sweep above), before hitting a real degradation floor. Marc's explicit design requirements: recompute the *real* theoretical minimum for short-read specifically (don't reuse the long-read pad40k value), add more intermediate levels than the 5 used for LR, and add one level below the new floor to positively confirm degradation.

**Do not reuse the long-read padding levels (fullpad/200k/100k/40k/5k) — they're calibrated to 15-20kb HiFi reads, not ~150bp short reads, and would almost certainly never reach the real short-read floor.**

Design:
1. **Measure the real short-read insert-size/fragment-length distribution first** — reuse `bwa mem`'s own `mem_pestat` output (already emitted in every SpecHLA log we've seen, e.g. "(25,50,75) percentile" lines) or run `samtools stats` on the sliced SR BAM. Get median/p95/p99/max, exactly the same rigor as the LR read-length measurement that set pad40k (median 16.5kb, p99 34kb, max 42.6kb).
2. **Recompute a new short-read "theoretical minimum"** from this real distribution — expect it in the low-hundreds-of-bp to low-kb range, nowhere near 34kb.
3. **More intermediate levels than the LR sweep's 5** — finer granularity between fullpad and the new floor (Marc's explicit ask).
4. **One level below the new floor**, to positively confirm degradation occurs rather than just infer it — same "risky control" role pad5k played for LR.
5. **New QC check not done for the LR sweep**: verify properly-paired fraction / mate-dropout rate at each candidate window. Flagged as a risk in `reference/AOU_DATA_ACCESS_NOTES.md` months ago and never actually tested — a narrow window can silently drop one mate of a pair, which the LR sweep didn't need to worry about (long reads aren't paired).
6. Reuse the same 4 gene-cluster coordinates already pinned via Ensembl REST (A | C+B | DRB1+DQA1+DQB1 | DPA1+DPB1) — no need to re-derive.
7. Run on person 2522883 (same as the LR sweep) for direct comparability, same single-variable-at-a-time methodology.

### Experiment C — SpecImmune gene-panel restriction (investigate feasibility, potentially bigger than B)

Goal: SpecImmune currently types ~30 loci (8 classical + 22 non-classical/pseudogenes we never use downstream). If restricting to just the 8 is buildable, this could be a bigger lever than padding — a real chunk of its runtime is fixed per-gene overhead (index-building, seen as dozens of `bwa index`/`minimap2` calls per gene in every log) that scales with gene *count*, not input size.

Previously ruled out (2026-07-08) as "not worth speculative effort" — reconsidered now given real padding payoff evidence.

Steps:
1. **Read `~/tools/SpecImmune/scripts/make_db.py` directly** (not just `--help`) to determine whether a custom, filtered `--HLA_fa` database (8 genes only) is actually buildable, and what input format it expects.
2. **Real open question to resolve, not assume:** SpecImmune builds `individual_ref/<GENE>/` per-person reference FASTAs for genes we never asked about (HLA-DRB5, HLA-T, MICA, etc., seen in our own logs) — this suggests genes are typed based on candidates found against a fixed full-panel DB. Need to determine whether restricting the DB itself actually stops non-classical genes from being processed, or whether some other internal logic always walks the full panel regardless (in which case this needs a deeper code change, not just a smaller DB file — don't assume it'll just work).
3. **If feasible:** build the restricted DB, re-run the same controlled single-variable timing methodology (same person, same FASTQ, only the DB changed) to quantify the actual speedup before adopting it.

### Experiment D — Fused ancestry-stratified 3-way comparison (this IS the n=25-100 scale-up)

**Decided 2026-07-09: this replaces BOTH the earlier "vertical/horizontal scaling for n=25/100" question and the separately-proposed "n=5-10-per-ancestry 3-way comparison" — they are the same experiment, not two.** The scaling-approach question (parallel VMs vs. single-VM optimization) folds into this experiment's execution plan rather than being a separate preliminary step.

Design:
- Compare AoU-native-SR / SpecHLA-SR / SpecImmune-LR across ~5-10 people per ancestry group, using AoU's genetic-ancestry TSV categories (AFR/AMR/EAS/EUR/MID/SAS/OTH — 7 categories, per DECISIONS.md's existing preference for genetic ancestry over self-reported race).
- **Candidate pool is already settled (confirmed by Marc, 2026-07-09):** an existing matched SR+LR cohort already covers adequate people across ancestries — most likely the previously-built `HLA_Pilot_Matched_SR_LR_100` cohort (9,358 candidates, originally selected via 6 self-reported races) or a newer resource. No further cohort-building or feasibility check needed before starting. Still apply the `/revio/` BAM-path filter (quirk #13) per candidate when actually picking individuals — that's a per-person data-quality check independent of which named cohort is used.
- **Must use the confidence-weighted concordance approach already built into `compare_hla_results.py`** (identity/reads/tied/ambiguous-list fields) — never a naive vote-based metric. This was a hard-won, expensive lesson (see the "confidence adjudication" entry above) — repeating it at 10-40x the people would be far more costly to catch.
- **Sequencing: pilot on ONE ancestry group first** (5-10 people), confirm the pipeline/metric works cleanly at that scale, then scale to the remaining ~6 groups.
- **Decided 2026-07-09: proceed to this pilot without waiting on the "supervisor conversation" open item** (DECISIONS.md). Compute cost is low enough (~$10, well under 40 hours serial even fully unoptimized) that gating a cheap pilot behind a conversation that can happen in parallel isn't worth the delay. The supervisor conversation remains valuable and stays open — it's just no longer a blocking prerequisite for starting D.
- **Use whichever optimized config validates from Experiments B/C by the time D starts** (padding level, aligner, gene-panel restriction) — but don't block D indefinitely waiting on B/C either; launch the pilot once a reasonable attempt at both has been made. Marc's expectation: if B and C both land, total compute drops from the current ~38-55hr estimate to "under a day," under $10 — a working target, not yet confirmed.
- **Explicitly designed to help resolve two carried-over open findings** from the original 3-person bake-off, not just produce a generic concordance number: is AoU-native's DQA1 discordance (3/3 people so far) ancestry-correlated or universal? Is the one confident SpecImmune DPA1 divergence ancestry-specific or general? Look for these patterns specifically.
- Rough compute estimate pre-optimization: ~56-80 people x ~41min/person (SpecHLA + SpecImmune at pad100k) ~ 38-55 hours serial.
