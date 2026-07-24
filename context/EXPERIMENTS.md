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

## 2026-07-10 — Experiment A complete: AoU-native distribution stats (n=535,658)

Ran `experiment_a_aou_stats.py` over the full `hla_genotypes.tsv` (535,658 research_ids) joined to `ancestry_preds.tsv` (535,658/535,658 joined, 100%).

**Missingness — effectively zero, not "assumed near-full."** All 8 classical genes: 0.0% missing (either or both alleles). 535,649/535,658 (99.998%) fully resolved across all 8 loci. The earlier "never actually measured" note in the roadmap is resolved: coverage is as complete as data gets.

**Resolution ceiling: capped at 3-field, universally, zero 4-field ever observed.** Every gene resolves to 3-field in 99.4-100% of calls, with a 0-0.6% 2-field remainder, and **not a single 4-field call across the entire cohort**. This answers the open "2-field vs 4-field" question in DECISIONS.md for the AoU-native arm specifically: it structurally cannot deliver 4-field resolution — that has to come from SpecHLA/SpecImmune if the downstream aims need it.

**Homozygosity per locus — DPA1 is a striking outlier.** A: 11.4%, B: 5.65%, C: 9.43%, DRB1: 7.53%, DQA1: 18.28%, DQB1: 14.4%, DPA1: **57.43%**, DPB1: 20.45%. DPA1 sits at ~3x the next-highest locus. New candidate mechanistic lead for the one confirmed SpecImmune-DPA1 divergence from the 3-person bake-off (2026-07-09 entry above): either DPA1 genuinely has low allelic diversity in this cohort (plausible, known in HLA population genetics), or AoU's ensemble has a reference-bias tendency toward over-calling homozygosity at this locus. Worth raising as a specific testable hypothesis, not just an unexplained discrepancy.

**DQA1 by genetic ancestry (first ancestry-resolved look at the standing open finding):** missingness stays 0.0% in every group (not an availability problem). Homozygosity: afr 18.08% (n=100,798), amr 19.05% (n=107,928), eas 18.25% (n=15,893), eur 18.03% (n=302,712), sas 18.33% (n=6,176) — all tight. **mid 24.17% (n=2,151)** — elevated, but n is much smaller than the other groups, so directional only. Does not yet resolve whether the AoU-native DQA1 discordance itself (vs. SpecHLA/SpecImmune) is ancestry-correlated — that needs Experiment D's actual 3-way comparison; this is only a distributional check of AoU-native's own self-consistency by ancestry.

**Presentation asset:** 4-panel visual built (stat tiles + resolution-depth stacked bar + homozygosity bar with DPA1 highlighted + DQA1-by-ancestry bar with mid highlighted) — see chat history 2026-07-10, not saved as a repo file (regenerate from the numbers above if needed for the deck).

Script: `experiment_a_aou_stats.py` (repo root). Outputs are aggregate-only (no participant-level data printed).

---

## 2026-07-10 (cont.) — Experiment A extension: locus x ancestry heatmap + PCA scatter

Ran `experiment_a_ext_ancestry_clustering.py` (full 535,658 for the heatmap, 50,000-person subsample for the scatter).

**DPA1's homozygosity "outlier" status is not locus-uniform — it's a strong ancestry gradient.** afr 27.0%, amr 60.1%, eas 39.0%, eur 67.7%, mid 64.6%, sas 50.3% — a ~2.5x range, tracking almost exactly the expected Out-of-Africa serial-bottleneck direction (African populations retain the deepest ancestral genetic diversity; non-African populations show reduced diversity from the bottleneck, which shows up as more homozygosity at a locus with few common alleles to begin with). This direction argues *against* simple reference-panel bias as the primary driver — GRCh38's EUR lean would be expected to push the opposite way, inflating apparent homozygosity in populations most divergent from the reference (afr), not least. **Reframes the earlier "confident SpecImmune-DPA1 divergence" finding (2026-07-09):** the underlying AoU-native homozygosity pattern at DPA1 looks like real population structure, not an obvious AoU-native calling artifact — though this doesn't resolve whether the specific SpecImmune-vs-AoU-native divergence itself is ancestry-driven; that still needs Experiment D.

A weaker version of the same afr-low/other-high gradient shows at DPB1 (13.5% → 19.8-23.9%), A (6.3% → 9-14%), and DRB1 (6.1-6.2% afr/amr → 8.3-10.9% elsewhere) — consistent with genome-wide diversity differences by ancestry, not something DP/DQ-region-specific.

**DQA1 (the standing AoU-outlier locus from the 3-person bake-off) does NOT show this gradient** — flat at 18.0-19.1% across afr/amr/eas/eur/sas, with only mid elevated (24.2%, n=195 in this subsample / 2,151 full cohort — small group, noisier). Reinforces the earlier hypothesis that DQA1's discordance is method/database-specific (AoU's ensemble likely leaning on HLA-HD alone at this locus), not a population-diversity signal the way DPA1's is.

**PCA scatter (n=50,000 subsample):** left panel (colored by ancestry) reproduces the expected textbook ancestry-PCA fan structure — afr and eur as separated poles, amr spread between them (consistent with real-world admixture), eas as a tight distinct cluster, sas/mid clustering near the eur-amr boundary. Good sanity check that the ancestry variable behaves as expected; not novel on its own. Right panel (colored by per-person cross-locus homozygosity count, 0-8) shows no obvious strong visual clustering at this resolution/alpha — the real signal is in the aggregated heatmap, not visible at the individual-point level.

**Caveat:** all of this is still AoU-native's *internal* self-consistency by ancestry — it doesn't yet tell us whether the SpecHLA/SpecImmune-vs-AoU-native discordances themselves are ancestry-correlated. That's still Experiment D's job.

Outputs (VM-local, not committed — view via Jupyter file browser): `~/pipeline_outputs/experiment_a_ext/locus_ancestry_homozygosity_heatmap.png`, `~/pipeline_outputs/experiment_a_ext/pca_ancestry_homozygosity_scatter.png`.

---

## 2026-07-10 (cont.) — Experiment B complete: SpecHLA (short-read) padding sweep

Person 2522883, same 4 gene-cluster coordinates as the LR sweep, single variable (padding only — SpecHLA hardcodes bwa internally, no aligner axis for this tool). Real fragment-size distribution measured first (`samtools stats` IS histogram, n=621,201 pairs): median 419bp, p95 586bp, p99 649bp, max 688bp — tight, as expected for standard short-insert Illumina WGS. Derived floor_pad=650bp (p99 rounded up) and subfloor_pad=325bp (floor/2, risky control) — both far smaller than the LR sweep's 34kb floor, confirming the roadmap's warning that reusing LR pad levels here would have been meaningless.

### Runtime + read-count scaling

| Pad level | Reads sliced | Properly paired | Runtime (real) |
|---|---|---|---|
| pad500000 | 1,101,569 | 98.47% | 18m33.8s |
| pad200000 | 561,771 | 98.33% | 11m40.3s |
| pad100000 | 303,398 | 98.08% | 8m6.9s |
| pad40000 | 156,671 | 97.83% | 5m38.2s |
| pad10000 | 87,236 | 97.56% | 4m35.2s |
| pad2000 | 70,719 | 97.74% | 4m22.3s |
| pad_floor_650 | 67,522 | 97.82% | 4m14.5s |
| pad_subfloor_325 | 66,832 | 97.81% | 4m14.4s |

(fullpad baseline reused from 2026-07-09, ~21m8s — a single contiguous 4Mb window, not directly read-count-comparable to these carved multi-region windows.)

**Runtime plateaus hard below ~2-10kb padding, around 4m14s** — a ~4.4x speedup from pad500000 to the subfloor. Read count plateaus even harder (67,522 → 66,832, effectively flat from pad2000 down to the risky subfloor). Short-read mirror of the LR sweep's fixed-overhead finding: SpecHLA has a substantial runtime floor independent of input size once the window is already narrow.

**Mate-pair dropout is real but small, and never propagates into call changes.** Properly-paired % declines mildly and non-monotonically from 98.47% (pad500000) to a minimum of 97.56% (pad10000), then recovers to ~97.8% at the narrowest levels — a real ~1-point range, not the cliff the roadmap worried about. Untested hypothesis: coarser windows capture more non-classical-gene flanking sequence (repetitive/pseudogene-dense), which may itself carry a different baseline pairing rate than the gene bodies — worth checking on a 2nd person.

**Zero call changes at any of the 8 genes, across all 9 configurations (fullpad down to 325bp).** Every gene's SpecHLA allele calls are identical from the full 4Mb window down to 325bp padding, including the historically hardest loci (DRB1, DPB1). Materially different from the LR sweep, which found real degradation at DRB1 and Gene A at narrow padding. Short fragments (~400-650bp) need only ~2kb of padding around the gene bodies to capture essentially all informative read pairs (read count is already >95% saturated by pad2000 vs pad500000) — padding beyond that adds flanking reads, not gene-body signal.

**Recommendation: pad2000 (or pad10000 for extra safety margin) looks like a strong, well-evidenced new default for SpecHLA specifically** — ~4.3-4.6min runtime vs ~18.5min+ at pad500000, no observed accuracy cost. **n=1 person — directional, confirm on a 2nd/3rd person before adopting**, same caveat as the LR sweep's pad100k call. **This is a SpecHLA-specific recommendation — do not conflate with the LR sweep's pad100k finding for SpecImmune; the two tools/read-types have entirely different floors (650bp vs 34kb) and were evaluated on different window mechanics.**

**Cosmetic-only note:** at pad_subfloor_325, gene C's two SpecHLA alleles print in swapped order (`C*07:02:01:01,C*01:02:01:01` vs every other config's `C*01:02:01:01,C*07:02:01:01`) — same unordered genotype, not a discordance, just hap1/hap2 label ordering (SpecHLA doesn't guarantee consistent ordering across independent runs).

Outputs: `~/pipeline_outputs/2522883/spechla_sweep/progress.log`, `comparisons.md` (9 matrices); `~/pipeline_outputs/2522883/comparison_log.csv` (master long-format log, now spans both the 2026-07-09 LR sweep and this SR sweep).

## 2026-07-10 (cont.) — correction to Experiment B: wrong QC metric, and a missing positive control

Marc pushed back on the "zero degradation, even below the floor" result as suspicious, correctly. Two real problems found on review, not just a shallow write-up:

**(1) The mate-dropout QC metric was wrong, not just imperfect.** `run_spechla_pad_sweep.sh` measured `samtools flagstat`'s "properly paired %" on the sliced BAM — but that flag is set during the *original* genome-wide alignment and is structurally insensitive to this script's own region-slicing. It could never have detected the thing it was built to check. The real signal was already sitting unused in the log: `samtools fastq`'s own "discarded N singletons" count. Recomputed as a rate (discarded / reads sliced): pad500000 1.21%, pad200000 1.29%, pad100000 1.55%, pad40000 1.82%, pad10000 2.19%, pad2000 2.04%, pad_floor_650 1.96%, pad_subfloor_325 2.02%. **There is a real mate-dropout signal — it roughly doubles as padding narrows** — the original write-up's "mate-dropout is real but small" conclusion happened to be directionally right, but for the wrong reason (a metric that couldn't see it, not one that saw a small effect). Fixed in `run_spechla_pad_sweep.sh` for future reruns; the fix is a general lesson (verify a QC metric actually responds to the manipulation being tested), not SpecHLA-specific.

**(2) The padding sweep only ever tested the flanking margin, never the gene body itself.** Gene clusters are large relative to the padding range tested (Gene A ~8.3kb; the DRB1/DQA1/DQB1 cluster ~90kb) — shrinking padding from 500,000bp to 325bp changes the margin around the gene, not the gene body where the diagnostic exons live. Once fragments are short enough (~650bp) to not need much margin, there's no mechanistic reason to expect further degradation from margin alone. So "zero degradation at 325bp padding" is plausible and not obviously a bug — but nothing in the original design could tell that apart from something silently broken (e.g. cached results being reused). **Fix: a genuine positive control, `run_spechla_truncation_sanity_check.sh`**, which truncates *into* each gene cluster around its midpoint (widths 5000/1000/650/300/100bp, deliberately below what a single fragment needs) rather than shrinking the margin. If this also shows zero degradation, that's real evidence of a pipeline bug; if it shows real degradation (missing calls, garbage identity, or SpecHLA erroring outright), that validates the original padding-sweep result rather than undermining it. New helper subcommand: `spechla_pad_helpers.py truncated-windows --width W`. Not yet run — queued as the next immediate step, ahead of Experiment C/D.

## 2026-07-10 (cont.) — truncation sanity check run: degradation confirmed, but a new success-check gap found

Ran `run_spechla_truncation_sanity_check.sh` on person 2522883 (widths 5000/1000/650/300/100bp, centered on each gene cluster's midpoint).

**Real, monotonic degradation confirmed — the positive control worked.** Genes-with-a-real-call collapses from 2/8 (trunc5000) to 1/8, still wrong (trunc1000/650/300) to 0/8 (trunc100). At trunc1000/650/300, Gene A's SpecHLA call changes to `A*33:03:01:01` in place of the established `A*31:01:02:01` — a genuinely different, wrong allele under starvation, not just a missing call. At trunc5000, DPB1 returns a novel homozygous call (`DPB1*04:01:01:13/DPB1*04:01:01:13`) never seen in any prior config across this whole project. This validates the earlier padding-sweep conclusion rather than undermining it: margin-only shrinkage (500,000bp → 325bp around already-large gene bodies) genuinely doesn't touch the genes, and the pipeline demonstrably *can* and *does* break when the window is actually inadequate.

**New problem found: the earlier "SpecHLA has a ~4m14s fixed runtime floor" claim (this file, same day) is wrong as stated.** These truncated runs — with similar or fewer reads than the padding sweep's narrowest configs — finished in 22-62 seconds, 4-11x faster than the claimed floor, despite `hla.result.txt` existing and the file-existence check passing. The likely mechanism: SpecHLA silently skips/aborts per-gene processing (phasing/consensus, presumably) when a gene has too few reads, finishes fast, exits 0, and writes a result file that's mostly `-`/empty per gene. **The "does the expected output file exist" check cannot detect this** — the file existed at every single one of these configs. Same category of lesson as the SpecImmune/sniffles bug (2026-07-09): exit 0 and file-existence are necessary but not sufficient; need a per-gene completeness signal too. Revises the earlier fixed-overhead claim: it's only fixed in a regime where every gene has enough reads to trigger full processing — it collapses to sub-minute once genes are starved and silently return empty.

**Open, unexplained:** DPA1 and DPB1 share the identical tiny window/read pool at trunc5000, yet DPA1 returned empty while DPB1 returned a full (and novel) call. No explanation yet — flagged, not resolved.

**Action item:** any future automated SpecHLA sweep should count non-`-` genes per config (or check `hla.result.txt`'s actual line count/content), not just confirm the file exists. Not yet fixed in `run_spechla_pad_sweep.sh`/`run_spechla_truncation_sanity_check.sh` — do this before the n=2/n=3 confirmation runs.

---

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

## 2026-07-10 (cont.) — Experiment C investigated locally: restriction confirmed feasible, two real upstream bugs found and fixed

Cloned SpecImmune from GitHub locally (not the VM) and read the actual source, per Marc's explicit instruction to investigate and iterate freely on a local copy before touching the VM, since this optimization is bottlenecked by SpecImmune's own code, not by anything VM-specific.

**Mechanism confirmed by reading the source, not guessing:** `make_db.py`'s `--HLA_fa`/`--HLA_exon_fa` flags accept a custom reference FASTA in place of the auto-downloaded IMGT file — `create_HLA_directories_and_save_sequences()` builds one directory + index per gene found in *whatever FASTA it's given*, nothing else. Every per-gene stage in the actual typing pipeline (`read_binning.py`, `long_read_typing.py`, `refine_typing.py`, `select_best_reference_alleleV2.py`, `visualization.py`) discovers its gene list dynamically via `get_folder_list(db_folder)` — a plain `os.listdir()` on the `--db` directory — confirmed by `db_objects.py`'s `My_db` class (`full_db`, `cal_gene_len()`, etc. all keyed off the same directory listing). The one file still using the hardcoded 37-gene list (`get_ref.py`) is dead code — its only call site in `main.py` is commented out. **A restricted DB genuinely restricts processing, with no code changes needed to the typing pipeline itself.** This reverses the 2026-07-08 "ruled out" assessment.

**Two real, previously-undiscovered bugs found and fixed (not specific to gene restriction, but the first one blocked it entirely):**
1. `make_db.py`'s `--HLA_fa`/`--HLA_exon_fa` branches referenced `local_release_version`/`local_g_group_annotation` before assigning them — `UnboundLocalError`, 100% reproducible on every call. **This is almost certainly why the 2026-07-08 assessment called it "not worth speculative effort"** — the documented lever was broken from the start, so nobody could have tried it successfully. Fixed by assigning the paths in that branch too (mirrors the default-download branch).
2. `alignment_modules.py`'s `cout_read_num()` opened a gzipped FASTQ with plain `open(path, "rb")` and counted raw compressed bytes as newline-delimited lines — a nonsense count on any platform, silently miscalibrating the subsample-vs-copy branching in `subsample_fastq()`. Fixed to use `gzip.open()`; also switched the shell-out from `zcat` to `gzip -dc` (more portable — macOS's BSD `zcat` failed outright on this exact invocation locally, producing empty subsampled files and cascading into a `TypeError` deep in `refine_typing.py`).

Both fixes are packaged as `patch_specimmune_for_gene_restriction.py` (idempotent, checks before touching a file) — same "patch lives on the VM, documented here, not upstream" pattern as the existing SpecHLA `less`→`cat` patch (ENVIRONMENT.md Fixes list).

**Real, measured local evidence (SpecImmune's own shipped `test/HLA/test_HLA_lite.fastq.gz` + known-truth `test.HLA.hap.alleles.txt`, both DBs built from the identical downloaded IMGT release for a clean comparison):**
- Restricted DB (8 genes): directory build confirmed exactly 8 gene subdirectories under `HLA/` and `HLA_CDS/` (no MICA/TAP1/HFE/DRB3-5/etc.) — **built in ~107s**.
- Full DB (~37 genes): confirmed all ~37 gene subdirectories — **built in ~225s**.
- **~2.1x faster at the DB-build stage alone**, before even running a typing comparison.

**Full local typing-run comparison hit a wall — a real, but orthogonal, macOS-specific issue, not a restriction problem.** After fixing both bugs above and installing minimap2/samtools/bwa/blast/seqtk locally (none were present), typing runs against both DBs hit a third, different failure: a malformed line in a per-gene depth file (`ValueError: invalid literal for int() with base 10: '0HLA-B*40:10:01:02'`), landing on a *different* random classical gene each rerun (HLA-B, then HLA-C) — consistent with a nondeterministic race under parallel (`-j 6`) execution, most likely a `samtools sort`/temp-file collision specific to this local macOS/Homebrew samtools build. Traced far enough to rule out a code-level shared-file bug (each gene's depth file has its own unique path) and to confirm it hits *both* the full-panel and restricted DB equally (so it's not something my change introduced). Not chased further: single-threaded (`-j 1`) avoids the race but is impractically slow locally (~75s per gene's minimap2 step alone on this test machine), and this exact tool/samtools combination has already completed dozens of real runs on the VM this project without ever hitting this — the VM environment is the right place to finish this validation, not further local debugging with diminishing returns.

## 2026-07-10 (cont.) — Experiment C, bwa leg complete: restricted DB was SLOWER and less accurate at 2/8 genes (n=1)

Ran `run_experiment_c_comparison.sh 2522883` (restricted DB, same FASTQ/padding/aligner as the pad100k_bwa baseline — DB is the only variable). **This disconfirms the "restriction is a free lunch" hypothesis, at least at n=1** — reporting it straight, not spun.

**Runtime: 26m3s, slower than the pad100k_bwa baseline's 21m11s** — a ~5min regression, opposite direction from the hypothesis.

**Read counts per gene inflated substantially in the restricted run**, most at A (9→38, +29) and DRB1 (71→120, +49); milder at DQB1 (+7), DPB1 (+5), DQA1 (+5); flat at B/C/DPA1. Working hypothesis: reads ambiguous between a classical gene and one of its removed pseudogene/paralog relatives (A has several close class-I pseudogene relatives — H, J, K, L, N, P, V, Y, T, U, W; DRB1 has DRB3/4/5) get partly absorbed by those non-classical competitors in the full panel. Remove the competitors and the same reads reassign to the nearest classical gene — more volume to process per gene, apparently outweighing the time saved by not indexing ~22 unused genes. Not proven, just the most plausible mechanism given the pattern.

**Real accuracy cost at 2 of 8 genes, not just noise:**
- **Gene A**: restricted called `A*24:667`/`A*32:101Q` — completely different from the baseline's `A*02:804`/`A*31:01:02:01`, and identity is blank (uncalculated/unclear) on both haplotypes where baseline had a real score on one.
- **DRB1**: restricted called `DRB1*11:01:01:05`/`DRB1*04:03:01:01` vs baseline's `DRB1*04:05:01:01`/`DRB1*15:01:01:01` — different alleles, again blank identity. SpecHLA's independent short-read call (`DRB1*15:02:01:01`/`DRB1*04:05:01:01`) matches the *baseline* SpecImmune call much better than the restricted one — the baseline looks more trustworthy here.
- **DPB1**: hap1 call unchanged; hap2 identity dropped (0.9997→0.9873) and flipped from untied to tied among 15 candidates — a real confidence loss, not a changed top pick.
- **B, C, DQA1, DQB1, DPA1 (5/8): unaffected** — same calls, same or near-identical identity.

**Mechanism confirmed working structurally, independent of the accuracy result:** `individual_ref/` under the restricted run contains exactly the 8 classical gene directories, nothing else — the restriction itself does what it's supposed to; the cost shows up in read reassignment, not leakage.

**n=1 — directional only.** Same-FASTQ minimap2 leg (`run_experiment_c_comparison_minimap2.sh`, reusing the pad100k_minimap2 baseline) queued next, run sequentially after this one (not concurrently, to keep the timing measurement clean). If minimap2 shows the same pattern (slower + A/DRB1 miscalls), that's a much stronger signal restriction genuinely costs something at these two loci specifically, not a one-off.

## 2026-07-10 (cont.) — Experiment C, minimap2 leg complete: same pattern, cross-aligner confirmed — restriction does not work as hoped

Ran `run_experiment_c_comparison_minimap2.sh 2522883` (sequentially, after the bwa job fully finished, to keep timing clean). **Confirms the bwa result was not a fluke.**

**Runtime, both aligners slower under restriction:**

| Aligner | Baseline (full DB) | Restricted DB | Delta |
|---|---|---|---|
| bwa | 21m11s | 26m3.4s | +4m52s (+23%) |
| minimap2 | 14m29s | 17m9.9s | +2m41s (+18.5%) |

**Read-count inflation is nearly identical in shape across both aligners** — dominated by Gene A (+29 bwa / +28 minimap2) and DRB1 (+49 / +52), everything else flat or mild (DQA1 +5/+10, DQB1 +7/+2, DPB1 +5/+2, B/C/DPA1 ~0). Same genes, same magnitude, regardless of aligner — this is the DB restriction driving it, not aligner behavior.

**Accuracy hit is cross-aligner reproducible, which is the important part:**
- **Gene A: both aligners converge on the exact same wrong alleles under restriction** — `A*24:667`/`A*32:101Q` in both bwa and minimap2, replacing each aligner's own clean baseline call (`A*24:02:01:01`/`A*31:01:02:01`-family, matching SpecHLA). Identical wrong answer regardless of aligner rules out aligner noise as the explanation.
- **DRB1: both aligners' second haplotype shifts from the `*15` allele family to the `*11` family** under restriction — not byte-identical alleles between aligners, but the same group-level shift in both.
- **B, C, DPA1: clean and unaffected in both aligners** — robust to restriction.
- **DQB1 (minimap2) / DPB1 (bwa): mild hap2 confidence erosion** (identity drop + newly tied) — present at one locus per aligner, not symmetric across both — a weaker, secondary effect.

**Mechanism now well-supported by the specific genes affected, not just plausible:** A and DRB1 are exactly the two classical genes with the most/closest pseudogene relatives in the removed panel (A: ~10 related class-I pseudogenes — E/F/G/H/J/K/L/N/P/S/T/U/V/W/Y; DRB1: DRB3/DRB4/DRB5, near-identical paralogs). Removing those decoy references, which previously absorbed ambiguous cross-mapping reads in the full panel, appears to redirect those reads onto the wrong allele within the remaining classical gene.

**Conclusion: gene-panel restriction to exactly these 8 genes, as currently implemented (simple FASTA filtering), should NOT be adopted.** Slower in both aligners (+18-23%) and produces real, reproducible miscalls at 2 of 8 loci (A, DRB1) — the opposite of the original hypothesis that motivated this experiment. n=1 person, but cross-aligner agreement on both the runtime direction and the specific genes/alleles affected makes this a substantially stronger signal than either run alone.

**Options going forward (not decided, for discussion):**
1. Drop gene-panel restriction as a lever entirely — treat this as a real, informative negative result, move to Experiment D.
2. Design a "smarter" restricted panel that keeps the specific pseudogene/paralog competitors that appear to matter for A and DRB1 (e.g., retain E/F/G/H/J/K/L/N/P/S/T/U/V/W/Y for A, DRB3/4/5 for DRB1) while still dropping the genes with no apparent interaction (MICA/MICB/TAP1/TAP2/HFE, DMA/DMB/DOA/DOB/DRA, DQA2/DPA2/DPB2) — could preserve some runtime benefit without the accuracy cost, but is real new design/validation work, not a quick follow-up.
3. Confirm this pattern on a 2nd/3rd person before finalizing either call — still n=1.

---

**Ready for the VM, in order:**
1. `patch_specimmune_for_gene_restriction.py` — apply both fixes to `~/tools/SpecImmune` (idempotent).
2. `build_restricted_specimmune_db.sh` — filters the *already-downloaded* IMGT source FASTAs sitting in the existing `~/tools/SpecImmune/db/HLA/` (no re-download) to the 8 classical genes, builds `~/tools/SpecImmune/db_classical8/`, and verifies only 8 gene directories exist.
3. `run_experiment_c_comparison.sh <person_id>` — reuses the exact same FASTQ as the existing pad100k-bwa baseline from the LR sweep (isolating the DB as the sole variable), runs SpecImmune against the restricted DB, checks per-gene completeness (same lesson as Experiment B — file existence isn't sufficient), lists what actually landed in `individual_ref/` (directly answers whether restriction leaks non-classical genes), and diffs against the existing baseline via `compare_hla_results.py`.

---

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

## 2026-07-11 — Experiment D pipeline built and staged for launch (results pending)

Three new scripts implement the fused ancestry-stratified 3-way comparison as a single unattended, resumable run. No results yet — this entry records the design + config so the eventual results entry has its provenance.

**Scope chosen: ~8-10 people × all 6 genetic-ancestry groups (~56), across multiple nights via resumability** (departs from the roadmap's literal "pilot one group first" — the resumable orchestrator + per-person isolation + ntfy person-1 ping make a full-cohort launch safe: a bad person can't sink the run, and person 1's result is checkable ~25 min in before trusting the batch).

**Config (applies the resolved B/C decisions, no new decision):** SpecHLA short-read at **pad10000** (Experiment B validated zero degradation across the pad2000-pad10000 band; 10k for margin), SpecImmune long-read at **pad100k / `--align_method_1 bwa`** (Experiment C's validated LR default). **Gene-panel restriction deliberately NOT used** — rejected in Experiment C (slower + miscalls at A/DRB1). Each person therefore doubles as an "n≥1 more person" confirmation of B's and C's n=1 padding optimizations; watch the concordance for any config-induced drift.

**Robustness design (the point of the build — "launch once and forget," incl. across nights):**
- `build_experiment_d_cohort.py` derives the cohort straight from the mounted v9 manifests (CRAM `manifest.csv` + LR `manifest.tsv` + ancestry `ancestry_preds.tsv`), applies the `/revio/` filter (quirk #13), inner-joins on person_id, picks N per ancestry deterministically (sorted, first-N — stable across re-derives). No notebook/BigQuery step. Prints per-ancestry availability (MID/SAS expected thin in the LR sub-cohort). Refuses to overwrite an existing cohort.tsv without `--force`.
- `run_experiment_d.sh` — no `set -e`; every person written to disk immediately; resumable at per-step granularity (SR fastq / LR fastq / SpecHLA / SpecImmune each skipped if already present-and-complete). "Complete" = per-gene completeness ≥ threshold (SpecHLA ≥ 8/16, SpecImmune ≥ 4/8), never exit code or file existence alone (quirks #17, #18). `ensure_mount()` remounts + ls-verifies if the FUSE layer dropped. `expd.done` marker skips finished people; `MAX_ATTEMPTS=4` per person backstops against infinite retry of a deterministically-broken person. All samtools/SpecHLA/SpecImmune calls are `pixi run -e <env>`-pinned so the run is independent of which shell launched it. Optional `ntfy.sh` pings per person + on completion.
- `analyze_experiment_d.py` — reads the accumulated per-person `comparison_log.csv` rows (run_label `experiment_d`), dedups (person,gene) to the latest timestamp, joins ancestry from cohort.tsv. Concordance is 2-field unordered genotype match with a `norm2()` that maps every method's format (`A*01:01:01:01`, `01:01`, `DRB1*15:01`) to a comparable 2-field key — verified that SpecImmune's deeper long-read resolution still matches short-read at 2-field (so it isn't penalized as a disagreement, per the DECISIONS.md "don't dismiss over-resolution" caution) and that uncallable slots drop out of denominators rather than counting as disagreements. NOT a naive vote (AoU+SpecHLA share short-read data). Emits, per locus × ancestry: full 3-way agreement, three pairwise agreements, and the two targeted patterns — AoU-native lone-outlier rate (the DQA1 question) and SpecImmune divergence-from-agreeing-SR-pair rate, raw + confidence-filtered (not-tied & identity ≥ 0.99, the DPA1 question). Markdown + 2 PNG heatmaps + counts CSV under `~/pipeline_outputs/experiment_d/analysis/`, aggregate-only (keep on VM per egress caveat).

Per-person budget ~25 min → ~56 people ≈ 23-33h wall. Launch sequence in STATUS.md.

## 2026-07-12 (cont.) — Experiment D complete: n=60 results, resolution-tier reframing of DRB1/DQA1, new open padding hypothesis

Full cohort finished cleanly: 60/60 done, 0 failed. Total compute: 20.6 hours across all 60 people.

**Runtime.** SpecHLA: mean 4.7min/person, tight across every ancestry group (4.2-5.3min) — very stable. SpecImmune: mean 15.9min/person, much wider spread (11.2-34.9min) — **AMR stands out at mean 21.7min, ~40-50% higher than every other group** (13.8-15.8min elsewhere), with 5 of 10 AMR people over 20 minutes. Plausible mechanism: AMR (admixed American) ancestry likely carries more divergent/heterozygous HLA alleles, forcing more work in SpecImmune's read-binning/reconstruction — not yet confirmed, just the working explanation.

**3-way concordance headline (`analyze_experiment_d.py`).** DRB1 is the worst locus in every single ancestry group (0-67% full 3-way agreement). DQA1 shows a real, large ancestry gradient in AoU-native's lone-outlier rate (20% in AFR/SAS up to 90% in EAS, 70% in EUR) — answering the carried-over question from the n=3 smoke test: yes, ancestry-correlated, and in the EAS/EUR direction (argues against simple reference-bias as the mechanism). DPA1 (the other carried-over question) comes back mostly null — SpecImmune's divergence from the agreeing SR pair stays low (0-20%) almost everywhere, with a mild AFR uptick (30% raw / 10% confident) not strong enough to call a real ancestry effect.

**Resolution-tier analysis (new script, `analyze_experiment_d_resolution.py`) — reframes both headline findings.** Treats SpecImmune-LR as ground truth and classifies every AoU-vs-truth and SpecHLA-vs-truth genotype call into exact / close (right allele family, wrong sub-type) / partial-miss (one allele from a different family) / full-miss (both alleles from a different family), instead of a binary match. Pooled across all 60 people (n=467 callable gene-rows):

- **DRB1 flips the earlier framing.** AoU-native's DRB1 "errors" are 62% exact, 34% close, only 4% partial-miss, 0% full-miss — almost entirely harmless near-misses. SpecHLA's DRB1 errors are 26% exact, 19% close, **49% partial-miss, 6% full-miss** — mostly real, wrong-family errors. DRB1's poor 3-way agreement is being driven overwhelmingly by SpecHLA's short-read calling specifically, not by the locus being symmetrically hard for everyone.
- **DQA1 flips the other way, reassuringly.** AoU-native's DQA1 discordance (the ancestry-correlated one above) is 42% exact, **58% close, 0% partial-miss, 0% full-miss** — every single "wrong" AoU call in this cohort is still in the right allele family. SpecHLA gets DQA1 right 90% of the time. So AoU's DQA1 problem is bounded imprecision, not unreliability, and SpecHLA is simply the better source for that locus.
- **DPB1 has zero "close" calls for either method** — when either caller misses at DPB1, it's always a jump to a completely different allele family, never a near-miss. A distinct failure mode from every other locus.
- **Overall:** SpecHLA has a slightly higher raw exact-match rate (75% vs AoU's 70%), earned almost entirely at the two DQ genes, but produces 5x more full-misses than AoU-native (5 vs 1), concentrated at DRB1. There is no single "better" method — it's gene-dependent, and now there's per-locus evidence for which method to trust where, rather than defaulting to one.

**This run also doubles as the "≥1 more person" confirmation Experiment B's SpecHLA-pad10000 recommendation was waiting on (n=60, not n=1) — but with an important caveat, not a clean pass.** Experiment B's original padding sweep (2026-07-10, person 2522883, n=1) found "zero call degradation" for SpecHLA down to 325bp padding, using a simple exact-2-field-match methodology that cannot distinguish a near-miss from a real error — the same blind spot the resolution-tier script above was built to fix. It also never tested DRB1's behavior across a padding range with SpecHLA specifically (the LR padding sweep found DRB1 degradation, but that was SpecImmune, a different tool/read-type, per DECISIONS.md's explicit warning not to conflate the two).

**New open hypothesis (Marc, 2026-07-12) — not yet tested, flagged for follow-up, not adopted:** SpecHLA's DRB1-specific breakdown at pad10000 (49% partial-miss + 6% full-miss) could be a padding artifact rather than a pure calling-algorithm weakness. pad10000 is a substantially narrower SR window than any previously-run production default, and DRB1 is exactly the locus with the closest short-read paralogs (DRB3/4/5) — the same mechanism already confirmed, independently, in Experiment C (2026-07-10 entries above): removing SpecImmune's decoy pseudogene/paralog references via gene-panel restriction caused reproducible miscalls at exactly Gene A and DRB1, by redirecting ambiguous paralog-crossing reads onto the wrong allele in the remaining panel. If narrowing SpecHLA's padding similarly discards flanking sequence that previously absorbed DRB1-paralog-crossing reads, that would produce exactly the "same family sometimes, wrong family other times" pattern seen here — but this has never actually been tested at the resolution-tier level. **Concrete next check, not yet done:** re-run SpecHLA at a wider padding level (e.g. fullpad or pad100000) for a handful of people already run at pad10000, resolution-tier-compare each against the same SpecImmune truth, and see whether DRB1's partial/full-miss rate drops as padding widens. If it does, pad10000 is not a safe default for DRB1 specifically, regardless of Experiment B's original (coarser) conclusion — this would need a locus-specific padding exception, not a blanket reversion.

Script: `scripts/analyze_experiment_d_resolution.py`. Aggregate counts only (no genotypes) — output pasted to chat by Marc, not committed.

## 2026-07-18/19 — AoU-native callset validation (origin, exhaustiveness, quality)

Supervisor-requested check on the AoU-native HLA callset itself (535,658 people), separate from
Experiment D's cross-tool comparison. Scripts: `experiment_a2_callset_quality.py` (diversity,
external frequency concordance, Hardy-Weinberg proxy) and `experiment_a3_allele_space_coverage.py`
(catalogued-vs-observed allele space, IMGT DB-version lower bound). Local plotting:
`local_plot_allele_ancestry.py`, `local_plot_coverage.py`.

Headline: callset passes every check run (real hyperpolymorphic diversity, frequencies
consistent with 3 published anchors, 99.96% IMGT catalogue match rate). AoU's own QC report
(GIAB-based) does not benchmark HLA specifically — this pilot is the only HLA-locus reliability
evidence found for this dataset. DRB1 is independently flagged as the least reliable locus by a
second method (Hardy-Weinberg excess homozygosity, survives ancestry stratification) that agrees
with Experiment D's cross-tool finding.

One real bug caught and fixed mid-analysis: initial ancestry-coloring logic ranked alleles by raw
carrier count, which is dominated by cohort composition (56.5% EUR) rather than biology — fixed
to rank by within-group frequency instead (see `local_plot_allele_ancestry.py`).

Full writeup, tables, and figures: **`reports/aou_callset_validation/README.md`** (delivery-ready,
sent to supervisors). Not duplicated here — this entry is a pointer, not a summary of the summary.

## 2026-07-12 — repo reorganized into folders (housekeeping, not an experiment)

Not a result — a structural note for anyone reading older entries above, which still reference bare filenames like `run_experiment_d.sh` or `compare_hla_results.py` at what used to be the repo root. As of this date: all context docs (this file, ENVIRONMENT/STATUS/DECISIONS/TASK_CONTEXT.md, SMOKE_TEST_PICKS.local.md) live under `context/`; all pipeline `.py`/`.sh` files live under `scripts/`; `reference/` and `pixi.toml` are unchanged. See the root `README.md` for the full layout and the context-update system it documents. Every script's internal cross-references were updated and syntax-checked (`bash -n`, `py_compile`) before this was committed. Done on the Mac only — deliberately **not** synced to the VM until the then-running Experiment D job (60-person cohort) finished, since overwriting `compare_hla_results.py`/`spechla_pad_helpers.py`'s location mid-run would have broken every subsequent person's processing.

---

> **Backfill note (2026-07-22 context-cleanup pass):** the four entries below were added to restore this log's top-to-bottom continuity. Experiment E, the field cascade, the LR census, and Experiment F had all been committed to git and (for the first three) written up as reports, but never got an EXPERIMENTS.md pointer — the append-discipline slipped across fast-moving multi-session work. The E / cascade / census entries are sourced from their committed report folders; the F entry is a stub (its results aren't yet recorded anywhere in `context/`).

## 2026-07-20 — Experiment E ("Operation DQ"): HLA×disease sanity check (celiac + narcolepsy)

First-principles end-to-end validation of the whole pipeline (HLA calling → BigQuery diagnosis join → phenotype label → model) via two known large-effect HLA-disease associations, run on the AoU-native calls (n=535,658). L1-penalized logistic regression, 5-fold out-of-fold CV, evaluated pooled and per-ancestry. Headline: both diseases recover the known biology with zero disease-specific hints; **celiac** shows real, modest out-of-sample signal (pooled AUROC 0.74) with a well-powered ancestry gap (EUR 0.681 vs AFR 0.542); **narcolepsy** predicts at chance (pooled AUROC 0.592) despite the model landing on a true causal-haplotype partner — a concrete demonstration that recovering the right biology and predicting well on unseen people are different questions. Full writeup + raw CV log: `reports/disease_sanity_check/README.md`. (Committed in git as "Experiment E: Operation DQ"; the LR-census script was later renamed to resolve a naming collision with it — commit `606a563`.) Open follow-on logged in DECISIONS.md (Methodological): the 2-field vs 3-field vs hierarchical resolution/generalization study.

## 2026-07-20 — Experiment D field-level cascade (follow-on analysis, no new pipeline runs)

Re-scoring of the existing Experiment D n=60 calls at cumulative field-by-field resolution (field 1 allele-group → field 2 protein → field 3 synonymous), AoU-native and SpecHLA each measured against SpecImmune-LR as truth, per gene. Quantifies where each method's agreement decays with depth (e.g. DQA1 100% at field 1 → 64% at field 2; AoU-native structurally caps at 3 fields, never reaching field 4 in anyone). Pure re-analysis of existing Experiment D output — no new runs. Report + heatmaps: `reports/experiment_d_field_cascade/` (git commit `cbe2fc1`).

## 2026-07-21 — Long-read data-format census: all 14,521 people have a real aligned BAM

Exhaustive per-column existence census of the full v9 lrWGS manifest (15,424 rows / 14,521 unique people, 29 primary data-bearing columns), superseding the earlier 2,763-person "floor" (see ENVIRONMENT.md quirk #13). Every person has a real `grch38_bam`; AoU backfilled the previously-missing revio BAMs between 07-11 and 07-21, confirmed by spot-checking quirk #13's own documented dangling case (person 1008366). The long-read validation cohort is no longer capped near 2,763 — all 14,521 are viable for the SpecImmune-LR pipeline. Script `lr_manifest_format_census.py` (checkpointed/resumable after a real 5h data-loss incident — quirk #22). Full per-platform data-shape table + caveats: `reports/lr_data_census/README.md`. Resolves DECISIONS.md's "where are the aligned BAMs" open question.

## 2026-07-22 — Experiment F: Immuannot (assembly-based caller) vs SpecImmune-LR, 60-person pilot

Tested the "assembly-based HLA typing" lead from DECISIONS.md: ran Immuannot (minimap2 vs
IPD-IMGT/HLA, per-haplotype assembly contigs) on the same Experiment D cohort. 57/60 completed,
~2x faster than SpecImmune-LR (median 6.9 vs 13.9 min/person, both measured on this VM). Field 2
(protein) concordance 65–94% for most genes, but only 40% at DRB1 — now the hardest locus by five
independent signals (this comparison, both callers' own confidence flags, plus the earlier
field-cascade/callset-validation findings). No ground truth here (cross-method agreement only);
Aleix's HLA-Resolve/Lai arm remains the truth-anchored comparison. A padding-tightening experiment
(`run_immuannot_pad_sweep.py`) found the production trim's generous whole-block approach is
necessary, not wasteful — narrower cuts lose real calls, especially at DQA1/DQB1 — not adopted.
Parallelizing the two haplotypes (`--parallel-haps`) gave no net speedup (halved per-hap threads
roughly cancels the concurrency gain) — not adopted either. Full writeup + figures:
`reports/immuannot_pilot/README.md`.

## 2026-07-24 — Experiment F follow-on: Immuannot-as-truth field cascade + confidence-matched truth comparison

Two follow-on analyses on Experiment F's cohort/calls, no new pipeline runs. (1) Re-ran the
field-cascade re-scoring (originally SpecImmune-truth only) with Immuannot substituted as
truth (`analyze_experiment_d_field_cascade_immuannot.py`) — DRB1's Field-2 concordance looks
better under Immuannot-truth than SpecImmune-truth (AoU 89% vs 77%), plausibly a selection
artifact (Immuannot's own confidence is worst at DRB1 of any gene) rather than real signal;
DQA1's AoU-native weakness is now cross-validated under both truth sources. (2) Built a single,
mathematically-grounded (not per-gene-tuned) confidence filter for both truth sources —
`si_identity >= 0.9995`/not-tied, `template_distance = 0`/no-warning — sample-size-matched by
count against the shared 480-slot denominator so neither truth source gets an inflated share
(`analyze_confidence_matched_truth.py`, `plot_confidence_distributions.py`). DRB1 is nearly
eliminated by the filter under both truth sources; DQA1's finding survives. Full writeup:
`reports/confidence_matched_truth/README.md`, `reports/experiment_d_field_cascade/README.md`.
Threshold derivation (sequencing-error-rate math, IMGT inter-allele spacing): `context/DECISIONS.md`.
