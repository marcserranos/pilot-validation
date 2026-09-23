# 41 — KIR scoping + pilot (WS5, Sprint S03 call #8 §7)

Cole (call #8, 2026-09-22): the production Immuannot run only ever trimmed assemblies to the
chr6 HLA window — **KIR genes (chr19) were never called**, even though they interact closely
with HLA (NK-cell receptors) and are expected to show more germline diversity and novel alleles
than HLA itself. Marc's action item: scope the region and rerun cost. This report is that scoping
pass plus a 20-person pilot, ahead of a full-cohort go/no-go decision.

Code: `scripts/hla_popgen/41_kir_pilot.py` (+ unit tests in
`scripts/hla_popgen/tests/test_kir_pilot.py`). No changes to `scripts/run_immuannot_person.py` —
its existing `--region`/`--pad` flags already generalize past chr6 HLA, so KIR only needed a
different `--region` value and a person-selection/aggregation wrapper.

## 1. Coordinates (GRCh38)

Immuannot's own reference data (IPD-KIR) and detection logic already support all 17 KIR genes —
see §2 — so the only thing missing was trimming assemblies to the right window.

Primary-source coordinates (GeneCards/Ensembl, GRCh38), not re-derived from a browser screenshot:

| Gene | Role | GRCh38 coordinates |
|---|---|---|
| LILRB1 | centromeric anchor, immediately outside the KIR array | chr19:54,616,309–54,638,022 |
| **KIR3DL3** | **centromeric-most framework gene** (array start) | chr19:54,724,442–54,736,632 |
| **KIR3DL2** | **telomeric-most framework gene** (array end) | chr19:54,850,443–54,867,207 |
| FCAR | telomeric anchor, immediately outside the KIR array | chr19:54,874,235–54,891,420 |

The core KIR tandem array (KIR3DL3 → KIR3DL2) spans **~143 kb**, matching the literature's
"~150 kb" figure for the gene cluster (NCBI Bookshelf, *The KIR Gene Cluster*) — and matching
Cole's own estimate on the call ("region may be comparable size to HLA (~150kb)").

**Trim window used: `chr19:54,600,000-54,920,000`** (320 kb) — `KIR_REGION_DEFAULT` in
`41_kir_pilot.py`. This is deliberately wider than the textbook array:
- ~110 kb margin on the centromeric side (past LILRB1's start), ~30 kb margin on the telomeric
  side (past FCAR's end).
- Justification for the extra margin, beyond HLA's precedent: the KIR region is structurally the
  most complex locus scoped so far in this project — large insertion/deletion CNV haplotypes are
  well documented (gene content varies haplotype-to-haplotype, unlike HLA's fixed gene set), so a
  tight boundary risks clipping a real gene near a per-person structural-variant edge. The anchors
  (LILRB1, FCAR) are themselves outside the variable region, so trimming past them is safe.
- This 320 kb is only the **outer candidate-contig-selection boundary** passed to
  `run_immuannot_person.py --region`. Its existing Tier-1 `.paf`-based trim
  (`regions_from_paf()`) then adds its own further **±100 kb pad** (`--pad`, unchanged default,
  already proven safe for the analogous HLA case — pad100k, zero degradation across 8 HLA genes,
  EXPERIMENTS.md) around whatever actually aligns inside that window — so the real structural
  margin per person is wider than 320 kb once their own contig alignment is factored in.

## 2. Immuannot KIR support

Read `~/tools/Immuannot`'s scripts and `reference/README_Immuannot.md` directly (not assumed from
the paper alone). Findings:
- Immuannot's refdata bundles **IPD-KIR** alongside IPD-IMGT/HLA and RefSeq (C4) — its own
  "Gene coverage" section lists all **17 KIR genes**: KIR2DL1, 2DL2, 2DL3, 2DL4, 2DL5A, 2DL5B,
  2DP1, 2DS1, 2DS2, 2DS3, 2DS4, 2DS5, 3DL1, 3DL2, 3DL3, 3DP1, 3DS1 — matching SCHEMA.md's `kir`
  `gene_class` row exactly.
- **KIR is annotated automatically once a KIR-region contig is provided — there is no separate
  gene-family flag to set.** Immuannot's detection strategy is alignment-driven: it maps every
  reference allele (HLA + KIR + C4) against the input contig(s) and reports whatever overlaps
  with >90% coverage. The pipeline's own shipped example (`example/test.fa.gz`) literally
  contains two contigs, "one for HLA region and the other for KIR region," run through the exact
  same `immuannot.sh` invocation — confirming region selection is entirely the caller's
  responsibility (i.e., ours, via the trim step), not a pipeline mode switch.
- The refdata's 1,530 KIR genomic alleles (recon'd from `~/tools/Immuannot_refdata` on the VM)
  is the allele panel novelty is measured against below.
- **Known bug, not fixed here (out of backward-compat scope):** `run_immuannot_person.py`'s
  Tier-3 self-align fallback (`ensure_chr6_ref()`) hardcodes `samtools faidx <ref> chr6` — not
  parameterized by the target region's chromosome. `41_kir_pilot.py` works around this without
  touching that file: `ensure_chr6_ref()` is a no-op once its cache path already exists and is
  non-empty, so `41_kir_pilot.py`'s own `ensure_kir_ref_cache()` pre-extracts **chr19** into a
  separate cache file and passes that path via the (misleadingly named but not chr6-specific)
  `--chr6-ref-cache` flag. Not exercised in this pilot (see §4 — the pilot's candidate pool had
  zero sequel2/Tier-3-only people), but the workaround is ready and unit-relevant if a full run
  needs it.

## 3. Pilot: 20 people, mixed ancestry + platform

Selection (`pick_pilot_cohort()`): drawn from the same unrelated-cohort definition as
`24_novelty_by_field.py`/`29_hla_ld_by_ancestry.py` (KIN_MIN 0.0442, third-degree-or-closer per
KING — logic copied, not re-derived), restricted to assembly-holding platforms
(revio/sequel2e/sequel2), round-robin stratified across (ancestry, platform) buckets. Person IDs
never left the VM (VM_CHANNEL.md hard rule) — only this aggregate mix and the summary numbers
below did.

**Selected mix (aggregate-only):** 2 people each for AFR/AMR/EAS × {revio, sequel2e} except
EUR/MID/SAS×sequel2e (1 each) and one UNASSIGNED/revio — 6 ancestry groups (AFR, AMR, EAS, EUR,
MID, SAS) plus 1 person with no confident ancestry call, split revio (13) / sequel2e (7). The
candidate pool (`cohort_membership.tsv`, the same table the production HLA run's people are drawn
from) had **zero sequel2 people** — that platform's Tier-3-only cohort was apparently excluded
upstream (consistent with sequel2's self-align tier being unvalidated at HLA-calling time too) —
so this pilot exercises **Tier 1/2 only** (`.paf` sub-range trim / whole-contig `.bam` fallback),
not Tier 3. Tier 3 remains untested for KIR, same caveat as it already carries for HLA.

**Result: 20/20 people (100%), 40/40 haplotypes (100%) produced at least one KIR gene call.**
All 40 haplotypes used the tight Tier-1 `.paf`-based sub-range trim (`trim_method=paf_region`,
zero whole-contig or self-align fallbacks needed) — mean trimmed contig size **1.02 MB/haplotype**
(vs. HLA's historically larger trimmed input), which is the main reason KIR runs faster than HLA
per person (see §4).

| Metric | Value |
|---|---|
| People run | 20/20 succeeded (100%) |
| Haplotypes with ≥1 KIR call | 40/40 (100%) |
| Mean KIR genes called / haplotype | 9.32 (of 17 possible) |
| Mean wall time / person (2 threads, this pilot's cap) | 192.3 s (~3.2 min) |
| Mean CPU time / person | 333.5 s (~1.7× the wall time — real ~2-thread parallelism) |
| Total pilot wall time (20 people, sequential) | 64.1 min |
| Total KIR allele calls | 383 |
| **Novel alleles (no exact IPD-KIR match)** | **60.8%** (see caveat below) |
| Undetermined calls | 0.0% |

**Framework genes** (KIR3DL3, KIR2DL4, KIR3DP1, KIR3DL2 — expected on ~all haplotypes):
KIR3DL3 97.5%, KIR3DL2 90.0%, KIR3DP1 90.0%, KIR2DL4 85.0% of the 40 successful haplotypes. Not
100% — flagged, not waved away: since every haplotype here used the tight Tier-1 trim (not the
looser whole-contig fallback), a trim-boundary clipping explanation is unlikely to be the whole
story; the ~3–15% framework-gene misses are more likely either genuine rare structural haplotypes
or an Immuannot-detection-sensitivity edge at KIR's higher sequence divergence. Worth a targeted
look (which specific haplotypes, do they cluster by ancestry or novelty rate) before a full run,
not a blocker to piloting further.

**Gene-content haplotype split:** 23/40 (57.5%) classified cA (fixed, inhibitory-dominated
content), 17/40 (42.5%) cB (carries ≥1 activating/variable-content gene) — a plausible split for
a mixed-ancestry cohort, well inside the range reported across human populations.

**Sanity checks (both passed):**
- **KIR2DL2 / KIR2DL3 mutual exclusivity:** 0.0% of haplotypes carrying either gene carry both —
  clean allelic-alternative behavior, as expected (these occupy the same locus).
- **KIR3DL1 / KIR3DS1 mutual exclusivity:** 2.8% co-occurrence (1 haplotype of 36 carrying
  either) — near-exclusive as expected, with a small real/measurement tail consistent with the
  literature (these are allelic alternatives but not as strictly exclusive as 2DL2/2DL3 in every
  published cohort).

**Novelty caveat, load-bearing for the 60.8% number:** an earlier version of this script's
novelty detector checked for a literal `":new"` suffix on the `consensus` field (matching
`IMMUANNOT_GTF_SPEC.md`'s HLA-style illustrative example, `"HLA-A*01:01:new"`) and silently
measured 0.0% novelty on all 383 real calls. Spot-checking one real `hap1.gtf.gz`'s raw consensus
strings found the true convention for KIR is different: HLA alleles are colon-field-delimited so
`":new"` fits, but **KIR alleles are not colon-delimited** (e.g. `KIR3DL2*00201`), so Immuannot
appends `new` directly with no colon (`KIR3DL2*00201new`). Fixed in `parse_hap_gtf()` (checks
`endswith("new")`, covers both conventions) and re-verified with a unit test
(`test_parse_hap_gtf`'s `test_novel_and_known_and_undetermined`) before trusting the number above.
**60.8% novel-allele fraction is high, but directionally exactly what Cole predicted** (KIR
"expected to show more germline diversity and novel alleles than HLA itself") — plausible given
IPD-KIR's reference panel (1,530 genomic alleles) is far smaller/less complete than IPD-IMGT/HLA's,
so a real novelty rate this much higher than HLA's is not surprising, but it has not been
cross-checked against a second truth source the way HLA's novelty numbers were (script 03's
recurrence-across-≥2-people QC) — that cross-check belongs in the full run's own QC pass, not this
scoping pilot.

Selection mix (aggregate-only, matches §3's summary): 2 people each for AFR/AMR/EAS ×
{revio, sequel2e}, 1 each for EUR/MID/SAS × sequel2e (2 each × revio), 1 UNASSIGNED/revio.

### 3a. Novelty decomposition (2026-09-23 follow-up)

60.8% "novel" does not mean 60.8% protein-changing. Decomposed each of the 383 calls into 4 tiers,
mirroring `24_novelty_by_field.py`'s known/f4_noncoding/f3_synonymous/f2_protein framework but
derived directly from Immuannot's own already-computed CDS-vs-refdata comparison
(`cds_distance`/`cds_mut` GTF fields) rather than rebuilding that script's `RefIndex`/
`classify_sequence()` machinery from raw contig sequence — `cds_distance` is itself the edit
distance of a CDS-vs-CDS alignment against **IPD-KIR's own `CDSseq/*.fa.gz` reference** (confirmed
present for all 17 KIR genes under `~/tools/Immuannot_refdata/CDSseq/`, the same directory
`RefIndex` globs), so this reuses the same reference data and search result, just without redoing
the alignment ourselves — a documented tradeoff (see `classify_novelty_tier()`'s docstring), not a
shortcut that changes what's being measured.

| Gene | n calls | known | novel, genomic-only (intron/UTR) | novel CDS, synonymous | **novel protein** |
|---|---|---|---|---|---|
| KIR2DL1 | 37 | 48.6% | 45.9% | 0.0% | 5.4% |
| KIR2DL2 | 12 | 0.0% | 91.7% | 8.3% | 0.0% |
| KIR2DL3 | 26 | 0.0% | 92.3% | 0.0% | 7.7% |
| KIR2DL4 (fw) | 36 | 58.3% | 38.9% | 0.0% | 2.8% |
| KIR2DL5A | 8 | 87.5% | 12.5% | 0.0% | 0.0% |
| KIR2DL5B | 11 | 36.4% | 45.5% | 0.0% | 18.2% |
| KIR2DP1 | 36 | 44.4% | 47.2% | 2.8% | 5.6% |
| KIR2DS1 | 9 | 55.6% | 44.4% | 0.0% | 0.0% |
| KIR2DS2 | 9 | 33.3% | 66.7% | 0.0% | 0.0% |
| KIR2DS3 | 9 | 33.3% | 66.7% | 0.0% | 0.0% |
| KIR2DS4 | 28 | 21.4% | 75.0% | 0.0% | 3.6% |
| KIR2DS5 | 9 | 33.3% | 66.7% | 0.0% | 0.0% |
| KIR3DL1 | 31 | 35.5% | 61.3% | 0.0% | 3.2% |
| KIR3DL2 (fw) | 37 | 18.9% | 70.3% | 0.0% | 10.8% |
| KIR3DL3 (fw) | 39 | 30.8% | 64.1% | 0.0% | 5.1% |
| KIR3DP1 (fw) | 38 | 73.7% | 23.7% | 0.0% | 2.6% |
| KIR3DS1 | 8 | 75.0% | 25.0% | 0.0% | 0.0% |

**Confirmed: novelty is dominated by non-coding/genomic-only differences**, exactly as
hypothesized — "novel, genomic-only" is the largest novel tier for every gene except KIR2DL5A,
protein-level novelty tops out at 18.2% (KIR2DL5B, n=11 — small-n, treat as noisy) and is ≤11% for
every other gene, usually single digits. "novel CDS, synonymous" is nearly always 0% — when the
CDS itself differs from every known CDS, it's almost always a real amino-acid change (protein
tier), not a silent one; genuine synonymous-only CDS novelty appears only in KIR2DL2 (8.3%) and
KIR2DP1 (2.8%). This reframes the headline finding: KIR alleles genuinely diverge from IPD-KIR far
more than HLA diverges from IPD-IMGT/HLA (consistent with IPD-KIR's much smaller reference panel,
1,530 genomic alleles), but that divergence is overwhelmingly in non-coding sequence, not protein
sequence — the same "bounded imprecision, not unreliable calling" pattern this project already
established for AoU-native's DQA1 calls (`context/DECISIONS.md`).

### 3b. Framework-gene miss classification (2026-09-23 follow-up)

All 15 framework-gene misses (across 40 haplotypes × 4 framework genes) classified by whether
genes flanking the missing one (centromeric AND/OR telomeric, per the canonical
KIR3DL3→KIR2DL4→KIR3DP1→KIR3DL2 array order) were also called on that same haplotype — see
`classify_framework_miss()`'s docstring for the exact rule and its limitations.

| Classification | Count | Meaning |
|---|---|---|
| **flanked** | 11 | Genes on **both** sides of the missing one are called — the contig clearly spans across its expected position. A real deletion or a miscall-as-a-different-gene is more plausible than an incomplete assembly. |
| **ambiguous_partial_flank** | 4 | Only one side has a flanking framework gene called — doesn't cleanly fit either explanation. |
| **edge_fragmented** | 0 | (Missing at an array end with nothing beyond it — none observed.) |

**Zero of the 15 misses look like simple assembly fragmentation.** 11/15 are clearly "flanked" —
the contig has real sequence on both sides of the gap, so contig-too-short is not the explanation;
these are either genuine deletions (KIR3DP1 and KIR2DL4 in particular have documented real-deletion
haplotypes in the literature, per the coordinator's own note) or the gene was called under a
different name (a miscall this classification cannot distinguish from a real deletion — would need
per-haplotype inspection of what other gene sits at the expected coordinates, out of scope for an
aggregate-only pilot pass). This strengthens, not weakens, the go/no-go case: the framework-gene
misses are not a trim/assembly artifact.

## 4. Scale-up estimate

The original full-cohort **HLA** production run (`run_production_orchestrator.py`,
`context/DECISIONS.md` 2026-08-04 entries) is the closest precedent:
- Ran on a **dedicated `n2-highcpu-96` Workbench VM instance** (96 vCPU/96GB, not the current
  4-vCPU `big_run` VM this pilot used), launched directly on that interactive Workbench VM —
  **not** dsub/Cloud Batch (considered as a future cost lever, never adopted; Cloud Life Sciences,
  the older managed-genomics-batch API, is fully shut down as of July 2025).
- Config: 24 people concurrently × 4 threads/person (=96 cores), `setsid nohup` + heartbeat to an
  external monitor, automatic two-phase run (main cohort, then the self-align-only tier).
- Measured/extrapolated: ~249 people/hour → **13,000–14,521 people in ~52–58 hours wall-clock**,
  **~$160–200** (N2 ~$0.0316/vCPU-hr × 96 vCPU ≈ $3.03/hr), 2TB disk, unpruned (~81 MB/person).

KIR-specific adjustments now that the pilot's real timing has landed:
- Pilot mean trimmed contig: **1.02 MB/haplotype** — smaller than HLA's typical trimmed input,
  and pilot per-person wall time (192.3 s at only 2 threads, this pilot's compute cap) is well
  under HLA's ~6.9 min/person median (Experiment F, 60-person Immuannot pilot, 4 threads) even
  before accounting for the extra 2 threads a full run would give each person.
- **Throughput extrapolation** (same "throughput scales ~linearly with total core count for this
  per-person-embarrassingly-parallel workload" finding the original HLA scale-up decision already
  relied on — `scaling_probe.py`, confirmed flat 32×1/16×2/8×4 throughput at a fixed 32-core
  budget, `context/DECISIONS.md` 2026-08-04): this pilot ran 20 people **sequentially on 2 cores**
  in 64.1 min wall-clock = **18.7 people/hour at 2 cores**. Scaling that by 48× (to the HLA run's
  96-core `n2-highcpu-96` machine) → **~897 people/hour → 13,000–14,521 people in ~14.5–16.2
  hours**, well under HLA's own 52–58 hours on the same machine — consistent with KIR's much
  smaller trimmed-contig size. **This is a bigger extrapolation factor (48×) than the original
  HLA estimate ever validated (that one only confirmed flat throughput up to 32 cores) — treat
  the 14.5–16.2 h figure as optimistic, not confirmed, until a larger pilot (e.g. 4 or 8 threads,
  concurrency >1) checks the scaling holds.**
- Estimated cost at that runtime: ~15 h × $3.03/hr (96-vCPU N2 rate) ≈ **$45** — well under the
  $300 cap, with more margin than the HLA run had.
- KIR runs on the **same already-produced haplotype contigs** — no new assembly/alignment step,
  just a second trim + Immuannot pass per person. A full KIR run is a genuinely separate
  orchestrator invocation (different `--region`, different output paths, `--out-suffix .kir` so
  it never collides with the canonical `immuannot_calls.tsv`), but reuses every other piece of
  infrastructure the HLA run already proved (mount check, PID lock, resumability, heartbeat).

**`run_production_orchestrator.py` now takes `--region`/`--pad`/`--out-suffix`** (added 2026-09-23,
backward compatible — defaults reproduce the prior HLA-only behavior byte-for-byte; unit-tested in
`scripts/production_orchestrator/tests/test_orchestrator_region_flags.py`, 6 tests). `--out-suffix`
isolates BOTH the per-worker fragment files AND the canonical `immuannot_calls<suffix>.tsv`/
`immuannot_timing<suffix>.tsv` merge target — the orchestrator now **refuses to start** (loud
`FATAL`, before touching the mount) if `--region` differs from the HLA default without a
non-empty `--out-suffix`, specifically to prevent a KIR run from silently merging into the real
HLA canonical files (the class of bug ENVIRONMENT.md quirk #29 already cost days to recover from
once). Writing that test also caught a real bug in the first cut of this change: a plain glob
pattern for the default (empty) suffix would ALSO match a `.kir.`-suffixed run's fragments (glob's
`*` matches across dots) — fixed with a regex-anchored fragment filter before this ever ran for
real.

**Exact command for the full run** (once pilot quality clears go/no-go — same VM sizing as the
HLA run, since nothing about KIR's I/O/CPU profile is expected to differ enough to justify a
different machine type until a larger pilot says otherwise):

```bash
# On a dedicated n2-highcpu-96 Workbench VM instance (not the shared 4-vCPU big_run VM):
cd ~/repos/pilot-validation
gcsfuse --billing-project wb-glacial-potato-8710 --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled
setsid nohup pixi run -e specimmune -- python3 scripts/production_orchestrator/run_production_orchestrator.py \
    --cohort ~/pipeline_outputs/immuannot_cohort_full.tsv \
    --concurrency 24 --threads-per-person 4 \
    --region chr19:54600000-54920000 --pad 100000 --out-suffix .kir \
    --vm-rate <live Workbench-UI-confirmed USD/hour for n2-highcpu-96> \
    --monitor-url http://46.225.123.54:8943 \
    > ~/pipeline_outputs/kir_full_run.log 2>&1 < /dev/null &
disown
```

**Still a real risk even with `--out-suffix`:** `run_immuannot_person.py` writes each person's
per-person intermediate output (`<outroot>/<pid>/immuannot_output/hap{1,2}.gtf.gz`) to a path keyed
only by `person_id`, not by region/suffix. Running this against the **same `--outroot`** as the
HLA run would overwrite each person's HLA `hap*.gtf.gz` with their KIR one (or vice versa) and
would also make the orchestrator's resumability scan (`scan_already_done()`, which checks for
literal `hap1.gtf.gz`/`hap2.gtf.gz` existence) think KIR people are "already done" from the HLA
run and skip them entirely, silently producing an empty KIR run. **Use a separate `--outroot`
(e.g. `~/pipeline_outputs_kir`) for the full KIR run** — not yet enforced in code (would need
`run_immuannot_person.py`'s own intermediate-path convention changed, out of this backward-compat-
only change's scope), so this is a documented operational rule, not a guardrail.

**Disk:** reuses the HLA run's 2TB disk — KIR outputs are per-person `hap{1,2}.gtf.gz` +
trimmed FASTA, same order of magnitude as HLA's ~81 MB/person, so no separate disk provisioning
needed if run as a second pass on the same VM/disk (recommended), or budget ~1.2 TB extra if run
as a fully separate instance.

## 4a. Extended pilot (launched 2026-09-23, running in background)

Launched a 170-person run (the original 20 + ~150 new people, same `--seed 41` — the selection
function is a deterministic prefix: requesting more people with the same seed reproduces the exact
same first 20 picks, so the 150 new people are genuinely additional, not overlapping, and the
already-done 20 skip instantly via the resumability check). Stratified the same way as the pilot
(ancestry × platform round-robin over AFR/AMR/EAS/EUR/SAS/MID), giving roughly 25+ people per
ancestry group once complete (MID capped by however many exist in the candidate pool). 2 cores,
`setsid nohup`, resumable (kill-and-relaunch-safe), running unattended for hours — this is expected
and also keeps the VM from idling out. Launched, not waited on.

```bash
# Launch (already running as of this report):
cd ~/repos/pilot-validation && setsid nohup pixi run -e specimmune -- python3 ~/s03/41_kir_pilot.py \
    --n-people 170 --seed 41 --threads 2 --outroot ~/s03/results/41/pipeline_outputs \
    --out-suffix .kir41 > ~/s03/results/41/run_extended.log 2>&1 < /dev/null &
disown

# Check progress any time (safe, read-only):
tail -20 ~/s03/results/41/run_extended.log

# Aggregate once complete (re-running the SAME command above is enough -- every person already
# done skips instantly via the hap1/hap2.gtf.gz resumability check, and it recomputes+overwrites
# ~/s03/results/41/pipeline_outputs/41_kir_pilot_summary.kir41.tsv from all 170 people's real
# output in well under a minute):
cd ~/repos/pilot-validation && pixi run -e specimmune -- python3 ~/s03/41_kir_pilot.py \
    --n-people 170 --seed 41 --threads 2 --outroot ~/s03/results/41/pipeline_outputs \
    --out-suffix .kir41 > ~/s03/results/41/run_extended_aggregate.log 2>&1
tail -5 ~/s03/results/41/run_extended_aggregate.log   # prints the final aggregate-only Quality dict
```

## 5. Go/no-go

**GO**, with two conditions before launching the full run:

1. **Add `--region`/`--pad`/`--out-suffix` flags to `run_production_orchestrator.py`** (it
   currently hardcodes the HLA window in two places — LOG.md's recon brief). Small, mechanical,
   mirrors `run_immuannot_person.py`'s own already-existing `--region` flag (backward-compatible
   default unchanged) — not done in this scoping pilot since it touches shared production
   infrastructure another workstream may also be mid-edit on.
2. **Spot-check the ~3–15% framework-gene misses** (§3) on a few more people (or during the full
   run's own QC pass) to rule out a systematic detection issue before trusting per-haplotype
   gene-content calls at scale — all 40 pilot haplotypes used the tight Tier-1 trim, so a
   trim-boundary explanation looks unlikely, but it hasn't been directly ruled out.

Why GO: 100% person and haplotype success rate on a real, ancestry- and platform-mixed 20-person
pilot; both structural sanity checks (KIR2DL2/2DL3 and KIR3DL1/3DS1 mutual exclusivity) passed
cleanly; runtime is fast (~3.2 min/person even under-resourced at 2 cores) and projected full-run
cost (~$45) is a small fraction of the $300 cap even under the more conservative reading of the
scaling extrapolation. The 60.8% novel-allele fraction is a genuinely exciting result for Cole's
stated interest (KIR novelty as a paper angle) but is *not yet cross-validated* against a second
source the way HLA's novelty numbers were — that validation is full-run/QC-pass work, not a
pre-condition for starting the run.

Why not a stronger caveat: Tier 3 (self-align, the untested-at-scale fallback) was not exercised
in this pilot (the candidate cohort had zero sequel2 people — see §3) — this pilot's 100% success
rate does not extend to that platform's people. If Tier 3 is needed for KIR at full-cohort scale
(likely, since sequel2 people exist in the broader manifest even if excluded from this
already-HLA-filtered candidate pool), it should get its own small confirmation pilot first, using
the chr19-reference-cache workaround already built and unit-tested here (`ensure_kir_ref_cache()`).

## Distilled

- KIR window: `chr19:54,600,000-54,920,000` (320 kb, anchored on LILRB1/FCAR outside the ~143 kb
  core array) — no code change to `run_immuannot_person.py`, just a different `--region` value.
- Immuannot calls KIR automatically once given a KIR-region contig; no gene-family flag exists or
  is needed. One known bug (Tier-3 self-align hardcodes chr6) worked around, not fixed, via a
  separate chr19 reference cache passed through the existing (misleadingly-named) flag.
- **Pilot: 20/20 people, 40/40 haplotypes succeeded (100%).** Mean 9.32 KIR genes/haplotype,
  correct KIR2DL2/2DL3 and near-correct KIR3DL1/3DS1 mutual exclusivity, framework genes present
  85–97.5% of haplotypes (flagged for a closer look, not blocking). **60.8% of allele calls are
  novel** (no exact IPD-KIR match) — a real number, but only after fixing a novelty-detection bug
  found by spot-checking raw output (KIR alleles append `new` with no colon, unlike HLA's
  `:new`) — the first version of this pilot silently reported 0% novelty.
- Runtime: 192.3 s/person mean at 2 cores (this pilot's compute cap); 1.02 MB mean trimmed contig
  (all Tier 1, no fallback needed). Extrapolated full-cohort estimate ~14.5–16.2 h / ~$45 on the
  same 96-core machine the HLA run used — faster/cheaper than HLA, but the 48× core-scaling
  extrapolation is a bigger leap than the original HLA estimate ever confirmed.
- **Recommendation: GO**, after (1) adding region/pad/suffix flags to the production orchestrator
  (currently hardcoded to HLA) and (2) a brief look at the framework-gene misses. Tier 3 untested
  here (pilot cohort had 0 sequel2 people) — needs its own small check before a full run relies on
  it.
