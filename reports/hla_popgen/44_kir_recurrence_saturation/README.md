# 44 — KIR recurrence / saturation, and KIR vs HLA catalogue coverage (v4b)

**Status: v4b, 2026-09-27 (commit `9ca7053`, figures re-rendered this pass).** This supersedes
every v1–v3 number ever written for this analysis. See "Version history" below before trusting
any number from an older copy of this README, a Slack message, or an earlier report draft.

## Question

Two questions Marc asked (S04 WS-A / WS-B):

- **(A)** How does allele discovery saturate as we add people, split by how often each distinct
  allele recurs (seen 1x / 2x / >2x / ≥20x), separately for "any sequence difference from the
  reference" (general novelty) and "translated-protein change" (protein novelty), overall and per
  ancestry (equal-N) — with sample coverage and Chao2 richness alongside?
- **(B)** Which reference catalogue — KIR (IPD-KIR 2.13.0) or HLA (IPD-IMGT/HLA 3.55.0) —
  represents its own region better, per gene where possible?

Both questions share one cohort (11,856 people with both a KIR and an HLA call, unrelated,
joint-QC'd) and one identity/novelty framework (below), computed once by
`44_kir_recurrence_saturation.py` and rendered by `45_kir_recurrence_figure.py` (question A) and
`46_kir_vs_hla_catalogue.py` (question B, in `../46_kir_vs_hla_catalogue/`).

## Method

**Identity levels** (four, per haplotype call, both species):

1. **genomic** — a hash of the true gene-span sequence, taken from `hap{1,2}.trimmed.fa` at the
   gene's own GTF coordinates (reverse-complemented on the `-` strand), compared against the
   IPD genomic reference (`gen.fa.gz`; "novel" = no reference record is equal to / contains / is
   contained in the observed span). **This is an UPPER BOUND on novelty, not a clean number** —
   see Caveat 1.
2. **cds** — a hash of the coding sequence only (`hap{1,2}/cds.fa.gz`), after the artifact gate
   below.
3. **protein** — a hash of the translated CDS, after the same artifact gate. `pct_novel_protein`'s
   denominator is `protein` s_obs (not `genomic`) — same granularity as its own numerator.
4. **any_novel** / **protein_novel** — the genomic-level and protein-level views restricted to
   calls whose identity has no match in the corresponding IPD reference set (a novelty flag, not a
   fifth granularity).

**Artifact gate** (CDS/protein tracks only, both species, same rule): a call is excluded from
`cds`/`protein`/`protein_novel` if it is `partial_cds` (doesn't span a full ORF),
`homopolymer_indel` (single-base indel inside a homopolymer run — very likely a long-read
sequencing artifact, not a real allele), or `frameshift_or_stop` (a premature stop or frame-shift
that the reference doesn't have). Genomic-level identity does **not** pass through this gate (it
can't — it's sequence, not an ORF), which is the main reason it's an upper bound. Artifacts removed
by species (from `diagnostics_identity.tsv`, summed over genes):

| species | total calls | clean | partial_cds | homopolymer_indel | frameshift_or_stop |
|---|---|---|---|---|---|
| HLA | 875,179 | 640,305 | 218,838 | 16,036 | ~0 |
| KIR | 222,771 | 151,925 | 64,874 | 5,932 | 709 |

**Pseudogenes**: KIR2DP1 and KIR3DP1 have no catalogued reference protein (every IPD-KIR entry for
them is itself frameshifted/stopped) — their `protein`/`protein_novel`/`completeness_protein` are
the literal string `"NA"` throughout, **never 0** (`kir_protein_catalogue_qc.tsv`). Rendered as an
explicit note in every figure that would otherwise show them, never as a zero-height bar or a
(0, 0) point.

**Recurrence classes**: `eq1` (seen in exactly 1 haplotype cohort-wide), `eq2` (exactly 2),
`gt2` (more than 2, includes `ge20`), `ge20` (≥20 — `gt2` and `ge20` are reported both because
`ge20` is Marc's specific disclosure-adjacent threshold, not because they're disjoint).

**Coverage / completeness** (kept strictly distinct, do not conflate — Caveat 2):
- **Good-Turing incidence coverage** (Chao & Jost 2012): P(the next carrier's allele was already
  seen). ~92–99% for both species at CDS/protein — most PEOPLE aren't surprising, which is not the
  same as most ALLELES being known.
- **Chao2 richness completeness** = S<sub>obs</sub> / Chao2 (1987): the fraction of the
  *estimated total allele richness* already observed. This is the number that answers "how well
  covered is this catalogue" — 20–60%, not 90%+.

**Two-proportion z-test** (pooled variance, two-sided, no scipy dependency,
`45_kir_recurrence_figure.two_proportion_ztest`): tests whether HLA's and KIR's singleton share
(`eq1 / s_obs`) differ, at each identity level. No multiple-testing correction is applied across
the four levels tested (genomic/cds/any_novel/protein_novel) — read each row as its own test.

Pooled-ALL counts used for every headline number and hypothesis test below are recovered from
`saturation_curves.tsv` at `n == n.max()` per gene, summed (`pooled_recurrence_from_curves()`),
**not** from `recurrence_classes.tsv`'s per-gene `<20`-masked cells (masking independently per gene
under-counts a naive sum — see the function's docstring). Verified exact:
`eq1 + eq2 + gt2 == s_obs` for every species × level checked (printed by `45`'s own run log).

## Results

### A — Discovery, saturation, and recurrence class

Pooled-ALL (n=11,856), from `coverage_chao2.tsv` / `saturation_curves.tsv`:

| level | species | S<sub>obs</sub> | Good–Turing coverage | Chao2 | completeness (S<sub>obs</sub>/Chao2) | singleton share (eq1/S<sub>obs</sub>) |
|---|---|---:|---:|---:|---:|---:|
| genomic (upper bound) | HLA | 17,188 | 92.3% | 97,392 | 17.6% | 78.5% |
| genomic (upper bound) | KIR | 38,024 | 73.0% | 197,651 | 19.2% | 79.4% |
| CDS | HLA | 1,363 | 99.6% | 3,515 | 38.8% | 50.4% |
| CDS | KIR | 2,100 | 81.7% | 5,048 | 41.6% | 56.1% |
| protein | HLA | 1,079 | 99.7% | 2,703 | 39.9% | 47.0% |
| protein | KIR | 1,444 | 92.8% | 3,687 | 39.2% | 56.2% |
| any-level novel | HLA | 14,591 | 39.7% | 118,184 | 12.3% | 85.5% |
| any-level novel | KIR | 37,503 | 61.2% | 198,813 | 18.9% | 80.2% |
| protein-level novel | HLA | 198 | 5.8% | 2,490 | 8.0% | 97.0% |
| protein-level novel | KIR | 1,063 | 74.3% | 3,378 | 31.5% | 69.1% |

Reading the curves (`fig_saturation_by_recurrence.png`, ALL ancestry, general vs protein novelty
side by side, both species): the "seen 1x" (singleton) line dominates every panel and never
visibly bends toward saturation at this cohort size — for both species, at both novelty levels.
The seen-2x/>2x/≥20x lines are an order of magnitude lower and grow more slowly. Nobody's allele
space is "explored" here in either species; what differs is the *slope*, not saturation.

**Per-ancestry (equal-N)** supplements
(`fig_saturation_by_recurrence_ancestry_{kir,hla}.png`, `fig_saturation_per_ancestry.png`,
`fig_coverage_completeness_ancestry_protein.png`): the same singleton-dominated, non-saturating
shape holds within every one of the 5 well-powered ancestries (AFR/AMR/EAS/EUR/SAS; MID excluded,
underpowered per the 39/44 precedent) for both species and both novelty levels — this is not a
pooling artifact. One genuinely new finding from breaking out coverage by ancestry: **HLA's
protein-level Chao2 completeness is higher than KIR's in every single one of the 5 ancestries**
(AFR 50.8% vs 34.3%; AMR 52.2% vs 38.2%; EAS 57.0% vs 48.3%; EUR 47.7% vs 28.6%; SAS 66.9% vs
40.5% — same ordering holds at CDS level too), even though the **pooled-ALL** numbers look nearly
tied (39.9% vs 39.2% at protein, 38.8% vs 41.6% at CDS, where KIR edges ahead). This is a real
aggregation effect (the pooled mixture of ancestries with different sample sizes and allele
richness can reverse an ordering that holds in every subgroup — the same shape as Simpson's
paradox), not a contradiction or a bug; both numbers are correct for what they measure
(within-ancestry vs. pooled-across-ancestry). **Report the per-ancestry numbers as the more
robust version of "which catalogue is more complete," and flag the pooled near-tie as
aggregation-sensitive** rather than picking one as "the" answer.

### Is novelty private or shared? (two-proportion z-test on singleton share)

| level | HLA singleton share | KIR singleton share | z | p | direction |
|---|---:|---:|---:|---:|---|
| genomic (upper bound) | 78.5% (13,491/17,188) | 79.4% (30,178/38,024) | −2.34 | 0.019 | KIR very slightly more private |
| CDS | 50.4% (687/1,363) | 56.1% (1,178/2,100) | −3.28 | 0.0010 | KIR more private |
| any-level novel | 85.5% (12,468/14,591) | 80.2% (30,090/37,503) | 13.83 | 1.8×10⁻⁴³ | HLA more private |
| protein-level novel | 97.0% (192/198) | 69.1% (735/1,063) | 8.15 | 3.7×10⁻¹⁶ | HLA more private |

**Do not overclaim from this table.** Four things temper it:
1. No multiple-testing correction is applied across these 4 rows.
2. The **direction flips** between genomic/CDS (KIR slightly more private / no strong difference)
   and any_novel/protein_novel (HLA clearly more private). "Any_novel"/"protein_novel" restrict to
   *already-novel* calls, a different (smaller, selected) population than "all calls at this
   level" — the two rows are not measuring the same underlying quantity, and a reader should not
   average them into one verdict.
3. At CDS/genomic (**the two least-selected, most trustworthy denominators**), the two species'
   singleton shares are close (50–56% vs 78–79%) and the effect sizes, while significant, are
   modest — this is not the dramatic "KIR novelty is systematically shared, HLA's is
   private" story one might expect from catalogue-coverage differences alone.
4. At protein_novel specifically, HLA's 198 novel proteins are a small, artifact-gate-cleaned set
   — a single re-annotation or one recurrent low-frequency real allele could shift the singleton
   share by several points. Treat p≈3.7×10⁻¹⁶ as "very unlikely to be pure chance," not as "the
   effect size is large or stable."

**Bottom line for A's private-vs-shared question**: both species' protein-level and CDS-level
novelty is **majority private** (singleton) — CDS: HLA 50%, KIR 56%; protein-level-novel: HLA 97%,
KIR 69%. KIR's novel *proteins* are recurrent more often than HLA's (consistent with KIR's known
gene-content/copy-number diversity producing the same derived protein independently in unrelated
haplotypes), but neither species' novelty is predominantly "shared" in any level tested.

## How to read each figure

- **`fig_saturation_by_recurrence.png`** (a–d): rows = species (HLA, KIR), columns = novelty level
  (any-level, protein-level), ALL ancestry pooled. Four lines per panel (seen 1x/2x/>2x/≥20x),
  direct end-of-line labels (no legend). Y = distinct alleles reaching that recurrence class, X = N
  unrelated people.
- **`fig_saturation_by_recurrence_ancestry_kir.png`** /
  **`_hla.png`**: same panel content, one row per ancestry (AFR/AMR/EAS/EUR/SAS), columns =
  any-level / protein-level, one figure per species. End-of-line labels shown only on the top
  (AFR) row to avoid repeating the same 4-entry key 5 times.
- **`fig_saturation_per_ancestry.png`** (a–b): genomic-level (upper-bound) discovery curve per
  ancestry, HLA vs KIR side by side — context for *why* the pooled-ALL slope is shallower than any
  single ancestry's slope (pooled N is capped by the smallest ancestry's own endpoint once you
  equalize N across ancestries).
- **`fig_coverage_completeness.png`** / **`_cds.png`** / **`_protein.png`** (a–b): pooled-ALL
  Good-Turing coverage (a) and Chao2 completeness (b), HLA vs KIR, at genomic (upper bound, titled
  explicitly), CDS, and protein identity respectively.
- **`fig_coverage_completeness_ancestry_protein.png`** (a–b): the same two metrics, protein level,
  broken out per ancestry (grouped bars, HLA vs KIR) — this is the figure behind the
  per-ancestry-reversal finding above.

## Version history (v1 → v4b) — do not use pre-v4b numbers

- **v1/v2**: allele "identity" was **name-based** (Immuannot's `<nearest template>...new`
  consensus label), which collapses distinct novel sequences that happen to share a nearest-
  template name. HLA's `protein`/`protein_novel` used the wrong CDS path and came back as a
  literal 0 for every classical gene (a real bug, not a biological finding). KIR protein-level
  coverage looked "complete" for 2DP1/2DL5B/2DL2 only because a filename-based pre-filter silently
  skipped bundled/mis-stemmed CDSseq reference files, and `pct_novel_protein` could read >100%
  because its denominator was genomic-level (not protein-level) s_obs.
- **v3**: fixed the protein-catalogue-coverage bugs above (pseudogene handling, no filename
  pre-filter, same-level `pct_novel_protein` denominator). **Still name-based identity** at the
  genomic level, and the KIR CDS/protein path still lacked HLA's artifact gate — this combination
  produced an internally-impossible number (6,300 distinct KIR "novel proteins" > 1,460 distinct
  KIR genomic alleles, impossible if protein were a true coarsening of a real sequence identity)
  that the orchestrator's own sanity check caught and rejected. **Every v3 number reported
  anywhere (including an earlier version of this README) is invalid** — see LOG entries around
  commit `751ee17` for the full incident writeup.
- **v4 / v4b (this version)**: sequence-hash identity at every level for both species (no more
  name-collapsing), the same artifact gate applied to both species' CDS/protein tracks, plus a
  **true genomic (span-level) identity** derived from `hap{1,2}.trimmed.fa` gene-span coordinates
  compared to IPD's `gen.fa.gz` (v4 originally set genomic == CDS, which silently dropped
  non-coding novelty — v4b adds the real span-level track Marc asked for). `check_identity_invariants()`
  enforces protein ≤ CDS ≤ genomic S<sub>obs</sub> and raises if violated. v4→v4b also fixed an
  HLA-specific bug (`_bare_gene()`: GTF gene_name `'HLA-A'` vs the bare `'A'` gene keys used
  elsewhere, which had made HLA's true-genomic track NA for every gene). This pass (2026-09-27)
  re-renders every 45/46 figure and both READMEs against the v4b tables and adds the CDS-level
  private-vs-shared test, per-ancestry catalogue-coverage check, and a third (CDS) panel to the
  KIR-vs-HLA completeness figure.

## Caveats

1. **Genomic (span-level) identity is an UPPER BOUND on novelty, not a clean measurement.** ~79%
   of KIR's genomic-level alleles and ~78% of HLA's are singletons, and the raw distinct-span
   counts are large relative to the number of catalogued alleles (KIR 38,024 distinct spans, HLA
   17,188 — `diagnostics_identity.tsv`'s `n_distinct_genomic`, summed over genes). Two things this
   pipeline cannot separate at the genomic level: **person-specific gene-span/UTR boundary
   differences** (the same allele, called with a slightly different 5′/3′ extent by the assembly)
   and **intronic assembly noise** that the CDS-based artifact gate structurally cannot see (it
   only inspects the coding sequence). **The headline KIR-vs-HLA comparison in this README and in
   `46`'s README is based on CDS and protein levels, not genomic**, for exactly this reason. A
   requested follow-up (not done in this pass) is a PAF-trimmed-core comparison — aligning each
   observed span against the reference and trimming to the aligned core before hashing, which
   would separate "different UTR extent" from "different core sequence." The exported TSVs from
   this v4b run do **not** include a per-gene span-length distribution (no such column in
   `diagnostics_identity.tsv`) — that would need a small addition to `44` itself, flagged as a
   follow-up, not fabricated here.
2. **Good-Turing incidence coverage ≠ Chao2 richness completeness.** 92–99% Good-Turing coverage
   does NOT mean "92–99% of the allele space is explored" — that claim is false. Chao2
   completeness (S<sub>obs</sub>/Chao2), 8–42% depending on species/level, is the number that
   answers "how much of the estimated true richness have we seen."
3. **KIR2DP1/KIR3DP1 are pseudogenes** with no catalogued reference protein — their
   protein-level fields are the literal string `"NA"`, never 0, in every table and figure.
4. **`<20` disclosure masking**: `recurrence_classes.tsv`'s per-gene eq1/eq2/gt2/ge20 cells are
   independently masked below 20; summing them naively under-counts the true pooled total (see
   Method). `coverage_chao2.tsv`'s per-gene s_obs is **not** masked (several KIR genes' s_obs at
   any_novel/protein_novel are below 20, e.g. KIR2DS3=51 at protein — check against
   `context/DECISIONS.md`'s open disclosure question before quoting an individual gene's
   raw s_obs outside the team; the pooled sums used throughout this README are all ≥100).
5. **No multiple-testing correction** across the 4 two-proportion z-tests reported above.
6. **44's own equal-N-slope headline table doesn't match a direct sum of `equal_n_slope.tsv`'s
   per-gene rows** (~18x off) — flagged, pre-existing (confirmed via `git show` predates v3), NOT
   fixed in this pass (out of scope for the WS-A/WS-B figures ask; a separate follow-up).

## Distilled

- Neither species' allele discovery saturates in this cohort at any identity level — singletons
  dominate and keep climbing, for both species, both novelty definitions, overall and within every
  ancestry.
- Genomic (span-level) identity is an upper bound (person-specific span/UTR + intronic noise
  inflate it, ~79%/78% singleton, 38k/17k distinct spans) — the KIR-vs-HLA headline uses CDS and
  protein levels instead.
- Pooled-ALL, CDS/protein completeness is near-tied (KIR 41.6%/39.2%, HLA 38.8%/39.9%) — but
  **within every one of the 5 well-powered ancestries, HLA is more complete than KIR** at both
  levels (a real aggregation effect, not a contradiction).
- Novelty is majority-private for both species at CDS (HLA 50%, KIR 56% singleton) and
  protein-novel level (HLA 97%, KIR 69%); direction flips at genomic/any-level — don't average
  across levels into one verdict.
- Good-Turing coverage (92–99%) and Chao2 completeness (8–42%) are different questions; only the
  latter answers "how much of the allele space is explored."
- KIR2DP1/KIR3DP1 protein = NA (pseudogenes), never 0.
- v1–v3 numbers are all invalid (name-based identity collapsed distinct alleles; KIR lacked HLA's
  artifact gate). Use only this v4b table.
