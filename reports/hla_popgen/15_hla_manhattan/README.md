# Per-site nucleotide diversity along HLA genes (`15_hla_manhattan.py`)

**What this measures:** at every base position of an HLA gene, how genetically variable the
population is at that exact position — and whether the variable positions cluster inside the
protein-coding sequence (CDS) or outside it (introns / UTRs).

**Headline result:** for every classical HLA gene tested, diversity is roughly **1.5–1.9× higher
inside the CDS than outside it**, concentrated at the exons encoding the peptide-binding groove.
This is the *opposite* of the textbook expectation for most genes, and it is the expected signature
of **balancing / diversifying selection** on antigen presentation.

Cohort: **12,261 people = ~24,000 haplotypes per gene** (All of Us long-read assemblies, annotated
by Immuannot).

---

## 1. What π (nucleotide diversity) actually means

π is the standard population-genetics measure of variability at a single position. Concretely:

> **π at a position = the probability that two haplotypes drawn at random from the population
> differ at that exact base.**

- π = 0 → everyone in the cohort carries the same base there (invariant / conserved).
- π = 0.5 → two random people differ there half the time.
- π ≈ 0.75 → the theoretical maximum when all four bases are equally common.

The formula used here is the standard **unbiased** estimator:

```
π_i = n/(n−1) × (1 − Σ_a p_a²)
```

where `p_a` is the frequency of allele `a` (a base, or a deletion) among the `n` haplotypes that
actually cover position `i`. The `n/(n−1)` term corrects for the fact that we observe a *sample*,
not the whole population — without it, π is systematically underestimated.

**Why π and not just "count of people who differ":** a raw count answers "how many differ from
*this particular reference*", which depends on which reference you picked. π is a property of the
population itself — it does not care which sequence you happened to call the reference. That makes
it comparable between genes and against published values. It is also the statistic used in the HLA
literature (sliding-window π along full-length genes, DnaSP-style), so our numbers can be checked
against theirs.

**Mean π for a region** (e.g. "mean π in CDS") is just the average of π over every position in that
region. So "CDS 0.0397 vs non-CDS 0.0266, ratio 1.49×" means: pick a random position inside the
coding sequence and a random position outside it; the coding position is ~1.5× more likely to
differ between two random people.

---

## 2. The core methodological problem, and how it is solved

Different HLA alleles are **different lengths**. In this cohort, HLA-DRB1 reference alleles range
from **10,850 bp to 16,110 bp** — over 5 kb of difference, almost entirely intron-length
polymorphism. So "position 5,000" is not the same biological place in two different alleles, and
per-position statistics cannot simply be pooled.

### What does NOT work (and why — this was a real failed attempt)

An earlier version of this analysis (`11_gene_diversity_track.py`) tried to get a shared coordinate
system by **restricting to haplotypes whose best-matching reference template was identical**. That
is self-defeating: a haplotype "matched template X" *precisely because it is nearly identical to X*.
Most such haplotypes have `template_distance = 0` — **zero variants by construction**. The result
was 8–15 variant positions across 90–260 haplotypes: the analysis had selected the least-diverse
possible subset and then plotted its emptiness.

More generally: measuring each sample against *its own* best-matching reference makes "difference"
mean "difference from a different thing for each sample", which is not a comparable quantity at all.

### What does work

Project **every** haplotype onto **one fixed canonical reference allele** per gene. This is the
field convention (it is what the IPD-IMGT/HLA alignment numbering does, and what MAFFT/DnaSP
sliding-window studies of full-length HLA genes do).

Crucially, this required **no re-alignment and no multiple-sequence-alignment build**, because of a
property of the existing data: `mm2.ipd.gen.paf.gz` is not a shortlist of plausible candidates —
it contains the **entire IPD/IMGT allele database aligned against every person's contig**. Verified
directly: all 20 probed people had exactly the same 4,810 distinct `HLA-A*` query names, 5,681 for
`HLA-B*`. So the canonical allele's alignment row is already sitting in every haplotype's PAF.

**Effect of the fix** (identical cohort, HLA-A):

| | variant positions found |
|---|---|
| old, template-stratified | 15 |
| new, canonical-projected | **1,733** |

And per-haplotype divergence from the reference went from ~0 for almost everyone to a realistic
**NM 0–226** distribution.

---

## 3. How a position is read off the alignment

For each haplotype we take the canonical allele's PAF row and walk its `cs` tag. In this
orientation **query = the canonical IPD allele, target = the person's assembled contig**:

| `cs` op | meaning here | effect |
|---|---|---|
| `:N` | N identical bases | person matches the reference for N positions |
| `*xy` | `x` = contig base, `y` = canonical base | the person carries base **`x`** at this position |
| `+seq` | bases in the canonical, absent from the contig | person has a **deletion** across those positions |
| `-seq` | bases in the contig, absent from the canonical | person has an **insertion** — *no canonical coordinate exists* |

Insertions therefore cannot be placed on the axis (there is no reference column for them). They are
counted in a **separate indel track** rather than being forced onto a position they do not have.
This is standard practice — inventing alignment columns would be worse than reporting them apart.

### Two subtleties that materially affect correctness

**Strand.** PAF rows for the same gene come back on *both* strands — measured almost exactly 50/50
(e.g. 12,191 `−` vs 11,976 `+` for HLA-A). For a `−` row minimap2 aligned the *reverse-complemented*
query, so the k-th query base consumed corresponds to canonical position `qend − 1 − k`, **not**
`qstart + k`. Handling this wrong silently mirrors half the dataset and smears the whole track.
(The predecessor script ignored strand entirely.) There is now an explicit regression test for it.

**Per-position coverage.** Forcing a fixed reference onto genuinely divergent haplotypes produces
*partial* alignments — observed query coverage ranged 0.53–1.00. So each position's denominator `n`
is **the number of haplotypes that actually cover that position**, not the cohort size. Positions
covered by fewer than 2 haplotypes get π = 0 (undefined, reported as zero to keep the axis
continuous).

---

## 4. Reference annotation

CDS and exon boundaries come from Immuannot's own reference metadata (`alleles.csv.gz`), in the
canonical allele's own coordinates — authoritative, not re-derived. Example record:

```
allele=HLA-A*01:01:01:01   geneRange=1..3503   nexon=8
CDS=301..373,504..773,1015..1290,1870..2145,2248..2364,2807..2839,2982..3029,3199..3203
UTR=1..300,3204..3503
```

**Peptide-binding groove exons** (highlighted in the figures): exons **2 and 3** for class I
(they encode the α1/α2 domains that form the groove), exon **2** for class II (β1 on
DRB1/DQB1/DPB1, α1 on DRA/DQA1/DPA1).

---

## 5. Choosing the canonical reference

Preference order: a well-known full-genomic reference allele if it is present for essentially all
haplotypes; otherwise auto-select **the longest allele present in ≥85% of probed haplotypes**.

The 85% floor matters. Requiring *strict* universality initially selected `HLA-C*16:85` (3,369 bp)
over full-length `C*07:02:01:01` (~4.3 kb) purely because the latter was missing from 1 of 12
probes — trading away ~20% of the gene to avoid a <10% dropout. Length now wins, and the resulting
dropout is reported explicitly as `n_missing_canonical` (**1.5–3.7%** of haplotypes across genes).

---

## 6. Results — full cohort (12,261 people)

| gene | canonical reference | length (bp) | haplotypes | variant sites | mean π CDS | mean π non-CDS | **ratio** |
|---|---|---|---|---|---|---|---|
| HLA-A | `A*01:01:01:01` | 3,503 | 24,167 | 1,733 | 0.0397 | 0.0266 | **1.49×** |
| HLA-B | `B*07:02:01:01` | 4,081 | 23,813 | 2,317 | 0.0372 | 0.0202 | **1.84×** |
| HLA-C | `C*17:01:01:30` | 4,325 | 23,823 | 2,326 | 0.0411 | 0.0233 | **1.77×** |
| HLA-DRB1 | `DRB1*13:02:01:01` | 13,941 | 23,622 | 4,534 | 0.0288 | 0.0151 | **1.91×** |

*(The wider 12-gene panel — adding DQA1/DQB1/DPA1/DPB1 plus the conserved controls DRA/E/F/G — is
reported in `panel_summary.md` alongside.)*

### Reading the figures

Each `*_manhattan.png` shows:
- **thin vertical bars** — π at each individual position; red = inside CDS, grey-lilac = outside
- **dark navy line** — 151 bp sliding-window mean, so the regional trend is visible above the noise
- **amber band** — the peptide-binding-groove exons
- **ribbon beneath the axis** — the gene model (exon boxes, groove exons in red)

---

## 7. Honest limitations

1. **Insertions relative to the canonical reference have no position** and are held in a separate
   indel track. A gap-padded IPD-IMGT/HLA multiple-sequence alignment (`ANHIG/IMGTHLA`,
   `A_gen.txt` etc.) would give those their own columns. That was judged unnecessary for this
   result and would require external download; it is the natural upgrade if insertion-level detail
   is ever needed.
2. **Reference-choice sensitivity is not yet quantified.** π itself is reference-independent in
   principle, but the *set of positions that exist* is not. Re-running a gene against a different
   canonical allele and confirming the profile is stable is a cheap, worthwhile robustness check
   that has not been run.
3. **1.5–3.7% of haplotypes are dropped** per gene because their PAF lacks the chosen canonical row.
4. **Multi-copy genes are excluded.** DRB3/4/5 have copy-number variation that makes copy
   assignment ambiguous; C4A/C4B go through a different Immuannot alignment path entirely.
5. **π here is computed over {reference, A, C, G, T, deletion} labels.** Deletions are treated as a
   distinct allele rather than missing data, which is a deliberate choice: a deletion is a real
   difference between two haplotypes.
6. **This is diversity, not selection *per se*.** Elevated π in the groove exons is the expected
   signature of balancing selection and matches the HLA literature, but a formal test of selection
   (e.g. dN/dS, Tajima's D) is a separate analysis. dN/dS in particular is blocked on per-codon
   synonymous classification, which failed against real Immuannot output for reasons documented in
   `scripts/hla_popgen/NEEDLE_VIEW_BRIEF.md`.

---

## 8. Reproducing

```bash
# full cohort, 4 classical genes, with figures
python3 scripts/hla_popgen/15_hla_manhattan.py --threads 6 --plot

# the wider panel (classical class I + II, plus conserved DRA/E/F/G controls)
python3 scripts/hla_popgen/15_hla_manhattan.py --panel --threads 6 --plot

# restyle figures from saved JSON without re-reading any PAF (~35 min saved)
python3 scripts/hla_popgen/15_hla_manhattan.py --panel --replot --out-dir ~/results/15_hla_manhattan
```

Runtime: **~37 minutes** for the full cohort. All genes share a single pass over each haplotype's
PAF, so adding genes is nearly free — the cost scales with haplotype count, not gene count.

Unit tests (including the strand regression guard):
```bash
python3 scripts/hla_popgen/tests/test_hla_manhattan.py
```
