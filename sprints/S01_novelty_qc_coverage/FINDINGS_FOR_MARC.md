# Sprint S01 — what I did, what changed, and what I think we should do

*Written for Marc, 2026-09-17, assuming no memory of the details. Every number here is reproducible
from committed files; each section says where.*

---

## 1. What I set out to do

You asked three things: pin down the novel-allele finding (you suspected the clustering had created
redundant groups), stress-test the QC that the paper leans on (0 phase switches but ~95%
concordance — where does the 5% come from?), and turn "how much of the allele space have we seen"
into something defensible, including "what would it take to see the rest".

I also re-derived the numbers behind the Figure 1 drafts, because the first two answers changed
them.

---

## 2. The five things that changed

### 2.1 You were right about the clustering — and the effect is large

Novel alleles were grouped by their **coding sequence only**, and no intron/UTR sequence was stored
anywhere. So every haplotype that shared a known coding sequence but differed *somewhere outside*
it was collapsed into a single "novel allele". The largest such group covered 7,713 people.

I recovered the actual genomic differences from alignment files that were already on disk
(`mm2.ipd.gen.paf.gz` stores, for every reference allele, exactly how a person's sequence differs).
For the classical genes:

- 708 old clusters actually contain **15,112 distinct genomic sequences**
- the biggest single old cluster holds **287** different sequences
- only **184** of those sequences are carried by 20 or more unrelated people
- **21% of them differ only by homopolymer indels** — the classic long-read error — so they are
  probably not real at all

*Consequence:* my 2026-09-15 callout ("a DPA1 allele in 1,531 people missing from IMGT") was wrong
and is retracted. Source: `reports/hla_popgen/25_noncoding_novelty_paf/`.

### 2.2 Most "protein-altering" novelty was artifacts, and the honest count is small

The old headline was "~3,000 novel alleles, 90.8% protein-altering". Recomputed per nomenclature
field, with each call checked against the reference database that Immuannot ships:

| step | calls |
|---|---|
| flagged novel by the caller | 280,695 |
| of which new protein (field 2) | 75,442 |
| …carrying no artifact flag | 22,122 |
| …coding sequence genuinely absent from IPD-IMGT, no frameshift/stop | 5,288 |
| **distinct novel proteins** | **1,404** |
| …seen in ≥2 unrelated people | **231** |

Two specific problems in the old pipeline explain the gap: 87% of "protein-altering" clusters carry
an artifact flag (mostly homopolymer indels producing fake frameshifts), and **25,282 calls labelled
novel have a coding sequence that is already catalogued** — a naming artifact, not biology.

Source: `reports/hla_popgen/24_novelty_by_field/`.

### 2.3 The novelty is not where we assumed — and this may be the better paper

Of the recurrent, clean novel proteins, **only 8 are in the eight classical HLA genes**. The rest
are in genes nobody types at population scale: TAP1/TAP2 (147), the class II accessory genes
DM/DO/DQA2/DQB2 (67), the DRB paralogs (49), MICA/MICB (36), and HLA-E/F/G (32).

So the defensible story is not "we found thousands of new HLA alleles". It is: **the classical HLA
genes are close to exhausted at protein level in a cohort this size, and the remaining undiscovered
MHC is the non-classical part that only long reads can type.** That is a cleaner claim, and it is
one nobody has made at this scale.

*Caveat I would check before leaning on it:* MIC/TAP novelty may partly reflect a thinner reference
for those genes — 22,657 of the 25,282 "already catalogued under another name" calls are MIC/TAP.

### 2.4 The QC now answers your question, and one leg of it was not evidence at all

**Where the 5% comes from.** Taking all 574 first-degree-or-closer pairs (chosen on All of Us's own
kinship, never on HLA data) and comparing the actual assembled sequences:

| what we see between relatives | share of comparable gene comparisons |
|---|---|
| identical sequence | **87.0%** |
| gene called in only one of the two | 6.5% |
| one relative looks homozygous (likely a collapsed haplotype) | 3.2% |
| differs by ≤3 bases in one place (ordinary sequencing error) | 2.3% |
| **genuinely different sequence, no benign explanation** | **0.8%** |
| gene alone on its own contig (fragmented assembly) | 0.2% |

So the "5%" was mostly presence/absence and ordinary error, not miscalled alleles.

**The phasing claim had to be rebuilt.** The old test excluded exactly the genes that disagreed
*before* counting switches, and never reported how many transitions it actually tested — so "0
switches in 545/545 pairs" could not distinguish "no switches" from "no test". The new version
reports the denominator (**3,021 testable transitions**) and measures its own power by injecting
synthetic switches: it detects them **100%** of the time, and still finds **zero real switches**.
That is now a real result.

**A number we should quote and did not have:** the five duplicate/identical-twin pairs are the only
true technical replicates. They differ at 121 bases out of 306,033 → **about Q34, one error per
2,500 bases**. That is worse than the Q46 implied indirectly, and it is the honest figure.

Two supporting signals: long-read calls show slightly more homozygosity than the same people's
short-read calls (1.10 vs 1.07 observed/expected), consistent with the small dropout above; and
short-read agreement is high for catalogued alleles (1.78 of 2 alleles matching) but low exactly
for novel-protein calls (0.93), which is expected but also a caution.

Source: `reports/hla_popgen/26_qc_relatives_v2/`.

### 2.5 "How much have we explored" — the answer flips depending on what you measure

The old answer ("5–20% discovered") came from estimating **how many alleles exist**, which is
notoriously unstable — and it was computed over artifact-contaminated clusters, with relatives left
in.

The stable quantity is **coverage**: the chance that the next person's haplotype carries an allele
we have already seen. For the eight classical genes, at protein level, coverage is
**99.0–99.8% in every ancestry group**. In plain terms: a catalogue built from this cohort already
types about 99 of every 100 haplotypes. Both statements are true at once — for HLA-B in
African-ancestry samples we see 105 distinct proteins, estimate ~177 exist, and still cover 99.2%
of haplotypes, because everything undiscovered is very rare.

**And the part I think is genuinely novel:** catalogues do not transfer across ancestries, and the
asymmetry is stark. Building catalogues of *equal size* (575 haplotypes) and using them to type
other groups, averaged over the classical genes:

- African-ancestry catalogue → types 97% of European haplotypes
- European-ancestry catalogue → types **83%** of African haplotypes
- East Asian catalogue → types **72%** of African haplotypes

This is a direct, quantitative reference-panel design argument: sampling African-ancestry genomes
buys more global coverage per genome than sampling European ones. It is the rigorous version of
your "if we sampled this many people from Africa…" question.

Source: `reports/hla_popgen/27_allele_space_coverage/`.

---

## 3. What is solid, and what I still doubt

**Solid:** the funnel in 2.2 (counts recomputed from raw calls with sequence-level checks); the
pooling in 2.1 (the cs-string parsing was hand-verified by an independent reviewer on both DNA
strands); the QC decomposition and switch power in 2.4; the coverage and transfer matrix in 2.5
(estimators implemented from the published formulas with 22 unit tests, including brute-force
checks).

**Still uncertain, and stated as such:**
1. Whether MIC/TAP novelty is biology or a thinner reference (test: reference completeness per gene).
2. Coverage assumes haplotypes are independent draws; relatives are removed, but structure within an
   "ancestry" is not.
3. The novel proteins have no orthogonal confirmation yet. The strongest available check is the
   RNA-seq data: a novel protein should show its variant in that person's transcripts. That needs
   Aleix's cohort and is not done.
4. Extrapolation beyond about twice the current sample size is unreliable; flagged in the outputs.
5. The Middle Eastern group (~575 haplotypes) is thin; its numbers carry wide intervals.

**Two mistakes worth knowing about**, both now regression-tested: a sequence file's header number
is a file-wide record counter, not a copy number (misreading it briefly produced an impossible 78%
discordance); and a run once silently executed a stale copy of a script, producing a result that
looked real. I now clear caches and verify which file ran.

**One old bug found in passing:** `21_hla_manhattan.py` had an off-by-one for deletions on one DNA
strand. I verified it only affects a field nothing reads — **the published per-site diversity
figures are unaffected** — and patched it anyway.

---

## 4. Figures

`reports/hla_popgen/28_figure1/figure1.png` (+ PDF) is a redrawn Figure 1 with four panels: the
audit trail from 280,695 flagged calls to 231 recurrent novel proteins; per-gene composition showing
novelty is gene-specific; the cross-ancestry transfer matrix; and coverage against catalogue size.
Its README holds a draft caption. Supporting figures for each analysis are in the numbered folders,
each with a plain-language README.

---

## 5. Decisions I need from you

1. **Reframe the paper around the non-classical MHC?** (§2.3) I think yes, with the classical genes
   as the validation that the method works.
2. **Which headline for novelty** — "231 novel proteins recurrent in unrelated people" (strict) or
   "1,404 distinct novel proteins" (inclusive)? I would lead with the strict one.
3. **Push the branch.** I could not push (`fig1-drafts-and-research-map`); publishing to the public
   repo is blocked for me. Everything is committed locally.
4. **The disclosure question stands:** `novel_alleles.tsv`, already committed to the public repo,
   lists participant counts below 20. All my new outputs suppress them. Your call, with the sponsor.
5. **Next sprint priority?** My ranking: (a) confirm novel proteins with RNA-seq reads, which also
   opens the join with Aleix; (b) the non-classical MHC atlas (TAP/MIC/DM/DO) as the paper's spine;
   (c) the remaining supervisor items (LD, structural variation, KIR).
