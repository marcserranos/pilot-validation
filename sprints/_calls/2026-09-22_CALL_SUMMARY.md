# Call summary — 2026-09-22, Progress update #8

*Post-call writeup. Source: full transcript pasted into the assistant session on 2026-09-22,
cross-referenced against `REPORT #8 - HLA_AoU - Marc Serrano and Aleix Ruf.pdf` (the slide deck
presented) and the repo's existing KIR/pipeline docs. Companion to
[`CALL_BRIEF_2026-09-22.md`](CALL_BRIEF_2026-09-22.md), which was the *pre*-call prep note; this is
what actually happened on the call.*

**Presenting:** Marc. **Supervisors:** Cole (drives most of the discussion), joined partway by
**David Bonet** (tabular-modeling expert) for the prediction/phenotype discussion. Aleix was not on
the call — his update was relayed via Slack and addressed indirectly.

---

## Transcription corrections

Two terms were garbled by the call's automatic transcription. Both are corrected throughout this
document; the raw transcript still has the original wording.

| In the raw transcript | Corrected to | Why |
|---|---|---|
| "current genes" / "curve genes" | **KIR genes** | Context is unambiguous — "natural killer [cell] interactions with the HLA," a region never yet called by this pipeline, worth adding to Figure 1. The project's own docs already name this gene family and flag it as intentionally excluded: `scripts/hla_popgen/SCHEMA.md:69` — *"the 17 KIR genes — expected to be entirely absent (chr19, outside the trim window)"* — and `context/EXPERIMENTS.md:720` / `sprints/CALL_BRIEF_2026-09-22.md` §Block E both call this "KIR: zero calls, as expected." |
| "Chromosome 14, I think" | **chromosome 19** | Same cross-check as above — KIR is chr19 in every existing repo reference (`SCHEMA.md`, `EXPERIMENTS.md`, `RUNBOOK.md`). Cole misspoke or the transcript mis-heard him; chromosome 14 is not corrected anywhere else in the call, and nothing in the discussion depends on the number itself, so this doesn't change any action item — just don't scope the rerun around chr14. |
| "VJ CAMERS" | **VJ k-mers** | "Camers" isn't a real term in this space; "k-mers" is the standard TCR/BCR repertoire feature (V/J-gene-usage plus CDR3 k-mer counts) and is phonetically what an ASR system garbles into "camers." No exact repo precedent for this term (Aleix's repertoire work lives on his branch, not merged here), but it's the only sensible reading given the surrounding sentence ("do like a normal regression... on the CAMERS as well"). |

**Resolved after the call:** "the onshore paper" Cole assigned as reading is confirmed (via the
bioRxiv link he sent that evening — see References below) to be **BenchRep-T** (Im, Cohen-Lavi,
Buendia, Kundaje, Boyd) — not "onshore" at all, just another ASR mangling. Its own headline result
— "simple baselines prove competitive, with tree-based models trained on V- and J-gene usage and
short sequence motifs approaching the classification performance of more complex methods" — is
independent confirmation that "VJ k-mers" is the right correction above (this is exactly the
feature type the paper benchmarks) and directly backs Cole's "start with the base model, don't get
cute" instruction.

---

## References Cole sent after the call (Slack, evening of 2026-09-22)

Five items, sent as a batch after the meeting ended — mapped here to the discussion topic each one
supports.

| Sent | Link / file | Identified as | Maps to |
|---|---|---|---|
| 7:34 PM | [sciencedirect.com/.../S0198885910005033](https://www.sciencedirect.com/science/article/pii/S0198885910005033) | Human Immunology article (pii prefix). **Could not confirm title/authors** — ScienceDirect returned HTTP 403 to automated fetch. Sent immediately before the PMC G1/G2 paper below, so likely an earlier/companion reference on HLA-DQ heterodimer compatibility — **open the link directly to confirm** before citing it anywhere. | §3 Linkage disequilibrium (DQ G1/G2 rule) |
| 7:44 PM | [nature.com/articles/s41586-026-10667-5](https://www.nature.com/articles/s41586-026-10667-5) | **"Analysis of 173,303 exomes and genomes in the Pakistan Genome Resource"**, *Nature*, 2026. A biobank with high familial relatedness, broadening the catalogue of human genetic variation. (Paywalled — could not fetch full text automatically, title/scope confirmed from the shared preview.) | §6 Figure 1 discussion — source of the Pakistan/endogamy tangent and very likely the paper with the saturation-curve figure ("Fig 3, panel E") Cole referenced for inspiration |
| 7:56 PM | [pmc.ncbi.nlm.nih.gov/articles/PMC9121842](https://pmc.ncbi.nlm.nih.gov/articles/PMC9121842/) | **Petersdorf EW, Bengtsson M, Horowitz M, McKallor C, Spellman SR, Spierings E, Gooley TA, Stevenson P. "HLA-DQ heterodimers in hematopoietic cell transplantation." *Blood* 2022;139(20):3009–3017.** This is **the** source of the G1/G2 rule — see the precise definition folded into §3 below. | §3 Linkage disequilibrium — resolves the "might be something in the literature" open question about *why* G1/G2 incompatibility exists: it's not just a population-genetics curiosity, it's a clinically established transplant-relapse risk factor |
| 8:04 PM | [biorxiv.org/content/10.64898/2026.06.09.727013v1](https://www.biorxiv.org/content/10.64898/2026.06.09.727013v1) | **Im C, Cohen-Lavi L, Buendia A, Kundaje A, Boyd SD. "BenchRep-T: A Systematic Evaluation of T-Cell Repertoire-Based Disease Diagnostics."** bioRxiv, 2026. Benchmarks 9 methods (statistical → deep learning) for TCR-based disease classification; finds tree-based models on V/J-gene usage + short motifs (k-mers) competitive with complex methods. | §9 Prediction workstream — this **is** "the onshore paper" (see correction above); the assigned reading/methodology template |
| 8:41 PM | `dq_ld.pdf` — now saved at [`reference/cole_dq_g1g2_target_figure.pdf`](../reference/cole_dq_g1g2_target_figure.pdf) | Cole's own recreation target, now confirmed by inspection (see full spec in §3 below). Cole: **"I would recreate this."** | §3 Linkage disequilibrium — direct, concrete spec for the replotted DQ heatmap (see action items) |
| 10:48–10:50 PM (later Slack message, same evening) | — (text only, no link) | Cole, looking at his own figure: *"I'm so interested in the phasing because I wonder if the light blue in the 'predicted incompatible' cells are actually phase errors? ... these are rare but occur at low frequency because there is still some recombination, but generally these are probably deleterious."* | §2 QC/phasing × §3 LD — a concrete cross-check request, folded into §3 below |

---

## The one-paragraph version

Every QC/LD/deletion/diversity result from last sprint held up and got positive reactions — Cole
called this "a pretty good start for what Figure 1 could look like" and is "super pumped." The big
new items are: (1) **KIR genes** (chr19) get added as a priority rerun of the long-read pipeline —
never called before because the pipeline was trimmed to the HLA region only; (2) the project
formally splits into two halves — **prediction/ML** (Marc + Aleix) and **stat-gen** (Cole) — both
anchored by Figure 1; (3) concrete, actionable feedback on every Figure 1 panel; (4) a specific
reading/methodology assignment for the disease-prediction workstream, with explicit guidance to
start simple (TCR/BCR repertoire only, baseline regression → XGBoost) and add HLA later.

---

### 1. Novel allele discovery (slides 4–5)

- **Slide 4** (novel protein alleles per gene, general vs. strict) and **slide 5** (the ≥20-carrier
  callout list — TAP1/TAP2/MICB/HLA-G/DQB2, colored by ancestry) were shown without pushback — this
  is the material from `32_novel_callouts` / `24_novelty_by_field`.
- Cole raised **AlphaMissense**-style rare-variant collapsing (treat rare protein-altering alleles
  as presumed-deleterious and collapse them for association testing) as an interesting technique,
  but **explicitly parked it** — "I wouldn't say it's focused right now" — and noted it fits
  short-read data better than long-read.
- Marc asked directly: use short reads to *validate* the long-read calls, or to *discover* new
  alleles too? **Cole's answer: validation only, and stay focused on this cohort** — "I think we
  should be careful... to keep it to focusing on this cohort." Soft no on expanding scope to
  short-read discovery right now; validation of the 38-allele list itself remains open (see Open
  Questions).

### 2. QC and relatedness (slides 8–11)

- Slide 8: the KING-style relatedness coefficient (φ) and its threshold table — what Marc called
  the "king sheep" score.
- Slide 9: 5 duplicate/twin pairs, 121/306,033 bases differing → **~1 error per 2,500 bases**, used
  as the assay's baseline error rate.
- Slides 10–11: physical phasing yield per gene pair (DQA1~DQB1 96.31%, DPA1~DPB1 96.96%,
  DRB1~DQB1 91.72%, B~C 86.88%, A~B 30.24%).
- **Cole's request:** an individual-level file marking whether each person's DQ/DR pairs are
  physically phased (same contig) vs. not — he wants to condition downstream analyses on phasing
  confidence rather than discard the ~3–10% unphased. Marc clarified the per-haplotype calls
  already exist as a table; what's missing is the phasing-confidence flag alongside them.
- Marc pushed back usefully: not-same-contig doesn't mean *not phased*, just *less confident* —
  worth keeping that nuance when this file gets built.

### 3. Linkage disequilibrium (slides 12–13)

- Slide 12 (D′ by ancestry, 5 gene pairs) and slide 13 (DQA1~DQB1 r² heatmaps per ancestry) were
  both very well received.
- Cole flagged specific cells as striking (DQA1\*04-type patterns differing between
  American/East Asian vs. European) and called it **"a great supplement figure... this should be a
  table too — a simple CSV with r² and D′."**
- Noticed the LD is **bimodal** — pairs are either near-perfect LD or near-zero, rarely in between —
  and explained *why*: DQA1/DQB1 alleles split into **group 1 (G1) / group 2 (G2)** protein groups;
  a G1+G2 pairing on the same haplotype is non-functional and gets purged by background selection
  (cis vs. trans inheritance). **He already ran this analysis in parallel and got the same answer
  Marc did.**
- **The precise rule** (Petersdorf et al. 2022, *Blood* — see References below, sent 7:56 PM):
  **Group 1 (G1)** = any **DQA1\*02/\*03/\*04/\*05/\*06** α-chain paired with any **DQB1\*02/\*03/\*04**
  β-chain. **Group 2 (G2)** = **DQA1\*01** α-chain paired with any **DQB1\*05/\*06** β-chain. A G1α
  with a G2β (or vice versa) forms a non-functional heterodimer in cis — the thing background
  selection purges. This isn't just a population-genetics curiosity: the same paper shows
  **G1G2/G2G2 genotypes carry significantly higher post-transplant relapse risk** than G1G1 in
  hematopoietic cell transplantation, and relapse risk scales with the number of G2 molecules
  present — resolving the "might be something in the literature" open question from the call itself.
- **Concrete ask:** replot the DQ heatmap reordered by G1/G2 group so the incompatible quadrant
  reads as visually near-zero — he called this a strong candidate panel. Noted DP has **no
  equivalent G1/G2 rule** (unresolved — see below).
- **Direct spec via `dq_ld.pdf`** (Slack, 8:41 PM — now saved at
  [`reference/cole_dq_g1g2_target_figure.pdf`](../reference/cole_dq_g1g2_target_figure.pdf)). Cole's
  recreation target, confirmed by inspection:
  - A **DQA1 (rows) × DQB1 (columns) heatmap**, one cell per specific allele pair (4-field
    resolution, e.g. `05:01:01:01`), colored by **signed phased D′** on a diverging scale from
    **−1.00 (blue) to +1.00 (red)** — explicitly labeled "Signed phased D′" on the colorbar. This is
    a different statistic from what's in slides 12–13: those used **multiallelic D′ (Hedrick 1987,
    unsigned)** and per-ancestry r², summarized per gene pair or per top-12-alleles; Cole's version
    is **pairwise, allele-specific, signed, and phased** — the sign is what makes "avoided together"
    (deep blue) visually distinct from "no signal" (white), which an unsigned/r² statistic can't do.
  - The alleles are **sorted into exactly two blocks per axis**, splitting at the G1/G2 boundary:
    DQA1 rows split into **DQA1\*01 (the G2 α-allele)** on top vs. **DQA1\*02/03/04/05/06 (G1
    α-alleles)** below; DQB1 columns split into **DQB1\*05/06 (G2 β-alleles)** on the left vs.
    **DQB1\*02/03/04 (G1 β-alleles)** on the right — producing a clean 2×2 grid with a black
    crosshair dividing it.
  - **Top-left quadrant (DQA1\*01 × DQB1\*05/06, i.e. G2×G2)** and **bottom-right quadrant
    (DQA1\*02–06 × DQB1\*02/03/04, i.e. G1×G1)** show real, patterned structure — clusters of strong
    positive D′ (near +1, specific allele pairs that travel together) amid mostly weak/negative
    background.
  - **Top-right quadrant (DQA1\*01 × DQB1\*02/03/04)** and **bottom-left quadrant
    (DQA1\*02–06 × DQB1\*05/06)** — the cross-group, biologically-incompatible pairings — are
    overwhelmingly **uniform deep blue (D′ ≈ −1)**, and Cole's own figure labels each of these two
    quadrants **"Predicted incompatible"** directly on the plot.
  - **Recreation spec:** compute pairwise signed, phased D′ (not r², not multiallelic/unsigned D′)
    between individual DQA1 and DQB1 alleles at high-enough resolution to match (looks like 4-field,
    i.e. `05:01:01:01`-style calls, at whatever frequency floor keeps cells estimable); sort/group
    rows and columns by the Petersdorf G1/G2 rule; render as a diverging heatmap (−1 to +1) with the
    2×2 grid and "Predicted incompatible" labels on the off-diagonal blocks.
- **Cole's follow-up (Slack, 10:48–10:50 PM, same evening), looking at his own figure:** *"I'm so
  interested in the phasing because I wonder if the light blue in the 'predicted incompatible' cells
  are actually phase errors? ... these are rare but occur at low frequency because there is still
  some recombination, but generally these are probably deleterious."* In other words: the
  "incompatible" quadrants aren't *perfectly* uniform −1 — a handful of cells are lighter blue
  (D′ closer to 0, i.e. the pairing does occur at low frequency), and Cole isn't sure whether that
  residual signal is (a) genuine rare recombination between the DQ genes producing a real,
  presumably-deleterious G1/G2 heterodimer, or (b) a **phasing artifact** — a switch error that
  falsely places a G1 and G2 allele on the same called haplotype when they're really in trans. This
  is directly testable with data this project already has: **cross-tabulate the light-blue
  (nonzero) cells in the incompatible quadrants against the per-individual physical-phasing-
  confidence flag from §2** (same-contig vs. not, from the DQA1~DQB1/DPA1~DPB1 phasing-yield work).
  If the light-blue occurrences cluster in the ~3–8% of individuals with *lower*-confidence phasing,
  that's evidence for phase error; if they're spread evenly regardless of phasing confidence, that
  supports Cole's "real, rare, and deleterious" reading. **New, concrete action item** (folds
  together the §2 phasing-confidence file Cole already asked for and this figure).
- Also flagged: this pattern could be partly confounded by DRB1's physical proximity rather than
  being independently causal; he's tried conditioning on DRB1 / higher-order LD without a clean
  answer yet.

### 4. Deletions (slides 14–15)

- Slide 14 (gene deletion % across all genes, green=positive control DRB3/4/5, grey=false-positive
  control, orange=everything else) and slide 15 (deletion frequency by ancestry for
  DRB5/Y/DRB4/DRB3/C4B/K/U/T) — **"Oh wow, deletions, okay, this is very interesting."**
- Cole connected this to a paper he found (published this past summer) arguing selective pressure
  favors **higher HLA copy number** for antigen-presentation diversity — DRB3/4/5 variability fits
  that story, and he wants a panel showing DRB3/4/5 copy number by ancestry (largely already
  delivered by slide 15).
- The **HLA-A deletion signal** got the most scrutiny — "that is a big result... how is that
  plausible?" He raised alternative explanations (misassembled inversions are "crazy" but real in
  genomics) and wants **short-read cross-validation** before trusting it.
- Marc proposed a statistical framing: check whether the cohort has **zero individuals with two
  independent HLA-A deletions** (biallelic loss would be near-lethal/implausible), and whether
  that's meaningful given the expected rate — left as a follow-up, not resolved live.
- **Concrete ask:** simplify the deletion figure's color scheme (drop the orange/green distinction,
  keep everything else) for the supplement.

### 5. Peptide groove diversity (slide 16)

- Slide 16 (DRB1 per-residue diversity track + groove/non-groove enrichment ratios, with pseudogene
  HLA-E as contrast) landed well — "it's really the class II alleles that have the most diversity."
- Cole floated a **future idea, not scoped**: relating peptide-groove distance (or embeddings of it)
  to distances between TCR/BCR motif clusters — speculated the SCEPTR embedding model may already
  implicitly pick up groove signal via its training. No action assigned; "worth thinking about."
- Marc asked what to build next on this thread; Cole's answer: keep it simple for now, "we can keep
  thinking about this... might be more obvious later."

### 6. Figure 1, panel by panel (slide 17)

| Panel | Content | Feedback |
|---|---|---|
| **a** — admixture barplot | Ancestry proportions, 6 groups | Good as-is; **drop the "≥98%" percentage annotations** — visually obvious already |
| **b** — HLA-B ternary (AFR/AMR/EUR corners) | Novel vs. catalogued alleles by ancestry composition | Confirmed as the right gene choice; **add 0–100 axis tick marks**; **filter to a stricter ancestry-probability threshold (95–98%)** — current admixed individuals (esp. AMR, and Europeans with African admixture) visibly drag the cloud toward the triangle's center |
| **c** — reference-catalogue incompleteness, by gene/ancestry | Bar chart | Cole prefers **a different, already-existing figure** (the ancestry × gene novelty-rate panel from `33_figure1_v3`, described as "the 800,000 one") — says it communicates the point better and could even become its own panel |
| **d** — novel proteins by gene ("not in classical genes") | Same data as slide 4 panel a | Clarify it's "proteins not in IMGT/HLA," not literally class I only; **candidate to merge with panel c**, or expand to show more genes instead of dropping the ancestry split (ancestry breakdown → supplement only) |

- Marc separately proposed reworking the novelty/coverage math (recalculate % of allele space
  explored per ancestry, using IMGT richness) and building **per-ancestry saturation/discovery
  curves** — inspired by an example figure Cole referenced from another paper (a
  homozygous-loss-of-function-variant paper, Figure 3 panel E), which also prompted a tangent on
  Pakistani/South Asian cohort under-saturation due to strong population structure/endogamy (PCA
  reportedly recovers caste structure). Expectation: African-ancestry curve will show the least
  saturation.

### 7. KIR genes — new priority

- Cole: the pipeline currently only calls genes in the HLA region; **KIR genes were never called**.
  He wants the long-read Immuannot pipeline **rerun restricted to the KIR region** (chr19 — see
  correction table above).
- Rationale: KIR interacts closely with HLA (receptors on NK cells), expected to show **more
  germline diversity and novel alleles** than HLA itself, and Cole wants it in — potentially with
  its own Figure 1 panel(s).
- **Action for Marc:** scope the region size and rerun cost for next week. Cole: cost isn't a
  concern ("you guys have been very reasonable"); region may be comparable size to HLA (~150kb) but
  structurally complex (more duplications/deletions).

### 8. Project scope and structure (strategic discussion)

- Cole is explicit about **competitive pressure** — worried other groups could publish something
  similar first, even if lower quality ("the worst thing is when someone does something first and
  not well... they absolutely steal the thunder"). Wants to target a high-tier journal, which argues
  for **focus over breadth** — not stuffing 10 sub-projects into one paper.
- **Project splits into two halves**, both anchored by Figure 1 as the introductory figure:
  1. **Prediction / ML** (disease classification from HLA + TCR/BCR repertoire) — owned by
     **Marc and Aleix**.
  2. **Stat-gen** (population genetics: LD, deletions, structure) — Cole will drive this side
     himself, citing higher risk of subtle errors.
- Cole's belief: **HLA alone is low-signal**; the TCR/BCR repertoire is where the real signal is,
  with HLA acting as a structuring/conditioning variable rather than a direct predictor.

### 9. Prediction workstream — concrete guidance

- **Data:** only **common, two-field HLA alleles** should go into the bucket for this — novel
  alleles are "interesting but dangerous," too rare to draw conclusions from yet. (Marc has already
  pushed most of what's needed; add the KIR calls once available.)
- **Start with TCR/BCR repertoire features alone**, not mixed with HLA — HLA and TCR features are
  highly correlated, so adding HLA early mostly just increases dimensionality without adding
  independent signal. Mix in HLA later, especially for cases where the repertoire isn't fully
  observed.
- **Baseline approach:** replicate **BenchRep-T** (Im, Cohen-Lavi, Buendia, Kundaje, Boyd, bioRxiv
  2026 — the paper Cole called "the onshore paper," sent 8:04 PM, see References) — it benchmarks
  9 methods (statistical → deep learning) for TCR-based disease classification and finds that
  **tree-based models on V/J-gene usage + short sequence motifs (k-mers) are competitive with much
  more complex methods**, with a strong technical-methods appendix on splitting/evaluation. This is
  independent published confirmation of Cole's "start simple" instruction, not just his opinion.
  **Action:** read it closely, replicate its exact train/test/cross-validation splitting
  methodology on the AoU dataset before improvising.
- David Bonet (tabular-modeling expert) will advise, but the explicit instruction is **"start with
  the base model... don't get cute"** — beat a simple baseline before reaching for fancier
  architectures. SCEPTR-style embeddings are fine to explore but not the starting point.
- **Philosophy:** the goal isn't pure predictive accuracy — a strong, interpretable predictor is
  valuable because it points to concrete, testable biology (e.g., a feature-interaction worth
  chasing mechanistically), not because of the accuracy number itself.

### 10. Phenotype selection guidance

- **Avoid cancer phenotypes** — Cole considers cleaning them "maybe an impossible task" (remission
  status, wildly variable treatment, immunotherapy regimens not reliably extractable from EHR).
- **Prefer chronic/autoimmune conditions** — explicitly named HIV and B-cell-related disorders as
  good candidates, presumably because they're more reliably captured and more stable over time in
  EHR data.

---

## Decisions made

1. Short-read data is for **validating** long-read calls, not for discovering new alleles — stay
   scoped to this cohort.
2. AlphaMissense-style rare-variant collapsing is a parked idea, not current work.
3. **KIR genes get added** — rerun the long-read pipeline restricted to chr19.
4. Project formally splits: **prediction (Marc + Aleix)** vs. **stat-gen (Cole)**, unified by
   Figure 1.
5. Prediction work starts from **TCR/BCR repertoire only**, common two-field HLA alleles only (no
   novel alleles), baseline models before fancy ones, and Cole's assigned paper as the
   methodological template.
6. Phenotype focus: autoimmune/chronic conditions (HIV, B-cell disorders) over cancer.
7. Figure 1: ternary plot gets stricter ancestry filtering + axis ticks; panel c likely replaced by
   the ancestry×gene novelty figure; panels c/d may merge; admixture panel loses its percentage
   annotation.

## Action items

**Marc (primary):**
- Scope size/cost of rerunning the long-read Immuannot pipeline over the KIR region (chr19) —
  report back next week.
- Push KIR calls to the bucket once run; already done for common two-field HLA alleles.
- ~~Get `dq_ld.pdf` out of Slack and into the repo~~ — **done**, saved at
  [`reference/cole_dq_g1g2_target_figure.pdf`](../reference/cole_dq_g1g2_target_figure.pdf).
- **Recreate Cole's DQ figure exactly**: pairwise **signed, phased D′** (not r², not
  multiallelic/unsigned D′) between individual DQA1 and DQB1 alleles, using the precise Petersdorf
  et al. 2022 rule (G1 = DQA1\*02/03/04/05/06 + DQB1\*02/03/04; G2 = DQA1\*01 + DQB1\*05/06) to sort
  rows/columns into the 2×2 grid, diverging colormap −1 to +1, "Predicted incompatible" labels on
  the two cross-group quadrants. Full spec in §3 above; reference figure at the path above.
- **New cross-check Cole asked for (10:48–10:50 PM Slack):** pull the individual-level
  physical-phasing-confidence flag from the §2 phasing work (same-contig vs. not, per DQA1~DQB1
  pair) and cross-tabulate it against which individuals carry the light-blue (nonzero D′, not −1)
  cells inside the "predicted incompatible" quadrants — tests whether that residual signal is
  **phase error** (concentrated in low-confidence-phasing individuals) or **genuine rare
  recombination** (evenly spread, Cole's working guess, and "probably deleterious" either way).
- Turn the LD figure into a clean supplement: CSV of r² and D′, with a distinct color for "not
  observed / below n=20" rather than folding it into the existing scale.
- Simplify the deletions figure's color scheme for the supplement.
- Rework Figure 1 per the panel-by-panel feedback above (axis ticks, drop annotation, stricter
  ancestry filter, panel c/d reconsideration).
- Recompute novelty/coverage percentages and build per-ancestry saturation curves (inspired by the
  Pakistan Genome Resource paper's saturation-curve figure — confirm which panel once the Nature
  paywall is cleared).
- Read **BenchRep-T** (Im et al., bioRxiv 2026.06.09.727013v1) closely; replicate its exact
  train/test/cross-validation splitting methodology on the AoU TCR/BCR data.
- Sync with Aleix on TCR/BCR data shape and kick off the prediction baseline.
- Follow up on whether "zero people with two independent HLA-A deletions" is statistically
  meaningful at this cohort size.
- Confirm the ScienceDirect link's actual title (automated fetch was blocked, 403) before citing it.

**Cole:**
- Already sent: the Petersdorf DQ G1/G2 paper (PMC9121842), the BenchRep-T preprint, the Pakistan
  Genome Resource paper, a ScienceDirect reference (title TBD), and `dq_ld.pdf` (Slack, 2026-09-22
  evening — see References section below).
- Continue leading the stat-gen half (LD, deletions, structural variation).
- Think further about peptide-groove ↔ TCR/BCR embedding relationship (no concrete plan yet).

## Open questions (unresolved)

- Is the DP-locus LD pattern independent biology, or confounded by DRB1 proximity / the DQ–DP
  recombination hotspot?
- Are the HLA-A / DQB1 / DRB1 apparent deletions real, or an assembly artifact (e.g., misassembled
  inversion)? Needs literature check + short-read validation.
- Whether/when to formally revisit short-read validation of the 38 novel callout alleles — flagged
  as high-value by Marc, not explicitly re-authorized in this call.
- ~~Biological mechanism behind the G1/G2 DQ incompatibility rule~~ — **resolved**: Petersdorf et al.
  2022 (*Blood*, PMC9121842) ties it to post-transplant relapse risk, i.e. functional
  antigen-presentation consequence, not just a population-genetics pattern. See §3.
- How to operationalize the peptide-groove ↔ TCR/BCR distance idea.
- **New:** are the light-blue (non−1) cells inside the "predicted incompatible" DQ quadrants phase
  errors or genuine rare recombination? Cole's working guess is the latter ("probably deleterious"),
  but he explicitly wants it checked against phasing confidence rather than assumed — see the
  cross-check action item above.
- The ScienceDirect reference (pii S0198885910005033) — title unconfirmed, automated fetch blocked
  (403); open it manually to see whether it's an earlier/companion G1-G2 reference or something
  else.
- Whether the Pakistan Genome Resource paper is in fact the source of the saturation-curve figure
  Cole referenced (Fig. 3, panel E) — plausible given the topic match (large biobank, high
  familial relatedness) but not confirmed page-by-page (Nature blocked automated fetch behind a
  login wall).
