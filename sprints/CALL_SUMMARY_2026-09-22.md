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

One more term is *probably* also a transcription error but is left uncorrected because there's no
way to confidently resolve it from context: **"the onshore paper"** — the TCR-repertoire-classification
benchmark paper Cole assigned as reading (VJ-usage + CDR3 features, XGBoost/logistic-regression
baselines, technical-methods appendix). Cole said he'd send it to Marc directly — get the real title
from that link rather than guessing here.

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
- **Concrete ask:** replot the DQ heatmap reordered by G1/G2 group so the incompatible quadrant
  reads as visually near-zero — he called this a strong candidate panel. Noted DP has **no
  equivalent G1/G2 rule** (unresolved — see below).
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
- **Baseline approach:** replicate a specific benchmark paper Cole is sending directly to Marc (the
  "onshore paper," title unconfirmed — see correction note above; benchmarks TCR-based disease
  classification using VJ-usage/k-mer features with linear/logistic regression and XGBoost
  baselines, and has a strong technical-methods appendix). **Action:** read it closely, replicate
  its exact train/test/cross-validation splitting methodology on the AoU dataset before
  improvising.
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
- Replot the DQA1~DQB1 LD heatmap grouped by G1/G2 protein groups.
- Turn the LD figure into a clean supplement: CSV of r² and D′, with a distinct color for "not
  observed / below n=20" rather than folding it into the existing scale.
- Simplify the deletions figure's color scheme for the supplement.
- Rework Figure 1 per the panel-by-panel feedback above (axis ticks, drop annotation, stricter
  ancestry filter, panel c/d reconsideration).
- Recompute novelty/coverage percentages and build per-ancestry saturation curves (inspired by the
  reference paper's Fig. 3E).
- Read Cole's assigned TCR-classification benchmark paper (get the real title/link from Cole);
  replicate its data-splitting methodology on the AoU TCR/BCR data.
- Sync with Aleix on TCR/BCR data shape and kick off the prediction baseline.
- Follow up on whether "zero people with two independent HLA-A deletions" is statistically
  meaningful at this cohort size.

**Cole:**
- Send the G1/G2 DQ-grouping reference paper, and the TCR-classification benchmark paper, directly
  to Marc.
- Continue leading the stat-gen half (LD, deletions, structural variation).
- Think further about peptide-groove ↔ TCR/BCR embedding relationship (no concrete plan yet).

## Open questions (unresolved)

- Is the DP-locus LD pattern independent biology, or confounded by DRB1 proximity / the DQ–DP
  recombination hotspot?
- Are the HLA-A / DQB1 / DRB1 apparent deletions real, or an assembly artifact (e.g., misassembled
  inversion)? Needs literature check + short-read validation.
- Whether/when to formally revisit short-read validation of the 38 novel callout alleles — flagged
  as high-value by Marc, not explicitly re-authorized in this call.
- Biological mechanism behind the G1/G2 DQ incompatibility rule (functional vs. antigen-presentation
  vs. something else) — "might be something in the literature."
- How to operationalize the peptide-groove ↔ TCR/BCR distance idea.
- Exact identity of "the onshore paper" — get the real title from Cole's link before citing it
  anywhere durable.
