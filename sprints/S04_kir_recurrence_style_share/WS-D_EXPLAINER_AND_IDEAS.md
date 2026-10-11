# WS-D — WS6 explained, experiment ideas ranked, naive ML/simple tests

*S04 WS-D deliverable. Written for Marc, plain language, no jargon without a one-line gloss.
Companion code: `scripts/hla_popgen/47_naive_ml_tests.py`, `48_kir_hla_ligand_cooccurrence.py`
(local, synthetic-tested; not yet run on the VM). Source material: `42_repertoire_baseline.py` +
its README, `ORCHESTRATOR_HANDOFF.md` §2 WS-D (the seed lists), the 2026-09-22 call notes.*

---

## Part A — WS6, in plain language

**What was asked:** can we predict whether someone has a disease from their immune "repertoire"
— the actual T-cell and antibody receptor sequences found in their blood?

**The data.** Aleix's pipeline (TRUST4) reads whole-blood RNA-seq and pulls out receptor
sequences: **TRB** (T-cell receptor beta chain) and **IGH** (antibody heavy chain). Each person
ends up with a list of a few hundred to a few thousand distinct receptor sequences, each one
built from a **V gene** and a **J gene** (two of the several dozen "building block" genes the
immune system shuffles together) plus a short variable loop in the middle (**CDR3**) that does
most of the actual target-recognition.

**The recipe (BenchRep-T)**, a published method we're replicating:
1. Count how often each V gene and each J gene is used, per person — this is "V/J usage."
2. Chop the CDR3 loop into overlapping 4-letter chunks ("**4-mers**," where the letters are amino
   acids) and count how often each chunk appears.
3. Feed all those counts into a **logistic regression** — the simplest kind of classifier, which
   is just "a weighted sum of the counts, squashed into a probability." **L1** means it's
   penalized for using too many features at once, so it tends to pick out a small, interpretable
   set of genes/chunks that matter, instead of using everything a little bit.
4. Test it with **cross-validation (CV)**: split people into 3 groups, train on 2, test on the
   held-out 1, rotate, repeat. This avoids the classifier just memorizing the training data.
5. Score it with **AUROC** ("area under the ROC curve"): **0.5 = coin flip, no better than
   guessing; 1.0 = perfect.** In practice, 0.55-0.65 is "weak but real," 0.7-0.85 is "clearly
   predictable," and disease classifiers in this kind of setting rarely beat 0.85-0.9.

**What the numbers mean, run by run:**
- **Ancestry from repertoire alone: AUROC up to 0.83 (AFR).** This is expected and it is
  important — not because we care about predicting ancestry, but because it's a **warning sign**.
  V and J gene usage differs by ancestry for genuine population-genetic reasons (some gene
  versions are more common in some ancestries). That means **any future "repertoire predicts
  disease X" result could just be "repertoire predicts ancestry, and ancestry correlates with
  disease X (e.g. through access to care, environment, or genuine ancestry-linked biology)."**
  The rule going forward: never trust a repertoire-only disease result unless it's also checked
  against a model that already includes ancestry as a covariate — if adding ancestry kills the
  signal, the "disease" signal was actually ancestry.
- **Sex from repertoire: AUROC ~0.55**, barely above chance. Real but weak.
- **No disease could be tested yet.** All 9 chronic/autoimmune conditions we tried (rheumatoid
  arthritis, lupus, HIV, etc.) are rare enough that even in a 1,500-person sample, none had 100
  people with the disease — the minimum needed to get a statistically trustworthy answer. This is
  a **power problem** (not enough sick people in the sample), not a sign the method doesn't work.
  The fix is scaling up to the full ~7,640-person cohort, which is still pending (it was queued
  behind the month-long KIR full-cohort run).

**Bottom line:** the pipeline works (it correctly found the ancestry and sex signals we expected
to find), but we haven't yet been able to ask the actual disease question at a big enough sample
size to get an answer either way.

---

## Part B — Ranked experiment ideas

Feasibility is scored against **data already on the VM**: HLA calls (`hla_calls_rich.tsv`), KIR
calls (per-person GTFs, 12,261 people), ancestry (`cohort_membership.tsv`), sequencing platform,
and the relatedness table. "VM cost" is a rough wall-clock/dollar guess for a 4-80 vCPU instance.

| Rank | Idea | Feasibility | Novelty | Relevance to Cole's paper | VM cost |
|---|---|---|---|---|---|
| **1** | **KIR-HLA ligand co-occurrence by ancestry** (script 48) | High — both loci already called per person, ancestry already assigned | Medium-high — long-read cohorts that call both loci per person at this scale are rare | High — direct "what does long-read give you that arrays don't" story for the paper | Low (~minutes; the per-person GTF parse dominates, already amortized by 43) |
| **2** | **HLA-C1/C2, Bw4/Bw6 epitope frequencies per ancestry** | High — cheap byproduct of #1, uses only HLA calls | Low on its own (standard population genetics) but feeds #1 directly | Medium — a clean, citable table for the methods/supplement | Low (seconds once HLA calls are loaded) |
| **3** | **Naive ML confound/QC battery** (script 47: ancestry, platform, cA/cB from HLA) | High — carriage matrices are simple to build from data already extracted | Low as science, high as **QC** — platform predictability would catch a batch effect before it contaminates a real result | High — this is exactly the due-diligence Cole/reviewers will ask "did you check for this?" about | Low (minutes; L1-LR + depth-3 tree on <=2,000 features is fast) |
| 4 | KIR novel-allele recurrence vs HLA (private vs shared) | High — comes directly out of WS-A (44/45), no new extraction needed | Medium — completes the "is KIR novelty real or private noise" question S03 left open | High — a headline comparison ("58.9% of KIR calls are novel, mostly non-coding — how much of that recurs across unrelated people vs HLA's much lower novelty rate?") | Already covered by WS-A's run |
| 5 | KIR2DL2/KIR2DL3 co-occurrence per-haplotype deep-dive (0.8% vs expected ~0) | Medium — needs manual per-haplotype GTF inspection to separate miscall / duplication / phasing error, not just aggregate counting | Medium — resolves a real S03 open question | Medium — a QC footnote more than a headline result | Medium (targeted, not full-cohort — a few hours of one person's attention plus a short VM read pass) |
| 6 | KIR framework-gene misses: real deletion vs miscall | Medium — the aggregate classifier (flanked/edge_fragmented/ambiguous) already exists in 41/43; going further needs per-haplotype coordinate inspection | Medium | Medium — assembly-completeness caveat for the paper's limitations section | Medium |
| 7 | WS6 full-cohort run (7,640 people + XGBoost) + HLA-carriage positive control | High feasibility once the platform-column issue is fixed, but **long-running** and was explicitly deprioritized behind KIR in S03 | High — this is the actual disease-prediction question, unanswered so far | High — but only once it can actually test a disease at n>=100 cases | **High** (the full run is the multi-hour/GPU-adjacent cost item on this list; the earlier one this deprioritized was queued behind a 31h KIR job) |

### Top 3 recommended, with rationale

1. **KIR-HLA ligand co-occurrence (idea #1, script 48 delivered this sprint).** Cheapest
   incremental VM cost of anything on this list (both inputs already exist), directly answerable
   with textbook population-genetics statistics (odds ratio + Fisher/permutation), and it's a
   distinctive result for the paper: showing genuine trans-chromosomal receptor-ligand structure
   from single-person long reads is something short-read genotyping+separate KIR typing studies
   have to do by imputation, not directly.
2. **Epitope frequencies per ancestry (idea #2)** — essentially free once #1's code runs, and
   gives Cole a clean supplementary table independent of whether the co-occurrence result itself
   pans out.
3. **Naive ML/QC battery (idea #3, script 47 delivered this sprint)** — not a science result, but
   the platform-predictability check in particular is the kind of check a reviewer will ask for
   before trusting any carriage-based association in the paper, and it is cheap enough to run
   every time the pipeline changes.

Ideas #4-6 are worth doing but are extensions of already-scoped work (WS-A) or need more manual
attention per unit of result than #1-3; #7 (WS6 full cohort) is the highest-value single item on
the list but is also the most expensive and was already explicitly deferred once — recommend
scheduling it as its own dedicated VM session rather than folding it into this sprint.

---

## Part C — Naive ML / simple tests: what each one reveals

For every test: **the "expected result" listed is what NOT finding it would mean**, and vice
versa — read the section, not just the headline number.

| Test | What it reveals | Expected result | How to read it |
|---|---|---|---|
| **HLA/KIR allele-carriage -> ancestry** (one-feature-per-allele LR or small tree) | Which specific alleles are ancestry-informative; a confound audit for every future carriage-association analysis | Strongly predictable (AUROC > 0.7-0.8 for several ancestries) — HLA is *the* textbook ancestry-informative locus | If NOT predictable, something is wrong with the ancestry labels or the carriage matrix, not a surprising biology finding |
| **PCA of HLA+KIR carriage, colored by ancestry and by platform** | Whether ancestry structure dominates the carriage data (expected) and whether sequencing platform also forms visible clusters (a batch-effect red flag) | Clear ancestry separation; **no** visible platform separation | Platform clustering here would mean the calls partly reflect *how* someone was sequenced, not just their genotype — would need investigating before trusting any cross-platform comparison |
| **Predict sequencing platform / assembly quality from the calls** | Same batch-effect question as PCA, but as a direct classifier instead of eyeballing a plot | AUROC near 0.5 (platform should be unpredictable from genuine HLA/KIR biology) | A high AUROC here is a bigger red flag than the PCA version — it's a quantified, statistically testable batch effect, not just a visual impression |
| **HLA -> KIR cA/cB** | Whether KIR gene-content (the "cA/cB" haplotype-family label) can be predicted from a person's unrelated-by-inheritance HLA genotype | Near chance (AUROC ~0.5) — HLA (chr6) and KIR (chr19) don't co-segregate | A clearly-above-chance result would be surprising and worth a second look — most likely explanation would be a shared ancestry confound (both loci vary by ancestry) rather than direct genetic linkage |
| **Nearest-neighbour rare-allele sharing** | Sanity check of the relatedness/unrelated-set filter: does a person's closest genetic neighbour (by overall allele-sharing) disproportionately share their rare/novel alleles? | Yes for people who *should* be related but weren't caught by the kinship filter (an early-warning check); background rate near-zero otherwise | A high sharing rate concentrated in a few pairs flags an incomplete relatedness filter, not a real recurrent-allele signal |
| **Permutation baselines (all of the above)** | Whether an observed AUROC/OR could plausibly arise from chance/label-shuffling alone | The real value sits clearly above the shuffled-null distribution's 95th percentile | If it doesn't, the raw AUROC number is not trustworthy even if it looks impressively far from 0.5 — small-sample class imbalance can inflate AUROC on its own |

Scripts 47 and 48 (this sprint) implement the first four rows plus the KIR-ligand version of the
"co-occurrence" idea and permutation baselines throughout. The nearest-neighbour rare-allele
sharing check is not yet implemented — flagged as a good next small script (re-uses the
relatedness table and 03/24's rare-allele identification, no new data extraction needed).

---

## Results (S04 WS-D VM run, 2026-09-25/26 -- supersedes "What's pending" below where they conflict)

Both scripts ran on the full cohort (11,845 unrelated people); figures added (`47b`, `48b`, see
`FIGURES_INDEX.md`). **47:** ancestry is strongly encoded in HLA+KIR carriage (AUROC 0.81-0.98,
6/6 ancestries), platform only weakly (0.573, both revio and sequel2e -- expected near-chance),
and HLA-derived KIR cB content not at all beyond a modest lift (0.550). All permutation p-values
sit at the 20-shuffle floor (~0.048) for every task including near-chance ones -- p is
uninformative here; read the observed AUROC against the null band directly (that's what
`fig_naive_ml_auroc.png` panel a does). No ancestry-adjusted or within-EUR platform number exists
in the current TSVs, so whether platform's 0.573 is itself ancestry-confounded is unresolved --
flagged below as a follow-up, not answered. **48:** no KIR-HLA receptor-ligand pair shows an
ancestry-consistent enrichment/depletion across the 30 (pair x ancestry) tests; the two nominal
single-ancestry signals (EAS 2DL2xC1 OR 0.62 depleted, EUR 2DL3xC1 OR 1.43 enriched) don't survive
Bonferroni (0.05/30). EAS has the most distinct epitope profile (lowest C2 38.0%, highest C1
89.5%). Both are QC/leads, not findings — see each script's README for full numbers and caveats.

---

## What's pending

*(Superseded by "Results" above for the two items it resolves — left here, corrected, rather than
deleted, since the "not yet run" and "needs a cross-check" wording below was written before the
2026-09-25/26 VM run and is stale as originally worded.)*

- ~~Both scripts (47, 48) are local-only and synthetic-tested so far — not yet run on real AoU
  data.~~ **Done** — both ran on the full cohort (11,845 people) 2026-09-25/26; see "Results"
  above and each script's own README.
- `--platform-col` (script 47) verification against `cohort_membership.tsv`'s actual columns —
  resolved during the VM run (no crash, ran cleanly); the column-name mismatch this note warns
  about hit the sibling WS6 script, not 47/48.
- ~~The HLA-C1/C2 and Bw4/Bw6 epitope-group lookup tables in script 48 ... need a VM-side
  cross-check against IPD-IMGT/HLA's own group assignments before the co-occurrence numbers are
  treated as final.~~ **Superseded**: script 48 was rewritten (commit `d84348a` and after) to
  derive C1/C2/Bw4/Bw6 directly from each allele's own IPD-IMGT/HLA CDSseq reference sequence
  (translated, leader-stripped, residues 77-83 read directly) as the *primary* path, with the
  lookup table demoted to a fallback + cross-check only (used for novel/unresolvable calls). The
  cross-check this note asked for is now a standing output
  (`ligand_seq_vs_lookup_crosscheck.tsv`, `ligand_lookup_qc.tsv`): 235 two-field groups compared,
  228 agree, 7 disagree — see the 48 README's Method/Result sections for the full framing.
- WS6 full-cohort run (idea #7) is still pending, as documented in `42_repertoire_baseline.py`'s
  own README.
