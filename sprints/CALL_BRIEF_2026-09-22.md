# Call brief — 2026-09-22

*What was asked, what exists, what to show, what is left. Written for a one-hour prep.
Sources: the 2026-09-17 call transcript, sprints S01 and S02.*

---

## The one-minute version

Every one of Cole's eight action items from 09-14 and all thirteen asks from the 09-17 call are
now either **done** or **partly done with the remainder scoped**. Nothing is untouched.

Three things changed the story since you last spoke to them:

1. **The novel-allele headline was wrong and is now right.** Not 3,000 alleles — **1,404 distinct
   novel proteins, 231 recurrent**. You quoted the old numbers in the call; correcting this is the
   first thing out of your mouth.
2. **The discovery is not in the classical genes.** Of 38 novel alleles common enough to name,
   **zero** are HLA-A/B/C/DR/DQ/DP. They are TAP1, TAP2, DQB2, MIC, DM/DO. This is a candidate
   reframe for the paper's first section.
3. **Cole's LD prediction was right, and his class I/class II prediction was wrong.** Both are
   properly tested now.

---

## BLOCK A — Novel allele discovery

**Asked:** how many novel alleles are real, how many are artifacts, which are common and missing
from the database ("a really great list we can call out in main text"), broken out by ancestry and
by gene. *(09-14 items 3 and 5; 09-17 asks A1c, A6.)*

**Done.**
- Rebuilt the whole novelty definition per nomenclature field — new protein / new synonymous CDS /
  new non-coding — with sequence-level truth checks against IPD-IMGT 3.55.0 (script 24).
- Found and quantified the clustering artifact: 708 old "alleles" were really **15,112 distinct
  sequences**; 21% of them differ only by homopolymer indels (script 25).
- Built the main-text list: **38 alleles with ≥20 unrelated carriers**, 29 proteins + 9 synonymous
  (script 32). Largest is a TAP1 protein in **426** unrelated people, ~139 per 1,000
  African-ancestry participants.
- Pre-answered the "TAP is just under-catalogued" objection with an independent measurement: share
  of a gene's haplotypes already absent from IPD-IMGT is **TAP1 8.0%, TAP2 6.9%, HLA-A/B/C
  0.1–0.2%**. The two rankings agree.

**Show:** `32_novel_callouts/fig_novel_callouts.png` + `main_text_candidates.tsv`,
`24_novelty_by_field/fig_novelty_funnel.png`.

**Left:** orthogonal confirmation. Cole parked short-read validation ("I wouldn't do that right
now") — **raise it again**, because the list now exists and every one of those people has
high-coverage short reads. That is the single highest-value open item in this block.

---

## BLOCK B — Quality check and phasing

**Asked:** *"We don't want there to be a switch between DPA1 and DPB1."* Estimate the phasing
error properly; Cole called this the critical question. *(09-14 item 1; 09-17 ask A2.)*

**Done.** The old test could not have failed — it excluded disagreeing genes *before* counting
switches and never reported a denominator. Rebuilt (script 26) so it can fail:
- **3,021 testable transitions**, **100%** detection power on injected synthetic switches, **zero**
  real switches.
- The 5% discordance decomposed: **87.0%** byte-identical, 6.5% gene called in one member only,
  3.2% dropout, 2.3% single-base, **0.8% unexplained**.
- Technical replicates (twins/duplicates) give the honest error rate: **Q34**, one error per 2,500
  bases.
- Physical phasing yield per interval (script 29): **DPA1–DPB1 on one contig in 96.96%** of
  assemblies carrying both; DQA1–DQB1 96.31%; HLA-A~HLA-B only **30.2%** (1.4 Mb is beyond what
  these assemblies usually span).

**Show:** `26_qc_relatives_v2/decomposition_by_gene.png`, `switch_power.png`,
`29_hla_ld/fig_phasing_yield.png`.

**Left:** nothing blocking. This block is finished and is the strongest part of the paper's
methods.

---

## BLOCK C — Population structure and Figure 1

**Asked:** four panels — (a) admixture, (b) a class I ternary "the best one you can make",
(c) novel rate by ancestry broken out by gene, (d) SFS with novel and known as two colours in one
plot. Plus a stricter admixture threshold (~98%) and a supplement figure dump. *(09-14 items 4 and
8; 09-17 asks A1, A9, A11.)*

**Done — all four panels are showable today.**
- (a) `06_figures_structure/lr_full/admixture_barcode.png` — exists, you already showed it.
- (b) `10_allele_ancestry_geometry/gene_B/ternary_B.png` — exists, and Cole said *"I quite like the
  ternary as it is."*
- (c) and (d) rebuilt at the corrected novelty definition: `33_figure1_v3/figure1_v3_panels_cd.png`.

Panel (c) is the one that improved most. Instead of one lumped rate it splits by what the novelty
*is*, with the **artifact rate drawn underneath as a flat negative control**. That control is the
argument: artifacts sit at **0.90–3.22%** across the five well-powered ancestries while clean
novelty at DRB1 runs **44.7% (EUR) → 69.6% (EAS)**. A gradient that steep against a control that
flat is reference incompleteness, not uneven assembly quality.

**Left (this is the main Figure 1 gap):**
1. Rerun scripts 06 and 10 at `--strict-threshold 0.98` (ask A9) **and commit the underlying
   frequency tables**, which they never did — that is why (a) and (b) cannot currently be restyled
   offline.
2. Regenerate the SFS at the corrected novelty definition; panel (d) currently reuses the older
   spectrum, so its *shape* is right but the absolute novel counts are over-estimates.
3. Assemble the supplement dump (per-gene SFS, the other ternaries, the UMAP/PCA panels).

*Needs the VM. Spec is written in `S02_supervisor_items/WS4_figure1_v3.md`.*

---

## BLOCK D — Selection and diversity

**Asked:** find citable methodology; compare diversity at peptide-binding vs non-binding residues
at the **amino-acid** level; is class II more differentiated between ancestries than class I?
*(09-14 item 7; 09-17 asks A7, A8, A10.)*

**Done.**
- Literature review with citations (Hughes & Nei 1988/1989; Solberg 2008; **Brandt 2018** for why
  raw Fst at HLA misleads). Honest flag: the canonical ARS residue tables **could not be verified**
  from primary sources, so we derived peptide-contact residues ourselves from six crystal
  structures — 207 residues at 4.5 Å, `reference/ars_peptide_contacts.tsv`.
- Per-residue amino-acid diversity along every protein (script 31). Groove exons enriched
  **2.2×–8.0×** in eight genes, all p ≤ 0.0015; **nothing** in the DRA and HLA-F controls.
- **Cole's class I/II hypothesis is not supported.** On Hedrick's standardised G′st the most
  differentiated gene is **HLA-B (0.518)**, then DPB1 (0.483), DRB1 (0.444), HLA-A (0.417). Class I
  and class II interleave. Raw Fst gives a completely different ranking driven by within-population
  heterozygosity — the Brandt 2018 artifact, which we can now demonstrate in our own data.

**The validation to lead with.** Ranking residues purely by diversity, with nothing told about
disease or KIR, the top positions are DRB1 β11/13/71/74 (the RA shared epitope), DQB1 **β57**
(T1D/celiac), HLA-B 77/80 (the Bw4/Bw6 KIR epitope) and 116. 8/10, 6/10 and 8/10 of each gene's
top ten land on canonical functional positions.

**Show:** `31_aa_diversity/fig_protein_track_DRB1.png`, `fig_groove_enrichment.png`,
`fig_allele_differentiation.png`.

**Left:**
1. The **transfer function is now written and verified (09-22)** — `contact_indices_for()` in
   `_ars_residues.py`, plus `reference/ars_chain_sequences.tsv`. It maps the 206 structural
   peptide-contact residues onto our own numbering by alignment, and is verified offset-invariant
   under 24/29/32-residue leaders, so **no signal-peptide length is assumed anywhere in the
   analysis**. What remains is one VM rerun of script 31 to produce the contact-vs-non-contact
   numbers; script 31 already calls the function.
2. The DQB1 −4 oddity was a *presentation-layer* artifact of the mature-numbering conversion, not
   an analysis error. Still verify leader lengths before any figure caption uses mature numbering.

---

## BLOCK E — Structural variation and KIR

**Asked:** *"Did you look at any deletions or duplications? Let's do that, it could be pretty
crazy."* And: did you look at KIR? *(09-14 item 6; 09-17 asks A4, A5.)*

**Done.** A gene is called deleted only from a **bridged absence** — one contig carrying genes on
both sides of it — so assembly fragmentation cannot masquerade as biology (script 30).
- **Positive control passes:** DRB3/4/5 presence matches the textbook DR51/52/53 expectation in
  **92–98%** of haplotypes across all 13 DRB1 groups, including the groups whose right answer is
  "no second DRB locus".
- **False-positive rate: 0.04–2.0%** at genes that are never deleted.
- Deletions: DRB5 83.1%, DRB4 69.6%, DRB3 49.3%, **C4B 19.6%, C4A 11.0%, MICA 4.2%**.
- C4 copy number phased per haplotype: 54.1% carry 1A/1B, 14.8% A-only, 10.5% B-only — the thing
  short reads do badly.
- **KIR: zero calls, and that is correct** — chr19, outside the chr6 trim window. Not an oversight;
  typing it needs a separate extraction from the original BAMs.

**Show:** `30_hla_sv/fig_deletion_rates.png`, `fig_c4_copy_number.png`, `drb_positive_control.tsv`.

**Left:**
1. **A live lead:** four class I pseudogenes — H, K, T, U — have deletion rates within 0.3 points
   of each other (12.15–12.42%) while every other pseudogene is under 1.4%. Either a **block
   deletion on ~12% of haplotypes** (a real finding) or a shared artifact. Deciding it means
   checking whether the *same* haplotypes are missing all four. Short, needs the VM.
2. KIR genotyping from the original BAMs — a separate mini-project, worth scoping with Cole.

---

## BLOCK F — Linkage disequilibrium

**Asked:** *"Compute the LD between alleles within each ancestry and see what looks different"* —
r², DQA1–DQB1 and DPA1–DPB1. Raised twice and repeated in his closing next-steps. *(09-14 item 2;
09-17 ask A3.)*

**Done** (script 29), with the trap handled: r² is bounded by allele frequencies and every
contingency-table statistic inflates with allele count, and ancestries differ in both. So the
cross-ancestry comparison is made only at a **common rarefied haplotype count**.

- **Controls behave:** B~C (90 kb apart) D′ 0.87–0.92; A~B (1.4 Mb) lowest at 0.58–0.74.
- **The result:** at equal n, African-ancestry haplotypes carry weaker class II linkage than
  European — DQA1~DQB1 **0.908 [0.890–0.925] vs 0.974 [0.964–0.983]**; DRB1~DQB1 **0.880 vs
  0.947**. Non-overlapping.
- **DPA1~DPB1 is the exception** and does not follow the pattern — which is exactly what Brandt
  2018 reports for the DP locus (directional rather than balancing selection). Independent
  agreement, worth saying out loud.

**Show:** `29_hla_ld/fig_ld_multiallelic.png`, `fig_r2_DQA1_DQB1.png`.

**Left:** scope the claim carefully — AFR vs EUR is established; "AFR lowest of all six" is not
(the AFR and EAS intervals touch at DRB1~DQB1).

---

## BLOCK G — HLA × TCR/BCR (Aleix)

**Asked:** Cole framed HLA as a **conditioning variable** for the repertoire work, not a prediction
target — which also means the disease-embedding negative is fine to report as a negative.
*(09-17 ask A13.)*

**Status: blocked on Aleix's join, by design.** 8,327 people have both long reads and RNA-seq.
Nothing is joined yet. Plan is in `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` §3.

**Left:** everything. This is the paper's second half and the biggest single unstarted block.

---

## What to show, in call order

| # | Figure | The one sentence |
|---|---|---|
| 1 | `24_novelty_by_field/fig_novelty_funnel.png` | The headline changed: 280,695 flagged calls → 1,404 novel proteins → 231 recurrent. |
| 2 | `32_novel_callouts/fig_novel_callouts.png` | Here are the 38 we can name — and none is a classical gene. |
| 3 | `33_figure1_v3/figure1_v3_panels_cd.png` | Reference bias by gene and ancestry, with the artifact rate as a flat control. |
| 4 | `10_allele_ancestry_geometry/gene_B/ternary_B.png` | The panel you already approved. |
| 5 | `06_figures_structure/lr_full/admixture_barcode.png` | Admixture, for Figure 1. |
| 6 | `29_hla_ld/fig_ld_multiallelic.png` | Your LD prediction was right, at equal sample size. |
| 7 | `30_hla_sv/fig_deletion_rates.png` | Deletions, with DR51/52/53 as the positive control. |
| 8 | `31_aa_diversity/fig_protein_track_DRB1.png` | Diversity lands on β11/13/71/74 — we never told it about disease. |
| 9 | `31_aa_diversity/fig_allele_differentiation.png` | Class II is *not* the most differentiated; HLA-B is. |

---

## Parallelisable work packages

Independent of each other — each can go to a separate agent or session.

| # | Package | Needs VM? | Size | Why it matters |
|---|---|---|---|---|
| **P1** | Rerun scripts 06 + 10 at strict 0.98, commit the frequency tables | yes | S | Unblocks Figure 1 panels a/b and answers ask A9 |
| **P2** | Regenerate the SFS at the corrected novelty definition | yes | S | Panel d's absolute counts are currently over-estimates |
| **P3** | H/K/T/U co-deletion: same haplotypes or not? | yes | S | Either a real ~12% block deletion or an artifact to retract |
| ~~P4~~ | ~~ARS contact transfer function~~ — **done 09-22**; rerunning script 31 to get the numbers still needs the VM | was: no | S | Sharpens Block D from exon-level to residue-level |
| ~~P5~~ | ~~Verify signal-peptide lengths~~ — **moot for the analysis**: transfer is alignment-based and verified offset-invariant under 24/29/32aa leaders. Still worth checking before any *mature-numbering figure caption* | no | XS | |
| **P6** | Short-read validation of the 38 callout alleles | yes | L | The orthogonal confirmation the whole novelty claim lacks |
| **P7** | Supplement figure dump | yes | M | Ask A11; mostly assembly of things that exist |
| **P8** | KIR genotyping from original BAMs | yes | L | Cole asked; currently impossible from these assemblies |
| **P9** | HLA × TCR/BCR join with Aleix | yes | XL | The paper's second half |

**Suggested parallel split if you want three agents running:** P1+P2+P7 (one agent, all Figure 1,
all VM); P3 (one agent, VM, short); the script-31 rerun that finishes Block D (VM, short — can
ride along with either). P6 and P9 need a decision from Cole first.

---

## Decisions to get from them

1. **Reframe the first section around the non-classical MHC?** The evidence says yes.
2. **Headline: 231 recurrent novel proteins, or 1,404 distinct?** Recommend 231.
3. **Restart short-read validation now?** He parked it; the list now exists, so re-ask.
4. **Is HLA-E the right negative control** for groove diversity? It has the largest raw ratio in
   the table (20×, non-significant) because its diversity is near zero.

## Housekeeping

- Branch `fig1-drafts-and-research-map` is **18 commits ahead of main and still unpushed** —
  pushing is blocked for the agent, so it needs you.
- The one disclosure issue blocking that push (a bare participant-pair count in S01's findings) is
  **fixed as of today**.
- The VM was showing "Reboot is required" on 09-17; Chrome/VM was unreachable from this session
  today, so anything marked "needs VM" above is waiting on that.
