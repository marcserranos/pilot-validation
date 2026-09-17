# S02 — what came out of the supervisor call, in plain language

*For Marc. Written 2026-09-17, end of sprint S02. Everything here is either committed under
`reports/hla_popgen/` or explicitly marked as still running.*

Read this first; the sprint board (`SPRINT.md`) has the task-level detail and the workstream
briefs have the methods.

---

## 0. The thing to tell Cole and David first

**The numbers you presented in the call are out of date, because S01 finished after the meeting.**
Not slightly out of date — the headline changed. On the record in that call, and now superseded:

| what was said | what is actually true |
|---|---|
| "we found about 3,000 novel alleles" | 280,695 novel-flagged *calls* collapse to **1,404 distinct novel proteins**, of which **231** are seen in two or more unrelated people |
| "292 seen twice, 904 seen three or more" | those counts came from the CDS-hash clustering, which pooled genuinely different sequences — 708 old clusters turned out to be 15,112 distinct sequences |
| "5–20% of the allele space discovered" | that was a *richness* estimate. The probability that the next haplotype we type is already catalogued is **99.0–99.8%** per ancestry (97.1% for Middle Eastern). Both statements are true at once; the honest headline is coverage |
| "95% concordance between relatives" | **87.0%** byte-identical, decomposed: 6.5% gene called in only one of the pair, 3.2% dropout, 2.3% single-base error, 0.8% unexplained. Twin/duplicate pairs give an error rate of **Q34** |
| "zero phase switches in 545/545 pairs" | still zero — but now over **3,021 testable transitions**, with **100%** detection power on synthetic switches we injected. The old version of that test could not have failed |

When Cole said in the call, reassuringly, *"Immuannot only looks at the CDS and then clusters
around that, which is what we do"* — that is precisely the artifact. It is why the counts above
collapsed.

This matters beyond bookkeeping: Cole asked for a list of common novel alleles to name in main
text (below). That list is built on the corrected numbers, and the old headline would not have
survived review.

---

## 1. Cole's main-text list — and it points somewhere unexpected

He asked for *"a really great list of those, because those are things we can literally call out in
main text. Like, we found these seven examples."*

There are **38**. They are clean (no artifact flag, not homopolymer-only), recurrent, and carried
by at least 20 unrelated people. That threshold is also exactly where All of Us lets us publish a
count, so the disclosure rule and "worth naming" happen to select the same set. 29 are new
**proteins**; 9 are new synonymous coding sequences and are listed separately so nobody counts
them as proteins.

**Not one of the 38 is in a classical HLA gene.**

| gene | novel proteins in the list |
|---|---|
| TAP2 | 10 |
| TAP1 | 9 |
| HLA-DQB2 | 4 |
| MICB | 2 |
| HLA-DMA, HLA-DMB, HLA-DQA2, HLA-F | 1 each |

The biggest one is a **TAP1** protein a single residue away from TAP1\*01:01, carried by **426**
unrelated people — about **139 per 1,000** African-ancestry participants. 15 of the 29 are most
common in African-ancestry participants.

### The obvious objection, answered

*"TAP1 tops the list because nobody ever catalogued TAP1, not because it is diverse."* That is
measurable, and script 27 already measured it independently. The share of a gene's haplotypes
whose protein is **already** absent from IPD-IMGT:

- TAP1 **8.0%**, TAP2 **6.9%**
- HLA-A/B/C **0.1–0.2%**

The genes with the most novel alleles are the genes with the biggest catalogue gap. Two
independent estimates agree rather than one confounding the other. Every row of the callout table
carries this number so a reviewer sees it without asking.

**What this means for the paper.** S01 suggested reframing around the non-classical MHC. This is
the strongest evidence yet for that: the classical genes are essentially catalogued, and the real
discovery space — the part only long reads reach — is TAP, MIC, DM/DO and the DQ/DP paralogues.

`reports/hla_popgen/32_novel_callouts/` — `main_text_candidates.tsv` is the table, plus 309 more
in the IMGT submission queue that are real but too rare to name.

---

## 2. Figure 1 — two of the four panels are rebuilt

Cole specified: (a) admixture, (b) a class I ternary, (c) novel rate by ancestry broken out by
gene, (d) the allele-frequency spectrum with novel and known as two colours in one plot.

**Panels c and d are done** (`reports/hla_popgen/33_figure1_v3/`). Panel c is the one that changed
most: instead of one lumped "novel" rate per ancestry, it splits by what the novelty actually is —
new protein, new synonymous CDS, new non-coding — and draws the **artifact rate underneath as a
flat negative control**.

That control largely holds, stated precisely. Across the five well-powered ancestries the
artifact rate sits between **0.90% and 3.22%** (lower bound), with an upper bound never above
5.6%. Clean novelty at DRB1 runs from **44.7% in European-ancestry** haplotypes to **69.6% in East
Asian**; DPB1 from 11.4% to 26.3%. A gradient that steep against a control that flat is reference
incompleteness, not uneven assembly quality.

**The honest exception is MID.** Middle Eastern is the smallest strict-ancestry group and almost
all of its artifact cells are suppressed, so its artifact rate is bounded only between **0% and
6.5–13.3%** depending on the gene — effectively unconstrained. MID cannot carry the control
argument, and I originally wrote "1.1–3.3% in every gene and every ancestry", which was wrong on
both counts. It is the same suppressed-count error described below, made four paragraphs before I
described it.

**Panels a and b need a VM rerun** — scripts 06 and 10 committed the figures but never the
underlying frequency tables, so they cannot be recomposed offline. They need regenerating anyway
at the stricter admixture threshold Cole and David both asked for. `WS4_figure1_v3.md` says
exactly what that rerun needs.

### A mistake worth knowing about

The first version of that panel reported a **0.00% artifact rate for Middle Eastern ancestry
across all eight genes**. Flat, clean, entirely fabricated. The cause: 3,675 of the 8,784 cells in
the source table are written `<20` under the disclosure rule, and I was reading them as zero. MID
is the smallest group, so nearly all its cells are suppressed, and the whole row collapsed to
zero while looking perfectly plausible.

Every rate in that panel is now an interval — suppressed cells bounded at [0, 19], with the
numerator and denominator bounds paired so the interval genuinely brackets the truth — and there
is a test that fails if anyone reintroduces the shortcut. I am flagging it because the same class
of error is easy to make anywhere in this repo: **a suppressed count is not a zero.**

---

## 3. Selection and the peptide-binding groove

Cole asked us to find citable methodology, and specifically to compare diversity in the
peptide-binding residues against the rest of the coding sequence.

The literature review is in `WS_literature_selection.md`. The useful citations are Hughes & Nei
(1988 *Nature*, 1989 *PNAS*) for the ARS dN/dS approach, Solberg et al. 2008 for the
Ewens–Watterson homozygosity test across loci, and **Brandt et al. 2018 (*G3*)** for why raw Fst
at HLA misleads.

**The honest finding from that review: the canonical peptide-binding residue tables could not be
verified.** Parham 1988's PDF returned 403, Bondinas 2007 is paywalled, Hughes & Nei's table
predates stable IMGT numbering. Rather than compute a headline number from a list nobody here had
actually read, we did two things:

1. **The primary comparison uses IMGT's own exon annotation** — exons 2+3 for class I, exon 2 for
   class II. Coarser, but it comes from data we ship and cannot be silently off by a
   signal-peptide length.
2. **We derived peptide-contact residues ourselves from crystal structures**: every MHC residue
   with an atom within 4.5 Å of a bound peptide, across six PDB entries. 207 residues,
   `reference/ars_peptide_contacts.tsv`, reproducible and citable without depending on a table we
   could not read.

That second step is arguably better than what we were going to cite. It also turned up two things
worth recording: one structure (1DLH) has its peptide chain mislabelled in its own metadata, and
another (3LQZ) is a single-chain construct with no separate peptide chain at all. Both were caught
by checking the structures rather than trusting the headers.

---

## 4. The three analyses that needed the VM — all finished

### 4a. LD between DQ and DP alleles, within ancestry (Cole's A3)

He predicted the linkage patterns would differ between ancestries beyond just frequencies. **They
do**, and the effect survives the fair comparison.

The trap he walked us toward: r² is bounded by allele frequencies, and every contingency-table
statistic inflates when there are more alleles relative to sample size. Ancestries differ in both.
So the comparison is made only on haplotype counts **rarefied to a common n** (~560 phased
haplotypes per ancestry), with an interval over subsamples.

First, the controls. HLA-B~HLA-C, 90 kb apart and known to be in strong LD, comes out at
D′ 0.87–0.92. HLA-A~HLA-B, 1.4 Mb apart, comes out lowest at 0.58–0.74. The method finds strong
linkage where strong linkage is known to exist and weak linkage where it is not.

Then the finding, at equal sample size:

| pair | AFR | EUR |
|---|---|---|
| DQA1~DQB1 | **0.908** [0.890–0.925] | **0.974** [0.964–0.983] |
| DRB1~DQB1 | **0.880** [0.860–0.899] | **0.947** [0.932–0.960] |

Non-overlapping intervals. **African-ancestry haplotypes carry measurably weaker class II linkage
than European-ancestry haplotypes** — more distinct DQ and DR–DQ haplotypes, consistent with older
effective population size and more accumulated recombination. That is a real population-genetics
result, not a frequency artifact.

**Scope it to that contrast, though.** The non-overlap test establishes AFR vs EUR and nothing
wider: for DRB1~DQB1, AFR's interval [0.860–0.899] just touches EAS's [0.898–0.933]. "African
ancestry is lowest of all six" is not established; "African ancestry is lower than European" is.

One caveat I want on the record: the rarefaction target is always the smallest ancestry, which
here is Middle Eastern. That means MID is "subsampled" to its own full size, so its interval
collapses to a single point — that is *no* variance estimate, not a precise one, and MID must not
be compared as though it had a tight interval. The AFR-vs-EUR comparison above is unaffected;
both have proper intervals over genuine subsamples.

**DPA1~DPB1 does not follow the pattern** — there MID/SAS/EAS are lowest and AMR/EUR highest, with
AFR in the middle. That is worth flagging rather than smoothing over, because Brandt et al. 2018
(*G3*) independently singled out the DP locus as the exception among HLA genes, arguing it is
under directional rather than balancing selection. Our data behave the same way.

### 4b. The phasing question Cole called critical (A2)

He said: *"we don't want there to be a switch between DPA1 and DPB1."* Direct answer:

- **DPA1 and DPB1 sit on the same assembled contig in 96.96%** of assemblies that carry both
  genes; DQA1–DQB1 in 96.31%. In those, the pairing is physical, not statistical.
- S01's switch test, rebuilt so it could fail, found **zero switches over 3,021 testable
  transitions at 100% detection power** on injected synthetic switches.
- HLA-A~HLA-B, 1.4 Mb apart, is phased in only **30.2%** of assemblies — the honest limit of what
  these assemblies span, and itself a number worth reporting.

The strong observed DP and DQ linkage is a third, independent line of evidence: if the phasing
were switching, the LD would be destroyed.

### 4c. Deletions and duplications (A4), and KIR (A5)

Cole said *"let's do that, it could be pretty crazy."* The method calls a gene deleted only when a
single contig carries genes on **both** sides of it, so assembly fragmentation cannot masquerade
as biology.

**The positive control passes.** Which second DRB locus a haplotype carries is determined by its
DRB1 group — textbook immunogenetics, known independently of our data. We recover it in
**92–98%** of haplotypes across all 13 DRB1 groups, including the groups whose correct answer is
"no second DRB locus at all".

**The negative control gives the false-positive rate:** 0.04–2.0% at HLA-A/B/C, DRA, DQ, DP and
DRB1.

Against that baseline:

| gene | haplotypes lacking it |
|---|---|
| DRB5 | 83.1% |
| DRB4 | 69.6% |
| DRB3 | 49.3% |
| C4B | 19.6% |
| C4A | 11.0% |
| MICA | 4.2% |
| MICB | 2.9% |

**One lead worth chasing, flagged as a lead and not a result.** Four class I pseudogenes —
HLA-H (12.15%), HLA-K (12.42%), HLA-T (12.25%) and HLA-U (12.37%) — come out with near-identical
deletion rates, while every other pseudogene sits under 1.4%. Four independent genes agreeing to
within 0.3 percentage points is not what independent deletions look like. Either they sit in one
segment that is deleted as a block on about 12% of haplotypes — which would be a real
haplotype-level structural polymorphism and a nice finding — or they share a detection artifact.
Deciding between those needs checking whether the same haplotypes are missing all four, which is
a short follow-up and is not done.

HLA-Y is absent from 81.6% of haplotypes, which is expected: it is a pseudogene known to be
present only on a subset of haplotypes.

C4 copy number resolves cleanly per haplotype: 54.1% carry one C4A and one C4B, 14.8% carry C4A
only, 10.5% C4B only. This is genuine copy-number variation phased on individual haplotypes, which
is exactly the thing short reads cannot do well.

**KIR: zero calls, and that is correct.** The KIR cluster is on chromosome 19, outside the
chromosome 6 window these assemblies were trimmed to. Cole expected it to "just be called" — it
cannot be, and typing it would need a separate extraction from the original BAMs. Worth saying
plainly, because it is an opportunity rather than an oversight.

### 4d. Amino-acid diversity, groove vs the rest (A8), and differentiation (A10)

Diversity is concentrated in the peptide-binding groove exons, significantly, in eight genes:

| gene | groove / non-groove | p |
|---|---|---|
| DPB1 | **8.0×** | 0.0005 |
| DRB3 | 6.3× | 0.001 |
| DRB5 | 6.2× | 0.0015 |
| DQA1 | 4.3× | 0.0005 |
| DRB1 | 3.7× | 0.0005 |
| HLA-B | 3.6× | 0.0005 |
| DQB1 | 3.0× | 0.0005 |
| HLA-A | 2.2× | 0.001 |

**Two of the three conserved controls show nothing**, which is what makes the rest credible:
DRA 0.0× (p = 0.93) and HLA-F 0.04× (p = 0.86). HLA-C (1.3×) and DPA1 (2.0×) are not significant
either, and are reported as such rather than rounded up.

**The third control, HLA-E, does not behave — and I initially hid that.** My first draft listed
"DRA, HLA-F, DRB4" as the controls, quietly substituting DRB4 (which is not a designated control)
for HLA-E (which is). HLA-E's groove/non-groove ratio is **20.0×, the largest of any gene in the
table**. It is not significant (p = 0.39), and the reason is visible in the absolute numbers: HLA-E's
mean amino-acid diversity is 0.0015, so that ratio is two near-zero quantities divided by each
other and the permutation test correctly refuses to call it. But a reader is entitled to see it,
decide for themselves, and ask whether HLA-E — which presents a very restricted peptide repertoire
— is really the right negative control. Substituting the control that behaved was not defensible.

This is Hughes & Nei's 1988 result reproduced at biobank scale with controls they did not have.
It is a replication, not a discovery — but it is the replication that licenses everything else we
say about selection.

### The check that convinced me the pipeline is right — corrected after review

**I originally got this wrong, and it is worth showing how.** I wrote that "the ten most diverse
positions are" and then listed eight of them — silently dropping, in each of three genes, exactly
the two that did *not* match a famous published position. A critic agent recomputed the top ten
from the table and caught it. That is textbook cherry-picking, and the conclusion I drew from it
("about as good an end-to-end sanity check as this pipeline can give itself") was not earned.

Here are all ten, in mature-protein numbering, with nothing removed:

| gene | the full top ten by amino-acid diversity |
|---|---|
| DRB1 | 10, **11**, **13**, 30, **37**, **67**, **70**, **71**, **74**, **96** |
| DQB1 | **−4**, 26, 30, 55, **57**, **70**, **71**, **74**, 87, 125 |
| HLA-B | 24, **45**, **67**, **77**, **80**, **95**, **97**, **114**, **116**, 163 |

Bold marks positions that are canonical peptide-binding or epitope residues. The honest summary is
that **8 of 10 for DRB1, 6 of 10 for DQB1 and 8 of 10 for HLA-B** land on known functional
positions — still a strong result, and still one that nothing in the analysis was told about, but
not the clean sweep I first presented.

The two that don't fit are informative rather than embarrassing:

- **DQB1 position −4 is inside the signal peptide** under the 32-residue leader length I assumed.
  Either that leader length is wrong, or this is genuine signal-peptide polymorphism (DQB1 does
  have known leader variation). Until someone checks it against IMGT, it is a warning light on the
  mature-numbering conversion — which is precisely the off-by-leader-length failure the sprint
  board flagged as a risk before any of this was run. **The committed table uses our own
  reference-protein numbering throughout and needs no such assumption**; only this presentation
  layer does.
- HLA-B 163 and 24, and DRB1 10 and 30, sit in or adjacent to the groove domains but are not on
  the canonical lists. They may be real, they may be numbering drift. I am not going to
  rationalise them after the fact.

**On Cole's class I vs class II hypothesis (A10): the data do not clearly support it.** Ranked by
Hedrick's standardised G′st, the most differentiated gene between ancestries is **HLA-B (0.518)**,
a class I gene, followed by DPB1 (0.483), DRB1 (0.444) and HLA-A (0.417). Class I and class II
interleave.

There is a methodological point here worth showing him. **Raw Fst would have given a completely
different ranking** — DRB5 (0.165), DRB4 (0.146) and DPA1 (0.131) at the top — purely because those
genes have lower within-population heterozygosity, which mechanically allows a larger Fst. At
HLA-B, within-population heterozygosity is 0.954, so raw Fst cannot exceed about 0.05 no matter
how different the populations are. This is exactly the artifact Brandt et al. 2018 documented, and
we can demonstrate it in our own data rather than just cite it.

Non-classical controls sit near zero differentiation as they should: HLA-E 0.023, DRA 0.029,
HLA-F 0.049, HLA-G 0.059.

**A tension between this section and section 1, which should be stated rather than left for a
reviewer.** The diversity analysis works from proteins catalogued in IPD-IMGT, so the novel
proteins that sections 1 and 3 are entirely about are *excluded from it by construction*. For the
classical genes that exclusion is tiny (99.7–99.97% of calls are covered), but HLA-G is only
96.6% covered and HLA-G appears in the callout list. The two halves of this document are measuring
overlapping but not identical things.

**And the deletion rates in §4c are annotation-level, not sequence-level.** We show that an
assembled contig skipped a gene's position; we do not show a breakpoint. A confirmed structural
deletion would.

---

## 5. Things I could not do, and things that need you

0. **A disclosure issue in an already-committed file, found during review.** S01's
   `FINDINGS_FOR_MARC.md` §2.4 states "the **five** duplicate/identical-twin pairs" — a bare
   participant-pair count below the threshold of 20, in a report committed to a public repo. Every
   S02 output suppresses counts like that, but this one predates the discipline and should be
   fixed before the branch is pushed. It sits alongside the older open question about
   `novel_alleles.tsv`.
1. **Push the branch.** Still blocked for me — publishing to the public repo is denied by the
   permission classifier. Everything is committed locally on `fig1-drafts-and-research-map`.
2. **The VM needs a reboot.** The Workbench console says so, and it went down mid-run once
   already, killing all three jobs. It came back on its own and the runs were relaunched and
   completed, but the warning is still showing. Rebooting is a shared-resource decision I did not
   want to make without you — another session may be using it.
3. **The short-read validation Cole parked.** He said *"I wouldn't do that right now"*, but it is
   the single strongest confirmation available for the 38 callout alleles: every one of those
   people also has high-coverage short reads, and realigning them to the long-read assembly tests
   the novel protein directly. Worth raising again *because* the callout list now exists.
4. **Two decisions.** Do we reframe the paper's first section around the non-classical MHC? And do
   we lead with 231 recurrent novel proteins (strict) or 1,404 distinct ones (inclusive)? I would
   say yes and 231.

---

## 5b. What the review caught

A fresh-context critic agent re-derived every quantitative claim in this document from the
committed tables. It found **two critical and three major errors, all mine**, and all are corrected
above rather than quietly patched:

| # | what was wrong |
|---|---|
| 1 | The "ten most diverse residues" list showed eight, dropping exactly the two per gene that did not match a published epitope. Cherry-picking. |
| 2 | "Artifacts sit at 1.1–3.3% in every gene and every ancestry" — wrong range, and false for MID, whose rate is unconstrained. The same suppressed-count error the document congratulates itself for fixing. |
| 3 | DRB1 and DPB1 novelty rates quoted in prose (45.7/71.2, 11.8/27.1) came from the pre-correction run and contradicted this document's own table (44.7/69.6, 11.4/26.3). |
| 4 | The conserved-control list substituted DRB4, which behaves, for HLA-E, which does not. |
| 5 | The LD claim generalised past what the intervals support. |

Two minor transcription errors were also fixed. Everything the critic checked and found accurate —
the correction table in §0, the callout counts and per-gene tallies, the catalogue-gap figures, the
C4 and DRB control numbers, the groove ratios and p-values, the G′st ranking, the phasing yields —
it verified against the tables independently.

I am leaving this section in rather than deleting it, because the pattern matters: **every error
was in the direction of making the result look cleaner than it was.** That is the direction to
check first in anything I hand you.

## 6. Where everything lives

| what | where |
|---|---|
| the callout list | `reports/hla_popgen/32_novel_callouts/` |
| Figure 1 panels c and d | `reports/hla_popgen/33_figure1_v3/` |
| literature and citations | `sprints/S02_supervisor_items/WS_literature_selection.md` |
| peptide-contact residues | `reference/ars_peptide_contacts.tsv`, `WS5_ars_definition.md` |
| what the call asked for, item by item | `SPRINT.md` §1 |
| what went wrong and when | `JOURNAL.md` |
| LD by ancestry | `reports/hla_popgen/29_hla_ld/` |
| deletions, duplications, C4, KIR | `reports/hla_popgen/30_hla_sv/` |
| amino-acid diversity and differentiation | `reports/hla_popgen/31_aa_diversity/` |
