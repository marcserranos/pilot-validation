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

That control is the panel's whole argument, and it holds. Artifacts sit at **1.1–3.3%** in every
gene and every ancestry. Clean novelty at DRB1 runs from **45.7% in European-ancestry** haplotypes
to **71.2% in East Asian**; DPB1 from 11.8% to 27.1%. A gradient that steep against a flat control
is reference incompleteness, not uneven assembly quality.

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

## 4. What is still running

Three analyses are executing on the VM as this is written, and their results are not in this
document:

- **LD between DQ and DP alleles within each ancestry** (Cole's A3, which he raised twice). The
  script reports the r² he asked for, plus a multi-allelic D′ and a bias-corrected Cramér's V, and
  makes the cross-ancestry comparison only on haplotype counts **rarefied to a common n** — because
  r² is bounded by allele frequencies and every contingency-table statistic inflates with allele
  count, and ancestries differ in both. DRB1–DQB1 and HLA-B–HLA-C are included as positive
  controls: if the method does not find strong LD there, the method is broken.
- **Deletions and duplications** (A4). A gene is called deleted only when one contig carries genes
  on *both* sides of it — a bridged absence — so assembly fragmentation cannot masquerade as
  biology. DRB3/4/5 against the textbook DR51/52/53 expectation is the positive control; HLA-A/B/C
  are the false-positive rate. C4 copy number and long/short composition come along for free.
- **Amino-acid diversity per residue** (A8/A10), groove vs non-groove, plus between-ancestry
  differentiation per gene reported as Hedrick's standardised G′st rather than raw Fst.

**KIR (A5) is already answered**: there are **zero** KIR calls, and that is correct — the KIR
cluster is on chromosome 19, outside the chromosome 6 window these assemblies were trimmed to.
Genotyping KIR would need a separate extraction from the original BAMs. It is a real opportunity,
not an oversight, and worth saying that way to Cole since he expected it to "just be called".

---

## 5. Things I could not do, and things that need you

1. **Push the branch.** Still blocked for me — publishing to the public repo is denied by the
   permission classifier. Everything is committed locally on `fig1-drafts-and-research-map`.
2. **Reboot the VM.** It went down mid-run (the Workbench console says a reboot is required). It
   came back on its own and the runs were relaunched, but if it needs a manual reboot again, that
   is a shared-resource decision I did not want to make without you — another session may be using
   it.
3. **The short-read validation Cole parked.** He said *"I wouldn't do that right now"*, but it is
   the single strongest confirmation available for the 38 callout alleles: every one of those
   people also has high-coverage short reads, and realigning them to the long-read assembly tests
   the novel protein directly. Worth raising again *because* the callout list now exists.
4. **Two decisions.** Do we reframe the paper's first section around the non-classical MHC? And do
   we lead with 231 recurrent novel proteins (strict) or 1,404 distinct ones (inclusive)? I would
   say yes and 231.

---

## 6. Where everything lives

| what | where |
|---|---|
| the callout list | `reports/hla_popgen/32_novel_callouts/` |
| Figure 1 panels c and d | `reports/hla_popgen/33_figure1_v3/` |
| literature and citations | `sprints/S02_supervisor_items/WS_literature_selection.md` |
| peptide-contact residues | `reference/ars_peptide_contacts.tsv`, `WS5_ars_definition.md` |
| what the call asked for, item by item | `SPRINT.md` §1 |
| what went wrong and when | `JOURNAL.md` |
| LD / structural variation / amino-acid diversity | `reports/hla_popgen/29_*`, `30_*`, `31_*` once the VM run lands |
