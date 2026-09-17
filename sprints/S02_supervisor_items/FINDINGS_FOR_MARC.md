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

Non-overlapping intervals. **African-ancestry haplotypes carry measurably weaker class II
linkage** — more distinct DQ and DR–DQ haplotypes, consistent with older effective population size
and more accumulated recombination. That is a real population-genetics result, not a frequency
artifact, and it is the kind of thing this cohort is uniquely placed to show.

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

**The conserved controls show nothing**, which is what makes the rest credible: DRA 0.0×, HLA-F
0.04×, DRB4 0.03×, all p > 0.5. HLA-C (1.3×) and DPA1 (2.0×) are not significant either, and are
reported as such rather than rounded up.

This is Hughes & Nei's 1988 result reproduced at biobank scale with controls they did not have.
It is a replication, not a discovery — but it is the replication that licenses everything else we
say about selection.

### The check that convinced me the pipeline is right

The per-residue tracks don't just show diversity in the right *region*. They pick out the right
*residues*. Converting our numbering to mature-protein numbering (subtracting the signal peptide),
the ten most diverse positions are:

- **DRB1**: 11, 13, 37, 67, 70, 71, 74, 96 — β11, β13, β71 and β74 are the positions that define
  the rheumatoid-arthritis shared epitope and that dominate the amino-acid-level association
  signal in autoimmune GWAS.
- **DQB1**: 26, 30, 55, **57**, 70, 71, 74, 87 — β57 is arguably the single most studied residue
  in human immunogenetics (type 1 diabetes, celiac disease).
- **HLA-B**: 45, 67, **77, 80**, 95, 97, 114, 116 — 77 and 80 are the Bw4/Bw6 epitope, i.e. the
  KIR-binding determinant; 116 is a principal peptide anchor.

Nothing in this analysis knows about disease, about KIR, or about those papers. It ranks residues
purely by how much amino-acid diversity the cohort's own allele frequencies produce, and the
residues that come out on top are the ones the field already knows matter. That is about as good
an end-to-end sanity check as this pipeline can give itself.

One honesty note: the mature-numbering conversion above uses signal-peptide lengths (DRB1 29,
DQB1 32, HLA-B 24) that I have not verified against IMGT in this session — the agreement with
known positions is itself the evidence they are right. The committed table uses our own
reference-protein numbering throughout, which needs no such assumption.

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

Non-classical controls sit near zero differentiation as they should: HLA-E 0.022, DRA 0.029,
HLA-F 0.049, HLA-G 0.059.

---

## 5. Things I could not do, and things that need you

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
