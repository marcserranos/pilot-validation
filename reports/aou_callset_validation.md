# AoU-Native HLA Callset — Exploratory Validation Report

> **Status: DRAFT.** Sections 1–2 populated from Experiment A (2026-07-10). Section 3 pending
> `experiment_a2_callset_quality.py` output. For supervisors. Aggregate stats only — no
> participant-level data.

**In one line:** AoU ships a pre-computed, near-complete, direct-calling HLA callset for
~535k short-read participants. It is exhaustive and internally consistent; its main structural
limits are a 3-field resolution ceiling and a heavy ancestry skew. Population-level reliability
checks [pending] to confirm the calls behave like real HLA data.

---

## 1. Origin — where these calls come from

| | |
|---|---|
| **Dataset** | `hla_genotypes.tsv`, v9 CDR (`C2025Q4R6`), Controlled Tier |
| **Input data** | Short-read WGS (srWGS) CRAMs |
| **Calling method** | Ensemble of **HLA-HD + Polysolver + OptiType** (v9 org PDF, p.17 — confirmed primary source) |
| **Output** | One row per participant; `[Gene]_1` / `[Gene]_2` columns for **30 genes** (8 classical + 22 non-classical/pseudogenes) |
| **Resolution** | 2–3 field (variable by evidence) |
| **No long-read equivalent** | AoU provides no lrWGS HLA callset — this is short-read only |

**AoU's own QC methodology — found (2026-07-18).** Now titled *"All of Us Genomics &
Multi-omics Quality Report"* for CDR v9 (`C2025Q4R6`, released to the Workbench 2026-06-26),
[publicly accessible, no Workbench login required](https://support.researchallofus.org/hc/en-us/articles/50655639562900-All-of-Us-Genomics-Multi-omics-Quality-Report).
145 pages. Ground truth: **8 GIAB/NIST reference samples** (HG-001/NA12878 ×2, HG-002/3/4
Ashkenazi trio, HG-005 Han/Chinese ancestry). Reported srWGS SNP/indel accuracy (Table 8, p.27-28,
against GIAB v4.2.1 high-confidence regions): **SNV sensitivity 0.984–0.9865 / precision
0.9991–0.9997; indel sensitivity 0.9708–0.9904 / precision 0.9961–0.9988**, across all 8 samples.

**Important negative finding: HLA is not locus-specifically benchmarked anywhere in this
report.** HLA calls are named only as an existing auxiliary data product; the "challenging
medically relevant genes" QC section (p.71-72) covers KCNE1/CBS/MAP2K3, not HLA/chr6 MHC. AoU's
own general-genome accuracy numbers above give confidence in the *underlying variant calls*, but
say nothing about the HLA-calling ensemble's own accuracy at this locus — **our n=60 3-way pilot
remains the only HLA-specific reliability evidence that exists for this callset**, AoU-produced
or otherwise. Worth stating explicitly to supervisors: this isn't a gap in our homework, it's a
real gap in the field's documentation of this specific dataset.

## 2. Exhaustiveness — coverage & completeness

- **535,658 individuals** with HLA calls (100% joinable to genetic-ancestry labels).
- **Missingness ~0.0% at every one of the 8 classical genes**; 99.998% of individuals fully
  resolved across all 8 loci. Coverage is as complete as data gets — not "assumed near-full,"
  measured.
- **Resolution ceiling: capped at 3-field, universally — zero 4-field calls across the entire
  cohort.** A real structural limit: if any downstream aim needs 4-field, this callset cannot
  supply it (our SpecHLA/SpecImmune tools would be required for that).

**Cohort composition by genetic ancestry** (heavily EUR-weighted — relevant to the diversity aims):

| Ancestry | n | % |
|---|---|---|
| EUR | 302,712 | 56.5% |
| AMR | 107,928 | 20.1% |
| AFR | 100,798 | 18.8% |
| EAS | 15,893 | 3.0% |
| SAS | 6,176 | 1.2% |
| MID | 2,151 | 0.4% |
| **Total** | **535,658** | 100% |

## 3. Quality / reliability — does it behave like real HLA data?  *(PENDING A2 output)*

Approach: we cannot per-sample-truth-check 535k people, so we test properties real
population-scale HLA data must satisfy. To be filled from `experiment_a2_callset_quality.py`:

- **Allelic diversity** per locus (distinct 2-field alleles) — confirms hyperpolymorphism, not
  calling collapse. *[table pending]*
- **External frequency concordance** — top alleles vs. published population frequencies (Allele
  Frequency Net Database). *Anchor check: does EUR A\*02:01 land near the known ~28–30%?* This is
  the single strongest reliability signal available without per-sample truth. *[pending]*
- **Hardy-Weinberg proxy** — observed vs. expected homozygosity per locus per ancestry. Excess
  homozygosity flags allele dropout / reference bias. *[table pending]*
- **Rare-allele burden** — share of rare (<1%) alleles per locus. *[pending]*

**Already-known reliability evidence (from our n=60 3-way pilot):**
- It's **per-locus trust**, not one verdict. AoU-native's discordances are overwhelmingly
  *near-misses* (right allele family), not wrong-family errors.
- **DQA1**: AoU shows an ancestry-correlated discordance, but every off-call stays in the right
  family — bounded imprecision, not unreliability.
- **DRB1**: AoU-native is *safer* than our own SpecHLA short-read calls here (SpecHLA produces
  real wrong-family errors at DRB1; AoU does not).

## 4. Bottom line & recommended next steps

- The callset is **exhaustive and, per our pilot, per-locus reliable** — a genuinely valuable,
  free resource for the frequency (Aim 1) work and as a validation baseline.
- **Two limits to design around:** the 3-field ceiling, and the EUR skew (EAS/MID/SAS are thin).
- **Next steps:** (1) retrieve the v9 Genomic Quality Report for AoU's own sensitivity/precision
  numbers; (2) run the external frequency-concordance check; (3) decide the strategic fork —
  build downstream directly on this callset where it's reliable, calling ourselves only where it
  isn't (e.g. DRB1 4-field needs).
