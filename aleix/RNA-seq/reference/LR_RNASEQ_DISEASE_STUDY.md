# The LR × RNA-seq overlap cohort — disease study design

**Written 2026-08-17, no VM access.** Everything here is either derived from primary
sources / this repo's verified notes, or is explicitly flagged as a number we must measure.
Three runnable scripts accompany this document and are ready to go.

---

## 1. What we are actually doing, and why it is worth doing

Two AoU sub-cohorts matter to this project:

- **Long-read WGS (lrWGS)** — 14,521 people. Long reads resolve HLA and other structurally
  complex regions that short reads garble. This is Marc's substrate.
- **RNA-seq** — 8,980 whole-blood samples. This is the repertoire substrate. Ours.

The people in **both** are the only people for whom we can eventually have *both* a
high-confidence HLA type from long reads *and* an immune repertoire. That intersection is
the only place the two workstreams physically meet.

**So the question this study answers is not "what diseases do these people have" for its
own sake. It is: _is this intersection a viable scientific cohort at all, and for which
diseases?_** If the overlap is 300 people and the most common immune-mediated condition in
it has 11 cases, that is a decisive, early, cheap finding that reshapes the plan. If it is
2,000 people with 400 asthma cases, that is a different project.

This is a feasibility study wearing the clothes of a descriptive one. Both deliverables are
real; the feasibility read is the one that changes decisions.

---

## 2. Background you need, from zero

**EHR (Electronic Health Record)** — the clinical record generated when a person receives
care: diagnoses, medications, procedures, lab results. AoU links participants' real EHRs
from their healthcare providers. It is *not* a research questionnaire; it is billing and
clinical documentation, with all the messiness that implies.

**OMOP CDM (Observational Medical Outcomes Partnership Common Data Model)** — the standard
schema AoU uses to make records from hundreds of different hospitals comparable. Every
hospital codes diabetes differently; OMOP maps them all onto one shared vocabulary. Tables
we care about:

| Table | Contents |
|---|---|
| `person` | One row per participant: year of birth, sex at birth |
| `condition_occurrence` | One row per diagnosis event. The core of this study. |
| `observation_period` | **When we were watching.** First/last date any record exists. |
| `concept` | The dictionary. Turns numeric concept ids into names and codes. |

**SNOMED CT** — the standardized clinical vocabulary OMOP maps conditions *to*.
**ICD-10-CM** — the billing code system hospitals actually record *in*. AoU keeps both:
the harmonized SNOMED concept, and the raw source code. We use SNOMED for naming and ICD-10
for cheap chapter-level grouping (first letter of the code = disease chapter).

**Phecodes** — a research-oriented regrouping of ICD codes into ~1,800 clinically
meaningful phenotypes, designed to remove the arbitrary granularity of billing codes (ICD
has dozens of near-identical diabetes codes; phecode has one). **Better than raw ICD for
serious phenotyping, but requires an external mapping file** and >75% ICD-10 coverage in
published evaluations. **Decision: start without it.** ICD-10 chapters plus direct SNOMED
concept counts answer the feasibility question fine, with zero dependencies. Add phecodes
only if we go on to real association work.

**Prevalence** — fraction of the cohort with a condition. **The trap:** raw prevalence is
confounded by how long each person's chart covers. Someone followed for 12 years accrues
more diagnoses than someone followed for 1, at identical true health. **Every number in
this study is normalized by EHR-years, or it is wrong.** This is why
`observation_period` is pulled as a first-class table and not an afterthought.

---

## 3. How big is the overlap? — the analysis I can do without the VM

We do not know the answer, but we can bracket it, which tells us what result would be
*surprising*.

Under **independent sampling** — if AoU picked the RNA-seq people and the long-read people
without regard to each other, from the same pool:

| Pool assumed | Expected overlap |
|---|---|
| All 535,662 srWGS participants | **~243** |
| The 245,388-person CDRv7 pool the original LR cohort was drawn from | **~531** |

So **~240–530 is the "nothing special happened" range.**

**But there is a concrete reason to expect more.** The documented selection criteria for the
original long-read cohort were: had srWGS and array data, **self-reported African American**,
sufficient high-molecular-weight DNA (>5 μg), **EHR available**, unrelated, no known rare
genetic disease. Several of those — enough biospecimen volume, EHR linkage, active
participation — are exactly the criteria that would also make someone eligible for a
multi-omics blood draw. Sub-cohorts selected on correlated criteria overlap more than chance.

**And one reason to expect a specific shape rather than just a size:** the v7 long-read
cohort was **deliberately African-American-enriched** (952 of 1,027 self-identified Black or
African American). Our own live check confirms the downstream effect is still visible in v9:
of long-read people with a resolvable BAM, the `sequel2` platform group is **940/991
(94.9%) AFR**. **If that enrichment carries into the overlap, this cohort will be
ancestry-skewed toward AFR — which is scientifically valuable (it is the group most
underserved by existing repertoire literature) but must be stated in every result, not
discovered later.**

**Two different overlap numbers, and they are not interchangeable:**

| Definition | Needed for | Expected scale |
|---|---|---|
| Person is in **both manifests** | **The disease study** — we need people and their charts, not their reads | ~240–530+ |
| Person also has a **real, existing long-read BAM** | Re-calling HLA from long reads | Much smaller. Only **2,763 of 14,521** long-read people have a verified existing GRCh38 BAM, so scale this down by ~5×: plausibly **~50–110** |

**This distinction matters and is easy to get wrong.** The disease study is *not* gated on
BAM availability. Run the overlap check without `--check-bams` first; the BAM count only
constrains the later HLA-recall step.

---

## 4. The three-step pipeline (all scripts written and ready)

```
  [1] check_lr_rnaseq_overlap.py     shell, on the VM, ~1 min      (mount only)
        │  reads: RNA-seq manifest, lrWGS manifest, ancestry TSV
        │  writes: VM-local cohort TSV + de-identified summary CSV
        ▼
  [2] query_overlap_phenotypes.py    JUPYTER NOTEBOOK, ~minutes    (BigQuery)
        │  pulls: demographics, EHR observation window, all conditions
        │  writes: three VM-local TSVs (individual-level clinical data)
        ▼
  [3] analyze_disease_burden.py      shell, seconds                (pandas only)
           writes: five de-identified, small-cell-suppressed CSVs into results/
```

**Note the environment split.** Steps 1 and 3 run in the ordinary pixi shell. **Step 2 must
run in a notebook** — the CDR is only reachable via BigQuery with a fully-qualified table
path, and `ENVIRONMENT.md` quirk #5 confirms it is not exposed via shell env vars. This is
the first time this workstream has needed a notebook; previous work was deliberately
shell-only.

None of these three is a batch job. Step 1 is a manifest join, step 2 is a database query,
step 3 is arithmetic. **This satisfies this week's rule that the VM is only used when an
experiment adds to results** — strictly speaking none of these is an experiment at all.

---

## 5. What the five outputs will tell us

1. **`lr_rnaseq_cohort_demographics.csv`** — size, ancestry mix, age, EHR-years by group.
   Answers: *is this cohort AFR-skewed as predicted, and how much chart do we have per
   person?*
2. **`lr_rnaseq_burden_summary.csv`** — median distinct conditions per EHR-year, by
   ancestry. Answers: *are these people sicker than average, and evenly so?* A heavily
   sick cohort would mean the repertoires we've been assembling are drawn from people with
   active disease — which changes their interpretation.
3. **`lr_rnaseq_top_conditions.csv`** — the ranking. Most prevalent diagnoses, with a
   `mean_codes_per_affected` column that separates a chronic followed diagnosis from a
   one-visit incidental code.
4. **`lr_rnaseq_icd_chapters.csv`** — burden by broad category (circulatory, endocrine,
   mental, neoplasm...). The one-slide summary of what this cohort's medical profile is.
5. **`lr_rnaseq_immune_feasibility.csv`** — **the decisive one.** Case counts across ten
   immune-mediated disease families (autoimmune rheumatologic / endocrine / GI / neuro /
   derm, allergy, chronic infection, acute infection, lymphoid malignancy,
   immunodeficiency & transplant), each with a blunt verdict: *enough for a real
   case/control model* (≥200), *descriptive only* (≥50), or *not viable* (<50).

Output 5 is the bridge between the two tasks. Every option in `POST_TRUST4_OPTIONS.md`
Tier 4 — the supervised repertoire-classification direction the supervisor papers point at —
requires labels and cases. **This table says which labels exist in sufficient number.**

---

## 6. Honest limits, stated before we see the numbers

- **This is a description, not an association study.** At a few hundred people, we can rank
  and quantify; we cannot claim disease A is *caused by* or even *reliably associated with*
  anything. Saying so up front is cheaper than being corrected later.
- **EHR data is care-seeking data.** A diagnosis appears when someone sought care and a
  clinician coded it. Under-diagnosis in under-served populations is real and systematic,
  and this cohort is likely AFR-enriched — so **apparent prevalence differences between
  ancestry groups may reflect healthcare access, not biology.** This is the single most
  important caveat and it must appear on any slide showing an ancestry breakdown.
- **Substring matching over-matches.** The immune feasibility table catches "osteoarthritis"
  under "arthritis". It is a screen. Anything promising gets a real concept-set definition.
- **Published evaluation of ICD→phecode reproducibility: 70–75% for chronic diseases,
  <10% for acute ones.** Acute-condition counts (pneumonia, influenza, UTI) are much less
  trustworthy than chronic ones. Weight conclusions accordingly.
- **Small cells are suppressed at n<20** throughout, per AoU's disclosure floor. This
  applies to slides too.
- **We still don't know how the 8,980 RNA-seq people were selected** (data report item N3).
  If that sub-cohort was itself disease-enriched, every prevalence number here inherits it.
  This is the largest unquantified uncertainty in the study and I could not resolve it from
  public sources.

---

## 7. What I need from you

**Nothing blocking until step 1 is run.** In order:

| # | Action | Where | Cost |
|---|---|---|---|
| **T1** | `git pull`, then run `check_lr_rnaseq_overlap.py` | VM shell | ~1 min |
| **T2** | Paste the printed output (all of it — it is already de-identified) | chat | — |
| **T3** | Confirm the current CDR dataset id from a Dataset Builder snippet | notebook | ~1 min |
| **T4** | Run `query_overlap_phenotypes.py` in a **notebook** | notebook | minutes |
| **T5** | Run `analyze_disease_burden.py`, paste output | VM shell | seconds |
| **T6** | *(optional, later)* rerun step 1 with `--check-bams` | VM shell | few min |

**T3 is the one thing I genuinely cannot guess.** The dataset id changes with every CDR
release; `ENVIRONMENT.md` records `wb-silky-artichoke-2408` / `C2025Q4R6` but that was
verified in July and I would rather you read it off a live snippet than have the script
fail on a stale constant.

**Decision I need from you, not from the data:** the phenotype files in step 2 are
individual-level clinical records — the most sensitive material this workstream has
handled. The scripts keep them VM-local by construction. **Confirm that posture is what you
and your supervisors want** before running step 2, rather than after.

---

## Sources
- AoU, *How the All of Us Genomic and multi-omic data are organized*, v9 (2026), pp. 4, 40 —
  cohort sizes.
- *Population-scale Long-read Sequencing in the All of Us Research Program* —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC12622093/> — long-read selection criteria,
  cohort composition, 95.9% EHR availability.
- AoU *Genomic Research Data Quality Report* (v7) — lrWGS selection wording, QC thresholds.
- *Understanding OMOP Basics*, AoU User Support; *allofus: an R package…* —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11631081/> — OMOP table structure,
  observation-period practice.
- *Mapping ICD-10 and ICD-10-CM Codes to Phecodes* —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6911227/> — coverage and reproducibility figures.
- This repo: `../../../context/ENVIRONMENT.md` quirks #5, #13, #14.
