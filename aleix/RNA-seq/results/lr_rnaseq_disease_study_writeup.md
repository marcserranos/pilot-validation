# LR × RNA-seq overlap cohort — disease study results — 2026-08-17

**Headline: the cohort is 8,327 people, ancestrally balanced, with complete EHR linkage,
and seven of ten immune-mediated disease families clear the threshold for real case/control
modelling. The repertoire direction is viable — but the cohort carries a severe
healthcare-access gradient across ancestry groups that will corrupt any naive
ancestry-stratified comparison.**

Data: `../scripts/run_disease_study.sh` (BigQuery pull + analysis), run against CDR
`wb-silky-artichoke-2408.C2025Q4R6`. Per-person clinical data is VM-local under
`~/pipeline_outputs/rnaseq/pheno/` and is not committed.

---

## 1. The cohort

| | |
|---|---|
| People with both long-read WGS and RNA-seq | **8,327** |
| Enrichment over independent-sampling expectation | **34.2×** |
| With any EHR record | **8,327 (100%)** |
| With a non-zero EHR observation window | **7,726 (92.8%)** |

**The 601-person gap matters and is easy to miss.** Every person has an
`observation_period` row, so "100% EHR coverage" is literally true — but 601 of them have a
window of *zero days*: a single date on which records begin and end. They cannot contribute
to any rate-based measure, so **the real denominator for everything below is 7,726, not
8,327.** Reporting the 100% figure without this caveat would overstate the usable cohort by
about 7%.

Ancestry composition (genetically inferred): EUR 2,475 · AMR 1,791 · AFR 1,369 ·
EAS 1,172 · SAS 1,092 · MID 427. **Five of six groups exceed 1,000 people and EUR is under
a third of the cohort** — for comparison, UK Biobank is ~94% European.

---

## 2. Disease burden — and the confound that dominates it

| Ancestry | n | Median distinct conditions | Median EHR-years | Conditions per EHR-year |
|---|---|---|---|---|
| EUR | 2,363 | **41** | **12.63** | **3.11** |
| MID | 407 | 26 | 9.59 | 2.37 |
| AMR | 1,606 | 27 | 10.52 | 2.18 |
| SAS | 1,062 | **16** | 8.28 | 1.78 |
| AFR | 1,176 | 20 | 11.69 | 1.68 |
| EAS | 1,111 | 17 | **8.54** | **1.57** |

**A 2.6× spread in recorded conditions (EUR 41 vs SAS 16) and a 2.0× spread in the
depth-normalized rate (EUR 3.11 vs EAS 1.57).**

**This is almost certainly not a 2.6× difference in how sick people are.** The two
strongest reads on it:

1. **EHR-years track the same ordering.** EUR has the longest recorded history (12.63
   years), EAS and SAS the shortest (8.54, 8.28). Longer engagement with the health system
   means more encounters, more codes, more distinct conditions — independent of health.
2. **Normalizing by EHR-years does not remove it.** If chart length were the whole story,
   the per-year rate would flatten. It doesn't — EUR stays highest at 3.11. So the residual
   is *encounter density*: how often someone sees a clinician within their covered window,
   which is a function of insurance, access, trust, and care-seeking behaviour, not
   pathology.

**This was flagged as the most important caveat in the study design before any data was
pulled, and the data confirmed it at a larger magnitude than expected. Consequence:
any ancestry-stratified disease comparison in this cohort measures healthcare access at
least as much as it measures disease.** It must appear on every slide that breaks disease
down by ancestry — not as a footnote.

**A new lead for the other workstream.** Disease burden ordering
(EUR > MID > AMR > SAS > AFR > EAS) is *not* the same as CDR3 recovery ordering
(AFR > AMR > SAS > EUR > MID > EAS) — EUR is highest in burden but fourth in recovery,
AFR the reverse. So burden is not trivially driving recovery. But EAS is lowest in **both**,
and there is a real mechanism worth testing: chronically ill people have activated immune
systems, more clonal expansion, and plausibly different repertoire recovery. **We now have
per-person disease counts and per-person CDR3 counts for the same people — this is a
directly checkable confound that did not exist before today.**

---

## 3. What this cohort is made of

Top conditions (denominator 7,726):

| Condition | n | % |
|---|---|---|
| Essential hypertension | 2,810 | 36.4 |
| Hyperlipidemia | 2,436 | 31.5 |
| Chest pain | 2,130 | 27.6 |
| Cough | 2,047 | 26.5 |
| Abdominal pain | 1,956 | 25.3 |
| Low back pain | 1,897 | 24.6 |
| Chronic pain | 1,858 | 24.0 |
| Obesity | 1,658 | 21.5 |
| Anxiety disorder | 1,443 | 18.7 |
| Type 2 diabetes (uncomplicated) | 1,307 | 16.9 |

**Note how many of the top entries are symptoms, not diseases** — chest pain, cough,
abdominal pain, back pain, dizziness. This is normal and expected for EHR data: billing
codes record *presentations*, not diagnoses. It is why the ICD-10 chapter table is led by
"Symptoms & abnormal findings" at 65.9%, and it is the reason a proper phenotype needs a
curated concept set rather than a single code.

Substantively, this is a **middle-aged-to-older, cardiometabolically loaded US clinical
population**: hypertension, hyperlipidemia, obesity, and diabetes in the top ten. Not a
healthy-volunteer cohort.

Top ICD-10 chapters by share of cohort: Symptoms & abnormal findings 65.9% ·
Musculoskeletal 57.7% · Endocrine/metabolic 54.9% · Digestive 49.3% · Genitourinary 46.5% ·
Nervous system 44.6% · Respiratory 44.4% · Circulatory 44.1%.

---

## 4. The table that decides the science

Immune-mediated disease families, denominator 7,726:

| Family | Cases | % | Distinct concepts | Verdict |
|---|---|---|---|---|
| Allergy & hypersensitivity | **2,297** | 29.7 | 92 | real case/control model |
| Acute / recurrent infection | **1,724** | 22.3 | 96 | real case/control model |
| Chronic infection | **562** | 7.3 | 46 | real case/control model |
| Autoimmune — rheumatologic | **444** | 5.7 | 43 | real case/control model |
| Immunodeficiency & immunosuppression | **344** | 4.5 | 42 | real case/control model |
| Autoimmune — endocrine | **327** | 4.2 | 33 | real case/control model |
| Autoimmune — dermatologic | **314** | 4.1 | 23 | real case/control model |
| Autoimmune — gastrointestinal | 183 | 2.4 | 41 | descriptive only |
| Lymphoid & haematologic malignancy | 151 | 2.0 | 74 | descriptive only |
| Autoimmune — neurologic | 76 | 1.0 | 11 | descriptive only |

**Seven of ten families clear the 200-case bar.** That substantially exceeds the
pre-measurement projection (which put only asthma/allergy and infection above the line).

**Discount these before believing them.** Matching is by substring against SNOMED concept
names, and it over-matches in known ways:
- "Allergy & hypersensitivity" catches every *drug hypersensitivity* code, which is
  extremely common and clinically trivial. The true asthma/atopy count is materially lower.
- "Immunodeficiency & immunosuppression" catches "transplant" in any context, including
  transplant *evaluation* and family history.
- "Autoimmune — rheumatologic" catches osteoarthritis via "arthritis", which is not
  immune-mediated and is very common in this age group.

**This is a screen, and its job was to tell us where to spend effort building real
phenotypes. It did that.** Even discounting aggressively, allergy/atopy, infection history,
and rheumatologic autoimmunity are clearly viable; the three underpowered families are
clearly not, and no amount of definitional care will rescue 76 neurologic cases.

---

## 5. What follows

**The repertoire direction is viable, on exactly the phenotypes where it should be.** The
argument for repertoire over genotype is that it records *what you were exposed to*, not
just what you inherited. The families that cleared the bar — allergy, acute infection,
chronic infection — are precisely the exposure-driven ones. That is a coherent result, not
a lucky one.

**Immediate, ordered:**
1. **Check whether disease burden predicts CDR3 recovery** (§2). New, cheap, and directly
   relevant to the open ancestry question. Data already on disk for both.
2. **Build a real concept set for one phenotype** — asthma is the obvious first, being both
   large and unambiguously immune-mediated — replacing substring matching.
3. **Scale TRUST4.** The binding constraint is no longer cohort size or confounds; it is
   that only ~154 of 8,327 people have repertoires assembled.

**Standing caveats for any presentation of this:**
- Descriptive, not association. No causal or even adjusted-association claim is supported.
- The healthcare-access gradient (§2) confounds every ancestry breakdown.
- ICD→phenotype reproducibility is 70–75% for chronic conditions and **under 10% for acute
  ones** — so the "acute / recurrent infection" row is the least trustworthy of the large
  families, despite its size.
- Small cells suppressed at n < 20, with derived percentages and medians blanked alongside
  them.

---

## Sources
- Live run, 2026-08-17: `../scripts/run_disease_study.sh` → `query_overlap_phenotypes.py`
  (BigQuery/OMOP) → `analyze_disease_burden.py`.
- Design and background: `../reference/LR_RNASEQ_DISEASE_STUDY.md`.
- ICD→phecode reproducibility figures: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6911227/>
