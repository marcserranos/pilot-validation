# 11 — W2: case/control labels for 33 chronic diseases

Script: `scripts/11_disease_labels.py` via `scripts/run_disease.sh`. Run 2026-10-12 (outputs
95f27af). Requested at the 2026-10-06 supervisor meeting: about 30 chronic conditions,
autoimmune and not classically immune, defined cleanly even at the cost of data. No figure.

## Definitions (`disease_definitions.csv`)

- **Case:** ≥ 2 recorded occurrences of the disease's ICD-10-CM or ICD-9-CM codes. A single
  occurrence is often a rule-out or billing artefact.
- **Control:**
  - no code from the disease's case or wider exclusion family (e.g. no diabetes code of any
    type for T2D);
  - EHR window ≥ 2 years;
  - at least one ICD-coded condition of any kind. A person with no ICD rows cannot be called
    disease-free.
- **Excluded:** exactly one occurrence, or an exclusion-family code.
- **Population:** `analysis_main` from report 09 (6,764 people). 5,557 have ICD rows; 92% of
  condition rows are ICD-coded.
- **Chronic only:** draw dates exist for few people, and a chronic diagnosis is permanent.

## Tiers by cases

| Tier | Rule | Diseases (cases) |
|---|---|---|
| classify + associate (W4, W5) | ≥ 150 | hyperlipidaemia 2,233; hypertension 2,166; obesity 1,186; T2D 1,059; allergic rhinitis 892; asthma 782; CAD 778; hypothyroidism 704; CKD 475; heart failure 373; AF 337; COPD 326; stroke 214; PAD 201; urticaria 188 |
| associate only (W4) | 50–149 | psoriasis 144; atopic dermatitis 143; RA 128; hepatitis B/C 128; thyrotoxicosis 127; lymphoid neoplasm 118; T1D 104; HIV 76; UC 60; SLE 59; Sjögren 58; Crohn 57 |
| describe | 20–49 | Parkinson 47; MS 43; alopecia areata 34; vitiligo 31 |
| suppressed | < 20 | coeliac, ankylosing spondylitis |

Classic autoimmune diseases are mostly under 150 cases in 6,764 people, so they can be
associated but not classified.

## Files

`summary.csv`, `disease_definitions.csv` (codes per disease), `case_counts.csv` (overall and by
ancestry; cells < 20 suppressed). Person-level labels stay on the VM
(`analysis/disease_labels.pkl`).
