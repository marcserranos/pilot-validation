# AoU small-cell / "<20" policy — what the rules actually say

**Purpose:** answer whether the project's blanket `<20 → "<20"` masking (applied today to
*any* count, including counts of alleles, not just counts of people) is required by All of
Us policy, or whether a looser rule is defensible for non-participant counts. This is
research and a proposal, not a decision — see `context/DECISIONS.md` (2026-09-28 entry) for
the open item this resolves, once resolved.

## Sources

1. **AoU Data and Statistics Dissemination Policy** (PDF, updated 5/12/2020), fetched in full:
   `researchallofus.org/wp-content/.../AoU_Policy_Data_and_Statistics_Dissemination_508.pdf`
   - "No **participant count** of 1 to 20 can be published or distributed directly (a count
     of 0 is permitted)."
   - "No data or statistics can be reported that allow a **participant count** of 1 to 20 to
     be derived from other reported cells... This includes... percentages or other
     mathematical formulas."
   - "data users... may obscure these values using scientifically accepted strategies,
     including collapsing data across cells, coarsening data, or cell suppression."
   - "If data users have a compelling reason... they may submit a request... for an
     exception. Exceptions will be rare."
2. **User Support Hub summary** — `support.researchallofus.org/hc/en-us/articles/22346276580372`
   — restates the same "participant count" language; no separate rule for variant/allele
   counts is stated anywhere in this policy.
3. **"How to comply..." article** — `support.researchallofus.org/.../360043016291` (search-cached
   text; direct fetch 403'd, likely login-gated) — gives the worked example: "report the
   total N... but the participant counts of the subgroups are suppressed" — i.e. the
   mechanism is always framed as *participants per subgroup*, not *properties per person*.
4. **Data User Code of Conduct (DUCC)** — `researchallofus.org/faq/data-user-code-of-conduct/`
   and PDF mirrors — "will NOT publish... any data or aggregate statistics corresponding to
   fewer than 20 **participants** unless expressly permitted." Same participant framing.
5. **All by All browser exception** (precedent for a formal exception request):
   multiple secondary sources (NIH Record, All by All docs) state the All by All
   genotype-phenotype browser "applied for an exception to the Data and Statistics
   Dissemination Policy to display genotype-phenotype associations with participant counts
   fewer than 20, and the program's Resource Access Board granted this exception on
   5/17/2024, in light of the browser's scientific utility and minimal risk." I could not
   fetch the All by All FAQ page directly (empty on WebFetch) to get AoU's own wording of
   the grant, so treat the RAB-approval fact as medium confidence, sourced secondarily.
   Important: even a *public, aggregate, variant-association* product needed an explicit
   exception to show n<20 — AoU's operating practice treats variant/genotype counts tied to
   carrier numbers as falling under the same rule, not as an automatic exemption.
6. **AoU genomics papers**: the 2023 Nature paper ("Genomic data in the All of Us Research
   Program") and the 2024 press coverage of "275 million new variants" report only
   **aggregate, cohort-wide** novel-variant tallies (e.g. "275,000,000 new variants," "3.9M
   with coding consequences") — never a specific variant's carrier count when that count is
   small. I did not find a published AoU paper that prints an individual rare variant's
   exact carrier count below 20; where per-variant frequencies appear (All by All, VAT),
   they are behind the Workbench/exception framework above, not in text.
7. **Context, not AoU rule — UK Biobank**: community guidance (`community.ukbiobank.ac.uk`,
   "Reporting small numbers...") gives thresholds of "a minimum of 5" participants per cell
   in tables/results and "100" for browser-served summary stats — a *stricter*, and also
   *participant-count-framed*, rule; no separate carve-out for variant/allele counts either.
8. **Context, not AoU rule — gnomAD**: public variant browsers commonly report singleton
   (AC=1) counts openly since gnomAD's ethics framework treats de-identified aggregate
   variant frequency as non-disclosive at population scale; some adjacent databases instead
   mask AC below a threshold (e.g. "<5"). Not binding on AoU data; different consent/IRB
   basis (gnomAD source cohorts largely pre-consented for open aggregate release; AoU
   participants did not).

## Answers

**(1) Is "number of distinct alleles/variants" (count of alleles, not people) covered by the `<20` rule?**
Medium-high confidence: **no, not directly** — the policy's operative term throughout is
"participant count." A statement like "we found 48 novel KIR alleles in this cohort" is a
count of *alleles*, not of *people*, and doesn't by itself reveal how many participants
carry any one of them. Caveat: if a downstream reader could combine that number with other
published cells to back out a small participant count for a specific stratum (policy clause
2, "derived... including... mathematical formulas"), it's covered again. Richness/QC-tally
counts of this kind are the strongest case for loosening.

**(2) Can we report the number of alleles seen in exactly 1 person (singleton *class* size)?**
Medium confidence, **conditionally yes**: "37 alleles are singletons" is again a count of
alleles, not an assertion "participant X carries allele Y" — no individual is identified and
no cell says "N=1 participants have trait Z." Risk: if the gene/locus has very few total
carriers overall (e.g., a gene where total N carriers <20), reporting "singleton count = 12"
combined with total-carrier N could let a reader derive that ≤12 people account for most of a
small stratum — apply clause 2 caution per-locus, not blanket.

**(3) Can an individual novel allele's sequence/name be published if <20 people carry it?**
High confidence: **no, without an exception.** This is exactly the shape the policy targets —
"allele Y, carried by N=3" ties a specific named entity to a small participant count, same as
the All by All precedent needing RAB sign-off for n<20 genotype-phenotype cells. Publish the
sequence/name only with the exact carrier count suppressed (`<20`) or omitted, unless an
exception is obtained.

**(4) Rates/frequencies that back-calculate to <20 people?**
High confidence: **no** — explicitly named in the policy text ("percentages or other
mathematical formulas that... would allow an individual to deduce a participant count of
less than 20"). This is the least negotiable of the five questions after (3).

**(5) Exception / review process?**
High confidence it exists, medium confidence on mechanics: submit a request to AoU's Resource
Access Board with a "compelling reason"; "exceptions will be rare." The All by All case is the
only documented precedent found, and its justification leaned on "scientific utility and
minimal risk" for a *public aggregate browser*, not a paper. For this project, treat a formal
RAB exception request as the correct channel if we ever want to publish an actual per-variant
n<20, rather than deciding it ourselves.

## Proposed rule for this project

Non-negotiable, no matter what: **any figure where the reported quantity is a count (or
share) of *participants* — carriers of an allele, cases in a QC category, people in an
ancestry×gene cell — stays masked at 1–19 → `"<20"`, censored (never printed as 0, never
back-calculable from an adjoining N or rate).** This is what the DUCC/policy text actually
says and is what a RAB exception is for, not a call this project should make unilaterally.

Loosened for counts that are properties of *alleles/variants/QC outcomes*, not of people,
because the current practice masks these needlessly:

| Count type | Allowed exact | Must be `<20` | Needs supervisor / RAB review |
|---|---|---|---|
| Participant count / carrier count per allele or cell (any stratum) | | ✓ | — |
| Rate or frequency with numerator or denominator 1–19 | | ✓ | — |
| A **named** novel allele's exact carrier count | | ✓ | Only via RAB exception |
| Total number of distinct novel alleles found (per gene or overall) | ✓ | | — |
| Per-gene allele richness (distinct-allele counts) | ✓ | | — |
| Count of alleles in a recurrence class (e.g. "N alleles seen in exactly 1 person") | ✓, if the underlying stratum has ≥20 total carriers | | ✓ if the locus/stratum total carrier N itself is small |
| QC call tallies (e.g. "14 samples failed depth filter") — no link to phenotype/trait | ✓ | | Flag if the QC category is itself a proxy for a sensitive trait |
| Any count published alongside a matching named/identifiable allele or participant subset | | ✓ | Treat as participant count |

Practical instruction for scripts/figures: keep the existing `<20`-censoring machinery
(reference/`43`-series disclosure helpers) exactly as-is for anything indexed by
participant/carrier; add a second, unmasked path only for allele-indexed/QC-indexed tallies,
and route the ambiguous middle case (recurrence-class counts on very-low-carriage genes) to a
supervisor flag rather than auto-publishing. Any single-allele carrier count <20 stays
suppressed regardless of the above, pending an actual RAB exception request — that is a
decision for Marc/Aleix/supervisors to initiate, not something this document authorizes.
