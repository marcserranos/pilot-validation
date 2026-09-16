# Sprint S01 — novelty re-definition, load-bearing QC, allele-space coverage

*Started 2026-09-16. Branch `fig1-drafts-and-research-map` (single branch for the sprint).
Orchestrator: Claude (Opus). Owner: Marc. Mode: autonomous, report back at the end.*

## 1. Why this sprint exists (how it fits the research)

The paper's first claims are: (1) long reads find a lot of HLA sequence absent from IPD-IMGT/HLA,
more so outside Europe; (2) the allele space is far from saturated; (3) the calls are trustworthy.
On 2026-09-15 a three-way split of novelty showed that claim (1) is mostly non-coding novelty, and
Marc suspects the novel-allele clustering itself inflates/merges groups. Claims (2) and (3) inherit
whatever is wrong with (1). So before any figure is finalized, this sprint re-derives the
foundation: **what exactly is novel, at what resolution, how sure are we, and how much of the space
have we seen.** Everything downstream (Figure 1, IMGT submission, disease/TCR work) sits on it.

## 2. Reflection before starting

**What we can doubt (ranked by how much of the paper rests on it)**
1. *Novelty definition.* Clusters are keyed on CDS hash. A `beyond_cds` cluster is "haplotypes whose
   CDS equals a known allele's CDS but whose genome didn't match any IMGT genomic sequence" — so
   "one novel allele in 1,531 people" may be many different intron variants pooled. Field-4
   novelty also depends on IMGT genomic completeness, which is itself ancestry-biased.
2. *The QC numbers.* "0 phase switches" and "~5% relative discordance" coexist. Either the switch
   test has little power (few testable transitions), or the discordance is mostly not phasing
   (allele dropout from collapsed haplotypes, consensus base errors, naming ties, misclassified
   relatives). Relatives were also *selected* by high sharing, so the rate is conditional.
3. *Saturation estimates.* Chao2 was LOW-CONFIDENCE everywhere; novel "alleles" include
   artifacts and pooled non-coding groups; relatives may inflate doubletons; richness is the
   wrong target for a heavy-tailed distribution — **coverage** (probability the next haplotype is
   already catalogued) is estimable far more reliably.
4. *Reference version.* Immuannot ships a fixed IPD-IMGT snapshot; some "novel" may be catalogued
   in newer releases.
5. *Ancestry labels.* Admixed people blur per-ancestry numbers; supervisors asked for stricter
   thresholds.

**What we can additionally test**
- Novelty at each nomenclature field: field 2 (new protein), field 3 (new synonymous CDS),
  field 4 (new non-coding), with protein-level clustering and position of changes (groove exons
  2/3 vs elsewhere).
- Direct read-level support (HiFi BAMs exist for every long-read person) for novel protein alleles.
- Identical-twin/duplicate pairs as true technical replicates; homozygosity excess as a
  dropout detector; per-base error implied by single-base discordances.
- Coverage-based extrapolation per ancestry and resolution; a cross-ancestry "catalogue built from
  population A covers X% of population B" matrix; sample-allocation scenarios ("n AFR + n EUR + …
  reaches X% protein-level coverage").

**Low-hanging fruit next to this problem**
- Novel-allele tables at protein level are exactly the IMGT submission candidates list.
- The coverage matrix is a reference-panel design argument (useful beyond HLA; KIR/IG loci next).
- Field-level novelty by ancestry is a cleaner reference-bias figure than raw novel-call rate.
- Relatives analysis can yield parent-child vs sibling split and MZ-twin error for free.

## 3. Workstreams and task board

Status: ☐ todo · ◐ in progress · ☑ done · ⏸ blocked (reason)

| WS | Goal | Brief | Status |
|---|---|---|---|
| WS0 | Scaffolding, commit drafts, STATUS pointer | this file | ◐ |
| WS1 | Novelty re-definition by field; artifact accounting; protein-level novel alleles; IMGT-known CDS check | `WS1_novelty_by_field.md` | ◐ audit running |
| WS2 | QC audit: source of relative discordance, switch-test power, twins, dropout, read support | `WS2_qc_audit.md` | ◐ audit running |
| WS3 | Allele-space coverage by resolution and ancestry; cross-ancestry matrix; sampling scenarios | `WS3_allele_space_coverage.md` | ☐ |
| WS4 | Figure 1 v2 from WS1–3 results | `WS4_figure1.md` | ☐ |
| WS5 | Literature positioning (coverage estimators, long-read HLA at scale, IMGT submission) | `WS5_literature.md` | ☐ |
| WS6 | Remaining supervisor items: LD (2), structural variation/KIR (6), selection lit (7) | `WS6_supervisor_items.md` | ☐ |
| — | HLA × TCR/BCR (needs Aleix's join) | plan in `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` §3 | ⏸ Aleix |

VM: ⏸ Workbench tab is at the login page (2026-09-16) — Marc must sign in. Local work proceeds.

## 4. Headline findings (distilled from briefs; newest first)

- (none yet this sprint) Carried in from 2026-09-15: after removing flagged artifacts, 93–99.7% of
  novel classical-gene haplotypes are non-coding-only; 7 recurrent CDS-changing novel alleles.
  **Status: being re-derived in WS1 — treat as provisional.**
