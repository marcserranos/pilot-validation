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
| WS1 | Novelty re-definition by field; artifact accounting; protein-level novel alleles; IMGT-known CDS check | `WS1_novelty_by_field.md` | ◐ audit done; scripts 24, 25 being implemented |
| WS2 | QC audit: source of relative discordance, switch-test power, twins, dropout, read support | `WS2_qc_audit.md` | ◐ audit done; script 26 being implemented |
| WS3 | Allele-space coverage by resolution and ancestry; cross-ancestry matrix; sampling scenarios | `WS3_allele_space_coverage.md` | ◐ `_coverage.py` library being implemented |
| WS4 | Figure 1 v2 from WS1–3 results | `WS4_figure1.md` | ☐ |
| WS5 | Literature positioning (coverage estimators, long-read HLA at scale, IMGT submission) | `WS5_literature.md` | ◐ research agent running |
| WS6 | Remaining supervisor items: LD (2), structural variation/KIR (6), selection lit (7) | `WS6_supervisor_items.md` | ☐ |
| — | HLA × TCR/BCR (needs Aleix's join) | plan in `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` §3 | ⏸ Aleix |

VM: ⏸ Workbench tab is at the login page (2026-09-16) — Marc must sign in. Local work proceeds.

## 4. Resume here (always current) — updated 2026-09-17

**State:** branch `fig1-drafts-and-research-map`, committed locally, **NOT pushed** (classifier
blocks publishing to the public repo; Marc must push).

**Done:** WS1 (scripts 24, 25 full cohort), WS2 (script 26 full cohort), `_coverage.py`, WS5
literature. All aggregate results are committed under `reports/hla_popgen/2{4,5,6}_*`.

**VM channel (ENVIRONMENT quirk #37):** JupyterLab tab →
`https://9394ec22-d949-4489-b46e-73790750472e.workbench-app-prod.verily.com/lab`;
helpers defined in the page: `window.bg(key,cmd)` / `window.res[key]` / `window.show(key,a,b)`,
`window.deploy(paths)` (laptop → VM, needs the local server:
`python3 <scratchpad>/serve.py`, GET serves scripts/+sprints/, POST /save writes into
reports/hla_popgen/), `window.pull([[vmPath, repoPath]...])` (VM → laptop, no file content passes
through the model). VM helper scripts: `~/run.sh <script.py> <outname> [args]` (detached, unbuffered),
`~/runtests.sh <module...>`. Worktree `~/repos/pv-s01`; python
`~/repos/pilot-validation/.pixi/envs/spechla/bin/python` (3.8/pandas 2.0.3).
**Trap:** after deploying a changed script, clear `__pycache__` and confirm the run used the new
file — a stale run once produced a completely wrong QC number.

**In flight:** script 27 (coverage) implementer agent; a verification agent on the suspected
off-by-one in the existing `21_hla_manhattan.py`.

**Next actions, in order:**
1. When 27 lands: deploy, test on VM, pilot (`--limit-genes 2`), then full run; pull results.
2. Compose Figure 1 v2 (new script 28) from: admixture panel (existing 06), class I ternary
   (existing 10), clean novelty rate by ancestry and field (24), coverage curves + cross-ancestry
   matrix (27). Then a supplement list.
3. Decide the 21_hla_manhattan question from the verification agent's verdict.
4. Write the final plain-language report for Marc (objective → finding → meaning → why novel),
   update `context/EXPERIMENTS.md` with pointer entries, and ask Marc to push the branch.

**Constraints:** subagents = sonnet only; expect interruption; commit after every step.

## 5. Headline findings (distilled from briefs; newest first)

**WS1 — what "novel" really is (script 24, full cohort: 12,233 people, 11,856 unrelated, 875,179
haplotype-gene calls, 2026-09-16):**

| step | count |
|---|---|
| calls the pipeline flagged novel (any field) | 280,695 |
| …novel only outside the CDS (field 4) | 172,272 |
| …novel synonymous CDS (field 3) | 30,795 |
| …novel protein (field 2) | 75,442 (71% carry an artifact flag; 22,122 clean) |
| of field-2/3 calls with a recovered sequence (104,532): CDS already in IPD-IMGT (naming artifact) | 25,282 |
| …frameshift or premature stop | 72,667 |
| …genuinely novel protein calls | 5,288 |
| **distinct novel proteins** | **1,404** (1,026 entirely clean) |
| **novel proteins recurrent in ≥2 unrelated people** | **231** |
| distinct novel synonymous CDS | 419 |

So the honest headline is **hundreds of new HLA proteins, not thousands**; ~25k novel-flagged calls
are a naming artifact (their CDS is catalogued); and the bulk of "protein-altering" novelty is
frameshift/stop calls dominated by homopolymer artifacts.

**WS2 — QC (script 26, full cohort):** relatives agree at 87.0% of 21,545 comparable gene
comparisons; only **0.8%** are unexplained sequence differences (6.5% is a gene called in one
relative only, 3.2% dropout candidates, 2.3% single-base). Duplicate/MZ replicates give the direct
error rate: **QV ≈ 34** (121 differing bases in 306,033). Phasing: **0 switches over 3,021 testable
transitions at 100% detection power** — unlike the old "545/545", this test could have failed.
Homozygosity excess vs the same people's short-read calls (1.10 vs 1.07) shows small real dropout.
Novel-protein calls are the one place short reads disagree (0.93 of 2 alleles matching vs 1.78 for
known alleles).

**WS1 — non-coding novelty (script 25, full cohort, classical genes):** the 708 former CDS-hash
clusters hide **15,112 distinct genomic sequences**; the largest single cluster pools **287** of
them. Only 184 signatures reach ≥20 unrelated carriers. **21.3% of non-coding novel haplotypes
differ only by homopolymer indels** (the long-read error mode). Marc's clustering hypothesis is
confirmed with numbers, and the 2026-09-15 "one allele in 1,531 people" callout is retracted.

**Carried in from 2026-09-15 (still true, now better quantified):** non-coding novelty dominates,
and the flagged-artifact rate is flat across ancestries (a built-in control).
