# Sprint S02 — supervisor items from the 2026-09-17 call

*Started 2026-09-17. Branch `fig1-drafts-and-research-map` (continuing S01's single branch).
Orchestrator: Claude (Opus). Owner: Marc. Mode: autonomous, report back at the end.*

Read `sprints/README.md` first for the L0–L3 layering. S01 (`../S01_novelty_qc_coverage/`) is the
foundation this sprint builds on; its `FINDINGS_FOR_MARC.md` supersedes several numbers that were
presented in the call.

---

## 1. What the call actually asked for

Attendees: Marc (presenting), Cole Shanks, David Bonet. Marc walked through novel-allele discovery,
allele frequencies / population structure, the disease-mapping negative, phasing, and the Manhattan
diversity tracks. Distilled asks, in Cole's priority order:

| # | Ask (paraphrased, Cole unless noted) | Where it lands |
|---|---|---|
| A1 | **Figure 1 panels.** Admixture bar plot; a class I ternary ("the best one you can make", A or B); novel-allele rate by ancestry **broken out by gene** ("where is the reference data bias worst?"); allele-frequency spectrum over **all** alleles with novel vs known as two colours, overlaid, not three separate plots. | WS4 |
| A2 | **Phasing error is the critical question.** "We don't want there to be a switch between DPA1 and DPB1" — the pair is what matters functionally. Estimate that error properly. | WS0 — **already answered in S01**; needs restating in the pair-specific form Cole asked for |
| A3 | **LD between DQA1–DQB1 and DPA1–DPB1 alleles, computed within each ancestry**, r². "I suspect the linkage patterns are very different between ancestries." Raised twice, explicitly listed in the closing next-steps. | WS1 |
| A4 | **Deletions and duplications.** "Let's do that. It could be pretty crazy." Notes the Immuannot/1000G paper reports dups and dels. | WS2 |
| A5 | **KIR** — "did you look at the KIR alleles at all? It should be called if you just ran it straight up." | WS2 (quick answer: it is not, by construction) |
| A6 | **A great list of novel alleles that are common in some ancestry and missing from the database** — "things we can literally call out in main text: we found these seven examples." | WS3 |
| A7 | **Selection / diversity literature review** — "so that we can cite somebody and at least their methodology". Specifically amino-acid-level diversity and peptide-binding vs non-peptide-binding residues within the CDS. | WS5 |
| A8 | **Amino-acid-level diversity, ARS vs non-ARS.** "What's actually more interesting is if we look at the CDS between the peptide-binding and the non-binding regions." | WS5 |
| A9 | Stricter admixture threshold (~98%) for structure plots, because the UMAP/PCA separation is confounded by admixture. | WS4 |
| A10 | Class II looks more diverse / more structured than class I — is that balancing selection differing between ancestries? | WS5 |
| A11 | Supplement: start a figure dump (SFS by gene, the extra ternaries, the UMAP/PCA panels). | WS4 |
| A12 | Short-read validation of novel alleles by realigning short reads to the long-read assembly. Cole: **"I wouldn't do that right now"** — parked deliberately. | ⏸ parked by Cole |
| A13 | HLA is not a prediction target; it is a **conditioning variable** for the TCR/BCR work. Confirms the disease-embedding negative is fine to report as a negative. | no work; framing note |

David's framing point: Figure 1 stays HLA-only; the TCR/BCR material is a later section, not a
combined overview figure. Cole agreed and wants a linear narrative.

## 2. What the call got wrong (must be corrected with the supervisors)

Marc presented S01's *pre-correction* numbers, because S01 finished after the meeting. On the
record in the call, and now superseded:

| said in the call | corrected (S01, `FINDINGS_FOR_MARC.md`) |
|---|---|
| "~3,000 novel alleles" | 280,695 novel-flagged **calls** → **1,404 distinct novel proteins**, **231** recurrent in ≥2 unrelated people |
| "292 seen twice or more, 904 seen three or more" | those were cluster counts under the old CDS-hash clustering, which pooled distinct genomic sequences (708 clusters = 15,112 real sequences) |
| "Immuannot only looks at the CDS and then clusters around that, which is what we do" — offered as reassurance | this is exactly the artifact: it is why the counts above collapsed |
| "5–20% of the allele space discovered" | that was a *richness* estimate. Protein-level **coverage** is 99.0–99.8% per ancestry (MID 97.1%). Both are true; the headline must be coverage, not richness |
| "95% concordance between relatives, probably ambiguity" | 87.0% byte-identical; decomposed as 6.5% gene called in one member only, 3.2% dropout, 2.3% single-base, 0.8% unexplained. Twin/duplicate QV ≈ 34 |
| "zero phase switches, 545/545 pairs" | still zero — but now over **3,021 testable transitions** with **100%** detection power on injected synthetic switches. The old test could not have failed |

**This is the single most important thing to tell Cole**, because A6 ("a list of common novel
alleles we can call out in main text") is built directly on the corrected numbers, and because the
old headline would not survive review.

## 3. Reflection before starting

**What we can doubt**
1. *LD at allele level.* r² is defined for biallelic markers. HLA alleles are multi-allelic, so the
   choice of statistic is the whole analysis: pairwise-biallelic r² per allele pair, a normalised
   multi-allelic D′ (Hedrick), or a normalised mutual information. Cole asked for "the R2 measure";
   we should give him that **plus** a multi-allelic summary, and say which is which.
2. *SV from Immuannot is presence/absence of an annotation, not a validated structural call.* A
   "deletion" can be assembly fragmentation (the gene fell off a contig end) and a "duplication"
   can be a mapping artifact in a segmental duplication. The negative control matters more than the
   positive: DRB3/4/5 copy number is *known* biology and must come out right, or we cannot trust
   anything else the method says.
3. *Ancestry-stratified anything* inherits the admixture confound Cole and David both flagged.
   Strict thresholds shrink n, particularly for AMR and MID.
4. *ARS residue lists* are usually published in mature-protein numbering; our codon positions come
   from IPD-IMGT CDS coordinates including the leader peptide. An off-by-leader-length error would
   silently invert the whole ARS-vs-non-ARS result. This needs an explicit alignment check against
   a known anchor residue.

**Low-hanging fruit**
- KIR is a one-line answer from `SCHEMA.md` + a grep of the calls table.
- The pairable-fraction number (genes split across contigs) is already computed for `hla_cis_pairs`
  and is itself a result about assembly fragmentation, which is half of A4.
- The phasing answer Cole wants (A2) exists; it needs re-cutting per gene pair.

## 4. Workstreams and task board

Status: ☐ todo · ◐ in progress · ☑ done · ⏸ blocked/parked

| WS | Goal | Brief | Status |
|---|---|---|---|
| WS0 | Scaffolding; restate S01's phasing QC in the DPA1–DPB1 / DQA1–DQB1 form Cole asked for (A2) | this file | ◐ |
| WS1 | LD between DQ and DP alleles within each ancestry (A3) | script 29 | ◐ written + tested + smoke-tested; needs the VM run |
| WS2 | Structural variation: gene deletions, duplications, DRB copy number; KIR answer (A4, A5) | script 30 | ◐ written + tested + smoke-tested; needs the VM run |
| WS3 | Main-text callout list of common novel alleles absent from IPD-IMGT (A6) | script 32 + `reports/hla_popgen/32_novel_callouts/README.md` | ☑ done |
| WS4 | Figure 1 v3 + supplement dump, strict-admixture variants (A1, A9, A11) | `WS4_figure1_v3.md` | ☐ |
| WS5 | Selection: ARS vs non-ARS amino-acid diversity; class I vs II differentiation; literature (A7, A8, A10) | `WS_literature_selection.md`, `WS5_ars_definition.md`, script 31 | ◐ literature + ARS definition done; script 31 needs the VM run |
| — | Short-read validation of novel alleles (A12) | — | ⏸ parked by Cole |
| — | HLA × TCR/BCR join | `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` §3 | ⏸ Aleix |

## 5. Resume here (always current)

*2026-09-17.* Scripts **29** (LD), **30** (structural variation), **31** (amino-acid diversity /
differentiation) and **32** (novel-allele callouts) are written, unit-tested (127 fixture tests
across the four) and smoke-tested against a synthetic cohort built in the session scratchpad.
**32 has already been run for real** — it needs only committed aggregates. **29, 30 and 31 still
need their full-cohort VM run**; that is the next action.

Literature (A7) and the structural ARS definition are done and committed.

VM: assume it restarted (ENVIRONMENT quirk #14); re-verify the auto-mount before trusting
anything that reads AoU data. **Do not touch `~/mnt/aou-controlled`** — quirk #35, it hangs the
process irrecoverably; all three scripts already default to the `~/workspace/` auto-mount. S01's
VM channel and traps are in `../S01_novelty_qc_coverage/SPRINT.md` §4 and quirks #35–38.

Command to run on the VM once the mount is verified:

```
for s in 29_hla_ld_by_ancestry 30_hla_structural_variation 31_aa_diversity_selection; do
  setsid nohup python3 -u scripts/hla_popgen/$s.py --out-dir ~/results/$s < /dev/null & disown
done
```

(Clear `__pycache__` after deploying changed files and verify the run used them — S01 lost a day
to a stale deployed script that produced plausible wrong numbers.)

## 6. Headline findings (distilled, newest first)

**WS3 — the main-text callout list (script 32, committed aggregates, 2026-09-17).**
38 clean, recurrent novel alleles are carried by **>= 20 unrelated people**; 29 are novel
proteins and 9 are novel synonymous CDS (reported separately, since they are *not* new proteins).
**None is in a classical HLA gene.** They are TAP1/TAP2 (26), DQB2/DQA2/DM/DO (10), MICB, HLA-G,
HLA-F. The largest: a TAP1 protein one residue from TAP1*01:01, in **426** unrelated people and
**~139 per 1,000** African-ancestry participants.

The obvious objection — "TAP is just under-catalogued" — is answered with script 27's independent
measurement, printed beside every callout: the share of a gene's haplotypes whose protein is
already absent from IPD-IMGT is **8.0% for TAP1 and 6.9% for TAP2**, against **0.1-0.2% for
HLA-A/B/C**. The genes with the most novel alleles are the genes with the largest catalogue gap,
so the two independent estimates agree. That is evidence for the reframing S01 proposed (the
paper's spine is the non-classical MHC), not against it.

A further 309 clean recurrent alleles sit below the disclosure threshold: the IMGT submission
queue, countable but not nameable.
