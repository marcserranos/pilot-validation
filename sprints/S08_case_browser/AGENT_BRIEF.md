# Agent brief — are the substitution-only "novel" HLA alleles real, and do they matter?

Read first: `CLAUDE.md`, `context/ENVIRONMENT.md` (VM quirks), `reports/hla_popgen/90_case_browser/README.md`, `sprints/S08_case_browser/LOG.md`.
Work on branch `s08-case-browser` (worktree `../pilot-validation-s08`). New scripts take numbers **91+** (S06 owns 70-75, 76-89 free).

## Mission
Take the promising non-indel novel HLA calls, QC them as hard as the data allows, find the ones recurring across people, say what each would
mean biologically, and test at **group level** whether carriers show the condition(s) that biology predicts. Systematic, logged, every claim traceable.

## What already exists (do not redo)
- Whole-cohort stats (`cohort_stats.json`): 8,015 novel calls; of 7,635 novel-protein only 449 are substitution-only, 8% carry a protein-altering substitution.
- Viewer cases + checks (local `case_browser/`, VM `~/case_browser/`): Immuannot's own CDS alignment (exact PAF replay), S06 long-read depth (partial coverage),
  short-read realignment per variant (`sr.jsonl`, zygosity-aware verdicts). Indels are ~90% unsupported by short reads: ignore them.
- **`case_browser/shortlist.tsv`** (`90_case_browser.py shortlist`): 337 candidates (novel, substitutions only, >=1 protein-altering). Tier A = every checked site supported (189),
  B = supported and none contradicted (5), C = the rest (128), plus 15 without short reads. HLA-F (145) and G (211) were capped at 60 each in the viewer; the full counts are in the stats.
- `case_id` -> person map exists **only on the VM** (`~/case_browser/case_map.private.tsv`, mode 600). Local files carry anonymous ids.

## Suggested stages (log each in LOG.md, one script + test per stage)
1. **Deeper QC** per candidate: strand balance / MAPQ / mismatches near the site in the short reads; paralog or multi-mapping reads; long-read support at the site (BAM or PAF);
   is the variant already in a newer IPD-IMGT release or a known allele (we only have IPD 3.55 snapshot); phase consistency with the other haplotype. Promote a candidate only if several independent checks agree.
2. **Recurrence**: group novel calls by (gene, set of CDS substitutions) from `hla_calls_rich.tsv` `cds_mut` (VM, one pass, aggregate output); counts by ancestry (`cohort_membership.tsv`, lowercase `ancestry_pred`),
   transmission in trios (57b / S06 trio tables), and shared haplotype background. Recurrent + supported + transmitted = real allele.
3. **Biology**: map each change onto the groove (`reference/ars_peptide_contacts.tsv`, `ars_chain_sequences.tsv`, `scripts/hla_popgen/_ars_residues.py`): pocket, conservative vs radical, same residue as a known disease-associated allele; literature for the nearest allele (web search is fine).
4. **Phenotype** (Workbench only): for each recurrent, well-supported allele choose a prior-justified condition list (start from `context/TASK_CONTEXT.md`: celiac, T1D, psoriasis, MS; add what the biology implies),
   then compare carriers with ancestry/age/sex-matched non-carriers. Report **counts and rates for groups only**.

## Hard rules
- **Participant-level data never leaves the Workbench and never enters the repo, a commit message or a report.** Individual-level phenotype is looked at on the VM and not written anywhere. Outputs are aggregates; counts 1-19 are written `<20` (never 0). A single carrier is never reported with a diagnosis.
- This repo is public. Local `case_browser/`, `*.private.tsv` are gitignored: keep it that way. Never broad-`git add`; add files by name. Do not touch other branches; another agent works on S06.
- **VM is shared and busy**: one core, `nice -n 10`, detached `nohup` with a checkpoint file, own folder under `~/case_browser/` or `~/s08/`, a throwaway `git clone --depth 1 --branch s08-case-browser` (never the shared `~/repos` checkout), remove it afterwards. The Jupyter websocket and the Chrome extension drop often: poll a log, do not hold a connection. `type` and `Return` are separate tool calls.
- The Claude Code classifier blocks agent-initiated downloads of participant-derived files; ask Marc to download, or keep the work on the VM.
- Ask Marc before: publishing anything, exporting any individual-level phenotype, using more than one core, or committing results with small counts (see `context/DECISIONS.md` open items).

## Known gaps and things nobody has verified
- **EHR / condition access is the unverified step.** The CDR lives in BigQuery (`C2025Q4R6`, project `wb-silky-artichoke-2408`, `ENVIRONMENT.md` quirk 5); the repo has no notes on the condition/EHR tables, who has EHR data, or the right query path. Probe once, read-only, and write what you find into `context/ENVIRONMENT.md`.
- Statistical power: 6% of novel-protein calls and 8,015 calls overall means most alleles will have few carriers; plan for pooled tests (by gene, groove position) and say plainly when it is underpowered. Multiple testing is a real risk with 4+ conditions x many alleles.
- S06 depth rows cover only part of the cohort; DRB3/DRB4 have no primary-assembly locus (no short-read check); verdicts are for the variants checked (up to 6 per case), not the whole gene.
- "Novel" here = Immuannot's `new` against IPD 3.55; relatedness (trios, kin) inflates apparent recurrence unless one counts unrelated people.

## Quick start
```
python3 scripts/hla_popgen/90_case_browser.py shortlist --cases case_browser/cases.json --out case_browser/shortlist.tsv   # local
python3 -m pytest scripts/hla_popgen/tests/test_90_*.py -q                                                                 # 33 tests
```
