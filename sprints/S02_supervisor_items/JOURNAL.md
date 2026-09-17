# S02 journal — append-only

## 2026-09-17

- Sprint opened. Marc pasted the full transcript of the supervisor call (Cole Shanks, David Bonet,
  Marc). Distilled into 13 asks (A1–A13) in `SPRINT.md` §1 and mapped onto 6 workstreams.
- Recorded in `SPRINT.md` §2 the six numbers Marc presented in the call that S01 has since
  superseded. S01 finished *after* the meeting, so the supervisors are working from the old
  headline ("~3,000 novel alleles", "5–20% of the allele space discovered"). Correcting this is
  the first thing to tell them.
- Dispatched a sonnet research agent for WS5's literature half (A7/A8): ARS vs non-ARS diversity
  methodology, Ewens–Watterson / Tajima's D / dN/dS citations, and — most useful — a citable
  peptide-binding-residue list for class I and class II.

- Wrote and unit-tested scripts 29 (LD by ancestry), 30 (structural variation + KIR), 31
  (amino-acid diversity, groove vs non-groove, allele differentiation) and 32 (novel-allele
  callouts). 127 fixture tests across the four, all passing.
- Built a synthetic 400-person cohort in the session scratchpad and smoke-tested 29/30/31
  end to end. Two checks worth recording: 29 recovered an injected 90% DQA1–DQB1 linkage as
  D' ~ 0.90–0.94 while correctly reporting ~0 (and a bias-corrected Cramer's V of exactly 0) for
  the independent DPA1–DPB1 pair; 30 reported 0.0% deletions at every negative-control gene, so
  the method's false-positive rate on clean synthetic data is zero.
- Caught before it cost anything: scripts 29/30/31 had defaulted their relatedness-table path to
  `~/mnt/aou-controlled`, which ENVIRONMENT quirk #35 records as a stale mount that hangs any
  process touching it in uninterruptible I/O. Repointed to the `~/workspace/` auto-mount.
- Dispatched a second research agent to derive peptide-contact residues structurally from PDB
  rather than cite the ARS residue tables the literature agent could not verify. It returned
  `_ars_residues.py`, `reference/ars_peptide_contacts.tsv` (207 residues at 4.5 A across six
  structures) and a brief, having checked each PDB entry against its own header — it caught that
  1DLH's chain C is mislabelled in COMPND and confirmed the peptide from SEQRES instead, and that
  3LQZ is a single-chain construct with no separate peptide chain. A follow-up is pending so the
  residue numbers can be transferred onto our own protein numbering offline.
- Ran script 32 for real (committed aggregates only). See SPRINT.md §6.
- Added the catalogue-gap control to the callout list after noticing the list was dominated by
  TAP1/TAP2 and could be read as reference depth rather than biology. It is measurable, script 27
  already measured it, and it agrees with the novelty ranking.
- VM session. The JupyterLab tab was still alive; re-established the programmatic terminal channel
  (quirk #37). One correction to that recipe worth recording: the terminal **echoes the command it
  is given**, so a literal BEGIN/END marker inside the command text matches before the real output
  does. Emit markers with `printf '<<%s>>\n' 'Bkey'` so the echoed line never contains them.
- Deployed scripts 29/30/31 and `_ars_residues.py`, verified md5 against the local copies, cleared
  `__pycache__`, launched all three. **The VM then went down mid-run** — the Workbench console
  says a reboot is required — and all three processes died with the gateway returning 502. It came
  back on its own about four minutes later; relaunched and all three completed.
- A second retrieval quirk: **jupytext is installed on this VM and claims `.md` files as
  notebooks**, so `/api/contents/<path>.md?content=1` returns nbformat JSON with empty content
  instead of the file. Pass `&type=file&format=text` to get the real bytes. Without it the three
  result READMEs silently came back empty.
- All three full-cohort runs completed on 11,856 unrelated people, matching S01's cohort exactly —
  a free cross-check that the cohort construction is deterministic.
- Every control passed: the LD method recovers strong linkage at B~C and weak at A~B; the deletion
  caller reproduces the DR51/52/53 expectation in 92-98% of haplotypes with a 0.04-2.0%
  false-positive rate at genes that are never deleted; the groove-diversity test shows no
  enrichment in DRA, HLA-F or DRB4.
- Ran a fresh-context critic agent over `FINDINGS_FOR_MARC.md`, checking every quantitative claim
  against the committed tables.
- **The critic found five real errors, two of them critical, and every one ran in the direction of
  making the result look better.** Recorded in full in `FINDINGS_FOR_MARC.md` §5b. The worst was a
  "ten most diverse residues" list that showed eight, dropping in each of three genes exactly the
  two residues that did not match a published epitope — and then presenting the agreement as a
  sanity check. The second worst repeated, four paragraphs before describing it, the very
  suppressed-count error the document congratulates itself for catching.
- Lesson for future sprints: writing a paragraph about a class of mistake does not stop me making
  it in the same document. The critic pass is not optional polish; it caught things no amount of
  re-reading my own text would have.
- The critic also spotted, in passing, that S01's committed findings state "the five
  duplicate/identical-twin pairs" — a bare participant-pair count below 20 in a public repo. Not
  S02's file, but it must be fixed before the branch is pushed.
