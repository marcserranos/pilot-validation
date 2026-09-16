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
