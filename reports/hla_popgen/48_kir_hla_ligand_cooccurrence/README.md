# 48 — KIR-HLA ligand co-occurrence by ancestry

*Status: script + synthetic tests done locally (S04 WS-D, 2026-09-25/26). Epitope-group assignment
is now sequence-derived (translated from IPD-IMGT/HLA CDSseq reference CDS, `--cds-dir`), with the
former two-field lookup table kept only as a fallback (missing CDS file, novel call, unmapped
allele) and a cross-check. Not yet run on real AoU data — needs a VM session with
`~/tools/Immuannot_refdata/CDSseq/` present to exercise the sequence path end-to-end (see
Caveats).*

## Question

Are functional KIR-HLA receptor-ligand pairs (KIR2DL1 x HLA-C2, KIR2DL2/2DL3 x HLA-C1, KIR3DL1
x HLA-Bw4, KIR3DS1 x HLA-Bw4) enriched or depleted relative to independence, within each
ancestry, in the AoU long-read cohort? Long reads give both loci (chr6 HLA, chr19 KIR) from the
same person, which is otherwise hard to get at scale.

## Method

- **Ligand groups, sequence-derived (primary) with a two-field lookup fallback + cross-check:**
  - Each called allele is mapped to an IPD-IMGT/HLA CDSseq reference record (`--cds-dir`, default
    `~/tools/Immuannot_refdata/CDSseq/<gene>.fa.gz`; headers like
    `>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp`) — exact match at the call's own field
    resolution, else the alphabetically-first reference allele sharing its two-field prefix
    (`mapping_level` = `exact` / `two_field_fallback`).
  - The CDS is translated respecting `frame=`, the 24-aa class-I leader peptide (HLA-A/B/C signal
    peptide) is stripped, and the mature protein's residues 77–83 are read directly.
  - **HLA-C1/C2:** C1 = Asn80 with Ser77; C2 = Lys80 (Colonna et al. 1993 PNAS; Winter & Long
    1997 J Immunol). Verified against known reference alleles: C*01:02/C*03:04/C*07:01 = C1
    (Asn80); C*02:02/C*04:01/C*05:01/C*06:02 = C2 (Lys80).
  - **HLA-Bw4/Bw6:** Bw6 reference pattern Ser77-Asn80-Leu81-Arg82-Gly83; Bw4 = Arg83 with Ile80
    or Thr80 (Gumperz et al. 1995 J Exp Med; Cella et al. 1994; Parham reviews, e.g. Parham 2005
    Nat Rev Immunol). Bw4-80I vs Bw4-80T is reported separately (`n_bw4_80I`/`n_bw4_80T` in
    `ligand_lookup_qc.tsv`) since KIR3DL1/3DS1 binding avidity tracks this dimorphism.
  - A residue pattern matching neither rule is `other` (counted, not silently folded into a
    group). A call the sequence path can't resolve — **novel allele** (protein at 77–83 not
    guaranteed to match any reference), missing CDS file for that gene, or no reference allele
    shares the two-field prefix — falls back to the two-field lookup table (`C1C2_TABLE`,
    `BW4_B_GROUPS`/`BW4_A_GROUPS`, same tables as before). Alleles unresolved by *both* paths are
    `unclassified`, counted and reported in `ligand_lookup_qc.tsv` (never silently folded into a
    group).
  - **Cross-check:** for every two-field allele group where the sequence path resolved,
    `ligand_seq_vs_lookup_crosscheck.tsv` records the sequence-derived label vs. the lookup-table
    label (allele two-field group names — catalogue facts, not person data — never paired with a
    carrier count). Aggregate agree/disagree counts are in `ligand_lookup_qc.tsv`.
  - B*15 Bw4 exceptions (B*15:13/15:16/15:17/15:24 and others) are handled correctly by the
    sequence path automatically; the lookup-fallback path still needs `--b15-as-bw4-list` since it
    has no access to the actual residues.
- **KIR receptors:** presence/absence per person (either haplotype) for KIR2DL1, KIR2DL2,
  KIR2DL3, KIR3DL1, KIR3DS1, via `41_kir_pilot.py`'s `parse_hap_gtf`.
- **Statistics:** per ancestry (n>=20 people), per pair: 2x2 contingency table (KIR present/absent
  x ligand present/absent), odds ratio (Haldane-Anscombe corrected on a zero cell), Woolf log-OR
  95% CI, Fisher's exact p, chi2 p, and a >=20-shuffle permutation null of the OR (shuffles KIR
  presence within ancestry, ligand held fixed).
- **Disclosure:** any of the four cell counts in [1,19] blanks the OR/CI/both p-values for that
  row; the masked cell counts themselves are still reported as `<20` (never blank, never 0 unless
  a true zero).
- **Epitope frequencies** per ancestry (`epitope_freq_by_ancestry.tsv`) are a cheap byproduct that
  also documents C1/C2/Bw4 carrier prevalence directly.

## Result

*(fill in after the VM run: `kir_hla_ligand_cooccurrence.tsv`, `epitope_freq_by_ancestry.tsv`,
`ligand_lookup_qc.tsv`)*

## How to read each output file

- `kir_hla_ligand_cooccurrence.tsv` — one row per (ancestry, KIR gene, ligand). `odds_ratio` > 1
  means the KIR gene and its ligand co-occur more than chance within that ancestry; `perm_p` is
  the permutation-based significance (preferred over `fisher_p`/`chi2_p` alone since it doesn't
  assume independence from ancestry-level allele-frequency structure beyond what's already
  stratified out by the per-ancestry split).
- `epitope_freq_by_ancestry.tsv` — % of people per ancestry carrying each epitope group.
- `ligand_lookup_qc.tsv` — coverage/source-mix counters: how many HLA-C/B/A calls resolved via
  sequence vs. lookup fallback (`n_{C,Bw}_source_{sequence,lookup_fallback}`), how many were a
  non-canonical residue pattern (`n_other_*_seq_pattern`) or unresolved by both paths
  (`n_unclassified_*`), Bw4 sub-type counts (`n_bw4_80I`/`n_bw4_80T`), and the sequence-vs-lookup
  agreement summary (`n_two_field_groups_seq_vs_lookup_{compared,agree,disagree}`). A high
  unclassified fraction or many disagreements means the epitope-group calls need a closer look
  before trusting the co-occurrence numbers.
- `ligand_seq_vs_lookup_crosscheck.tsv` — one row per two-field allele group where the sequence
  path resolved: `seq_label` vs. `lookup_label` and whether they `agree`. Allele group names only
  (catalogue facts), never paired with a carrier count.

## Caveats — read before trusting any number here

- **Sequence-derived assignment depends on `~/tools/Immuannot_refdata/CDSseq/<gene>.fa.gz` being
  present on the run host.** If a gene's CDS file is missing, every call for that gene silently
  (but countably, via `ligand_lookup_qc.tsv`) falls back to the two-field lookup table — check
  `n_{C,Bw}_source_sequence` is non-zero before trusting the run used real sequence, not just the
  fallback.
  - Not yet run on the VM with the real CDSseq directory present — local testing used only
    synthetic FASTA fixtures (`scripts/hla_popgen/tests/test_48_kir_hla_ligand_cooccurrence.py`),
    so the parsing/translation logic is verified but not yet exercised against the real IPD-IMGT/
    HLA 3.55.0 release files.
- **Novel alleles (Table 1's spliced `new` field) cannot be sequence-resolved** — their protein at
  77-83 isn't guaranteed to match the reference two-field group's — and fall back to the lookup
  table on their two-field prefix, counted in the source-mix QC.
- **B*15 Bw4 exceptions are handled correctly by the sequence path** (it reads the actual
  residues); the lookup-fallback path still needs `--b15-as-bw4-list` for any B*15 call that ends
  up unresolved by sequence.
- **KIR3DS1 has no confirmed direct HLA-binding ligand** — its inclusion here tests the published
  epidemiologic Bw4-80I association (Martin et al. 2002), not a biochemically confirmed
  receptor-ligand pair; interpret differently from the other four rows.
- Cannot distinguish real linkage-disequilibrium/selection from an ancestry confound or a miscall
  artifact on its own — a significant row is a lead to look into, not a finding to publish as-is.

## Distilled

- Not yet run on real data. Epitope-group assignment is now sequence-derived (CDS translation +
  leader-peptide stripping + residue 77-83 rules), with the former lookup table kept as a fallback
  + cross-check. Script + 44 synthetic unit tests pass locally
  (`scripts/hla_popgen/tests/test_48_kir_hla_ligand_cooccurrence.py`), including translation/frame
  handling, leader stripping, the C1/C2 and Bw4/Bw6 residue rules pinned to known reference
  alleles, exact/two-field-fallback/unresolved CDS mapping, novel-call handling, seq-vs-lookup
  agreement/disagreement, disclosure masking (a 1-19 cell blanks OR/CI/p but never blanks the
  masked count itself, and a true zero cell is never treated as disclosive), and a
  planted-enrichment detection test.
- Next step: a VM run with the real `~/tools/Immuannot_refdata/CDSseq/` directory present, so the
  Result section and the seq-vs-lookup cross-check can be filled in with real numbers.
