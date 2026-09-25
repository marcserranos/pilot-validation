# 48 — KIR-HLA ligand co-occurrence by ancestry

*Status: script + synthetic tests done locally (S04 WS-D, 2026-09-25/26). Not yet run on real
AoU data — needs a VM session, and the epitope lookup tables below need a VM-side cross-check
before numbers are trusted (see Caveats).*

## Question

Are functional KIR-HLA receptor-ligand pairs (KIR2DL1 x HLA-C2, KIR2DL2/2DL3 x HLA-C1, KIR3DL1
x HLA-Bw4, KIR3DS1 x HLA-Bw4) enriched or depleted relative to independence, within each
ancestry, in the AoU long-read cohort? Long reads give both loci (chr6 HLA, chr19 KIR) from the
same person, which is otherwise hard to get at scale.

## Method

- **Ligand groups**, from a documented two-field allele lookup (NOT a direct position-80
  translation — see Caveats):
  - HLA-C1/C2: `C1C2_TABLE` in the script, citing the standard Asn80(C2)/Lys80(C1) dimorphism
    (Colonna et al. 1993; Winter & Long 1997) at two-field-group resolution.
  - HLA-Bw4/Bw6: `BW4_B_GROUPS` (B*13/27/37/38/44/47/49/51/52/53/57/58/59/63/77) and
    `BW4_A_GROUPS` (A*23/24/32), the standard KIR-ligand-study allele-group list. B*15 defaults
    to Bw6 with a documented, real exception subset not yet applied (`--b15-as-bw4-list`).
  - Alleles outside these tables are `unclassified`, counted (never silently folded into a
    group), and reported in `ligand_lookup_qc.tsv`.
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
- `ligand_lookup_qc.tsv` — how many HLA-C/B/A calls fell outside the lookup tables
  (`n_unclassified_*`); a high fraction here means the epitope-group calls are unreliable for
  that run and should be fixed before trusting the co-occurrence numbers.

## Caveats — read before trusting any number here

- **C1/C2 and Bw4/Bw6 group assignment is a two-field-name lookup table, not a translated
  position-80 amino acid.** It needs a VM-side cross-check against IPD-IMGT/HLA's own
  ligand-group reference (or a position-80 translation from cohort CDS sequences) before these
  numbers are final. The lookup tables are documented in the script's module docstring with
  citations; treat `ligand_lookup_qc.tsv`'s unclassified fraction as the coverage check.
  Alleles not in the tables are excluded from both the numerator and denominator of the affected
  epitope group (not counted as absent).
- **B*15 Bw4 exceptions are not applied by default** — all B*15 alleles are treated as Bw6 unless
  `--b15-as-bw4-list` supplies the documented exception 2-field alleles (a real, known-but-not-
  reproduced-here IPD-IMGT/HLA exception list).
- **Bw4-80I/80T avidity sub-stratification is out of scope** for this delivery (needs the same
  position-80 translation gap noted above).
- **KIR3DS1 has no confirmed direct HLA-binding ligand** — its inclusion here tests the published
  epidemiologic Bw4-80I association (Martin et al. 2002), not a biochemically confirmed
  receptor-ligand pair; interpret differently from the other four rows.
- Cannot distinguish real linkage-disequilibrium/selection from an ancestry confound or a miscall
  artifact on its own — a significant row is a lead to look into, not a finding to publish as-is.

## Distilled

- Not yet run on real data. Script + 14 synthetic unit tests pass locally
  (`scripts/hla_popgen/tests/test_48_kir_hla_ligand_cooccurrence.py`), including disclosure
  masking (a 1-19 cell blanks OR/CI/p but never blanks the masked count itself, and a true zero
  cell is never treated as disclosive) and a planted-enrichment detection test.
- Next step: a VM run, plus verifying the C1/C2/Bw4 lookup tables against IPD-IMGT/HLA before the
  Result section is filled in.
