# 48 — KIR-HLA ligand co-occurrence by ancestry

*Status: run on the full cohort 2026-09-25/26 (S04 WS-D VM session), after fixing a
`GENE_CDS_FILENAME` mismatch (commit `d84348a`): the VM's CDSseq files are named
`HLA-A.fa.gz`/`HLA-B.fa.gz`/`HLA-C.fa.gz`, not `A.fa.gz`/etc. as first assumed, which had been
silently falling back to the lookup table for every call. Verified fixed on the VM: md5 of the
deployed script matches the commit, and the QC counters below show sequence resolution
dominating. Figure added 2026-09-26: `fig_kir_hla_ligand_forest.png` (+ `.pdf`), from
`48b_ligand_figure.py`.*

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

Full cohort: 11,845 people, 30 (ancestry x KIR-gene x ligand) pairs tested (5 pairs x 6
ancestries). `n_unclassified_C=1,093`, `n_unclassified_Bw=1,451` (calls resolved by neither
sequence nor lookup — mostly novel/unmatched alleles, counted not guessed).

**Sequence vs. lookup-fallback source mix** (`ligand_lookup_qc.tsv`) — sequence resolution
dominates, confirming the `GENE_CDS_FILENAME` fix took effect: HLA-C 20,671 sequence-resolved vs.
2,533 lookup-fallback; HLA-B/A (Bw) 41,910 sequence-resolved vs. 4,481 lookup-fallback. Bw4
sub-type split: 9,298 Bw4-80I (higher avidity) vs. 3,016 Bw4-80T. Two-field-group cross-check:
235 groups compared, 228 agree, 7 disagree between the sequence-derived and lookup-table labels.

**Epitope carrier frequency by ancestry** (`epitope_freq_by_ancestry.tsv`):

| ancestry | n | C1 | C2 | Bw4 |
|---|---|---|---|---|
| AFR | 3,009 | 68.8% | 71.9% | 74.0% |
| AMR | 2,656 | 79.0% | 62.3% | 71.7% |
| EAS | 1,460 | 89.5% | 38.0% | 75.5% |
| EUR | 2,974 | 79.1% | 62.8% | 74.6% |
| MID | 487 | 62.0% | 76.6% | 80.3% |
| SAS | 1,236 | 76.1% | 65.5% | 83.3% |

EAS stands out with the lowest C2 frequency (38.0%) and highest C1 (89.5%) — the expected
ancestry-stratified pattern for these epitope groups.

**Receptor-ligand odds ratios by ancestry** (`kir_hla_ligand_cooccurrence.tsv`; `perm_p` from
1,000 shuffles; blank OR/CI/p = a 2x2 cell was <20 and masked):

| ancestry | pair | OR | 95% CI | perm_p |
|---|---|---|---|---|
| AFR | 2DL1xC2 | 1.50 | 0.99-2.29 | 0.098 |
| AFR | 2DL2xC1 | 1.07 | 0.92-1.25 | 0.427 |
| AFR | 2DL3xC1 | 1.12 | 0.90-1.39 | 0.339 |
| AFR | 3DL1xBw4 | 1.38 | 0.83-2.28 | 0.245 |
| AFR | 3DS1xBw4 | 1.00 | 0.80-1.26 | 1.000 |
| AMR | 2DL1xC2 | 1.18 | 0.81-1.72 | 0.432 |
| AMR | 2DL2xC1 | 0.99 | 0.82-1.19 | 0.916 |
| AMR | 2DL3xC1 | 0.76 | 0.55-1.05 | 0.088 |
| AMR | 3DL1xBw4 | 1.32 | 0.96-1.80 | 0.109 |
| AMR | 3DS1xBw4 | 0.99 | 0.83-1.18 | 0.895 |
| EAS | 2DL1xC2 | — (masked, <20 cell) | — | — |
| EAS | 2DL2xC1 | 0.62 | 0.44-0.88 | 0.011 |
| EAS | 2DL3xC1 | — (masked, <20 cell) | — | — |
| EAS | 3DL1xBw4 | 0.97 | 0.60-1.57 | 0.909 |
| EAS | 3DS1xBw4 | 1.01 | 0.79-1.29 | 1.000 |
| EUR | 2DL1xC2 | 1.08 | 0.78-1.48 | 0.654 |
| EUR | 2DL2xC1 | 0.87 | 0.73-1.04 | 0.142 |
| EUR | 2DL3xC1 | 1.43 | 1.12-1.83 | 0.013 |
| EUR | 3DL1xBw4 | 1.01 | 0.73-1.38 | 1.000 |
| EUR | 3DS1xBw4 | 1.12 | 0.95-1.34 | 0.199 |
| MID/SAS | (all 5 pairs) | mostly 0.8-1.2, none significant | — | — (2 pairs masked in MID) |

No pair shows a strong, ancestry-consistent enrichment or depletion. The two nominal signals —
EAS 2DL2xC1 depleted (OR 0.62, perm_p 0.011) and EUR 2DL3xC1 enriched (OR 1.43, perm_p 0.013) —
are each seen in only one ancestry, not replicated across ancestries, and not Bonferroni-corrected
for the 30 pairs tested (0.05/30 ≈ 0.0017, which neither clears); per the script's own framing,
these are leads to look into (possible real LD-with-a-third-factor, selection, or assay artifact),
not a finding to publish as-is.

Full per-ancestry, per-pair table (masked cells included) is in
`kir_hla_ligand_cooccurrence.tsv`.

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

## Figure: `fig_kir_hla_ligand_forest.png` (+ `.pdf`, from `48b_ligand_figure.py`)

- **Panel a** — one forest-plot row per (KIR gene x ligand pair, ancestry), ordered by biology:
  the four inhibitory pairs first (2DL1xC2, 2DL2xC1, 2DL3xC1, 3DL1xBw4), then the one
  activating/epidemiologic pair (3DS1xBw4-80I — tests the Martin et al. 2002 association, not a
  confirmed receptor-ligand bond; see Caveats). Point = odds ratio, whiskers = 95% CI, x-axis is
  log-scaled, dashed vertical line at OR=1 (no association). Y-tick ancestry labels are
  colour-coded to match their point (no legend needed). A `*` marks the two nominal (uncorrected
  p<0.05) signals called out in Result — EAS 2DL2xC1 and EUR 2DL3xC1 — neither survives the
  Bonferroni threshold noted in the panel's own caption (0.05/30 tests ≈ 0.0017).
- **Any row with a hatched grey bar labeled "censored (<20 in a cell)"** is a masked row: one of
  its four 2x2 contingency cells has 1-19 people, so the OR/CI/p-values are blanked upstream in
  `kir_hla_ligand_cooccurrence.tsv` per the disclosure rule. This is drawn as an explicit, visually
  distinct placeholder — never a point, never OR=1, never silently dropped from the row list.
- **Panel b** — a compact carrier-frequency-by-ancestry dot plot for the three epitope groups
  (C1, C2, Bw4), same ancestry colour coding as panel a, for quick visual context alongside the
  ligand-side ORs (e.g. EAS's low C2 frequency visible directly next to its 2DL2xC1 forest row).

## Caveats — read before trusting any number here

- **Sequence-derived assignment depends on `~/tools/Immuannot_refdata/CDSseq/<gene>.fa.gz` being
  present on the run host.** If a gene's CDS file is missing, every call for that gene silently
  (but countably, via `ligand_lookup_qc.tsv`) falls back to the two-field lookup table — check
  `n_{C,Bw}_source_sequence` is non-zero before trusting the run used real sequence, not just the
  fallback.
  - Now run against the real `~/tools/Immuannot_refdata/CDSseq/` directory (IPD-IMGT/HLA 3.55.0,
    per ENVIRONMENT.md quirk #38) — file names on the VM are `HLA-A.fa.gz`/`HLA-B.fa.gz`/
    `HLA-C.fa.gz` (not `A.fa.gz`/etc., the original assumption; fixed in commit `d84348a`).
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

## Plain language, for Marc

KIR receptors on immune (NK) cells recognize specific HLA "flavors" (C1/C2/Bw4) as ligands. We
checked whether people who carry a given KIR receptor gene are more or less likely than chance to
also carry its matching HLA flavor, separately in each ancestry group, using odds ratios (>1 means
"co-occur more than expected"). Short answer: no. None of the five receptor-ligand pairs shows a
consistent pattern across ancestries. Two ancestries had a borderline result (East Asian people
carrying KIR2DL2 slightly less often alongside its ligand than expected; European people carrying
KIR2DL3 slightly more often) but neither holds up once we account for testing 30 things at once —
these are worth a second look, not something to report as a real finding yet.

## Distilled

- Run on the full cohort (11,845 people) after fixing a `GENE_CDS_FILENAME` mismatch that had
  been silently forcing every call through the lookup-table fallback (commit `d84348a`). Sequence
  resolution now dominates (20,671/23,204 HLA-C calls, 41,910/46,391 HLA-B/A calls resolved by
  sequence, not lookup).
- No receptor-ligand pair shows an ancestry-consistent enrichment/depletion; the two nominal
  single-ancestry signals (EAS 2DL2xC1 depleted, EUR 2DL3xC1 enriched) don't survive a multiple-
  testing correction for the 30 pairs tested — leads, not findings.
- Figure: `fig_kir_hla_ligand_forest.png` (forest plot + epitope frequency dot plot). Next step: a
  closer look at the 7 two-field groups where the sequence-derived and lookup-table labels
  disagree (`ligand_seq_vs_lookup_crosscheck.tsv`).
