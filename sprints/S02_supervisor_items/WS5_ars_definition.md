# WS5: peptide-contact (ARS) residues, recomputed from crystal structures

Scope: the classical published ARS/PBR residue tables (Hughes & Nei 1988 Table 1; Parham et al.
1988 PNAS; Bondinas et al. 2007) could not be verified from primary sources in this project (see
`WS_literature_selection.md` — Parham's PDF 403'd, Bondinas is paywalled, Hughes & Nei's Table 1
predates a stable IMGT/HLA numbering convention). Rather than cite a number nobody here re-derived,
this note defines peptide-contact residues directly: any MHC heavy atom within a distance cutoff of
any bound-peptide heavy atom, measured on solved crystal structures. Code: `_ars_residues.py`;
tests: `tests/test_ars_residues.py`; data: `reference/ars_peptide_contacts.tsv` (207 rows at the
chosen 4.5 Å cutoff, one row per gene/pdb/residue).

## Structures used

All fetched from `https://files.rcsb.org/download/<ID>.pdb` (not committed — PDB files are not
participant data, but they are also not this project's own derived output, so the brief records
the download command instead of the file). Each header/COMPND was read, not assumed:

| Gene | PDB | Allele (from header) | MHC chain(s) | Peptide chain | Notes |
|---|---|---|---|---|---|
| A | 1HHK | HLA-A\*02:01 | A | C | HTLV-1 Tax 11-19 peptide. Matches the suggested ID. |
| B | 1A1M | HLA-B\*53:01 | A | C | HIV-2 Gag peptide. Suggested 1A1M or 1A9E; used 1A1M (1A9E, B\*35:01, was fetched and verified but not needed once 1A1M worked cleanly). |
| C | 4NT6 | HLA-C\*08:01 | A | C | Influenza matrix-protein-1 peptide. Matches the suggested ID. |
| DRA | 1DLH | DRA / DRB1\*01:01 | A | C | Alpha chain. |
| DRB1 | 1DLH | DRA / DRB1\*01:01 | B | C | Beta chain, same structure. COMPND mislabels chain C as "enterotoxin type B precursor" but SEQRES for chain C is `PKYVKQNTLKLAT`, the influenza HA306-318 peptide the TITLE and the literature both name — verified from the sequence itself, not trusted from COMPND. |
| DQA1 | 1JK8 | DQA1\*03:01 (HLA-DQ8) | A | C | Insulin B-chain peptide. Suggested 1S9V or 1JK8; used 1JK8 because it is a single copy in the asymmetric unit (1S9V, DQ2/DQA1\*05:01+DQB1\*02:01 with a gliadin peptide, has two copies — verified and fetched but not needed). |
| DQB1 | 1JK8 | DQB1\*03:02 (HLA-DQ8) | B | C | Beta chain, same structure. |
| DPA1 | 3LQZ | HLA-DP2 | A | Z (derived, see below) | Matches the suggested ID. |
| DPB1 | 3LQZ | HLA-DP2 | B | Z (derived, see below) | Beta chain, same structure. |

**3LQZ is a covalent single-chain construct**, not three independent chains: chain B's SEQRES is
the peptide+flexible-linker followed directly by the mature beta-chain sequence, with author resnum
running -22..-8 for the peptide/linker then jumping to 3 for the mature beta domain (the linker
itself has no resolved coordinates). There is no separate peptide chain to point the script at, so
the PDB text was preprocessed once (not committed, reproducible): every chain-B `ATOM` record with
resnum < 0 was relabeled to a synthetic chain `Z` before running the same distance code unmodified.
Flagged here rather than silently absorbed into the numbers.

## Cutoff and its sensitivity

Default cutoff is 4.5 Å (any heavy atom to any heavy atom), chosen as the conventional
literature default for "peptide contact," not tuned to hit a target number. To show the exact
value is not load-bearing, the same 9 gene/structure pairs were run at 4.0 Å and 5.0 Å too
(committed only at 4.5 Å; the other two are reproducible from the same script and PDB files):

| Gene | 4.0 Å | 4.5 Å | 5.0 Å |
|---|---|---|---|
| A | 26 | 30 | 32 |
| B | 25 | 28 | 31 |
| C | 26 | 30 | 31 |
| DRA | 19 | 22 | 24 |
| DRB1 | 20 | 21 | 22 |
| DQA1 | 16 | 16 | 21 |
| DQB1 | 21 | 22 | 25 |
| DPA1 | 16 | 17 | 20 |
| DPB1 | 13 | 20 | 22 |

Counts move by roughly ±15-20% per 0.5 Å step, as expected for a distance shell around an
irregular surface. The class I vs class II gap (class I consistently high-20s/low-30s, class II
single-chain counts high-teens/low-20s) is stable across all three cutoffs. DPB1's jump from 13 to
20 between 4.0 and 4.5 Å is the largest single move and merits a second look before leaning on that
count specifically; it does not change the qualitative class I/class II comparison.

## Comparison to Parham (1988)

Parham et al. 1988 (PNAS 85:4005) reports **20 highly variable amino acid positions** clustered in
the class I peptide-binding groove, out of 91 positions surveyed across the then-known class I
sequences (see `WS_literature_selection.md` §2 — the primary PDF still could not be fetched, so
this is the number as reported in that literature-review note, not re-verified against the table
itself). Our 4.5 Å class I contact counts (28-30 residues per gene, A/B/C) are of the same order of
magnitude but **larger and not directly comparable**: Parham's list is a *sequence-variability*
criterion (which positions differ across alleles), ours is a *structural-contact* criterion (which
positions touch the peptide in one particular allele/peptide pair). A residue can be a peptide
contact and be invariant across alleles (contacts the peptide backbone, not side-chain-specific),
and a residue can be highly variable without directly contacting the peptide (TCR-facing, see
below). The two lists should overlap substantially but are not expected to be identical, and this
project cannot confirm the degree of overlap without the literal Parham position numbers.

## What this approach cannot capture

- **One structure is one allele's groove with one bound peptide.** Contacts depend on the peptide
  sequence and register; a different 9-13mer could engage a meaningfully different residue subset,
  especially near the peptide termini and class II's open-ended register. The rows in
  `reference/ars_peptide_contacts.tsv` describe *these specific PDB entries*, not "the ARS of
  HLA-A" as a peptide-independent property.
- **A single allele per gene stands in for the whole gene.** HLA-A\*02:01, B\*53:01, C\*08:01,
  DRB1\*01:01, DQA1\*03:01/DQB1\*03:02, and one DP2 haplotype are each one point in a highly
  polymorphic space; groove geometry and peptide register can shift somewhat between alleles of
  the same gene.
- **TCR-facing residues are excluded by construction but are still under selection.** The
  classical "ARS" concept (Bjorkman/Saper, and the Hughes & Nei dN/dS framework built on it)
  bundles peptide-contact and TCR-contact residues together, because both point up and out of the
  groove and both are hypothesized to be under balancing selection (peptide contacts via
  differential pathogen-peptide binding, TCR contacts via differential T-cell recognition). This
  script only measures the peptide-contact half. A dN/dS or diversity analysis that wants the full
  classical ARS, rather than "peptide contacts only," would need to add TCR-facing groove residues
  from a TCR:pMHC costructure — deliberately not attempted here since that is a different
  co-crystal and a separate verification problem.
- **Author residue numbering, not a codon-alignment position.** `author_resnum` follows each PDB's
  own numbering convention (not guaranteed consistent in signal-peptide handling across
  depositors). `align_and_transfer` in `_ars_residues.py` exists precisely so that using these
  residues against *our own* translated sequences never requires trusting shared resnums — it
  transfers by sequence alignment instead.

## Reproducing

```bash
for id in 1HHK 1A1M 4NT6 1DLH 1JK8 3LQZ; do
  curl -sf "https://files.rcsb.org/download/${id}.pdb" -o "${id}.pdb"
done
# 3LQZ needs the chain-B peptide/linker (author resnum < 0) relabeled to a synthetic chain
# before calling the script -- see "Structures used" above.
python3 scripts/hla_popgen/_ars_residues.py --pdb 1HHK.pdb --gene A --pdb-id 1HHK \
    --mhc-chains A --peptide-chain C --cutoff 4.5 \
    --source "RCSB 1HHK; this project's own contact recomputation" \
    --out reference/ars_peptide_contacts.tsv
# ...repeat per gene/chain pair listed in the table above.
```

Tests: `python3 scripts/hla_popgen/tests/test_ars_residues.py`.
