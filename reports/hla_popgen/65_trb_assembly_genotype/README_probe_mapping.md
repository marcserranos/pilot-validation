# TRBV probe mapping: GRCh38 primary vs chr7_KI270803v1_alt

`trb_probe_mapping_grch38_primary_vs_alt.tsv` places the IMGT TRBV *01 V-REGION probes on both GRCh38 versions of the TRB locus (minimap2 2.28 `-x sr`, best hit only). It uses public reference data only, computed on the VM on 2026-10-12.

What it shows: the ALT differs from primary by two separate insertions.

- **Insertion A (~21.6 kb)**:
  - Holds TRBV3-2, TRBV4-3 and the second TRBV6-2/6-3 copy.
  - Sits between TRBV4-2 and TRBV6-2 on primary (about chr7:142,345,986–142,349,378).
  - Every downstream gene shifts by +21.7 kb.
- **Insertion B (about 40 kb further downstream)**:
  - Holds TRBV7-8.
  - Primary has only TRBV7 paralog hits for it.

The 65 probe build (PR #21) uses only insertion A's block for the breakpoint. Its result: chr7:142,346,200–142,352,104, ins_len ≈ 21,666 bp.
