# Figure 1 (v2) — what long reads add to the HLA catalogue, and for whom

Built by `scripts/hla_popgen/28_figure1_compose.py` from committed aggregate results (scripts 24,
25, 27). Runs on a laptop; no participant data.

**Caption draft.** *Long-read HLA calling in 12,233 All of Us participants.*
**(a)** Audit trail from the caller's novelty flags to defensible novel proteins: of 280,695 calls
flagged novel, 75,442 are protein-level (field 2), 22,122 of those carry no artifact flag, 5,288
have a coding sequence absent from IPD-IMGT/HLA with no frameshift or premature stop, giving 1,404
distinct novel proteins, of which 231 recur in two or more unrelated people.
**(b)** Composition of haplotype-gene calls per gene: the eight classical genes (left of the dotted
line) are dominated by catalogued alleles and by novelty confined to introns/UTRs (DRB1: 60%),
whereas protein-level novelty is concentrated in the non-classical genes (TAP1, TAP2, MICA).
Hatching marks the share of each class carrying an artifact flag (chiefly homopolymer indels).
**(c)** Cross-ancestry transfer: catalogues of equal size (575 haplotypes) built from one ancestry
group and used to type another, averaged over the eight classical genes at protein resolution.
Transfer is strongly asymmetric — an African-ancestry catalogue types 97% of European haplotypes,
while a European one of the same size types 83% of African haplotypes.
**(d)** Sample coverage against catalogue size for HLA-B: the fraction of haplotypes a catalogue of
that size would already contain. Dotted segments are extrapolated beyond the observed sample.

**Reading note.** Panels a–b answer "what did we find"; c–d answer "what would it take to find the
rest". Coverage (panel d) saturates near 99% while richness does not — most undiscovered alleles are
very rare, so they matter for catalogue completeness but rarely for typing a random person.

Supplement candidates (already committed): admixture barcode and class I ternary
(`06_figures_structure/`, `10_allele_ancestry_geometry/`), non-coding pooling
(`25_noncoding_novelty_paf/fig1_*`), QC decomposition and switch power (`26_qc_relatives_v2/`),
per-gene coverage and design scenarios (`27_allele_space_coverage/`).
