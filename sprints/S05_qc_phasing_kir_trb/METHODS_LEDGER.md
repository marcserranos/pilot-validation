# S05 methods ledger

The raw material for the final report. One section per analysis. Each section answers, in
plain words: question; why it matters; inputs (exact files and versions); method step by step;
every parameter and why that value; what was run (script, commit, command); results (numbers);
interpretation; what it does not show; caveats. Written as the work happens, not after.

## 58 — Regional long-read depth (MHC, KIR, TRB)

**Question.** Is each person's long-read coverage over the MHC (chr6), the KIR region (chr19) and
the TRB locus (chr7) as good as over the rest of their genome? And, in the TRB locus, does depth
show whole-gene germline deletions?

**Why it matters.** Every HLA and KIR call comes from that person's phased assembly. If the reads
over the region are thin, the assembly breaks or collapses there, and genes look missing or get
phased wrongly. A regional coverage filter was a supervisor ask (2026-10-06): flag people whose
MHC depth is far below their genome-wide depth (e.g. ~8–10× over the MHC at ~17–20× overall). For TRB, Aleix found TRBV4-3 is never used in ~27% of
people's repertoires; if that is a germline deletion, the reads aligned to GRCh38 over TRBV4-3
should drop to about half (one copy deleted) or zero (both).

**Inputs.** Per person: the GRCh38 (noalt) long-read BAM from `v9/wgs/long_read/manifest.tsv`
(column `grch38_bam`, read through the controlled-data automount). Cohort and ancestry:
`~/pipeline_outputs/cohort_membership.tsv`. Gene coordinates: GENCODE basic annotation GTF
(version recorded in `regions.bed` provenance once run).

**Method.**
1. Build one BED of regions: five 2 Mb control windows on chr2/3/5/10/12 (away from the MHC, KIR,
   TRB, immunoglobulin loci, centromeres and the big inversion polymorphisms); the MHC window
   chr6:29.5–33.5 Mb (the exact window the HLA pipeline used); every HLA-region gene in our 65-gene
   list; every KIR gene GENCODE places on chr19; every TRB V/D/J/C gene on chr7; and padded KIR and
   TRB locus windows.
2. For each person, run `samtools bedcov -c` twice: counting all aligned reads (MAPQ ≥ 0) and only
   confidently placed reads (MAPQ ≥ 20). bedcov gives the summed depth per region; dividing by the
   region length gives mean depth.
3. Genome-wide proxy = the median mean depth of the five control windows (median, so one odd
   window does not move it). Relative depth = region mean depth / that proxy. 1.0 = covered like
   the rest of the genome; for a two-copy gene, 0.5 ≈ one copy, 0 ≈ none.
4. TRB copy number per gene: 0 copies if relative depth < 0.25, 1 if 0.25–0.75, 2 above. These
   cut points sit halfway between the expected values (0, 0.5, 1); the full histogram is
   published next to them so the separation can be checked by eye.
5. A TRB gene is only called if reads there are uniquely placed: median (MAPQ≥20 depth / MAPQ≥0
   depth) across people ≥ 0.8. Near-identical paralogs (e.g. TRBV6-2/6-3, TRBV12-3/12-4) fail
   this; their reads are shared between copies, so depth cannot count copies of one of them.
6. Per ancestry: genotype counts, deletion-allele frequency ((2·n0 + n1) / 2n), and an exact
   Hardy–Weinberg test (Wigginton 2005). A real germline biallelic deletion should fit HWE.

**Run.** Pending (VM). Test suite: `tests/test_58_regional_depth.py`, 5/5 on synthetic data.

**Results / interpretation / caveats.** Pending the VM run. Known caveats in advance: KIR genes
that GRCh38 carries only on ALT contigs (most activating KIR) cannot be measured on a noalt BAM;
depth from different platforms (HiFi vs ONT) is comparable as a ratio, not in absolute terms.
