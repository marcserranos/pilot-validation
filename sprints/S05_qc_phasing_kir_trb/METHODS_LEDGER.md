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

## 65 — TRB insertion/deletion genotype from phased assemblies (TRBV4-3, TRBV3-2, TRBV6-3)

**Question.** For each person, which of their two TRB haplotypes carries the ~21 kb segment that
holds TRBV4-3, TRBV3-2 and TRBV6-3? From that: how often each ancestry carries 0, 1 or 2 copies,
the frequency of the insertion allele, and whether the three genes really sit together on one
segment.

**Why it matters.** Aleix found TRBV4-3 is never used in ~27% of people's repertoires (SAS 43.5%,
EUR 16.7%). A known, common germline deletion removes TRBV4-3 together with TRBV3-2 and one copy
of TRBV6-2/6-3 (`READING.md` §B). If someone carries no copy, they cannot use the gene at all, so
the deletion is the first explanation to test before reading the usage gap as immunology. 58
cannot answer this: GRCh38's main chr7 is the haplotype *without* the segment, the segment exists
only on an ALT contig, and 58 reads depth from noalt BAMs. Reads from TRBV4-3 have nowhere to map
there, and GENCODE's primary-assembly gene list (58's regions) does not even contain TRBV4-3. The
phased assemblies do not depend on the reference: if a haplotype carries the segment, its contig
contains the genes.

**Inputs.** Per person: the two hifiasm haplotype assemblies (`assembly_hap1_fa`,
`assembly_hap2_fa`, `.fa.gz`) and their alignments to GRCh38 (`assembly_hap1_aln2_hg38_paf`,
`assembly_hap2_aln2_hg38_paf`) from `v9/wgs/long_read/manifest.tsv`, read through the
controlled-data automount (the same columns 57 uses). Cohort and ancestry:
`~/pipeline_outputs/cohort_membership.tsv`. Gene sequences ("probes"): IMGT/GENE-DB V-REGION
nucleotide sequences for every human TRBV allele (release recorded when built), plus the TRBC1
genomic sequence from Ensembl (GRCh38); or a FASTA given with `--probe-fasta`. Reference check:
the GRCh38 ALT contig that carries the segment (expected `chr7_KI270803v1_alt`, to be confirmed on
the VM from UCSC `altLocations`) and the primary chr7 TRB window.

**Method.**
1. Build the probes once (`--build-probes`): every TRBV allele from IMGT, with IMGT gap dots
   removed, and TRBC1 as genomic sequence (its IMGT C-REGION is spliced, so it would never align
   end to end to genomic DNA). All TRBV genes are kept so that each piece of a contig can be
   credited to the gene that fits it best; only ten genes are reported: TRBV4-3, TRBV3-2, TRBV6-3,
   TRBV7-8 (insertion candidates) and TRBV4-2, TRBV6-2, TRBV3-1, TRBV2, TRBV12-4, TRBC1 (controls
   that every haplotype should carry).
2. Merge genes this method cannot tell apart: if any allele of one gene is ≥ 99% identical to any
   allele of another (`--merge-ident 0.99`), the two become one class. TRBV6-2 and TRBV6-3 share
   the same V-REGION sequence in IMGT, so they become the class `TRBV6-2/TRBV6-3`. For that class
   the script counts copies per haplotype (1 on a deletion haplotype, 2 on an insertion haplotype)
   and adds the feature `TRBV6-2/TRBV6-3[copy2]` ("this haplotype has a second copy"), which stands
   in for TRBV6-3.
3. Check the probes against the reference (`--alt-fasta`, `--ref-fasta`): TRBV4-3 and TRBV3-2 must
   be found on the ALT contig and not on the primary window, and the 6-2/6-3 class must have two
   copies on the ALT and one on the primary. If not, stop (wrong probes or wrong nomenclature).
4. For each haplotype, keep only the contigs that carry TRB: those whose primary alignments in the
   person's own PAF cover ≥ 1 kb of chr7:142,250,000–142,860,000 (the IMGT TRB locus is
   chr7:142,299,011–142,813,287, padded ~50 kb). Each is cut to the outermost projected positions
   ± 50 kb (`--pad`). Outermost, not merged pieces, so the 21 kb segment, which has no place on the
   reference and may be left unaligned in the PAF, stays inside the cut. Without a PAF, contigs are
   found by exact 25-mers taken from the reported genes. The assembly is read once as a stream and
   never indexed on the mount (the mount is read-only; `samtools faidx` there fails, which is why
   57 copied files locally).
5. Map every probe to those contig pieces with minimap2 (in-process through `mappy`, or the
   `minimap2` binary), preset `sr`, base-level alignment (`-c`) and secondary hits kept
   (`--secondary=yes -N 50 -p 0.5`). Why `sr` and not `asm5`: the probes are ~290 bp; `asm5`'s
   scoring (mismatch −19, minimum score 200) cannot report a 300 bp probe with five or more
   mismatches (98.3% identity) and drops weaker secondary hits. Both are needed: the 97% cut must
   be what decides, and the paralog hit must be *seen* in order to be rejected. `sr` (dense 21-mer
   seeds, mismatch −8, secondary ratio 0.5) reports a 95% copy as a secondary hit (checked with
   mappy 2.31 on synthetic sequence). Identity = matching bases / alignment length; coverage =
   aligned part of the probe / probe length. A pure-Python k-mer fallback exists for testing.
6. Group overlapping hits into loci (overlap ≥ 50% of the shorter). Give each locus to the gene
   class whose best probe fits it best: probes covering ≥ 90% first, then highest identity. A locus
   counts as that gene only if identity ≥ 0.97 (`--min-ident`) and coverage ≥ 0.9 (`--min-cov`).
   0.97: alleles of one TRBV gene differ by at most a few bases (> 98.5% identical), while
   neighbouring genes such as TRBV4-2/4-3 differ by more; 0.97 sits between the two. 0.9: allows a
   probe whose ends are clipped by an indel or a contig break while still requiring the gene to be
   essentially whole. If two classes tie within 0.001 (`--tie-eps`), the locus is "ambiguous" and
   counts for neither (reported per haplotype).
7. Haplotype completeness: a haplotype is complete if one contig carries both TRBV2 (5′ end of the
   locus) and TRBC1 (3′ end); the polymorphic segment lies between them, so on a complete
   haplotype "not found" means "not there". Otherwise: `no_trb_contig`, `anchor_missing` or
   `anchors_split` (locus broken across contigs).
8. Per person and gene: copies = present on hap1 + present on hap2; status `ok` only when both
   haplotypes are complete. All genotype counts use `ok` people only.
9. Per ancestry, ALL and platform: genotype counts 0/1/2 copies, insertion-allele frequency
   ((n1 + 2·n2) / 2n), share with 0 copies, and the exact Hardy–Weinberg test (Wigginton 2005, the
   function from 58). Completeness rate per stratum. Co-presence: for each pair of TRBV4-3,
   TRBV3-2 and TRBV6-2/6-3[copy2], on the same complete haplotype, the 2×2 count, concordance, r²
   and D′; plus the joint +/− pattern of the three. If they are one insertion, concordance ≈ 1 and
   only two patterns (all + / all −) should appear.
10. Disclosure as 58: PUBLIC counts of 1–19 people (or haplotypes) are shown as `<20`, and a
    frequency, percentage or HWE p that depends on such a cell is left blank. Per-person tables
    exist only in `~/s05/internal/65/` on the VM. `--unrelated` repeats the tables on the greedy
    unrelated set (kin ≥ 0.0442, the `_s05_common` rule).

**Run.** Pending (VM). Commands: `vm_parts/65.sh`. Probe build:
`$PY 65_trb_assembly_genotype.py --build-probes --out-dir ~/s05/results/65 --alt-fasta ~/ref/trb/trb_alt.fa --ref-fasta ~/ref/trb/trb_primary.fa`;
smoke: `--limit 20 --workers 4 --out-dir ~/s05/results/65_smoke --internal-dir ~/s05/internal/65_smoke`;
full: `setsid nohup $PY -u 65_trb_assembly_genotype.py --probes ~/s05/results/65/probes.fa --workers 12 > ~/s05/logs/65.log 2>&1 < /dev/null & disown`
(`$PY` = `~/repos/pilot-validation/.pixi/envs/spechla/bin/python3`, run from
`~/repos/pilot-validation-main/scripts/hla_popgen`). Test suite:
`tests/test_65_trb_assembly_genotype.py`, 13/13 on synthetic assemblies (planted insertion at 45%
SAS / 60% EUR, a 95% TRBV4-2 paralog, a probe-less 95% TRBV4-3 decoy on every haplotype, a
fragmented haplotype, a haplotype without the locus, a person without a PAF, a missing assembly).

**Results / interpretation.** Pending the VM run. What to look at first: `reference_check.tsv`
(the probes behave as assumed), completeness by platform, TRBV4-3 insertion-allele frequency and
HWE per ancestry, and whether co-presence concordance is ≈ 1.

**What it does not show.** Whether a gene is *expressed* or *used*: that is script 60's join with
Aleix's repertoire calls. Small deletions or gene conversion inside an otherwise present gene
(identity ≥ 0.97 with ≥ 90% coverage counts as present). Copy number of genes outside the ten
reported.

**Caveats.**
- Nomenclature. IMGT names the duplicated pair TRBV6-2/TRBV6-3; some sources name the segment's
  genes differently, and GENCODE (Ensembl) gene names on the ALT contig and on chr7 do not always
  match IMGT. The literature describes the event as removing "TRBV6-2" or "TRBV6-3" depending on
  which copy is called which; here the pair is one class and only the copy count is used. TRBV7-8
  is included because some nomenclatures place it in or next to the polymorphic segment; if
  co-presence shows it on every haplotype, it is not part of the segment. Report IMGT names and
  the release date of the probe set.
- 58 cannot see this segment (noalt reference); 58's TRBV6-2 depth is a mix of both copies and is
  flagged unassessable there. 65 and 58 measure different things and should not be compared gene
  by gene for these four genes.
- A complete haplotype can still be a collapsed assembly; two near-identical copies (6-2/6-3) are
  the most likely to collapse, which would read as a deletion of the copy. TRBV4-3 and TRBV3-2 are
  unique sequence and are the primary readout; `[copy2]` is the secondary one.
- Assemblies may exist only for some platforms (57 only looks for them on Revio and Sequel IIe rows);
  Sequel II people (mostly AFR) and ONT people may be missing or less complete, so completeness is
  reported per platform and frequencies are on complete-both people only.
- Relatives are included unless `--unrelated` is used; for allele frequency and HWE, prefer the
  unrelated run.
- The IMGT query code (`7.2`), the ALT contig name and whether `mappy` installs into the python
  3.8 env were not verifiable from the cloud container (IMGT/Ensembl/UCSC blocked); each is
  checked by a step in `vm_parts/65.sh`.
