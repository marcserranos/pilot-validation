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

## 65 — TRB insertion/deletion genotype (TRBV4-3, TRBV3-2, TRBV7-8, TRBV6-3): long reads, validated on assemblies

**Question.** How many copies of the ~21 kb segment that holds TRBV4-3, TRBV3-2, TRBV7-8 and
TRBV6-3 (0, 1 or 2) does each person carry? From that: the frequency of the insertion allele in
each ancestry, and whether the genes really sit together on one segment.

**Why it matters.** Aleix found TRBV4-3 is never used in ~27% of people's repertoires (SAS 43.5%,
EUR 16.7%). A known, common germline deletion removes TRBV4-3 together with TRBV3-2 and one copy
of TRBV6-2/6-3 (`READING.md` §B). If someone carries no copy, they cannot use the gene at all, so
the deletion is the first explanation to test before reading the usage gap as immunology. 58
cannot answer this: GRCh38's main chr7 is the haplotype *without* the segment, the segment exists
only on the ALT contig `chr7_KI270803v1_alt`, and 58 reads depth from noalt BAMs. Reads from
TRBV4-3 have nowhere to map there, and GENCODE's primary-assembly gene list (58's regions) does not
even contain TRBV4-3.

**Two routes.** *Reads* (default, fast): count long reads that carry the segment against reads that
span the empty insertion site on chr7. *Assemblies* (validation): look for the genes in each
person's two phased haplotype assemblies. The assembly route was built first. On the VM it turned
out to be I/O-bound: 20–112 s per haplotype to stream a ~1 GB `.fa.gz` from the mount, which is
2–4 days for the cohort on 4 vCPUs. So reads mode runs on everyone, assembly mode runs on a
stratified 2,000, and 65b measures how often the two agree.

**Inputs.** Per person, from `v9/wgs/long_read/manifest.tsv` through the controlled-data automount:
the GRCh38 long-read BAM (`grch38_bam`, or `grch38_haplotagged_bam` when it carries HP tags), and,
for assembly mode, the two hifiasm haplotype assemblies (`assembly_hap{1,2}_fa`) and their GRCh38
alignments (`assembly_hap{1,2}_aln2_hg38_paf`). Cohort and ancestry:
`~/pipeline_outputs/cohort_membership.tsv`; relatedness: the AoU `samples_relatedness.tsv`. Gene
sequences ("probes"): IMGT/GENE-DB V-REGION nucleotide sequences for every human TRBV allele
(GENElect query 7.2; 240 alleles on 2026-10-11), plus TRBC1 genomic sequence from Ensembl (GRCh38).
Reference: UCSC hg38 `chr7_KI270803v1_alt` (the ALT, placed at chr7:142,038,120–143,088,503) and
UCSC hg38 `chr7.fa.gz` cut to that same span.

**Method — probes and breakpoint (once).**
1. Probes: every TRBV allele from IMGT, gap dots removed, and TRBC1 as genomic sequence (its IMGT
   C-REGION is spliced, so it would never align end to end to genomic DNA). All TRBV genes are kept
   so that each piece of sequence can be credited to the gene that fits it best. Only these genes are
   reported: TRBV4-3, TRBV3-2, TRBV6-3, TRBV7-8 (segment) and TRBV4-2, TRBV6-2, TRBV3-1, TRBV2,
   TRBV12-4, TRBC1 (controls every haplotype carries).
2. Genes this method cannot tell apart are merged: if any allele of one gene is ≥ 99% identical to
   any allele of another (`--merge-ident 0.99`), they become one class. TRBV6-2 and TRBV6-3 share
   their V-REGION sequence in IMGT and become `TRBV6-2/TRBV6-3`.
3. Reference check (`--alt-fasta`, `--ref-fasta`, `--ref-region`): TRBV4-3, TRBV3-2 and TRBV7-8 must
   be found on the ALT and not on primary, and the 6-2/6-3 class must have two copies on the ALT and
   one on primary.
4. Breakpoint: 1 kb pieces of the ALT, one every 250 bp, over 40 kb on each side of the segment's
   genes, are aligned to the primary span (minimap2 asm5). Pieces with a single forward hit (≥ 98%
   identity over ≥ 90% of the piece) give a shift between ALT and primary coordinates. The left side
   has one shift and the right side another; the difference is the segment length. The breakpoint
   is the interval between the last hit that follows the left shift and the first that follows the
   right shift. An interval rather than a point absorbs microhomology or a duplicated unit at the
   junction. Written to `probe_meta.tsv`; `--breakpoint` overrides it.

**Method — reads mode (everyone).**
5. For each person, `samtools view -F 0x904` (primary mapped reads only) of the GRCh38 BAM over the
   primary TRB span ± 50 kb, chr7:142,249,011–142,863,287 (`--reads-region`).
6. Reference-allele read: aligned without a break across [breakpoint left − 1,000, breakpoint right
   + 1,000] (`--flank 1000`), with no insertion, deletion or skip ≥ 50 bp (`--max-indel`) inside that
   interval. Only a haplotype without the segment produces such a read. A read from a haplotype with
   the segment either carries it as one huge insertion (if the read is longer than the segment) or is
   clipped at the breakpoint. 1,000 bp on each side gives the aligner enough anchored sequence that a
   clipped read cannot pass. 50 bp tolerates ordinary indels and HiFi homopolymer errors.
7. Insertion-allele read: the probes are aligned to every other read (one minimap2 `sr` index over the
   reads, secondary hits kept, `-N 5000`). The read counts as insertion if a TRBV4-3, TRBV3-2 or
   TRBV7-8 locus on it is called under the same best-hit rule as assembly mode (step 11). Counts are
   kept per gene and for `insertion_any`. TRBV6-3 cannot be used: it is identical to TRBV6-2, which
   both alleles carry.
8. Allele fraction AF = n_ins / (n_ins + n_ref). Copies: 0 if AF < 0.2, 1 if 0.2 ≤ AF < 0.8, 2 if
   AF ≥ 0.8 (`--af-cuts`). These cuts sit far from the expected 0, ~0.5 and 1, so a het call needs at
   least ~1 read in 5 from each allele. A call needs ≥ 6 informative reads (`--min-reads`, the
   agreed default), else `low_reads`. Note the cost of that floor: if a het's reads were drawn 50:50,
   the chance of AF falling outside 0.2–0.8 is 22% with 6 reads, 7% with 8–10, 4% with 12 and 2%
   with 15 (binomial). So if the median informative read count allows it, rebuild the tables with
   `--aggregate-only --min-reads 12`. PUBLIC tables recompute status and copies from the stored read
   counts. The INTERNAL per-person table (which 60 reads) keeps the floor used at run time.
   The published AF histogram shows whether the three modes really sit where the cuts assume.
9. With HP tags (haplotagged BAM): per-haplotype calls. A haplotype is "complete" with ≥ 3
   informative reads (`--min-hap-reads`) and carries the gene if its own AF ≥ 0.5. These calls feed
   the co-presence tables.

**Method — assembly mode (validation subset).**
10. TRB contigs of each haplotype = contigs whose primary PAF alignments cover ≥ 1 kb of
    chr7:142,250,000–142,860,000. Each is cut to the outermost projected positions ± 50 kb, so a
    segment that the PAF leaves unaligned stays inside. The `.fa.gz` is read once as a stream, never
    indexed on the read-only mount.
11. Probes are aligned to those pieces (minimap2 `sr`, `-c`, secondaries kept). Why `sr` and not
    `asm5`: a ~300 bp probe with ≥ 5 mismatches (98.3%) is below asm5's minimum score, and the
    paralog hit must be seen in order to be rejected. Overlapping hits form a locus. The locus goes
    to the class whose best probe covers ≥ 90% and has the highest identity. It is called if
    identity ≥ 0.97 (`--min-ident`; alleles of one gene are > 98.5% identical, neighbours such as
    TRBV4-2/4-3 less) and coverage ≥ 0.9 (`--min-cov`). Ties within 0.001 are "ambiguous" and count
    for neither class.
12. Completeness (relaxed after the VM pilot, where the locus was usually split over 2–3 contigs and
    only 2 of 9 haplotypes passed the one-contig rule): presence is called on *any* of the haplotype's
    TRB contigs. A haplotype is complete if TRBV2 (5′ end) and TRBC1 (3′ end) are both called among
    them. The strict rule (both on one contig) is kept as `complete_single_contig`. For the
    6-2/6-3 class, the number of copies per haplotype is reported, and `TRBV6-2/TRBV6-3[copy2]`
    (≥ 2 copies) stands in for TRBV6-3.
13. copies = present on hap1 + present on hap2; status `ok` only if both haplotypes are complete.
14. Subset (`--sample-stratified 2000`): people with both assembly paths, unrelated (greedy removal
    at kin ≥ 0.0442, the `_s05_common` rule), equal numbers per ancestry. A smaller group gives all
    it has and its unused share is split among the others. Seed 65. The list is INTERNAL
    (`--people-file`).

**Method — tables (both modes) and 65b.**
15. Per feature × ancestry, ALL and platform: people measured, people called, genotype counts 0/1/2,
    insertion-allele frequency (n1 + 2·n2) / 2n, % with 0 copies, and the exact Hardy–Weinberg test
    (Wigginton 2005, 58's function). Haplotype completeness and reasons. Co-presence of the segment
    genes on the same complete haplotype: 2×2 counts, concordance, r², D′, and the joint +/− pattern.
    If they are one insertion, concordance ≈ 1 and only all-plus or all-minus patterns appear.
    Reads mode adds the AF histogram (0.05-wide bins).
16. 65b (`--compare reads.tsv assembly.tsv`): for people called in both, the % with the same copies
    (Wilson CI) and the 3×3 confusion matrix, per feature and stratum.
17. Disclosure as 58: in PUBLIC tables, counts of 1–19 people or haplotypes become `<20`, and any
    frequency, percentage or HWE p that depends on such a cell is left blank. Per-person tables exist
    only under `~/s05/internal/65_*` on the VM. An internal directory records its mode, and a run in
    another mode or with an older column layout refuses to write there.

**Run.** Done on the VM (2026-10-11, worktree `~/repos/pv-65`): probe build and reference check,
plus an assembly-mode pilot (9 haplotypes, which exposed the speed and completeness problems).
Pending: probe rebuild with the breakpoint, reads smoke (20), assembly on the stratified 2,000,
reads on everyone, 65b. Exact commands are in `vm_parts/65.sh`. Full run:
`setsid nohup $PY -u 65_trb_assembly_genotype.py --mode reads --aligner minimap2 --minimap2 $ENVBIN/minimap2 --samtools $ENVBIN/samtools --workers 12 > ~/s05/logs/65_reads.log 2>&1 < /dev/null & disown`.
Tests: `tests/test_65_trb_assembly_genotype.py`, 23/23. Synthetic assemblies cover a planted
insertion at 45% SAS / 60% EUR, a 95% TRBV4-2 paralog, a probe-less 95% TRBV4-3 decoy, a
fragmented haplotype, a missing locus, no PAF and a missing assembly. Synthetic reads cover a fixed
primary/ALT pair, reads placed as an aligner would (clipped or with one big insertion), 90% HP-tagged,
and a missing BAM. Also tested: breakpoint recovery within 300 bp and segment length within 50 bp,
stratified sampling, and 65b.

**Results.** Reference check (VM, 10-11): PASS. TRBV4-3, TRBV3-2 and TRBV7-8 are on
`chr7_KI270803v1_alt` only; TRBV6-2/6-3 has 2 copies on the ALT and 1 on primary. This also settles
that TRBV7-8 is in the segment in IMGT naming. Genotype results are pending.

**What it does not show.** Whether a gene is expressed or used (that is 60's join with Aleix's
repertoire calls). Small deletions or gene conversion inside an otherwise present gene. Copy number
of genes outside those reported.

**Caveats.**
- Reads-mode allele fraction is not exactly 0.5 for a het. Insertion reads are reads that contain a
  segment gene; reference reads are reads spanning the site ± 1 kb. The two depend differently on
  read length and on where each gene sits in the segment, and reads lying wholly inside the 21 kb
  segment are unmapped or misplaced on noalt GRCh38 and mostly missed. The cuts are wide for this
  reason; the AF histogram and the 65b concordance are the checks. ONT reads (lower identity) may
  miss the 0.97 cut more often, so look at results per platform.
- Without HP tags there are no per-haplotype calls, and so no reads-mode co-presence (reason
  `no_hp_tag`); use `--bam-column grch38_haplotagged_bam` if the plain BAM has none.
- Assembly mode: relaxed completeness accepts a locus split over contigs. If a contig break falls
  inside the segment, part of it can be missing; `complete_single_contig` gives the strict subset to
  check against. A collapsed assembly can merge the 6-2/6-3 copies, so `[copy2]` is the weaker
  readout; TRBV4-3 and TRBV3-2 are unique sequence.
- Nomenclature. IMGT's TRBV6-2/6-3 pair is one class here; the literature names "the deleted copy"
  either way. GENCODE/Ensembl names on the ALT and on chr7 do not always match IMGT. Report IMGT
  names and the probe-set date.
- 58 cannot see this segment (noalt reference). Its TRBV6-2 depth mixes both copies and is flagged
  there as not assessable. 58 and 65 should not be compared gene by gene for these genes.
- Relatives are included unless `--unrelated` is used; for allele frequency and HWE, prefer the
  unrelated tables.

## 59 — Does low regional depth explain the calls we distrust? (depth × outcome join, coverage cutoff)

**Question.** People whose long-read coverage over the MHC or the KIR region is thin may get
broken or collapsed assemblies there. Does low regional depth (from 58) explain (a) incomplete KIR
haplotypes and missing KIR framework genes, (b) structurally anomalous DQ contigs and DQA1~DQB1
pairs that cross the G1/G2 groups, (c) absences of HLA genes that are on every normal haplotype,
and (d) MHC contig breaks? If it does, at what depth should a person be filtered?

**Why it matters.** A coverage filter is only worth applying if it removes the technical calls
and not many people. 50, 51 and 57b already showed that most suspicious HLA/KIR calls are
assembly artefacts (contig breaks, open-ended KIR contigs, anomalous DQ contigs). This analysis
asks whether those artefacts concentrate in thinly covered people, and turns the answer into one
cutoff chosen by a stated rule, rather than by eye. If depth does not explain an outcome, the
filter should not be sold as fixing it.

**Inputs.**
- 58's INTERNAL `~/s05/internal/58/per_person_depth.tsv` (columns `58.PER_PERSON_COLS`: person_id,
  region, mapq, mean_depth, n_reads, control_depth, rel_depth). Only the windows `window_MHC` and
  `window_KIR` are used, at MAPQ ≥ 0 and MAPQ ≥ 20.
- Table 1 `~/pipeline_outputs/hla_calls_rich.tsv` (HLA calls per person × hap × contig × gene).
- KIR annotation `~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz` (the files 51 read).
- `~/pipeline_outputs/cohort_membership.tsv` (ancestry, platform) and the AoU relatedness table
  `v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv`.

**Method.**
1. Depth per person: rel_depth over `window_MHC` and `window_KIR` at both MAPQ floors. Windows only:
   single genes are too noisy for a person-level filter (58 smoke test, n = 20: gene-level p5 of
   0.4–0.6 in ordinary people).
2. Two depth scales, both reported:
   - `abs` = rel_depth as 58 defines it (region depth / the person's genome-wide depth);
   - `norm` = rel_depth / the cohort median rel_depth of the same region, MAPQ and platform.
   The 58 smoke test showed regions run below genome depth even in ordinary people (MHC window
   median 0.825, p5–p95 0.76–0.91; KIR window 0.86, 0.74–0.97), so on `abs` a fixed cutoff mixes
   "this person is thin" with "this region is always thinner". On `norm`, 1.0 = typical for the
   region on that platform. Medians come from everyone with depth, before the kinship filter; a
   platform with < 20 people uses the all-platform median.
3. Cohort: people in the depth table, unrelated only (greedy removal at kin ≥ 0.0442, the rule used
   across hla_popgen, via `_s05_common.greedy_unrelated`), so the regression treats people as
   independent.
4. Outcomes, one 1/0 flag per person ("on any haplotype"), blank when not assessable:
   - KIR, using 51's own functions (`person_records`, `classify_hap`, `FRAMEWORK`), imported, not
     re-written: `kir_incomplete` (a haplotype that is not 51-complete: one contig carrying all KIR
     calls with both KIR3DL3 and KIR3DL2); `kir_fw_missing` (a framework gene 3DL3/2DL4/3DP1/3DL2
     missing); `kir_fw_missing_complete` (the same on complete haplotypes only). 51 found the
     latter are mostly contracted haplotypes, i.e. biology, so it serves as a negative control and
     is never used to pick the cutoff. Assessable = both haplotype GTFs exist.
   - DQ, using 50's `dq_group` (Petersdorf G1/G2 from field 1): `dq_split` (DQA1 and DQB1 on
     different contigs); `dq_anomalous` (one contig, DQA1–DQB1 start-to-start spacing ≥ 100 kb;
     50e found 98% of same-contig haplotypes at 16.6–24.2 kb and the ≥ 100 kb ones are the
     structurally anomalous contigs); `dq_cross_group_any` (a G1/G2 cross pairing, any contig
     layout); `dq_cross_group_physical` (the same on one contig only: the 57b physical-evidence
     rule, under which only a same-contig pair counts as a physical haplotype).
   - HLA absences: `hla_absent` = one of A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1 absent from a
     haplotype; `hla_hap_empty` = a haplotype with no HLA call. DRB3/4/5 are excluded because
     their absence is real copy-number variation set by the DRB1 lineage.
   - MHC contig breaks, using 50's `hap_gene_contigs` (50b's per-haplotype contig count):
     `mhc_contigs_mean` (continuous, mean over the person's haplotypes) and `mhc_fragmented`
     (a haplotype with ≥ 7 contigs; 50b's cohort mean is 4.3 on Revio and 4.9 on Sequel IIe, so 7
     is roughly 1.5× the mean).
5. Each outcome is judged against its own region's depth: KIR outcomes against `window_KIR`, the
   rest against `window_MHC`.
6. Rate by depth bin (`abs` edges 0.4…1.0, `norm` edges 0.5…1.1), per ALL, ancestry and platform.
7. Cutoff scan: for each cutoff c (`abs` 0.30–0.90, `norm` 0.50–0.95, step 0.05), people with depth
   < c are removed; the table gives people removed overall and per ancestry, outcome rate among
   kept vs removed, share of the outcome's events removed, and J = share of events removed − share
   of people removed.
8. Logistic regression per outcome × scale × MAPQ: outcome ~ depth (capped at 1.5) + platform +
   ancestry; OR per 0.1 lower depth with Wald 95% CI. statsmodels when available, otherwise a
   Newton–Raphson fit that matches statsmodels in the tests (the VM `spechla` env is Python 3.8
   without statsmodels; the engine used is written in each row). Guards: no model with < 5 events
   or non-events; covariate levels with < 5 events or non-events are pooled into "other" (a level
   with no events is perfectly separated and stops the fit); PUBLIC shows no OR when events or
   non-events are 1–19. `mhc_contigs_mean` gets the matching OLS slope.
9. Proposed cutoff, applied to exact numbers, on the `norm` scale by default (`--proposal-scale`),
   both scales reported: scan c upward and take the **smallest c that removes ≥ 50% of the
   outcome's events (X) while removing ≤ 10% of people (Y) overall and in every ancestry with
   ≥ 100 people (M)**. If no c qualifies, the outcome is reported as not explained by depth, with
   how many events the largest admissible c would remove and which ancestry hits the people limit
   first. Per region and MAPQ, the region proposal is the largest proposed c over that region's
   outcomes, excluding the negative control; removal only grows with c, so it still meets the
   people limit everywhere. Smallest-first keeps the filter as light as the rule allows. X, Y and M
   are CLI parameters and the rule text is written into every output row.
10. Disclosure: every participant count goes through `_disclosure` (1–19 → `<20`), every rate built
    on such a cell is blank, and where a row publishes a total and one part of it, the other part is
    blanked when it is 1–19 (complementary suppression within a row).

**Run.** Pending (VM; needs 58's full run first). Test suite: `tests/test_59_depth_outcome_join.py`,
13/13 on synthetic data (planted step at normalised depth 0.55 recovered within one grid step for
every outcome; per-person outcomes match the plant exactly through 50/51's rules; platform offset
removed by the `norm` scale; negative-control role; numpy fit matches statsmodels; masking; lock;
resume; INTERNAL refuses repo paths). Exact VM commands: `vm_parts/59.sh`.

```bash
cd ~/repos/pilot-validation
setsid nohup pixi run -e spechla -- python3 -u scripts/hla_popgen/59_depth_outcome_join.py \
    --out-dir ~/s05/results/59 --internal-dir ~/s05/internal/59 --workers 8 \
    > ~/s05/logs/59.log 2>&1 < /dev/null & disown
```

**Results / interpretation.** Pending the VM run.

**Caveats (known in advance).**
- Association, not cause: a person with a hard-to-assemble MHC haplotype (for example a large
  structural variant) may also have lower depth there because reads map poorly to GRCh38. Depth and
  the outcome then share a cause.
- `norm` is relative to the cohort: if a whole platform were thin, it would look normal on `norm`
  (the `abs` scale and the per-platform medians in `run_summary.tsv` show this).
- Only the windows are used. A deletion of one gene with normal window depth is not caught, by design.
- `hla_absent` mixes true absence with breaks. 52's atlas and 30's bridged-absence logic separate
  them further; this analysis only asks whether depth explains the total.
- The ≥ 7-contig fragmentation threshold and the 100 kb DQ spacing are fixed choices (CLI
  parameters), taken from 50's distributions, not tuned here.
- Complementary suppression is applied within each row, not across rows (for example ALL minus
  the other ancestries). This matches 58's standard. Cross-row disclosure stays with the open
  small-cell question in `context/DECISIONS.md`.
- statsmodels is absent from the VM env, so the fallback fit is used unless it is installed.

**Exclusion list.** 59 also writes INTERNAL `exclude_people.tsv` (one row per person per window, region proposal on the normalised scale at MAPQ 0; `exclude` = below the cutoff in any selected window). 61 drops excluded people; 62 takes its coverage flag from the window_MHC rows.

## 60 — TRBV germline copy number vs "never used" in the repertoire

**Question.** Some TRBV genes are never used in a sizeable share of people's T-cell repertoires
(TRBV4-3 is the clearest case in the W1 V-usage analysis). Is "never used" explained by having no
germline copy of the gene? How well does copy number 0 predict it, by ancestry? Does usage scale
with copy number (1 vs 2 copies)? And how much "never used" is left among people who do carry the
gene?

**Why it matters.** Differences in V usage between ancestries are a candidate signal of HLA or
selection effects. A germline deletion whose frequency differs by ancestry explains the same
pattern with no HLA involved, so it must be measured and removed before any usage difference is
interpreted. The deletion/insertion polymorphism spanning TRBV4-3, TRBV3-2 and TRBV6-2/6-3 is
known prior art (`READING.md` §B: repertoire inference in Omer 2022, long-read assembly in the
Rodriguez/Watson preprint, short-read genotyping in a 2025 preprint). What has not been done is
ancestry-stratified genotyping from long reads in a large, diverse cohort, linked to each person's
own repertoire. The residual (never used despite a copy) points to non-functional alleles or
regulation.

**Inputs.**
- Copy number, one of two sources (`--cn-table`, `--cn-source auto|58|65`):
  - 58's INTERNAL `~/s05/internal/58/trb_cn_per_person.tsv` (person_id, gene, rel_depth, copies,
    assessable). Depth over TRB genes on GRCh38 primary chr7; copies 0/1/2 at rel_depth cuts 0.25 /
    0.75; genes whose reads are not uniquely placed (median MAPQ20/MAPQ0 depth < 0.8) flagged not
    assessable.
  - 65's INTERNAL `~/s05/internal/65/per_person_gene.tsv` (person_id, feature, hap1/hap2_present,
    copies, status, method): presence per haplotype called from each person's assemblies.
- Usage: Aleix's per-person TRBV usage table. Path and format were not known when the script was
  written. The W1 code caches person × gene clonotype counts as `analysis/vj_usage.pkl`
  (`{"V": counts, "J": counts}`, index research_id). The script accepts that pickle, a pickled
  table, or a TSV/CSV/parquet in long or wide form. Without `--usage` it searches
  `~/workspace/aleix_rnaseq/` (depth ≤ 4) for `*usage*` and stops with a list of candidates unless
  exactly one file matches.
- `~/pipeline_outputs/cohort_membership.tsv` (ancestry) and the AoU relatedness table.

**Method.**
1. Copy number.
   - 58 source: GRCh38 primary chr7 lacks TRBV4-3, TRBV3-2, TRBV6-3 and TRBV7-8. They sit on the
     ~21 kb insertion, which GRCh38 carries only on an ALT contig, so 58 cannot measure them
     directly. A proxy is added and labelled as such (`cn_method = proxy_58_paralog_excess`). Reads
     from the insertion's copies map at MAPQ 0 onto their primary paralogs TRBV4-2, TRBV6-2 and
     TRBV3-1 and lift their depth: no insertion ≈ 1.0, one copy ≈ 1.5, two copies ≈ 2.0.
     Per person, proxy_rel = the median rel_depth over those three paralogs (the median damps
     single-gene noise). Insertion present = proxy_rel ≥ 1.5. Proxy copies of each insertion gene
     = 0 / 1 / 2 at cuts 1.25 / 1.75, halfway between the expected values. All four insertion genes
     get the same proxy call, so their concordance is identical by construction and is labelled so.
   - 65 source: copies = hap1_present + hap2_present. Rows with status ≠ ok are dropped: an
     incomplete TRB locus on a haplotype makes copies only a lower bound. Where 65 merges
     near-identical genes into one class (e.g. TRBV6-2/6-3), usage is summed over the member genes
     so both sides refer to the same genes. Derived `[copy2]` features (65's insertion proxy) are
     used only in the haplotype table.
2. Usage. Gene names lose allele suffixes. A gene is "never used" by a person when their count is 0.
   With counts, a zero is only judged where a chance zero is unlikely:
   P(chance zero) = (1 − p_g)^N < 0.01, where p_g = the gene's pooled usage over all people and
   N = the person's total TRBV count. This is the same estimator as the W1 zero-usage analysis.
   Pooled p_g is conservative, since people lacking the gene only lower it. With fractions instead
   of counts, totals are unknown and the filter is skipped (logged in `run_summary.tsv`).
3. Cohort: people in both tables, unrelated only (greedy, kin ≥ 0.0442), joined to ancestry. The
   usage id is assumed to be the AoU person_id.
4. Per gene × ancestry (and ALL):
   - contingency: never used × copies 0/1/2;
   - CN0 as a test for "never used": sensitivity = share of never-used people with CN0 (how much is
     explained); specificity; PPV = share of CN0 people who never use it (should be ~100% for a
     real deletion); residual = never-used rate among CN ≥ 1, shown next to the chance-zero rate
     expected there;
   - dosage: usage fraction quantiles by copies; Kendall tau-b of usage fraction on copies, over
     all people and over CN ≥ 1 only. The second asks whether 2 copies give more usage than 1,
     beyond presence/absence.
5. Haplotype genes (TRBV4-3, TRBV3-2, TRBV6-3, TRBV7-8, TRBV6-2, TRBV12-4): pairwise copy-number
   concordance (3 × 3 counts, % identical, Cohen's kappa) and never-used concordance
   (P(B never used | A never used)). If the deletion calls are one haplotype, the insertion genes
   agree almost perfectly and TRBV12-4 (a separate locus) does not. With the 58 proxy this is
   true by construction; with 65 it is a real test.
6. Disclosure: participant counts through `_disclosure` (1–19 → `<20`). A rate is blank when its
   numerator, denominator or complement is 1–19. When one half of a never-used/used split is 1–19,
   the other half is blanked.

**Run.** Pending (VM; needs 58's or 65's INTERNAL table and the usage table's path). Test suite:
`tests/test_60_trb_germline_usage.py`, 10/10 on synthetic data (planted insertion haplotype by
ancestry, dosage, 4% residual; proxy copies recovered > 97%; PPV > 97%; residual within 1.6
points; 65 layout with merged classes and incomplete rows; pickle and long-TSV formats give
identical results; not-found message lists candidates; masking; INTERNAL refuses repo paths).
Exact VM commands: `vm_parts/60.sh`.

```bash
cd ~/repos/pilot-validation
pixi run -e spechla -- python3 -u scripts/hla_popgen/60_trb_germline_usage.py \
    --cn-table ~/s05/internal/58/trb_cn_per_person.tsv --usage <path> \
    --out-dir ~/s05/results/60 --internal-dir ~/s05/internal/60 > ~/s05/logs/60.log 2>&1
```

**Results / interpretation.** Pending the VM run.

**Caveats (known in advance).**
- The 58 proxy is indirect: any duplication or mapping artefact that raises depth over TRBV4-2,
  TRBV6-2 or TRBV3-1 reads as an insertion. Under the proxy, the insertion genes agree by
  construction. 65's direct assembly calls are the real test; the proxy is a fallback and a
  cross-check.
- TRBV6-2 and TRBV6-3 share their V-region sequence. Neither depth (58 flags TRBV6-2 as not
  assessable when reads are not uniquely placed) nor 65's assembly calls can count them
  separately, and repertoire alignment may not separate their usage either.
- "Never used" comes from RNA-seq repertoires of limited depth. The chance-zero filter removes
  low-depth people for rare genes, so rare genes may be unassessable.
- Usage counts are unique clonotypes (the W1 unit), not reads.
- The usage table's id is assumed to be the AoU person_id. If Aleix's table uses another id,
  the join is empty and `run_summary.tsv` shows 0 people in both tables.
- Complementary suppression is within a row only (as in 58); cross-row disclosure stays with the
  open small-cell question in `context/DECISIONS.md`.

## 61 — KIR gene content and alleles by ancestry (Figure 1 panel)

**Question.** How do KIR gene content, cA/cB haplotype share, the number of genes per haplotype
and the commonest alleles differ between ancestry groups, when only fully assembled KIR
haplotypes are counted?

**Why it matters.** This is the KIR panel of Figure 1. 43 reported gene presence on all
haplotypes, but 25% of KIR haplotypes are fragmented or open-ended (51). A fragmented haplotype
loses genes, so gene frequencies from all haplotypes are biased downwards, and the bias could
differ by ancestry or by haplotype content. The proposed S05 rule (DECISIONS.md) is to make
KIR haplotype-level statements only on complete haplotypes.

**Inputs.** VM mode: per-person Immuannot KIR GTFs
`~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz` (the 43 run: 12,261 people,
Tier-3 `sequel2` excluded), `~/pipeline_outputs/cohort_membership.tsv` (predicted ancestry), the
AoU v9 relatedness table, and optionally the 59 exclusion list. Local draft: the committed PUBLIC
tables `reports/hla_popgen/51_kir_structure_qc/gene_content_haplotypes.tsv`,
`cB_share_by_completeness.tsv` and `reports/hla_popgen/43_kir_full_cohort/kir_allele_freq.tsv`,
`kir_gene_by_ancestry.tsv`.

**Method.**
1. Parse each haplotype's KIR calls with 44's contig-aware parser (`parse_hap_gtf_full`, which
   reuses 41's novelty tiers).
2. Classify each haplotype with 51's rule (`classify_hap`): complete = all KIR calls on one
   contig, and KIR3DL3 and KIR3DL2 both present.
3. Keep unrelated people (greedy removal at kinship ≥ 0.0442, the rule used across hla_popgen),
   minus anyone flagged `exclude = True` in 59's INTERNAL `exclude_people.tsv`
   (`--exclude-people`; by default below 59's proposal in the MHC or the KIR window).
4. Per ancestry, on complete haplotypes only:
   - (a) % of haplotypes carrying each of the 17 genes, plus person-level carrier frequency
     among people whose two haplotypes are both complete;
   - (b) top alleles per gene;
   - (c) cB share (any of 2DL2, 2DL5A/B, 2DS1/2/3/5, 3DS1 on the haplotype);
   - (d) the distribution of distinct genes per haplotype.
   - Also, per ancestry, the share of haplotypes in each completeness class, to check whether
     the rule itself is ancestry-neutral.
5. Allele resolution for (b):
   - default is the protein field (first three digits). Calls that keep a known protein (known,
     genomic-only novel, synonymous novel) take their nearest allele's protein field; protein-novel
     calls form one pool.
   - A second table gives full names of exact (known) matches only, which is 43's convention.
6. Disclosure:
   - counts 1–19 are shown as `<20` and their rates blanked;
   - an allele is named only when ≥ 20 haplotypes carry it;
   - complementary suppression: when a masked cell could be recovered as the total minus the shown
     cells, the smallest shown cell is also masked (`<20*`) or folded into `other`.
7. Local draft mode: rebuilds the same tables from the committed aggregates.
   - (a) and (d) come from the 27 content strings with ≥ 20 complete haplotypes. The lower
     bound sums their %; the upper bound adds the share those strings do not cover; the plotted
     estimate renormalises over the covered share.
   - (c) is copied from 51.
   - (b) sums 43's pooled, known-call, all-haplotype allele counts to the protein field.

**Parameters.**
- Kinship 0.0442: same unrelated definition as 24/29/43/51.
- n < 20 masking: project disclosure rule.
- Top 5 alleles per gene × ancestry in the tables; 3 per gene, 6 genes (2DL1, 2DL4, 3DL1, 3DL2,
  3DL3, 2DS4) in the figure: the genes with the most known calls.
- Protein-field default: full-resolution KIR names are dominated by non-coding novelty (59% of
  calls are novel, 43), so a 7-digit table would mostly describe the known minority.

**Run.**
- Local draft 2026-10-11:
  `python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode local`, outputs in
  `reports/hla_popgen/61_kir_by_ancestry_panel/`.
- VM: pending (`vm_parts/61.sh`).
- Tests: `tests/test_61_kir_by_ancestry_panel.py`, 12/12 on synthetic GTFs, including reading 59's exclusion file. They check the
  completeness restriction, removal of relatives and of excluded people, masking, the absence of
  person ids, the local draft and the figure.

**Results (local draft; to be replaced by the VM run).**
- Complete haplotypes, unrelated, with ancestry: 17,868 (AFR 4,374; AMR 3,984; EAS 2,310;
  EUR 4,500; MID 752; SAS 1,948).
- Content strings with ≥ 20 haplotypes cover 91.2–96.0% of complete haplotypes per ancestry
  (MID lowest).
- cB share, complete haplotypes, 95% CI:

  | ancestry | cB share |
  |---|---|
  | EAS | 32.0 (30.2–34.0) |
  | AFR | 38.7 (37.2–40.1) |
  | AMR | 42.7 (41.2–44.3) |
  | EUR | 44.5 (43.1–46.0) |
  | MID | 52.0 (48.4–55.5) |
  | SAS | 57.1 (54.9–59.3) |

  Each value is 0.4–1.6 points lower than on all haplotypes, and the ancestry order does not
  change.
- Gene content, estimates (lower–upper bound):
  - KIR3DS1 AFR 6.4 (6.1–10.5) against EUR 19.6 and SAS 27.0;
  - KIR2DL5B EAS 4.0 (3.8–9.3) against AFR 23.5;
  - KIR2DS2 EAS 13.8 (13.0–18.5) against SAS 39.1;
  - KIR2DL3 EAS 86.2;
  - KIR3DL1 AFR 89.7, SAS 67.7.
- Mean distinct genes per complete haplotype: 9.35 (EAS) to 9.97 (SAS). The 9-gene class is
  70.7% of EAS and 47.2% of SAS haplotypes.
- Top alleles, pooled, known calls, protein field:

  | gene | top allele | share |
  |---|---|---|
  | KIR2DL1 | *003 | 51.1% |
  | KIR2DS4 | *001 | 56.5% |
  | KIR2DL4 | *001 | 34.7% |
  | KIR3DL2 | *002 | 30.8% |
  | KIR3DL1 | *015 | 25.4% |
  | KIR3DL3 | *003 | 21.2% |

**Interpretation.**
- The ancestry pattern of 43 holds on complete haplotypes:
  - EAS is the most cA-rich group, SAS and MID the most cB-rich;
  - AFR is low for the telomeric 3DS1–2DL5A–2DS5–2DS1 block.
- Framework and common A-content genes are 2–9 points more frequent on complete haplotypes than
  on all haplotypes. This is the size of the fragmentation bias in 43's presence values.

**What it does not show.**
- The draft gives no per-ancestry allele frequencies; 43 published only pooled ones.
- In the draft, (a) and (d) are estimates. They assume the rare content strings look like the
  common ones, which pulls B-content genes down by up to ~3 points. Genes that always travel
  together on common strings (3DL1/2DS4, 3DS1/2DL5A, 2DP1/2DL1, 2DL4/3DP1) get identical
  estimates.
- Copy number beyond presence (duplicated 3DL1/S1 haplotypes, 51) is not shown.

**Caveats.**
- If completeness depends on ancestry or content, complete-haplotype frequencies are biased. The
  VM table `completeness_by_ancestry.tsv` measures this.
- The 991 `sequel2` people (mostly AFR) were never KIR-called.
- The coverage filter (59) is not applied yet.
- Ancestry is predicted, not self-reported.
- The protein field of a genomic-only novel call is inherited from its nearest named allele. It
  is correct only if Immuannot's nearest template shares that allele's protein, which is not
  checked here.

## 62 — Phased class II haplotype units for fine-mapping

**Question.** For each person and each assembled haplotype, which DQA1~DQB1 and DRB1~DQA1~DQB1
haplotype does it carry? How sure are we that the alleles really sit together on one
chromosome?

**Why it matters.**
- Fine-mapping and association work on DQ (the G1/G2 heterodimer question) needs phased units,
  not unphased genotypes.
- Our long-read assemblies give physical phasing for free, but only when the genes share one
  contig.
- The trio benchmark (57b) showed that same-contig phasing is reliable for class II and wrong
  about 40% of the time for A~B. Each unit therefore has to carry its own evidence level, plus
  the QC flags a downstream user needs to filter on.

**Inputs.**
- HLA Table 1 `~/pipeline_outputs/hla_calls_rich.tsv`: `person_id`, `hap`, `contig`, `gene`,
  `copy_index`, `consensus`, `gene_start`.
- `~/pipeline_outputs/cohort_membership.tsv` (platform, predicted ancestry).
- The AoU v9 relatedness table.
- 58's INTERNAL `~/s05/internal/58/per_person_depth.tsv` (optional).

**Method.**
1. Every (person, haplotype) in Table 1 gets two rows, one per unit (DQA1~DQB1,
   DRB1~DQA1~DQB1).
2. Allele strings: all copies of each gene, full consensus and 2-field (`_s05_common.two_field`;
   empty when field 1 or 2 is novel). The unit labels use copy 1.
3. Physical evidence: copy-1 contigs of the unit's genes, plus a pair rule from 57b trio truth.
   - The rule says, per gene pair, whether same-contig phasing counts:
     - allowed: DQA1~DQB1, DRB1~DQA1, DRB1~DQB1, DPA1~DPB1, B~C (trio mismatch < 7.7%);
     - flagged: DQB1~DPB1, 21%;
     - not allowed: A~B and A~C, ~40%;
     - any pair not listed: not allowed.
   - Unit level:
     - `physical`: one contig, all pairs allowed;
     - `physical_flagged`: one contig, a flagged pair;
     - `not_allowed_by_rule`: one contig, a pair that is not allowed;
     - `split_contigs`: two or more contigs;
     - `incomplete`: a gene has no call.
   - For the two units built here, every pair is allowed, so in practice the level is decided by
     contig sharing. The other levels exist so the same code serves DP or class I units.
4. DQ G1/G2 class of the copy-1 pair: Petersdorf 2022 on field 1, the rule of 37/50e.
5. DQ structural flag from the 50e follow-up:
   - `anomalous_spacing`: same contig and |gene_start(DQB1) − gene_start(DQA1)| ≥ 100 kb;
   - `anomalous_multicopy`: a second DQA1 or DQB1 copy;
   - `canonical`: otherwise, on one contig;
   - `not_assessable`: genes missing or split.
6. Coverage flag, two sources:
   - preferred: 59's INTERNAL `exclude_people.tsv` (`--exclude-people`). Its `window_MHC` rows give
     each person's MHC depth on 59's scale and 59's proposed MHC cutoff. `low_coverage` below the
     cutoff, `pass` otherwise (everyone passes if 59 proposed no MHC cutoff), `not_measured` if the
     person is not in the file; status `set_by_59`. 59's `exclude` column is not used, because it
     also reflects the KIR window, which says nothing about class II units;
   - fallback (no 59 file): 58's relative depth over `window_MHC` at MAPQ ≥ 0, divided by the
     cohort median on the person's platform (59's normalisation; platforms with < 20 people use the
     all-platform median). `low_coverage` below 0.6 on that scale; status `provisional_until_59`.
   - Each row records the scale, the value compared, the cutoff and the status.
7. Unrelated flag:
   - greedy removal at kinship ≥ 0.0442 among HLA-called people;
   - `not_assessed` if the relatedness table is missing, so it is never silently True.
8. Outputs:
   - INTERNAL `haplotype_units.tsv.gz`, with an `# INTERNAL` first line, plus the dictionary and
     an md5/sha256 manifest. Only under `~/s05/internal/62`; the script refuses repo and
     `reports/` paths.
   - PUBLIC:
     - 2-field unit frequency per ancestry (unrelated, `physical`, coverage not low; units named
       only with ≥ 20 haplotypes);
     - a QC summary per evidence level, DQ class, DQ structure, coverage flag, platform and
       unrelated flag;
     - complementary suppression in both tables.

**Parameters and why.**
- Pair statuses: proposed S05 rule from 57b (179 trios), configurable with `--pair-rules`.
- 100 kb spacing: same-contig DQA1–DQB1 spacing is 16.6–24.2 kb on 98% of haplotypes (50e
  follow-up), so 100 kb sits far outside the normal range. The spacing definition (distance
  between `gene_start` values) is an assumption; the raw `dq_spacing_kb` is written so the cut
  can be moved.
- Coverage: 59's proposal when available. The fallback cutoff 0.6 on the normalised scale is a
  placeholder. Normalising matters because the MHC window runs ~17% under genome depth in typical
  people (58 smoke test median 0.83), so an absolute cutoff would mean different things per region
  and platform.
- 2-field as the default unit resolution: 57b trio Mendelian error is about 4% at 2 fields and
  12.6% at full resolution.

**Run.**
- VM: pending (`vm_parts/62.sh`). Best after 58 and then 59, so the coverage flag is set by 59.
- Upload to the in-perimeter collaborator bucket is not done by the script. The `gsutil rsync`
  commands are in the script docstring, the dictionary and `vm_parts/62.sh`, for Marc to run after
  he has reviewed the MANIFEST.
- Data dictionary (no data): `reports/hla_popgen/62_cole_haplotype_units/README.md`, generated by
  `--dictionary-only`. A test checks that it matches the script.
- Tests: `tests/test_62_cole_haplotype_units.py`, 13/13 on synthetic Table 1. They cover:
  - every evidence level and configurable rules;
  - DQ classes and both anomaly flags;
  - the coverage flag from 59's file (MHC rows only, KIR-only exclusions ignored, no-cutoff case),
    the normalised 58 fallback, and no depth data at all;
  - the unrelated flag and `not_assessed`;
  - masking and complementary suppression;
  - no person ids in PUBLIC outputs;
  - refusal of repo and `reports/` paths for INTERNAL output.

**Results.** Pending the VM run.

**What it does not show.**
- `physical` means the genes share one contig. It does not mean the contig is free of phasing
  error; 57b bounds that at < 7.7% for class II pairs.
- `split_contigs` haplotypes keep their hap label, but 50e showed the label is not reliable
  across a contig break.

**Caveats.**
- The 991 `sequel2` people (mostly AFR) have no HLA calls.
- The pair rule and the coverage cutoff are proposals, not ratified decisions.
- Platform changes novelty rates (50c) and should be a covariate in any use of the units.

## 63 — HLA and KIR calls for the `sequel2` people (Tier-3 trim fix)

**Question.** Can the 991 people on the `sequel2` platform be HLA- and KIR-typed from their phased
assemblies, with the same caller and settings as the rest of the cohort, even though they have no
assembly-to-GRCh38 alignment files? If so, what do their calls add to the cohort tables?

**Why it matters.** `sequel2` is AoU's original long-read pilot cohort, and it is mostly AFR. None
of these people have HLA calls (50's waterfall: 0 of 991 HLA-called), and none have KIR calls (43,
61). So AFR is under-sampled in every HLA and KIR result to date, and any ancestry-purity rule
(DECISIONS.md, S05 open item) will cut AFR further. Adding them would make the AFR estimates more
precise and would let the platform effect be checked on a third platform.

**Diagnosis (full account in `SEQUEL2_DIAGNOSIS.md`).**
- Every `sequel2` person falls in trim tier `self_align_needed`. Their manifest rows have the hap
  FASTAs, but not `assembly_hap{1,2}_aln2_hg38_{paf,bam}`, which Tiers 1 and 2 need.
- The HLA production run and the KIR run (43) both excluded this tier
  (`--skip-trim-tier self_align_needed`). The Tier-3 self-align path had failed its only real test
  (EXPERIMENTS.md 2026-08-05, n = 2): 0 genes, for two coupled reasons.
  - `samtools faidx` cannot index the source FASTA, which appears to be plain gzip rather than BGZF.
  - Self-aligning a whole assembly to chr6 alone gave contigs from other chromosomes no true locus.
    With no filtering and a min/max span, the "trimmed" input was 1.3–1.9 GB per haplotype instead
    of about 4 MB.
- For KIR there is a third reason. The self-align reference was hardcoded to chr6, so a chr19
  region could never match.

**Inputs.**
- Cohort: `~/pipeline_outputs/immuannot_cohort_full.tsv`, read only. People are selected by
  `trim_tier == self_align_needed`, and the script logs their platform mix.
- Per person: `assembly_hap1_fa` and `assembly_hap2_fa` from `v9/wgs/long_read/manifest.tsv`, read
  through the controlled-data automount (`~/workspace/vwb-aou-datasets-controlled-v9`).
- Reference: `~/ref/Homo_sapiens_assembly38.fasta` with its `.fai`. Three things are derived from it
  into `~/s05/internal/63/ref/`:
  - two region slices;
  - a primary-only copy (chr1–22, X, Y, M);
  - that copy's minimap2 `asm5` index.
- Caller: Immuannot (`~/tools/Immuannot`), with refdata IPD-IMGT/HLA 3.55.0 (ENVIRONMENT.md quirk
  #38). This is the same caller and refdata as the production run.

**Method.** For each person, each locus and each haplotype, `run_immuannot_person.py` runs with
`--enable-self-align-fallback`:
1. **Stage 1 (candidate contigs).** minimap2 `-x asm5` aligns the whole haplotype assembly to a
   slice of GRCh38: the target window ±1 Mb. That is chr6:28.5–34.5 Mb for HLA and
   chr19:53.6–55.92 Mb for KIR. A contig is a candidate if it has at least one *primary* alignment
   block of ≥10 kb. If the candidates add up to more than 600 MB, the haplotype is refused, because
   that pattern is the whole-assembly pull-in seen in 2026-08.
2. **Stage 2 (confirmation against the whole genome).** Only the candidate contigs are extracted.
   They are read straight from the gzip file, in one pass, with no index, so plain gzip works. They
   are then aligned with minimap2 `-x asm5` to the whole *primary* GRCh38. In that alignment every
   contig competes against all chromosomes, so MAPQ means something again. ALT, HLA and decoy
   contigs are left out so that a true MHC or KIR contig is not split between chr6 and its ALT
   copies. Only blocks that meet all of these are used:
   - primary;
   - ≥10 kb long;
   - MAPQ ≥20;
   - PAF matches/block length ≥0.5;
   - overlapping the target window (HLA chr6:29.5–33.5 Mb; KIR chr19:54.60–54.92 Mb, the exact
     windows of the production and 43 runs).
3. **Trim.** On each contig, the kept blocks are grouped when they are ≤250 kb apart along the
   contig. The group with the most aligned bases is padded by 100 kb (the production `--pad`) and
   extracted.
   - If the trimmed input would be larger than 25 MB, the haplotype is refused rather than typed.
   - Expected sizes are about 4 MB for HLA and about 1 MB for KIR.
4. **Typing.** `immuannot.sh` runs on the trimmed FASTA, with the same settings as production. A
   person counts as done when both haplotype GTFs contain at least one call (the production
   orchestrator's rule).
5. **Merge.** Per-person fragments are combined into
   `~/s05/internal/63/<locus>/immuannot_calls.sequel2*.tsv`. The existing canonical tables are read
   but never written. The existing table plus the new rows is written to
   `~/s05/internal/63/merged/immuannot_calls{,.kir}.with_sequel2.tsv`.

**Parameters and why.**
- **Flank ±1 Mb (stage 1).** It covers the window plus the 100 kb pad with margin. The flank only
  affects which contigs are proposed; stage 2 decides.
- **Block ≥10 kb.** Interspersed repeats (L1 is ~6 kb at full length) do not reach it, while a true
  MHC or KIR contig has blocks of tens to hundreds of kb.
- **MAPQ ≥20.** This is the conventional "uniquely placed" floor, also used by 58. It is only
  meaningful because stage 2 aligns to the whole genome.
- **Matches/block ≥0.5.** Without base-level alignment (`-c`), PAF column 10 is a seed-coverage
  count, not an identity. 0.5 removes diverged paralog hits and keeps near-identical true hits.
  This threshold is loose on purpose; MAPQ and length do most of the filtering.
- **Cluster gap 250 kb.** It is larger than any KIR gene-content insertion relative to GRCh38's
  haplotype, so a B-haplotype contig stays one piece. A stray block 1 Mb or more away along the
  contig is dropped.
- **Size guard 25 MB.** About six times the normal HLA input and 25 times the normal KIR input. It
  is far below the 1.3 GB failure mode.

**Run.** Code and tests:
- `scripts/run_immuannot_person.py` ("Tier 3 v2"; Tiers 1 and 2 unchanged) and
  `scripts/hla_popgen/63_sequel2_recall.py`.
- `tests/test_63_sequel2_recall.py`: 16/16 pass on synthetic data. The synthetic data includes:
  - a plain-gzip haplotype with no PAF or BAM;
  - a chr1 contig carrying a 12 kb MHC segmental duplication (a stage-1 candidate, rejected in
    stage 2);
  - a repeat contig;
  - a stray block that would have stretched the MHC contig to 1.9 Mb.

  The tests also check that the default `regions_from_paf()` output is unchanged, and that a chr19
  region with a chr6 reference is refused.

VM commands are in `vm_parts/63.sh`:
1. `--prepare-refs`
2. `--smoke` (2 people, both loci)
3. the full run

Pending.

**Runtime estimate (4 vCPU, 31 GB; to be replaced by the smoke run's `trim_summary.tsv`).** The
estimate is built from these past measurements:
- The 2026-08-05 chr6 self-align took 157–215 s per haplotype at 4 threads. This is an upper bound
  for stage 1, because the slice reference is about 28 times smaller.
- The production HLA run handled ~249 people/h at 24×4 threads, which is ≈5.8 min of typing per
  person at 4 threads.
- The 43 KIR run handled ~390 people/h at 20×4, which is ≈3.1 min per person.

| Step | Estimate per person (both loci) |
|---|---|
| Stage 1 (2 loci × 2 haplotypes) | ≈12 min |
| Candidate extraction (gzip stream) | ≈4–8 min |
| Stage 2 (index load ≈0.5–1 min, plus alignment of tens to hundreds of Mb) | ≈8–16 min |
| Typing | ≈9 min |
| **Total** | **≈30–50 min, central ≈40** |

For 991 people at one worker × 4 threads, that is **≈500–800 h, central ≈650 h (~27 days, ~$240
at $0.37/h)**.

The work is CPU-bound, so more workers on 4 vCPUs do not help. Each concurrent stage 2 holds a
~5–8 GB index, so on 31 GB RAM `--workers` should stay ≤2. A bigger app would be faster: about
10 workers × 4 threads on 80 vCPU, limited by RAM, would take about 2.5–3 days. Resizing needs
Marc's authorisation.

**Results / interpretation.** Pending the VM run. To report when it is done:
- **Coverage:** people done per locus, haplotypes per trim method, and the median trimmed size
  (`run_summary.tsv`, `trim_summary.tsv`, `failure_reasons.tsv`).
- **Check against other platforms:**
  - genes per haplotype compared with `sequel2e`/`revio` of the same ancestry;
  - novelty rate by field compared with 24.
- **Merged tables:** the merged tables feed 01 (rich tables) and then 50, 61 and 62. That needs a
  re-run of 01 with `--outroot ~/s05/internal/63/hla/pipeline_outputs`, which is a separate step.

**What it does not show.**
- It does not show that the Tier-3 trims select exactly the contigs Tier 1 would have chosen. No
  `sequel2` person has a PAF to compare with. A concordance check is possible: run Tier 3 on ~20
  `revio` people (forced, into a separate outroot) and compare their calls with the production
  calls. This is suggested before the results are pooled with the other platforms.

**Caveats.**
- The plain-gzip explanation is still a strong hypothesis. `htsfile` has not confirmed it yet
  (`SEQUEL2_DIAGNOSIS.md` §5 U3).
- The `sequel2` assemblies come from an older pipeline (the v7 pilot), so their contiguity and
  base quality may differ. Platform must stay a covariate.
- The approximate matches/block ratio is a proxy, not an identity.
- The `<20` masking applies to every public count (the small-cell open item in DECISIONS.md).

## 64 — Novel HLA alleles: current IPD-IMGT/HLA release and read-level support

**Question.** Are the novel HLA calls (protein-novel first, synonymous next) real alleles, or
assembly and basecalling artifacts? Two checks per distinct novel sequence: (a) is it now in a
newer IPD-IMGT/HLA release than the 3.55.0 Immuannot typed against? (b) do the person's own HiFi
reads (and, optionally, short reads) carry the called sequence, or the nearest known allele?

**Why it matters.** 50d showed the classical protein-novel funnel collapsing: for the 8 classical
genes, 4,421 distinct protein-novel sequences → 174 with ≥ 2 unrelated carriers → 8 not
homopolymer-indel-only → 0 clean at ≥ 20 carriers. Sequel IIe calls more novel alleles than Revio
in the same ancestry. The open item in `context/DECISIONS.md` proposes demoting classical
protein-novel counts to a resource or supplement until a current release and read-level support
have been checked. Recurrence alone cannot settle it: a systematic HiFi error (homopolymer length)
recurs across people, and a real rare allele is a singleton. Reads settle it per call: a real
allele is carried by the reads of that haplotype, an assembly error is not.

**Inputs.**
- HLA Table 1 `~/pipeline_outputs/hla_calls_rich.tsv`: `person_id`, `hap`, `contig`, `gene`,
  `copy_index`, `consensus`, `novelty_depth`, `gene_start`, `gene_end`, `cds_mut`,
  `template_warning`, `template_allele`, `alleles`.
- Per haplotype, Immuannot's own files under `~/pipeline_outputs/people/<pid>/immuannot_output/`:
  `<hap>.trimmed.fa` (the contigs Immuannot read), `<hap>.gtf.gz` (CDS rows),
  `<hap>/cds.fa.gz` (Immuannot's observed CDS, used as a cross-check).
- `~/pipeline_outputs/cohort_membership.tsv` (platform, predicted ancestry); the AoU v9
  relatedness table.
- IPD-IMGT/HLA, current release, from the ANHIG/IMGTHLA GitHub repository at git ref `Latest`
  (`--ipd-release`): `hla_nuc.fasta`, `hla_gen.fasta`, `Allelelist.txt` (version line),
  `Allelelist_history.txt` (allele name per release). Version recorded in `ipd_release.tsv`.
- Long reads: `grch38_bam` from `v9/wgs/long_read/manifest.tsv` through the controlled-data
  automount. Gene coordinates on GRCh38: 58's `regions.bed` (GENCODE).
- Optional short reads: `cram_uri` from `v9/wgs/cram/manifest.csv`, decoded with
  `~/ref/Homo_sapiens_assembly38.fasta`.

**Method.**
1. Calls: Table 1 rows with novelty depth 2 (protein) or 3 (synonymous) (`--depths`; depth 4,
   beyond-CDS, is optional), copy 1, not `undetermined`, in the greedy unrelated set (kin ≥ 0.0442,
   the `_s05_common` rule).
2. Sequence key: `(gene, cds_mut)`, the 50d key, so the funnel is directly comparable. Depth 4 has
   no `cds_mut` (its CDS equals a known allele); its key is the hash of the genomic gene span and
   exists only for tested calls.
3. Which calls are tested: per key at most 5 carrier haplotypes (`--max-per-key`), picked in turn
   from each platform × ancestry stratum in a fixed hash order. The same calls are picked on every
   rerun. Every key therefore gets at least one test, and a recurrent key is tested across
   platforms and ancestries rather than five times in the largest stratum.
4. Called sequence. The CDS is rebuilt from the GTF `CDS` rows of that contig, gene and copy,
   sliced from `<hap>.trimmed.fa` and reverse-complemented on the minus strand. It is compared with
   the record in `cds.fa.gz`; calls where the two disagree are skipped and counted
   (`cds_check`), because the tested sequence must be the called one. A terminal stop codon is
   removed on both sides (the GTF CDS rows exclude it; IPD's CDS includes it).
5. Re-match (a). Exact CDS match against the current release. "Added since 3.55.0" = the allele's
   accession has no name in the `3550` column of `Allelelist_history.txt`. Fallback: accessions in
   Immuannot's refdata, which is the 3.55.0 snapshot. Classes per key, first that applies:
   `known_in_baseline` (exact match to an allele already in 3.55.0: the novel label was wrong),
   `protein_known_in_baseline` (depth 2 only: the protein is a 3.55.0 protein), `now_known`
   (exact match to an allele added since 3.55.0), `now_known_protein` (depth 2: the protein of an
   added allele), `now_known_partial` (an added partial-sequence allele lies inside the called CDS
   and not inside the nearest known CDS, so it covers at least one novel difference),
   `still_novel`. Depth 4 uses genomic sequences with 44's containment rule.
6. Two candidate haplotypes. Wc = the called contig over the gene ± 2 kb, in coding orientation.
   The nearest known allele is the reference allele named in `cds_mut` (Immuannot's own nearest).
   Its CDS comes from the current release, with Immuannot refdata as fallback. Called and nearest
   CDS are aligned (banded global alignment, unit costs), and each difference is placed on Wc
   through the exon coordinates. Tn = Wc with those differences reverted: the same haplotype, but
   carrying the nearest known allele. The two differ only at the novel positions. Calls with more
   than 40 differences, or with an indel crossing an exon junction, are skipped and counted.
7. Reads (b). `samtools view -F 0x904` (primary alignments only) over the gene ± 1 kb on the GRCh38
   BAM, at most 80 reads. Each read is aligned to Wc with minimap2 (`mappy`, preset `map-hifi`).
   Differences closer than 50 bp are grouped into windows (25 bp flank; a homopolymer run is
   always included whole). For each window the read fully spans, its edit distance to the Wc
   segment and to the Tn segment is computed. The read votes for whichever is closer; equal
   distance is no vote.
8. Which haplotype a read comes from. In a heterozygous person, reads from the other haplotype
   carry the other allele and would look like support for the nearest allele. When the other
   haplotype has a Table 1 call for the same gene, its window Wo is built the same way and every
   read is also aligned to it. Over the part of the read both alignments cover, edits are counted
   against Wc and against Wo, leaving out the tested windows themselves, so the tested sites never
   decide the haplotype. A difference of ≥ 1 assigns the read to this or the other haplotype;
   otherwise it is unknown. Other-haplotype reads are dropped.
9. Counting and verdict.
   - Three-target mode (≥ 3 voting reads assigned to this haplotype): a real allele is expected at
     a fraction near 1, an artifact near 0.
   - Two-target mode otherwise (this + unknown reads): about 0.5 for a heterozygous real allele,
     about 1 for a homozygous one, about 0 for an artifact.
   - `supported`: ≥ 3 reads for the called allele and fraction ≥ 0.7 (three-target) or ≥ 0.3
     (two-target).
   - `contradicted`: ≥ 3 reads for the nearest allele and fraction ≤ 0.15.
   - `ambiguous`: everything else, with the reason recorded (`low_reads`, `mixed`,
     `no_window_spanned`).
   - The thresholds sit between the expected values, so a single read error cannot move a call
     across them.
10. Homopolymer indels: an indel of one repeated base next to a run, with the longer side ≥ 3
    (25's sequence-context rule). For each spanning read, the run length is read out between
    6-bp exact flanks and counted as the called length, the nearest length or another length.
    This shows directly whether the reads agree with the assembly's run length. The 50d proxy on
    `cds_mut` is kept as the split for public tables, so they stay comparable with 50d; the
    sequence-context flag is reported per call.
11. Short reads (c, `--sr`): the same windows and votes with srWGS reads (preset `sr`). A read votes
    only on windows it spans; ≥ 3 voting reads are needed. Haplotype assignment rarely applies to
    150 bp reads, so the two-target thresholds are used.
12. Key-level read class: `supported` if any tested carrier is supported and none contradicted;
    `contradicted` in the reverse case; `mixed` if both occur; `ambiguous` or `untested` otherwise.
13. Disclosure: PUBLIC by default. Call, carrier and read counts of 1–19 are shown as `<20`, and a
    percentage over such a cell is blank. Distinct-key counts are exact (allele properties,
    `_disclosure.ALLELE_DISTINCT`). `now_known_alleles.tsv` names public IPD alleles with carriers
    masked and `review_flag` where carriers < 20. Per-call evidence (person ids) stays in
    `~/s05/internal/64/` on the VM.

**Run.** Pending (VM). Commands: `vm_parts/64.sh`. IPD download:
`$PY 64_novel_allele_validation.py --fetch-ipd --ipd-release Latest --ipd-dir ~/s05/ref/ipd_latest`;
re-match only, classical genes:
`--no-reads --genes A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 --out-dir ~/s05/results/64_ipd --internal-dir ~/s05/internal/64_ipd`;
smoke with reads:
`--limit 30 --workers 4 --out-dir ~/s05/results/64_smoke --internal-dir ~/s05/internal/64_smoke`;
full:
`setsid nohup $PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest --workers 12 --samtools $ENVBIN/samtools > ~/s05/logs/64.log 2>&1 < /dev/null & disown`
(`$PY` = `~/repos/pilot-validation/.pixi/envs/spechla/bin/python3`, run from
`~/repos/pilot-validation-main/scripts/hla_popgen`). Test suite:
`tests/test_64_novel_allele_validation.py`, 21/21 on synthetic data. Planted cases:
- a true exon-3 SNV allele carried by the reads (supported, three-target mode);
- a one-base homopolymer loss in the assembly with reads carrying the known run (contradicted;
  run-length readout favours the nearest length);
- a minus-strand allele present in a fake newer release and absent from its `3550` column
  (`now_known`, supported);
- a carrier without a BAM, a relative of a carrier (excluded), and a third carrier beyond
  `--max-per-key` (not tested).

**Results / interpretation.** Pending the VM run. What to look at first:
- `ipd_release.tsv`: the version and the number of alleles added since 3.55.0 per gene.
- `funnel_with_evidence.tsv` for the classical genes: how many of the 174 recurrent keys are
  `now_known`, read-supported or contradicted.
- `summary_by_stratum.tsv`: the `hp_only` versus `other` split, and Revio versus Sequel IIe.
- `hp_runlength_by_gene.tsv`.

**What it does not show.**
- Whether a supported allele is functional or expressed.
- Whether a read-supported allele is a somatic or cell-line change rather than germline (there
  is no second tissue).
- Calls that were not tested beyond the per-key cap. Their key verdict is inferred from the tested
  carriers.
- Depth-4 keys are counted only over tested calls, so depth-4 key totals are not the cohort's
  distinct non-coding sequences.

**Caveats.**
- The reads are the same HiFi reads the assembly was built from. A systematic basecalling error at
  a long homopolymer can be present in the reads too, so read support at a homopolymer indel is
  weaker evidence than at a substitution. The run-length readout and the short-read check (`--sr`,
  independent chemistry) address this. For homopolymer-only keys, short-read agreement is the
  stronger test.
- The nearest allele is Immuannot's (3.55.0). If an allele added later is closer, the test still
  compares the call with the 3.55.0 nearest. The re-match in (a) reports the newer allele
  separately.
- Reads come from the GRCh38 gene region, so reads from haplotypes too divergent to map there are
  missed. DRB3/4/5 (absent from the primary assembly) use the span of the DR genes GENCODE places
  on chr6, which is an approximation.
- Haplotype assignment needs at least one heterozygous site outside the tested windows inside the
  read. In long runs of homozygosity, reads stay `unknown` and two-target thresholds apply.
- The IPD download URLs (raw versus Git-LFS media), the `Allelelist_history.txt` layout (column
  `3550`), whether `mappy` is installed in the python 3.8 env, and the agreement of GTF-rebuilt CDS
  with `cds.fa.gz` were not verifiable from the cloud container (no network to GitHub, no real
  data). Each one is checked by a step in `vm_parts/64.sh` before the full run.
- Small-cell publication of `now_known_alleles.tsv` (named alleles with masked carrier counts) is
  flagged for a human decision; see the open disclosure item in `context/DECISIONS.md`.

## Review corrections (2026-10-12, external review; integrated in the S05 report)

- **59:** the cutoff rule (≥ 50% events at ≤ 10% people) cannot pass for an outcome with prevalence p > 20% (ceiling 10/p).
  The result is read from `logistic_or.tsv` and `outcome_by_depth_bin.tsv`: associated, partly reverse causation; no filter.
- **64:** `read_support` hp votes now use only reads assigned to the called haplotype (`hap == "this"`); before, unassigned
  reads (which include the other haplotype) were pooled, which pushes "known length" towards 50% even for a real het allele.
  `hp_runlength_by_gene.tsv` from the 10-12 run is confounded until the rerun. The 345 CDS disagreements: source undetermined.
- **65:** reads-mode het AF bias comes from the read rules (reference reads must span ~7.9 kb; insertion reads need a ~0.3 kb probe
  hit), so expected het AF ≈ (L+0.3)/((L+0.3)+(L−7.9)) ≈ 0.64 at L ≈ 18 kb. Added `--genotyper mixture` (3-component binomial EM per
  feature × platform, free f and w, call at posterior ≥ 0.9, writes `af_mixture_fit.tsv`).
- **56b/57b:** the relative filter biases the error down; the two checks share people; the per-call error is about half the
  pair rate; "< 7.7%" is a masking ceiling.
- **Disclosure:** sums across tables (65 genotype vs AF bins; 58 CN vs depth bins; older 50–56 tables) were re-masked.

## 65b / 65c, 67, 68 (2026-10-13)

- **65b/65c:** assembly calls on a stratified 2,000 (≈330 per ancestry), stopped at 913 on Marc's instruction; order was random,
  so the done set is a random subset. Concordance on people called by both methods; ppv_c / sens_c per copy number; HWE on
  people called by all three callers. Public table: rates rounded to whole percent and blanked where the implied number of
  disagreeing people is 1–19 (two genotypers on the same people otherwise leak small counts by subtraction).
- **67:** depth per call = 58's per-gene relative depth (MAPQ 0) ÷ the median for that gene × platform. Pair records use the
  56b sharing filter; error = no shared 2-field allele. Trio records use 57b parent of origin. Logistic: error ~ depth (clipped
  at 2.0) + confidence score. Switch topology: consecutive resolved genes on one contig, switch = origin changes, zones by hg38
  midpoint. Prune: worst x% of unrelated people by each metric vs 20 random removals; pair/trio errors re-scored with the same
  cut applied to everyone. Rates blanked when numerator or denominator is 1–19; ancestry shares of removed people blanked under 30.
- **68:** cohort quantiles of numeric columns in `rnaseq_metadata.tsv` and the RNA-SeQC 2 metrics file (ID columns excluded),
  and per-sample TPM shares by gene family streamed from the gene TPM GCT. Globin read share = 1 − Non-Globin / Mapped Reads
  (an earlier row divided by Total Reads, which includes unmapped reads, and was dropped).
