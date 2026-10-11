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
