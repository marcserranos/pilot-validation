# AoU / Verily Workbench Genomic Data Access — Research Notes

**Purpose:** practical reference for locating srWGS/lrWGS CRAM/BAM paths in the All of Us (AoU) Researcher Workbench, region-slicing them to the HLA locus, and converting to FASTQ for SpecHLA/SpecImmune. Written from web research only — **nothing in this document has been executed against the real Workbench environment.** Every command is an untested template to adapt and verify live.

Confidence key: **[HIGH]** official AoU docs or upstream tool docs directly state this. **[MED]** inferred by combining multiple credible sources, or from docs that may be stale. **[LOW]** best-guess / needs live verification, includes plausible reasoning but no direct source.

---

## 1. How AoU/Verily Workbench exposes WGS CRAM/BAM paths

**[HIGH] Mechanism: a per-sample manifest file, not a single well-known BigQuery table.** Multiple AoU support articles ("How the All of Us Genomic data are organized," various CDR versions — v5/v6/v7/v8 archived articles exist, meaning this article gets re-published per CDR release) describe a **manifest CSV** with one row per sample containing columns equivalent to `person_id`, `cram_uri`, `cram_index_uri` (naming may vary slightly by CDR version). This manifest is generated for you (e.g., via the Genomics extraction / Dataset Builder "Genomics" domain, or provided as a static file in the workspace) rather than being a table you query directly with arbitrary SQL joins against `T_ENT_person`.

**[MED] Bucket path pattern (short-read WGS, at least as of CDRv6-era docs):**
```
gs://fc-aou-datasets-controlled/pooled/wgs/cram/v6_base/<research_id>.cram
gs://fc-aou-datasets-controlled/pooled/wgs/cram/v6_base/<research_id>.cram.crai
```
This is very likely **stale** — bucket names are versioned by CDR release (the `v6_base` segment already tells you this), and your workspace's CDR is `C2025Q4R6` per your session notes, not v6. Treat this exact string as an *example of the pattern*, not the current path. The correct bucket for your workspace should come from:
- The **manifest file itself** (safest — it will contain the literal, current, correct paths for `C2025Q4R6`).
- The **"Controlled CDR Directory"** support article (`support.researchallofus.org` — I could not fetch its body directly, got HTTP 403 from the automated fetch tool; a search snippet confirms it exists and is described as the authoritative place for "bucket locations for accessing the data in analysis notebooks" per CDR version — **you should open this one yourself in-browser**, it's likely gated/JS-rendered).

**[MED] How you actually get the manifest, practically:** AoU's own workflow (per "Featured Workspaces" — specifically one literally named **"Working with All of Us genomic data (CRAM_processing and IGV)"**, and a related **"CRAM_Processing CT"** tutorial workspace) is: use the Cohort/Dataset Builder to select participants and the "Genomics" data domain, which exports a manifest CSV to your workspace bucket, then read that CSV in your notebook to get `person_id → cram_uri` mappings. This is almost certainly reachable in your Verily Workbench UI too, but the *UI labels* ("Cohort Builder," "Dataset Builder," "Genomics" tab) may have been renamed/reorganized in the Verily rebrand — **flag as staleness risk**: the mechanism (a manifest file mapping ID→path) is probably unchanged, but the click-path to generate it may differ from what these tutorials show.

**[MED] The JupyterLab "Snippets" menu** (mentioned in your session context) is the fastest way to check this live — AoU's snippet menu has historically included a "Genomic - CRAM" or similar boilerplate snippet that generates/reads exactly this manifest. Since the user is checking this in parallel, treat that as the primary source of truth over anything in this document.

**[HIGH, important operational finding] Direct remote (`gs://`) access from *viewer tools* is explicitly blocked by design, and AoU recommends localizing files first.** The AoU support article "Using the Integrative Genomics Viewer (IGV)" states that although IGV nominally supports direct GCS access, **"this feature unfortunately cannot currently be used in the Workbench due to data exfiltration controls,"** and the recommended workaround is to `gsutil cp` the CRAM + `.crai` into the current VM/environment first, then work locally. I could not fetch the full article body (403'd), so I don't have the exact accompanying command text — this is from a search-result snippet summarizing it, so treat the precise wording as **[MED]** even though the underlying fact (exfiltration controls exist and affect genomic file access patterns) is corroborated separately by AoU's "Researcher Workbench" security documentation, which describes a **perimeter/VPC-Service-Controls-style policy** restricting "copy, transfer and retrieval of data" across the workspace's security perimeter boundary **[HIGH — this general perimeter-control fact is stated directly in AoU workbench security docs]**.

  - **Important nuance for your samtools-remote-slicing plan [LOW, my own inference, needs live test]:** IGV being blocked doesn't necessarily mean `samtools view gs://...` from *inside* the Workbench VM is blocked the same way. The perimeter controls described are about data leaving the security perimeter (e.g., to a researcher's laptop, or to an external egress endpoint) — the Workbench VM itself, running under the pet service account, sits *inside* that perimeter, and reading from a controlled-tier bucket that's also inside the same perimeter/project family is a fundamentally different network path than IGV-desktop trying to reach out from a browser or unmanaged client. So there is a real chance `samtools view` region-slicing from the VM directly works fine even though IGV's approach doesn't. But do not assume this — the IGV article's exfiltration-controls warning is evidence this general class of "direct GCS streaming access" has tripped AoU's controls before, so this must be your **first live smoke test**, with a `gsutil cp`-then-slice-locally fallback ready to go immediately if the remote streaming approach errors out or hangs.

---

## 2. Long-read WGS sub-cohort specifics

**[HIGH] Source: "Population-scale Long-read Sequencing in the All of Us Research Program" (medRxiv 2025.10.02.25336942 / PMC12622093), a 2025 AoU-affiliated paper.**

- **Cohort composition:** 1,027 self-identified Black/African-American participants sequenced with PacBio HiFi at ~8x coverage, plus a 50-participant subset re-sequenced with **ONT R10.4** at 35x for validation. (Note: an older/different snippet elsewhere mentioned "1,773 New — Total 2,800 samples" for lrWGS in a different CDR-version context — **these numbers are inconsistent across sources/time**, consistent with a growing pilot; don't hardcode either number, the actual current count in `T_ENT_person.has_lr_whole_genome_variant` is your ground truth.)
- **Platform/chemistry:** PacBio Sequel IIe, Circular Consensus Sequencing (CCS) — i.e., **HiFi reads**, not raw/uncorrected long reads.
- **Aligner:** **pbmm2 v1.12.0**, preset `CCS`. (pbmm2 is PacBio's minimap2 wrapper.)
- **Reference genome(s):** aligned to **both GRCh38_no_alt and CHM13v2.0 (T2T)** — i.e., there may be *two* aligned callsets per long-read sample. **This matters a lot for your pipeline**: make sure whatever BAM path the manifest gives you for lrWGS is explicitly the GRCh38-aligned one, not T2T — SpecHLA/SpecImmune's HLA databases assume GRCh38/hg38 coordinates.
- **File format: delivered as BAM, not CRAM**, per this paper — this is a real, structural difference from short-read data (which is CRAM). This directly answers your question in section 3/5: **the "CRAM needs a reference FASTA" gotcha applies to your short-read data but likely does NOT apply to your long-read BAMs** (BAM is not reference-compressed).
- **Access mechanism:** the paper states the combined dataset ("CDRv7 AoU Long-Read dataset") "is already available to approved researchers via the AoU Researcher Workbench," but does **not** detail whether the access pattern (manifest structure, bucket naming) differs from short-read WGS. **[LOW/uncertain]** — my best guess is it follows the same manifest-CSV pattern but in a separate lrWGS-specific manifest/column (which lines up with your `T_ENT_person.has_lr_whole_genome_variant` boolean existing as a distinct flag from `has_whole_genome_variant` — i.e., AoU's own schema already treats these as separate data products, so expect a separate manifest export step for lrWGS, not a unified one).
- There's an older/earlier AoU-related paper, **"Utility of long-read sequencing for All of Us" (Nature Communications 2024 / PMC10822842)**, describing an even earlier technical pilot (HapMap + 2 AoU control samples) — this predates the production 1,027-person cohort and is a *methods validation* paper, not the actual delivered-data source. Don't confuse the two.

---

## 3. `samtools`/`htslib` remote region-slicing of `gs://` CRAM/BAM

**[HIGH] Basic GCS support exists in htslib/samtools since v1.4** via libcurl-based network file access, and your VM's samtools 1.21 build already advertises `gs+http`, `gs+https`, `gs` handlers per your own session notes — consistent with this being enabled.

**[MED] Auth: default GCS credentials should generally suffice on a Workbench VM** (per your notes, the VM runs under a pet service account with presumed default read access to CDR genomic storage) — but two manual-auth patterns are documented elsewhere (e.g., Terra, a sibling Broad/Verily platform with a very similar architecture) for cases where default credential discovery doesn't kick in automatically for htslib's GCS plugin:

```bash
# Option A: token file + HTS_AUTH_LOCATION
gcloud auth print-access-token > /home/jupyter/token.txt
export HTS_AUTH_LOCATION="/home/jupyter/token.txt"

# Option B: inline application-default token env var
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
```
Both are **[MED confidence, sourced from a Terra support community post]**, not AoU-specific — Terra and AoU Researcher Workbench/Verily Workbench share Broad/Verily lineage and similar GCP-backed notebook environments, so this is a reasonable first thing to try if plain `samtools view gs://...` fails with an auth error, but it is not confirmed AoU-specific guidance.

**[HIGH] Region-slicing syntax (BAM — straightforward):**
```bash
# Untested template — adapt path and region
samtools view -b "gs://<bucket>/<path>/<sample>.bam" chr6:29500000-33500000 \
  > sample_hla_region.bam
```

**[HIGH — this is the load-bearing gotcha] CRAM is reference-compressed and needs the reference FASTA reachable, and htslib will NOT auto-fetch a `gs://` reference from the CRAM header's `UR:` tag.** Per htslib/samtools docs and multiple GitHub issues: the CRAM `@SQ` header's `UR:` field is checked for a *local* reference (including `file://`), but **non-local URIs — explicitly including `http://`/`ftp://`, and by extension `gs://`, are not used automatically**. samtools' reference-resolution search order is: (1) `-T`/`--reference` command-line flag → (2) `REF_CACHE` env var → (3) `REF_PATH` env var → (4) fall back to fetching by MD5 from a public EBI reference server (works, but adds network dependency + latency, and is a bad idea to rely on for a pipeline). **Practical implication: always pass `-T <local_or_gs_reference.fasta>` explicitly for CRAM.**

The AoU srWGS reference FASTA is publicly documented as:
```
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
```
**[MED]** — sourced from a search snippet citing AoU support docs; this is a **public** Broad Institute reference bucket (not controlled-tier), so it should be reachable regardless of your workspace's own perimeter controls, and is the standard GATK-best-practices hg38 FASTA. You'll also want its `.fai` index alongside it (should already exist at that path; confirm live). Note: htslib can reportedly read a reference FASTA directly from `gs://` when passed via `-T` even though it won't auto-discover one from a CRAM header — this is a meaningfully different code path (explicit vs. implicit lookup), consistent with the "UR tag not auto-fetched, but -T works" distinction above. Verify live; if `-T gs://...fasta` doesn't work directly, download that reference once to local/persistent disk (it's a static public file, cheap to keep around) and point `-T` at the local copy instead.

**[HIGH] Combined CRAM region-slice command template:**
```bash
# Untested template
REF="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
samtools view -b -T "$REF" "gs://<controlled-tier-bucket>/<path>/<sample>.cram" \
  chr6:29500000-33500000 \
  > sample_hla_region.bam
```

**[MED] Known gotchas checklist:**
- CRAM `.crai` index must be reachable alongside the CRAM at the same `gs://` path (standard `<file>.cram.crai` sibling convention) — if the manifest only gives you the CRAM URI, confirm the `.crai` companion exists at the same prefix, or use the separate `cram_index_uri` manifest column if present.
- Contig naming: confirm whether the CRAM/BAM header uses `chr6` or bare `6` — AoU's srWGS is aligned to `Homo_sapiens_assembly38.fasta`, which uses **`chr`-prefixed** contig names (`chr1..chr22,chrX,chrY,chrM`) in the standard GATK/Broad hg38 distribution, so `chr6:...` should be correct, but confirm with `samtools view -H` on one file before assuming.
- If remote streaming is slow/unreliable or blocked by perimeter controls (see section 1), fallback is `gsutil cp` the CRAM+CRAI (or BAM) locally into the VM's persistent disk, then slice locally — this is literally AoU's own documented recommendation for IGV use, so it's a safe, precedented fallback, just slower and uses more local disk per sample.
- Long-read BAMs (section 2) don't need `-T` at all since BAM isn't reference-compressed — one less moving part for that half of the pilot.

---

## 4. HLA region coordinates and padding recommendation

**[HIGH] Coordinate reference points (GRCh38/hg38), from cross-referenced sources (GWASLab HLA region docs + general MHC literature):**

| Region definition | hg38 coordinates | Notes |
|---|---|---|
| Classical HLA/MHC region | chr6:29,602,238–33,409,896 | Gene boundaries GABBR1→KIFC1 |
| Extended MHC (xMHC) | chr6:25,726,063–33,400,644 | Gene boundaries HIST1H2AA→RPL12P1 |
| Common rounded convention ("xMHC mode") | chr6:25,000,000–34,000,000 | Used by tools like GWASLab as a safe LD-exclusion window |
| Common rounded convention ("HLA/MHC mode", tighter) | chr6:29,500,000–33,500,000 | Also GWASLab; closely matches classical boundaries + small pad |

**[MED] Recommended starting region for this pilot: `chr6:29,500,000-33,500,000` (hg38).** Rationale:
- This is essentially the classical MHC region (29.60–33.41 Mb) with a small explicit buffer already baked in (~100kb–370kb), which is the same rounded window multiple independent tools/papers converge on independently (GWASLab uses this exact range as a named preset), suggesting it's a reasonably battle-tested default rather than an arbitrary guess.
- It's tight enough to keep the region-slice small/fast for a smoke test (~4 Mb vs. the ~8.4 Mb xMHC window), while still safely covering all classical class I/II/III genes plus flanking sequence, which matters because: (a) alignment near region boundaries is generally less reliable (reads spanning the boundary may be soft-clipped/misplaced), and (b) HLA is unusually polymorphic/structurally variable, so some individuals' actual gene boundaries can shift slightly relative to the reference — a pure exact-boundary slice risks clipping true positive signal at the edges.
- If the smoke test shows read dropout or edge artifacts, the fallback is the wider **xMHC convention, chr6:25,000,000–34,000,000** — I could not find SpecHLA/SpecImmune's own docs specifying a recommended slice window (their READMEs describe extracting "HLA-related reads" via `ExtractHLAread.sh` / `ExtractReads.sh` but don't publish exact coordinates for that extraction — see below), so this recommendation is triangulated from general HLA/MHC literature, not from the tools' own documentation. Treat the exact number as **[MED]**, a good starting point, not gospel — the sprint doc's existing framing of padding as an open metaparameter to characterize empirically is the right posture.

**[MED] What SpecHLA/SpecImmune's own extraction scripts actually do (partial info only):**
- **SpecHLA**: `ExtractHLAread.sh` / `spechla-extract-hla-reads` takes a sorted+indexed BAM/CRAM, a reference build flag (`-r hg19|hg38`), sample ID, and output dir. Its README notes "We use the script of Kourami with minor revision for this step" — Kourami is a known HLA-typing tool with its own published HLA-region reference/bed file; the exact coordinates weren't visible in what I could fetch of SpecHLA's docs. **If you want the precise coordinates SpecHLA itself uses, check `script/ExtractHLAread.sh` directly in the repo** (I did not have terminal access to fully inspect the script body — this is a concrete "check this file" action item, not something I've confirmed).
- **SpecImmune**: `ExtractReads.sh -s <sample_id> -i <input_bam_or_cram> -g HLA -o <output_dir> [-r <reference>]`. Its README explicitly requires input mapped to **"whole hg38... contain[ing] alternative contigs and alleles"** — i.e., it wants a BAM/CRAM aligned against an hg38 reference build that *includes* the ALT/HLA-allele contigs (not just the primary assembly), specifically to retain reads that align better to alternate HLA haplotypes than to the primary chr6 sequence. **This is an important compatibility check for your pipeline**: if AoU's srWGS/lrWGS alignment was done against a no-ALT or primary-only reference (plausible — DRAGEN-based pipelines sometimes exclude ALT contigs, and the long-read pilot paper explicitly mentions "GRCh38_no_alt" as one of its two reference builds), SpecImmune's own stated assumption may not hold for reads that would otherwise fall on ALT/HLA contigs, potentially reducing sensitivity for divergent alleles. Given this: **when you do the CRAM/BAM→FASTQ extraction, favor the wider padded region-slice over a scalpel-precise one**, since you can't rely on ALT-contig realignment to rescue divergent-allele reads the way SpecImmune's own docs assume.

---

## 5. CRAM/BAM → FASTQ conversion

**[HIGH] Both SpecHLA and SpecImmune want FASTQ, not BAM/CRAM, as their main typing input** — the BAM/CRAM extraction scripts (`ExtractHLAread.sh`, `ExtractReads.sh`) are pre-typing filtering steps; you still need to go BAM/CRAM → FASTQ after extracting the HLA-region reads, before calling `SpecHLA.sh` or `SpecImmune/scripts/main.py`.

**[HIGH] Core `samtools fastq` mechanics** (from `samtools-fastq` manual + `bam_fastq.c` source):
- `-1 out_R1.fq -2 out_R2.fq` writes READ1-flagged and READ2-flagged reads to separate files (paired-end short-read convention).
- `-0 /dev/null` catches reads where both/neither R1/R2 flags are set (typically unpaired/singleton-like edge cases) — discard these for standard paired short-read typing.
- `-s /dev/null` — the singleton file; if given, only *true pairs* (both R1 and R2 present for a QNAME) go to `-1`/`-2`, and anything missing its mate gets diverted here instead — discard for standard paired workflows, but consider keeping if you're worried about pair-loss from the region slice (a read pair can span the slice boundary — mate outside the sliced region — see gotcha below).
- **Default flag filtering excludes secondary and supplementary alignments already** (`-F 0x900` is the default), so you generally don't need to add extra `-F`/`-G` flags for that specific concern — but it's worth being explicit/defensive about it given how important it is for typing accuracy (secondary/supplementary alignments duplicating or fragmenting reads would corrupt allele-balance signal).
- `-N` forces `/1` `/2` suffixes on read names — largely a legacy-compatibility flag; check whether SpecHLA/SpecImmune's downstream aligner cares (older Novoalign-based short-read paths inside SpecHLA might; modern aligners usually don't need it).

**[HIGH] Short-read paired-end template:**
```bash
# Untested template — assumes sample_hla_region.bam already region-sliced (section 3)
# and coordinate-sorted (or at minimum, name-sorted; samtools fastq works best name-sorted)
samtools sort -n -o sample_hla_region.namesorted.bam sample_hla_region.bam

samtools fastq \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -0 /dev/null \
  -s sample_singletons.fastq.gz \
  -F 0x900 \
  sample_hla_region.namesorted.bam
```
Note the `-n` name-sort step: `samtools fastq` needs mates adjacent to pair them correctly; a coordinate-sorted BAM slice (which is what you get straight out of a region `view`) will scatter mates that map far apart, and functionally *all* mates of an HLA-region read are likely to also fall within/near the region given short fragment sizes — but sort-by-name first regardless, it's cheap insurance and is the documented-correct order of operations.

**[MED] Long-read (PacBio HiFi / ONT) template — no real pairing concept, single-end style:**
```bash
# Untested template
samtools fastq -F 0x900 sample_hla_region_longread.bam > sample.fastq
```
Long reads don't have the R1/R2 pairing model, so this simplifies to one output file. Since long-read data here is delivered as BAM (section 2), there's no CRAM reference-fetch complication for this half of the pipeline.

**[MED] Region-slice boundary gotcha for short-read pairs:** if you slice tightly to the classical MHC region (29.5–33.5 Mb) using `samtools view <region>` on a *coordinate-sorted* BAM/CRAM, samtools' default region-view behavior includes a read if **either** it or its mate overlaps the region window in typical use, but this is more subtle for CRAM/BAM `view` with a single-region query — mate-rescue across the region boundary is not guaranteed the way it would be with a mate-aware tool. This is exactly the scenario the wider recommended padding (section 4) helps mitigate: a generously padded slice window reduces the odds that a read's mate falls just outside the sliced region and gets silently dropped, which would otherwise quietly reduce effective coverage right at the region edges — worth specifically checking for in smoke-test QC (e.g., compare read counts / properly-paired fraction between a tight slice and a wide slice for the same sample).

---

## What to try first, in order

1. **Check the JupyterLab Snippets menu** (user is doing this in parallel) — if it has a genomic-manifest snippet, that supersedes everything in section 1 of this doc; use its exact table/column names and bucket paths over any guesses here.
2. **Locate and read the actual manifest CSV** for your `C2025Q4R6` cohort (via Dataset Builder → Genomics domain export, or wherever the Snippets point you) to get real `person_id → cram_uri`/`cram_index_uri` (short-read) and the lrWGS equivalent. Do not hand-guess bucket paths from this document — the `v6_base` example path is very likely wrong for your CDR.
3. **On one participant's srWGS CRAM**, first just confirm read access and header contents locally/remotely:
   ```bash
   samtools view -H "<cram_uri_from_manifest>"
   ```
   This alone will tell you: (a) whether remote `gs://` header access works at all without auth setup, (b) the exact contig naming convention (`chr6` vs `6`), (c) whether the `UR:` reference tag is populated and what it points to.
4. **Attempt the region-slice directly on that same CRAM** with explicit `-T` pointing at the public Broad hg38 FASTA (section 3). If this errors or hangs, fall back immediately to `gsutil cp` the CRAM+CRAI locally (AoU's own documented workaround for the IGV case) and retry the slice locally — don't burn time debugging remote streaming if the perimeter/exfiltration controls are actually in play here as they are for IGV.
5. **Once one srWGS slice works, repeat for one lrWGS BAM** for the *same* participant (confirming it's the GRCh38-aligned BAM, not the CHM13 one, per section 2) — this should be simpler since BAM needs no reference flag.
6. **Run `samtools fastq` on both slices** (section 5) and sanity-check read counts/pairing before handing off to SpecHLA/SpecImmune's own extraction scripts (`ExtractHLAread.sh`/`ExtractReads.sh`) — or consider skipping straight to those scripts, since they may duplicate/replace your manual slice-then-fastq steps. Worth reading `ExtractHLAread.sh` and `ExtractReads.sh` source directly to see whether they expect a pre-sliced region or the whole CRAM/BAM (their usage strings suggest whole-file input mapped to hg38, which would make your region-slice step complementary/optional-but-faster rather than required).
7. **Only after all of the above works for 1 person**, scale to the full 2-3 person smoke test, then to ~100.

---

## Sources referenced (informational, not exhaustive — see inline confidence markers above)
- support.researchallofus.org: "How the All of Us Genomic data are organized" (CDRv8 and prior archived versions), "Using the Integrative Genomics Viewer (IGV)", "Controlled CDR Directory" (fetch blocked by 403; used via search snippets only), "Accessing Genomic Data in the All of Us Controlled Tier"
- PMC12622093 / medRxiv 2025.10.02.25336942 — "Population-scale Long-read Sequencing in the All of Us Research Program"
- PMC10822842 / Nature Communications 2024 — "Utility of long-read sequencing for All of Us" (earlier technical pilot, not the production dataset)
- github.com/deepomicslab/SpecHLA and github.com/deepomicslab/SpecImmune READMEs
- samtools/htslib GitHub issues on CRAM reference resolution (`UR:` tag, `REF_CACHE`/`REF_PATH`, `-T`), samtools-fastq/samtools-fasta man pages, `bam_fastq.c`
- Terra (support.terra.bio) community post on samtools GCS auth patterns (`HTS_AUTH_LOCATION`, `GCS_OAUTH_TOKEN`) — Terra-specific, used here as an architectural analogy, not AoU-confirmed
- cloufield.github.io/gwaslab HLA region documentation (MHC/xMHC coordinate conventions)

**Caveat on staleness:** several `support.researchallofus.org` articles are explicitly versioned/archived per CDR release (v5/v6/v7/v8), confirming AoU's own docs churn with each data release — and none of what's cited here was verified against the *Verily Workbench* rebrand UI specifically, since all indexed docs found were written in the AoU Researcher Workbench-era terminology. Treat UI navigation paths as the most likely thing to have changed; treat underlying GCS bucket/manifest *mechanisms* and htslib/samtools *behavior* as most likely to still be accurate, since those are lower in the stack than the UI rebrand.
