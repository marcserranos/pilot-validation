# aleix/ — HLA-Resolve workstream

Isolated workspace for Aleix's tasks, kept separate from Marc's Experiments A–D (`../context/`, `../scripts/`). The science, environment quirks, and prior results still live in `../context/` — read those for *why*; this folder is *my* new work.

> Data rule unchanged from the root repo: **no AoU participant data here.** Public reference material (the Lai truth table, tool docs) is fine — it is not participant data.

---

## The experiment, stated once

**Experiment D + a 4th arm.** Marc's Experiment D compared 3 callers on 60 AoU people at the 8 classical HLA genes: AoU-native (SR), SpecHLA (SR), SpecImmune (LR). I add **HLA-Resolve** (a 2nd long-read caller, PacBio HiFi) as a 4th arm, on the same cohort — and, critically, I anchor the whole thing to a **real ground truth**, which Marc never used (he treated SpecImmune-LR as "presumed truth", `../context/DECISIONS.md`).

Goal: not "pick a winner" but a **per-tool scorecard** — accuracy, ancestry-fairness, resolution, completeness, biological plausibility, compute cost — to decide which long-read caller is best *for the project* (trustworthy + unbiased across ancestry + affordable at 400k scale).

## Two phases

| | Phase 1 — Calibrate | Phase 2 — The experiment |
|---|---|---|
| Where | External, public data (off-AoU) | Inside AoU Workbench |
| Samples | HPRC/GIAB truth set (dial 3 → 44) | The 60 AoU people (10 × 6 ancestries) |
| Truth? | **Yes** — Lai et al. 2023, 4-field | No per-participant truth exists in AoU |
| Answers | Is each LR tool actually *correct*? Ancestry-fair accuracy? | Which is best *for us*: bias, completeness, cost |

Phase 1 is a **trust anchor**, not the verdict. Both callers ride through Phase 2 regardless of Phase 1 — a tool that ties on truth may win on ancestry-fairness or cost.

## Phase 1 truth set — `reference/Lai_Supplementary-6.xlsx`

Ships in the HLA-Resolve repo (`docs/`), copied here. **44 GIAB/HPRC samples**, 2 alleles/gene at **4-field**, 8 classical genes + DRB3/4/5. HLA-Resolve only *reported* on 3 (HG002, HG01258, HG03579); truth covers all 44, and HPRC hosts **public PacBio HiFi** for ~all of them (`s3://human-pangenomics`, no creds; HG002 via GIAB FTP). HPRC is ancestry-diverse → Phase 1 gives **ancestry-stratified accuracy against real truth**, the one axis the AoU-60 can't. Many are also 1000G samples with public **short-read** → can truth-test the SR arms too.

## HLA-Resolve — github.com/matthewglasenapp/hla_resolve

PacBio HiFi (ONT in dev). Input: single-sample FASTQ or unmapped BAM. Same 8 classical genes, IPD-IMGT/HLA. Install: conda `environment.yml` → `pip install -e .` → `hla_resolve setup`. Run: `hla_resolve --input_file X.bam --sample_name X --platform pacbio --scheme wgs --output_dir out --threads N`. AoU LR is WGS revio HiFi → **use `--scheme wgs`** (not the demo's `hybrid_capture`). Repo ships `demo/HG002.hifi_reads.fastq.gz` for a smoke test.

## Inherited config (settled by Marc — do NOT re-derive)

- SpecHLA (SR) → **pad10000**
- SpecImmune (LR) → **pad100k, `--align_method_1 minimap2`**
- No gene-panel restriction (Exp C: slower + miscalls at A/DRB1)
- Slicing window `chr6:29,500,000–33,500,000`; reuse `../scripts/slice_and_fastq.sh` pattern

## Status / next

- [x] Truth set located + inspected (44 samples, 4-field) → `reference/Lai_Supplementary-6.xlsx`
- [x] **Machine decided: Aleix's own Windows laptop, via WSL2.** Phase 1 is public data — no AoU compliance issue, no need for the shared Workbench. SpecImmune's own README confirms **"supports Linux and Windows WSL systems"** — officially supported, not a hack. Rebuilding SpecImmune from this repo's `pixi.toml` + `../context/ENVIRONMENT.md` runbook; not reusing Marc's build (separate per-user Workbench environments, no shared disk).
- [x] **Package manager: pixi, not conda.** Supervisor preference; pixi consumes the same bioconda/conda-forge + pip deps both tools' READMEs list (`environment.yml` for HLA-Resolve, `environment.yml` for SpecImmune) via `[pypi-dependencies]`, same pattern Marc's root `pixi.toml` already uses.
- [x] **First batch decided: HG002, HG01258, HG03579** (not a disjoint 6). Checked SpecImmune's own `evaluation/hprc_year1_level1_aws_locs.csv` (its real published validation set, 38 HPRC samples) — all 3 of HLA-Resolve's reported samples are in it, and all 3 are in the Lai truth. A genuine 3-way intersection: both tools were separately validated on these exact people, against the same ground truth. Scale toward the rest of the 44 later for the ancestry-fairness read.
- [x] **WSL2 + pixi 0.73.0 working.** Config gotchas hit and fixed, all recorded in `ENVIRONMENT_LOCAL.md`.
- [x] **SpecImmune built + HLA DB constructed** (62 dirs in `db/HLA`, 66 in `db/HLA_CDS`).
- [x] **HLA-Resolve installed** via `aleix/pixi.toml` (its `environment.yml` translated 1:1, plus `pip`).
- [x] **Smoke test PASSED — 16/16 alleles exact at 4-field vs Lai truth on HG002.** Demo file, 88,955 reads, 21m25s. Proves the install works; does NOT prove real-world accuracy (it's the tool's own curated hybrid-capture demo). See `results/hg002_demo_vs_truth.md`.
- [x] **Workbench (Phase 2 machine) fully provisioned** — see "Workbench state" below.
- [x] **Remote slicing solved — the single most important practical finding.** See "Data access" below.
- [ ] Slice HG002 → FASTQ, run HLA-Resolve on real (non-demo) data, score vs Lai
- [ ] SpecImmune smoke test on laptop (`test/HLA/test_HLA_lite.fastq.gz` + shipped `test.HLA.hap.alleles.txt`) — **never yet run; it is half the head-to-head and is the tool that exits 0 on failure**
- [ ] Locate HG03579's aligned BAM (HG002 + HG01258 found; see below)
- [ ] Scorer: normalize to 2/3/4-field, compare vs Lai as **unordered** pairs
- [ ] Repeat for all 3 samples, both callers

---

## Data access — remote slicing (solved 2026-07-18)

The raw files named in HLA-Resolve's README are **unaligned** and enormous: HG002 43.6GB, HG01258 **175.9GB**, HG03579 34.1GB (~254GB total). Unaligned means no region-slicing, so both callers would have to align every read genome-wide before filtering to chr6 — extrapolating from the demo, **~5h for HG002 and ~20h+ for HG01258, per tool.** Not viable on a laptop.

**Fix: use the *aligned, indexed* BAMs published alongside them.** With a `.bai`, samtools issues HTTP range requests and fetches only the HLA window — a few hundred MB instead of 254GB, and both tools then only ever see HLA reads.

| Sample | Aligned BAM | Size | Status |
|---|---|---|---|
| HG002 | `.../PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam` (GIAB FTP) | 70GB | ✅ verified working |
| HG01258 | `HG01258_aligned_GRCh38_winnowmap.sorted.bam` (`s3://human-pangenomics`, https endpoint) | 97GB | found, exact S3 key still needed |
| HG03579 | — | — | **not yet located** |

Verified on HG002: `samtools view -c <url> chr6:29500000-33500000` → **12,986 reads** in seconds. At ~17kb/read over the 4Mb window that's **~55× coverage**, consistent with the file's advertised 48× genome-wide. Plenty for confident typing.

Ignore HG01258's `.maternal`/`.paternal.GRCh38_no_alt.bam` (~2GB) — those are *assembly contigs* aligned to the reference, not reads.

**Slice recipe** (`-F 0x900` drops secondary/supplementary so no read is duplicated in the FASTQ):
```bash
samtools view -b -F 0x900 -o <out>.chr6.bam "<url>" chr6:29500000-33500000
samtools fastq <out>.chr6.bam | gzip > <out>.chr6.fastq.gz
```

**Open caveat, must be tested not assumed:** HG002's BAM uses `GRCh38_GIABv3_no_alt_analysis_set_...` — **no ALT contigs** (218 `@SQ`, none HLA). Good news: no HLA reads hiding on ALT contigs, so the chr6 slice is complete *with respect to what mapped*. Bad news: this is exactly Marc's open `grch38_noalt` question (`../context/DECISIONS.md`) — reads from highly divergent HLA alleles may have failed to map at all. `pbmm2` was run with `--unmapped`, so such reads are retained in the BAM's unmapped pool, not lost. **Test empirically:** HLA-Resolve scored ~100% at 1–3 field on HG002 using the full raw file; if the chr6-only slice reproduces 16/16 vs truth, slicing is provably lossless here. If it degrades, add unmapped reads.

---

## IMGT/HLA version tracking

**Both machines: IPD-IMGT/HLA 3.65.0** (dated 2026-07-14, pulled from ANHIG `Latest`). SpecImmune DB verified identical on both: 62 dirs in `db/HLA`, 66 in `db/HLA_CDS`, 46,647 sequences.

Marc logged a **DB-version confound** in `../context/EXPERIMENTS.md` — his SpecImmune (3.64.0) vs SpecHLA (3.38.0) produced apparent discordances that were really newer-release *renaming*. Three deltas to keep in view:
- **vs Marc's 3.64.0** — his stored Experiment D results aren't strictly comparable to ours without accounting for renames.
- **vs HLA-Resolve's own bundled IMGT** — version *not yet determined*. This one directly affects the Phase 1 head-to-head: if the two callers ship different releases, some "disagreements" are naming, not calling.
- **vs the Lai truth** (2023/24, older release) — alleles renamed or split since would score as false misses. The 16/16 smoke test suggests this isn't badly broken, but that's n=1.

Note `Latest` drifts — a rebuild months from now could silently shift allele names. Version is pinned in writing here for that reason.

---

## Workbench state (Phase 2 machine)

Instance: **AoU Jupyter** app, own instance (not Marc's `HLAcalling_pilot_v0_r` — collaborators get separate per-user environments, no shared disk). 500GB disk, autostop 4h idle, `us-central1-a`.

Provisioned and verified:
- Network egress **works** despite the workspace's perimeter warning (github / conda.anaconda.org / raw.githubusercontent all reachable) — Phase 2 builds are viable here.
- **No sudo** (`jupyter is not in the sudoers file`, Marc's quirk #6) — but `unzip` (`/usr/bin/unzip`) and `makeblastdb` (`/opt/workbench-tools/binaries/bin/`) already exist **and are visible inside the pixi env** (verified with `pixi run -e specimmune -- which`). This is why Marc could build the DB without sudo; it costs us nothing.
- Disk topology differs from the laptop: `/`, `$HOME`, `/tmp`, `/opt` are **one 492G overlay** (428G free) — no separate small tmpfs, so the DeepVariant `/tmp` failure can't recur the same way. But it's the same topology behind Marc's quirk #19: when his disk filled, `/tmp` filled too and pixi activation itself started failing. ~5× his headroom, same discipline required.
- pixi 0.73.0; repo cloned on branch `aleix/hla-resolve-phase1`.
- `specimmune` env + SpecImmune source + HLA DB built (62/66, IMGT 3.65.0).
- HLA-Resolve installed via `aleix/pixi.toml`; `which` resolves inside `.pixi/envs/`, not `/opt/conda`; `hla_resolve setup` completed (DeepVariant SIF pulled). **The container-in-container risk did not materialize** — apptainer works inside the Dockerized Jupyter app despite no sudo.

Still to do on the Workbench: **SpecHLA** (`spechla` env + `bash index.sh`), **gcsfuse mount** to the AoU bucket (Marc's quirk #11), and **Marc's 60-person `cohort.tsv`** so Phase 2 runs on the same people.

**Sequencing:** the Workbench is *infrastructure readiness only*. **Phase 1's result gates whether Phase 2 is worth running** — if HLA-Resolve can't beat SpecImmune against real truth, adding it to the AoU cohort answers a far less interesting question.
