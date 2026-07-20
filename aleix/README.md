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
- [ ] Pull HG002 (GIAB FTP) + HG01258, HG03579 (HPRC S3, no creds) HiFi BAMs; region-slice to chr6 HLA
- [ ] Run HLA-Resolve (`--scheme wgs`) + SpecImmune (minimap2, pad100k) on all 3
- [ ] Scorer: normalize calls to 2/3/4-field, compare vs Lai, per gene per tool
