# aleix/RNA-seq/ — immune repertoire workstream

Isolated workspace for the RNA-seq / TCR-BCR repertoire direction. HLA calling is handed off
to Marc; this folder is the new work. Kept separate from `../` (the HLA-Resolve workstream)
and from `../../context/` (Marc's Experiments A–D) — same isolation rule as before.

> **Data rule unchanged from the root repo: no AoU participant data here.** Public reference
> material (TRUST4's own small db files, tool docs) is fine — it is not participant data.

Confidence key, used throughout this folder same as `../reference/AOU_DATA_ACCESS_NOTES.md`:
**[HIGH]** confirmed live or from primary-source docs. **[MED]** inferred from credible
sources, not yet verified live. **[LOW]** best guess, needs live verification.

---

## The direction, stated once

AoU RNA-seq BAMs contain TCR/BCR (T-cell / B-cell receptor) reads that **STAR discards**
during alignment, because V(D)J-recombined receptor sequences don't exist in any reference
genome. That means the immune repertoire is completely absent from every expression file
AoU ships — nobody has mined it. **TRUST4** recovers it directly from the BAM via de novo
assembly.

HLA type is predictable from TCR repertoire alone (AUC 0.95, beta chain, published result) —
so our own HLA calls (AoU-native + long-read, from the HLA-Resolve workstream) become the
**labels** for this direction, rather than being discarded. Full reasoning, the three
supervisor papers, and the two-directions comparison (expression vs. repertoire) are in the
`RNAseq_directions*.pptx` decks in `slides/`.

**RNA-seq manifest — found 2026-08-10:**
`gs://vwb-aou-datasets-controlled/v9/multiomics/rnaseq/manifest.tsv`
(columns: `sampleid`, `research_id`, `markduplicates_bam_file_path`,
`markduplicates_bam_index_path`). Sibling dir `v9/multiomics/` also holds `proteomics/`.
The `rnaseq/` dir itself has `manifest.tsv`, `rnaseq_metadata.tsv`, `eqtl/`, `rnaseqc2/`,
`rsem/`, `sqtl/` — matches the deliverables table in the primary-source PDF. BAMs are
mark-duplicates-applied (`.md.bam`), physically under
`pooled/multiomics/v9_base/rnaseq/bam/`.

**Feasibility — confirmed 2026-08-10.** TRUST4 recovers a real repertoire from an AoU
whole-blood RNA-seq BAM: 2,013 CDR3s from person 1000291, all 7 chain types represented,
~9 min end-to-end. Full writeup: `results/1000291_trust4_smoke_test.md`. Required a
`--abnormalUnmapFlag` fix (now in `scripts/run_trust4_sample.sh`) because AoU's STAR run
uses `--outSAMunmapped Within`, which TRUST4 doesn't handle by default.

## How we work together on this (the git loop)

1. Scripts and docs are authored **here**, in this repo, on the local machine (with Claude).
2. Commit + push from here.
3. On the Workbench VM: `git pull` inside `~/repos/pilot-validation` to get the latest.
4. Run scripts from `~/repos/pilot-validation/aleix/RNA-seq/scripts/` on the VM.
5. Anything the VM produces that's worth keeping (numbers, short summaries, not raw
   BAMs/FASTQs) gets written back into `results/*.md` and pushed from the VM, or pasted back
   here — either way it ends up in git, not just on a VM that can be deleted.

No more pasting full command blocks back and forth — the VM always has the current script via
`git pull`.

## Layout

- **`scripts/`** — everything runnable. Read a script before running it; nothing here is
  proven yet (see confidence tags inside each file).
- **`reference/`** — TRUST4's small db files once fetched (`hg38_bcrtcr.fa`,
  `human_IMGT+C.fa` — a few MB each, not participant data, safe to vendor once we have them).
- **`results/`** — `results/*` is gitignored except `.md`/`.csv` — same policy as
  `../results/`. Bulky outputs (BAMs, assembled contigs) stay on the VM; written-up numbers
  come back here.
- **`slides/`** — the progress decks for this direction.

## Machine

New, separate VM from the HLA one (`HLAcalling_pilot_v0_m`) — own siloed app instance,
matching the existing two-collaborator pattern. TRUST4 is CPU/IO-bound, no GPU, far lighter
than the long-read HLA stack (no DeepVariant/sawfish/pbsv/sniffles). Recommended spec:
4-8 vCPU, 16-25 GB RAM, ~100 GB disk, Ubuntu 24.04 — see chat log for full reasoning.

## Status checklist

- [x] VM created (`n1-highmem-8`, 8 vCPU / 52 GB RAM, 100 GB disk)
- [x] repo cloned on VM, gcsfuse mounted (same recipe as the HLA workstream)
- [x] pixi env built (`pixi install` against this folder's `pixi.toml`)
- [x] RNA-seq manifest located — `v9/multiomics/rnaseq/manifest.tsv`
- [x] TRUST4 reference files fetched, vendored into `reference/`
- [x] one participant's RNA BAM path resolved (1000291)
- [x] TRUST4 run on that one sample — needed `--abnormalUnmapFlag`, now fixed in the script
- [x] CDR3 count recorded in `results/` — **2,013, feasibility confirmed**
- [ ] pick a real ~25-person pilot cohort via `ancestry_preds.tsv` (not another arbitrary first-row pick)
- [ ] check whether recovery is even across ancestry groups
- [ ] attach HLA labels (AoU-native + our own long-read calls) once repertoires exist for the pilot
