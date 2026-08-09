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

**Step 0, unresolved:** the RNA-seq manifest (`research_id -> RNA BAM path`) has not been
located yet. Every other AoU data type in this project (srWGS CRAM, lrWGS BAM, HLA calls,
ancestry) resolves through a `v9/...manifest.*` file — the RNA-seq one almost certainly
follows the same pattern but the exact path is unconfirmed. Finding it is the actual first
task on the machine. `scripts/find_rnaseq_manifest.sh` is a discovery helper, not an answer.

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

- [ ] VM created
- [ ] repo cloned on VM, gcsfuse mounted (same recipe as the HLA workstream)
- [ ] pixi env built (`pixi install` against this folder's `pixi.toml`)
- [ ] RNA-seq manifest located
- [ ] TRUST4 reference files fetched, vendored into `reference/`
- [ ] one participant's RNA BAM path resolved
- [ ] TRUST4 run on that one sample
- [ ] CDR3 count recorded in `results/` -- the feasibility answer
