# reports/hla_popgen/_briefs/

Narrative, read-on-demand docs for the `hla_popgen` pipeline that aren't a per-script `NN_*/`
result writeup and aren't pipeline code — a feature brief and an operator runbook. They belong in
the `reports/` tier per the repo's four-tier split (`scripts/` is "used, not read"; `reports/` is
"read on demand"), not at the `scripts/hla_popgen/` top level next to the live `README.md` and
`SCHEMA.md`. The leading underscore distinguishes this from the numbered `NN_description/` result
folders it sits alongside.

Moved here from `scripts/hla_popgen/` in the 2026-09-27 WS-F reorg (see `context/EXPERIMENTS.md`'s
matching entry for the old→new path). Bare-filename mentions of these two docs in pipeline script
docstrings/comments and in `context/STATUS.md` were updated to the new path.

- **`NEEDLE_VIEW_BRIEF.md`** — feature brief for the CDS-diversity "needle view" track
  (`scripts/hla_popgen/20_gene_diversity_track.py`): design rationale, the `cs`-string walking
  logic, and live-VM validation notes.
- **`RUNBOOK.md`** — operator runbook for running the `hla_popgen` pipeline on the Workbench VM:
  mount steps, the `people/` on-disk layout move (2026-09-04) and why it matters for the Jupyter
  file browser, and per-script gotchas referenced throughout `scripts/hla_popgen/*.py`.
