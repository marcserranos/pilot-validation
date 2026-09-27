# scripts/legacy/

Pre-numbered-convention code, kept for provenance rather than deleted. Everything here predates
the `scripts/hla_popgen/NN_description.py` pipeline convention (see the repo `README.md`) and is
superseded by it — nothing under `scripts/legacy/` is part of the live pipeline, and nothing on
the VM's production path calls into it.

- **`experiments/`** — the 37 one-off Experiment A/B/C/D-era scripts (`analyze_*.py`, `build_*.py`,
  `experiment_*.py`, `run_experiment_*.sh`, `compare_hla_results.py`, the `*_pad_sweep.sh/py`
  scripts, etc.). Moved here in the 2026-09-27 WS-F reorg from a loose `scripts/` top level (see
  `context/EXPERIMENTS.md`'s matching entry for the full old→new path table, and
  `sprints/S04_kir_recurrence_style_share/WS-F_REORG_PLAN.md` for the reference audit that
  justified moving them as a block). Several of these scripts are still cited by path from
  `reports/*/README.md` (as the "how this result was produced" record) and from
  `context/ENVIRONMENT.md`'s SR/LR runbook — those references were updated to the new path, so the
  commands there still work as written; the scripts themselves are otherwise not maintained.

Not here despite being candidates: `scripts/run_immuannot_person.py`, `scripts/bootstrap_vm.sh`,
`scripts/setup_immuannot.sh`, `scripts/spechla_pad_helpers.py` — these have real VM-production
wiring (hardcoded paths/`sys.path` imports from the orchestrator and fresh-VM bootstrap) and were
deliberately left at the `scripts/` top level. See the WS-F reorg plan's §3.1 for the dependency
evidence.
