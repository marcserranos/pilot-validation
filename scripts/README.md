# scripts/ — index

| Path | What | Status |
|---|---|---|
| [`hla_popgen/`](hla_popgen/) | The live numbered pipeline (`NN_description.py`), tests, research notes. Start at its README. | live |
| [`production_analysis/`](production_analysis/), [`production_orchestrator/`](production_orchestrator/), [`monitoring/`](monitoring/) | Production-run analysis, orchestrator and monitoring from the July–August full long-read run. `hla_popgen/01`, `06` and `tests/make_fixtures.py` still import from `production_*`. | partly live |
| [`legacy/`](legacy/) | 37 one-off Experiment A–D scripts, kept for provenance. | frozen |
| `bootstrap_vm.sh`, `setup_immuannot.sh` | VM bootstrap (called from the preflight runbook). Kept at this level: runbooks cite these paths. | live |
| `run_immuannot_person.py` | Per-person Immuannot runner (used by `hla_popgen/63`, cited in ENVIRONMENT.md). | live |
| `spechla_pad_helpers.py` | Helper for the `scripts/legacy/` pad-sweep scripts. | legacy |
