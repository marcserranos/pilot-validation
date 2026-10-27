# reports/ — index

| Folder | What | Era |
|---|---|---|
| [`pre_sprint/`](pre_sprint/) | Seven early reports: `aou_callset_validation`, `confidence_matched_truth`, `disease_sanity_check`, `experiment_d_field_cascade`, `full_immuannot_lr_calling`, `immuannot_pilot`, `lr_data_census`. Written by `scripts/legacy/` and `scripts/production_*`. | July–August 2026 |
| [`hla_popgen/`](hla_popgen/) | One folder per numbered pipeline script (`scripts/hla_popgen/NN_*.py` → `reports/hla_popgen/NN_*/`). Its README has a by-era table (pre-sprint, S01–S09). | September 2026 → |

Gitignored and private (not in the repo): `walkthrough_S03_S04/`, `walkthrough_S05/` (built on private call notes).
Sprint-level finished deliverables live next to their sprint (e.g. `sprints/S05_qc_phasing_kir_trb/S05_REPORT.html`).

Why `hla_popgen/` was not split by sprint: ~40 scripts hard-code `reports/hla_popgen/NN_*` as their output path
(`_viz_common.DEFAULT_REPORT_ROOT` is the single choke point) and the VM runs rely on it.
