# VM scratch-cleanup dry run (S04, read-only)

Generated 2026-09-26 on `AoU_Jupyter_ComputeEngine_20260805_big_run`. **Nothing was deleted.**
IDs masked throughout (no person directory names or 6+-digit IDs listed; per-person dir counts
are given as plain integers, never names). Hard-excluded paths (never proposed for deletion,
per CLAUDE.md / ORCHESTRATOR_HANDOFF.md): `~/pipeline_outputs/*.tsv`, `~/pipeline_outputs_kir/`,
`~/s04/results/`, `~/repos/`, `~/tools/`, `~/mnt/`, `~/workspace/`.

Raw data: `sprints/S04_kir_recurrence_style_share/cleanup_dryrun.tsv`.

## Candidates

| Path | Size | Files | Action | Reason | Referenced by |
|---|---|---|---|---|---|
| `~/s03/results/41/pipeline_outputs` | 1.2G | 343 person-dirs | **DELETE** | KIR-pilot scoping scratch copy, superseded by `~/pipeline_outputs_kir`; only appears as an example `--outroot` in 41's docstring, not an active default | none (docstring example only) |
| `~/s03/_patch2.py`, `_patch3.py`, `_patch4.py`, `_patch5.py`, `_patch_local.py` | ~20K | 5 | **DELETE** | one-off VM hotfix scratch scripts from S03 | none |
| `~/s03/_chunks` | 116K | few | **DELETE** | chunked-file-transfer deploy scratch (PUT-splitting artifacts) | none |
| `~/s03/_parts` | 224K | few | **DELETE** | same as `_chunks` | none |
| `~/s03/__pycache__` | 660K | - | **DELETE** | regenerable bytecode cache | none |
| `~/s03/kir_full` | 1.6M | many small | **DELETE** | S03 full-KIR-run scratch (logs, status files, `cohort_bench200.tsv`/`cohort_kir_main.tsv`/`cohort_test2.tsv` benchmark cohorts); aggregate already committed via 43 | none |
| `~/s03/results/38` | 368K | - | **DELETE** | S03 diagnostic scratch (frag_check/hist_export/refigure/deletion_validation), already pulled to repo | none |
| `~/s03/results/39_full` | 12M | - | **DELETE** | 39's full-cohort intermediate export; 44 re-derives via its own `m39()` call, does not read these files | none |
| `~/s03/results/40` | 44K | - | **DELETE** | already pulled to repo | none |
| `~/s03/results/42` | 296K | - | **DELETE** | already pulled to repo (42_repertoire_baseline) | none |
| `~/s03/results/43` | 44K | - | **DELETE** | already pulled to repo (43_kir_full_aggregate committed outputs) | none |
| `~/s03/results/37` | 2.4M | - | **KEEP** | explicit instruction: keep 37 (DQ/DP phasing-confidence artifact) | `49_cole_share_package.py`'s `--phasing-confidence` default path lives here (the file itself is currently absent, but the dir is kept) |
| `~/pipeline_outputs_kir_test` | 13M | 44 | **DELETE** | test/scratch KIR run, name says test, not part of the real 12,261-person cohort | none |
| `~/pipeline_outputs/*.corrupted-backup-*` | 544K+22M (~22.5M) | 2 | **DELETE** | historical corruption-incident backups (2026-08-10, 2026-08-11) of `immuannot_calls.tsv`, superseded by the current clean file (the `*.tsv` itself is never touched) | none |
| `~/s04/share_release_2026-09` | 11M | 6 | **DELETE** | Cole package now uploaded + verified in the bucket (`release_2026-09-25`); this local copy is redundant | none |
| `~/.cache/pip` | 347M | - | **DELETE** | pip download/wheel cache, regenerable | none |
| `~/.cache/rattler` | 3.5G | - | **DELETE** | pixi/rattler package cache, regenerable — single largest reclaim item | none |
| `~/allele_geometry_scratch` | 3.2M | 5 | **DELETE** | unreferenced viz-iteration scratch (10_allele_ancestry_geometry work), no script imports it | none |
| `~/results` (top-level, pre-sprint) | 32M | - | **FLAG** | pre-S01 pilot/scratch dirs (`15_hla_manhattan`, `24_novelty_by_field`, `25_noncoding_novelty_paf`, `26_qc_relatives_v2`, `27_allele_space_coverage`, `29_hla_ld_by_ancestry`, `30_hla_structural_variation`, `31_aa_diversity_selection`, `36_figure1`, plus loose `allele_frequency`/`clustering`/`completeness`/`confidence`/`drb1_capstone` dirs); predates the `reports/` convention. **Not explicitly requested — needs Marc's confirmation before deleting.** | unknown, not grepped exhaustively |
| `~/pipeline_outputs/rnaseq/ws6_scratch` | 4.0K | 0 (empty) | **KEEP** | negligible size, empty placeholder | `42_repertoire_baseline.py` |
| `~/repos/pilot-validation/scripts/production_orchestrator/run_production_orchestrator_s03kir.py` | 52K | 1 | **KEEP** | untracked file on the VM repo branch (`needle-view-cds-diversity-density`); keep until that branch merges | none (untracked, VM-repo-only) |

## Total reclaimable (confirmed DELETE rows only, excludes the FLAGged `~/results`)

**~5.1 GB**, dominated by `~/.cache/rattler` (3.5G) and `~/s03/results/41/pipeline_outputs` (1.2G).
Adding the flagged `~/results` (32M, needs Marc's OK) brings it to ~5.15 GB.

## To run after Marc's OK

Explicit paths only, no globs over IDs. Two backup files' exact timestamped names given in full;
everything else is a fixed directory/file path.

```bash
rm -rf ~/s03/results/41/pipeline_outputs
rm -f ~/s03/_patch2.py ~/s03/_patch3.py ~/s03/_patch4.py ~/s03/_patch5.py ~/s03/_patch_local.py
rm -rf ~/s03/_chunks
rm -rf ~/s03/_parts
rm -rf ~/s03/__pycache__
rm -rf ~/s03/kir_full
rm -rf ~/s03/results/38
rm -rf ~/s03/results/39_full
rm -rf ~/s03/results/40
rm -rf ~/s03/results/42
rm -rf ~/s03/results/43
rm -rf ~/pipeline_outputs_kir_test
rm -f ~/pipeline_outputs/immuannot_calls.tsv.corrupted-backup-20260810-150440
rm -f ~/pipeline_outputs/immuannot_calls.tsv.corrupted-backup-20260811-215920
rm -rf ~/s04/share_release_2026-09
rm -rf ~/.cache/pip
rm -rf ~/.cache/rattler
rm -rf ~/allele_geometry_scratch
```

`~/results` (top-level, 32M) is deliberately **not** in this list — flagged for Marc's separate
confirmation, not pre-approved here.

## Never proposed (hard exclusions, confirmed untouched by this dry run)

`~/pipeline_outputs/*.tsv`, `~/pipeline_outputs_kir/`, `~/s04/results/`, `~/repos/`, `~/tools/`,
`~/mnt/`, `~/workspace/`.
