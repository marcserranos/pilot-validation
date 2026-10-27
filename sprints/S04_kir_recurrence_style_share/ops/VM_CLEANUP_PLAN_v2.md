# VM cleanup plan v2: organise first, delete only evident trash

This supersedes the delete list in VM_CLEANUP_DRYRUN.md (kept as the inventory). Marc's steer (09-27): *"lower threshold for deletion but a higher organization or cleanup into folders… if something is evidently trash or corrupted or has no possible use or reuse, of course that can be deleted."*

## Target VM home layout

| Path | Role | Rule |
|---|---|---|
| `~/pipeline_outputs/` | HLA production (per-person tree + `*.tsv` tables) | read-only, never touched |
| `~/pipeline_outputs_kir/` | KIR production (12,261 people) | read-only |
| `~/s03/`, `~/s04/` | sprint workdirs: deployed scripts at top level (on `PYTHONPATH`), `results/NN*/` per script | scripts stay put (commands depend on the paths) |
| `~/archive/` | anything with possible reuse that is not live: dated, with a README | new |
| `~/tools/`, `~/repos/`, `~/mnt/`, `~/workspace/` | tooling, VM repo (other branch), mounts | untouched |
| `~/VM_LAYOUT.md` | one-page map of the above | new |

## Delete (evident trash: corrupted, regenerable, or an exact duplicate), ~0.4 GB

| Path | Why |
|---|---|
| `~/pipeline_outputs/immuannot_calls.tsv.corrupted-backup-20260810-150440`, `…-20260811-215920` | corrupted-incident backups, superseded by the clean table |
| `~/s03/_patch2.py … _patch5.py, _patch_local.py`, `~/s03/_chunks`, `~/s03/_parts` | one-off hotfix / deploy-transfer fragments |
| `~/s03/__pycache__`, `~/s04/__pycache__` | regenerable bytecode |
| `~/s04/share_release_2026-09/` | byte-identical copy of the verified bucket release (re-buildable with 49) |
| `~/.cache/pip` | regenerable download cache |

## Archive (move into `~/archive/`, nothing lost)

| From | To |
|---|---|
| `~/results/` (pre-S01 pilot outputs) | `~/archive/pre_s01_results/` |
| `~/allele_geometry_scratch/` | `~/archive/scratch/allele_geometry/` |
| `~/pipeline_outputs_kir_test/` | `~/archive/kir_test_run/` |

## Keep in place (organised already, or reused)

`~/s03/results/*` (per-script folders; 37 is the phasing-confidence location 49 reads, 41 holds the KIR-pilot outputs), `~/s03/kir_full/` (S03 full-run logs = audit trail), `~/.cache/rattler` (3.5 GB pixi cache; avoids re-downloading on env rebuild), the untracked `run_production_orchestrator_s03kir.py` in the VM repo (until the branch merges).

## Commands (Marc runs; the classifier blocks agent-side deletes)

```bash
cd ~ && mkdir -p ~/archive/scratch
mv ~/results ~/archive/pre_s01_results
mv ~/allele_geometry_scratch ~/archive/scratch/allele_geometry
mv ~/pipeline_outputs_kir_test ~/archive/kir_test_run
rm -f ~/pipeline_outputs/immuannot_calls.tsv.corrupted-backup-20260810-150440 ~/pipeline_outputs/immuannot_calls.tsv.corrupted-backup-20260811-215920
rm -f ~/s03/_patch2.py ~/s03/_patch3.py ~/s03/_patch4.py ~/s03/_patch5.py ~/s03/_patch_local.py
rm -rf ~/s03/_chunks ~/s03/_parts ~/s03/__pycache__ ~/s04/__pycache__ ~/s04/share_release_2026-09 ~/.cache/pip
md5sum ~/pipeline_outputs/*.tsv | diff - ~/s04/hla_tsv_md5_pre.txt && echo HLA_TABLES_UNCHANGED
```
Then paste `~/VM_LAYOUT.md` (the next block) into the VM:
```bash
cat > ~/VM_LAYOUT.md <<'MD'
# VM home layout (updated 2026-09-27, S04)
- pipeline_outputs/      HLA production: per-person immuannot tree + *.tsv tables. READ-ONLY.
- pipeline_outputs_kir/  KIR production (chr19), 12,261 people. READ-ONLY.
- s03/, s04/             sprint workdirs: deployed scripts at top level (on PYTHONPATH), results/NN*/ per script.
- archive/               non-live material kept for possible reuse: pre_s01_results/, kir_test_run/, scratch/.
- tools/                 immuannot + reference data (IPD-IMGT/HLA 3.55.0, IPD-KIR 2.13.0).
- repos/pilot-validation VM clone (branch owned by other work; do not checkout/commit here).
- mnt/, workspace/       bucket mounts.
Share bucket: gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/ (Cole package).
MD
```
