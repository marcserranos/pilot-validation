# Local environment (Aleix's laptop, WSL2) — quirks & fixes

> Same role as Marc's `../context/ENVIRONMENT.md`, but for the **local WSL2 box** where Phase 1 runs.
> **Edit:** append a quirk the moment you hit one; never rewrite an existing entry.
> Marc's quirks are for the AoU Workbench VM and mostly do NOT apply here (notably: **we have sudo, he doesn't**).

## Layout

- `~/repos/pilot-validation` — the repo (cloned into WSL's native filesystem, **not** `/mnt/c/...`, which is 5–20x slower for file-heavy work)
- `~/repos/pilot-validation/aleix/pixi.toml` — HLA-Resolve environment (separate from Marc's root manifest)
- `~/tools/SpecImmune` — upstream source + `db/` (HLA database, built locally)
- `~/tools/hla_resolve` — upstream source
- `~/apptainer_tmp` — apptainer/singularity scratch (see quirk 3)

Machine: 16 cores, 15.7GB physical RAM (WSL2 capped — see quirk 4), ~935GB free disk.

## Quirks & fixes

1. **`unzip` and `makeblastdb` are missing by default; SpecImmune's `make_db.py` needs both.**
   `unzip` fails *silently* — it prints `sh: 1: unzip: not found` to stderr but the script continues and only dies later with a confusing `FileNotFoundError: ./db/HLA/hla_gen.fasta`, because the download succeeded and only the extraction failed. `makeblastdb` is worse: it prints "not found" and the script **completes without error**, reporting success while the BLAST database was never built.
   Fixed locally with `sudo apt install -y unzip ncbi-blast+`.
   **This will NOT work on the AoU Workbench VM (no sudo, Marc's quirk #6)** — there, add `unzip` and `blast` to the pixi manifest instead. Same lesson as Marc's #17/#18: a clean-looking log is not a success signal.

2. **`pip` is not included with conda-forge's `python`.** `pixi.toml` must list `pip` explicitly or `python -m pip` fails with `No module named pip`.

3. **`/tmp` is a 3.9GB tmpfs, independent of the 935GB disk.** Broke apptainer's DeepVariant container pull with `no space left on device` while `df -h ~` showed 935G free. Fixed by pointing apptainer at real disk — four exports now persisted in `~/.bashrc`:
   `APPTAINER_TMPDIR`, `APPTAINER_CACHEDIR`, `SINGULARITY_TMPDIR`, `SINGULARITY_CACHEDIR` → `~/apptainer_tmp`.
   **Still an open risk:** DeepVariant also stages `make_examples` `.tfrecord` intermediates in `/tmp` (seen at `/tmp/tmpng1bp0cj/`). It fit for the small demo; **whole-genome samples may not.** Slice to chr6 before running, and redirect DeepVariant's temp dir if it recurs.

4. **WSL2 caps RAM at ~50% of physical by default (~7.8GB here) — this hard-crashed the entire WSL VM**, not just the process, during HLA-Resolve's whole-genome alignment step (terminal dropped straight back to PowerShell). Fixed via `%UserProfile%\.wslconfig`:
   ```
   [wsl2]
   memory=12GB
   swap=16GB
   ```
   then `wsl --shutdown` and relaunch. Now shows 11Gi + 16Gi swap. Swap is disk-backed and effectively free given available disk — cheap insurance against a repeat.

5. **`git clone` over HTTPS can fail with `GnuTLS handshake failed`** (a WSL2 network-layer quirk). Fixed globally:
   `git config --global http.version HTTP/1.1` and `git config --global http.postBuffer 1048576000`.

6. **A conda `(base)` env auto-activates in this shell.** It silently hijacked `pip install -e .`, installing HLA-Resolve into `~/miniconda3` instead of the pixi env — the tool then couldn't see pixi's `pysam`. **Always invoke via `pixi run --manifest-path <path> -- python -m pip ...`**, never a bare `pip`, and verify with `pixi run ... -- which hla_resolve` (must resolve inside `.pixi/envs/`, not `miniconda3`).

## Measured timings

| Run | Data | Time |
|---|---|---|
| HLA-Resolve, demo | 88,955 reads, hybrid-capture, `--threads 12` | 21m25s (DeepVariant `make_examples` alone: 11m18s wall / 87m54s CPU) |
