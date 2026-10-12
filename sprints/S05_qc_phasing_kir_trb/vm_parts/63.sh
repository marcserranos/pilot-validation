# S05 63 — HLA + KIR calls for the sequel2 (Tier-3, self_align_needed) people. VM commands.
# Paste ONE block at a time and check its output before the next (ENVIRONMENT.md quirk #2).
# Not meant to run as a whole script. Never touch ~/mnt/aou-controlled (hangs, quirk #35).
# Writes ONLY to ~/s05/internal/63* (INTERNAL, never pulled) and ~/s05/results/63* (PUBLIC aggregates).
# ~/pipeline_outputs* is read, never written (63 refuses any output path inside it).
# Diagnosis + read-only pre-checks: sprints/S05_qc_phasing_kir_trb/SEQUEL2_DIAGNOSIS.md §5.

# --- 0. Session log, env, mount, code (assume the VM restarted, quirk #14) ----------------------
exec > >(tee -a ~/.claude_session.log) 2>&1
PY=~/repos/pilot-validation/.pixi/envs/specimmune/bin/python3
ENVBIN=$(dirname $PY)
cd ~/repos/pilot-validation-main && git pull --ff-only && git log -1 --format='%h %s'   # note in LOG.md
timeout 10 ls ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv
ls -la ~/pipeline_outputs/immuannot_cohort_full.tsv ~/pipeline_outputs/immuannot_calls.tsv \
       ~/pipeline_outputs_kir/immuannot_calls.kir.tsv
ls -la ~/ref/Homo_sapiens_assembly38.fasta ~/ref/Homo_sapiens_assembly38.fasta.fai
ls ~/tools/Immuannot ~/tools/Immuannot_refdata | head
$ENVBIN/minimap2 --version; $ENVBIN/samtools --version | head -1; $PY -c "import pandas; print(pandas.__version__)"
nproc; free -g | head -2; df -h ~ | tail -1
mkdir -p ~/s05/results ~/s05/internal ~/s05/logs
test -e ~/s05/internal/63/63.lock && echo "63 lock present: check 'ps -p \$(cat ~/s05/internal/63/63.lock)'" || echo "no 63 lock"

# --- 0b. Read-only diagnosis checks U1-U7 (paste the block from SEQUEL2_DIAGNOSIS.md §5) --------
#     Must see before going on: tier3 people with calls = 0; sequel2 hap FASTAs plain_gzip (or bgzf:
#     either works); bucket root = gs://vwb-aou-datasets-controlled/... ; .fai present.
#     If the reference is missing (bootstrap_vm.sh does not fetch it):
# mkdir -p ~/ref && gcloud storage cp \
#   gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
#   gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai ~/ref/

# --- 1. Hashes of the production tables (to prove later they were not touched) ------------------
md5sum ~/pipeline_outputs/immuannot_calls.tsv ~/pipeline_outputs_kir/immuannot_calls.kir.tsv \
       ~/pipeline_outputs/hla_calls_rich.tsv > ~/s05/logs/63_prod_md5_pre.txt; cat ~/s05/logs/63_prod_md5_pre.txt

# --- 2. References (once; ~5-15 min, ~5-8 GB RAM for the index; ~3.5 GB disk) -------------------
$PY scripts/hla_popgen/63_sequel2_recall.py --prepare-refs 2>&1 | tail -8
ls -la ~/s05/internal/63/ref/
head -1 ~/s05/internal/63/ref/hla_stage1.fa ~/s05/internal/63/ref/kir_stage1.fa   # >chr6:28500001-34500000 / >chr19:53600001-55920000
grep -c '^>' ~/s05/internal/63/ref/hg38_primary.fa                               # expect 25

# --- 3. Smoke test: 2 people, HLA + KIR, own dirs (quirk #22). ~1.5 h on 4 vCPU (estimate) -------
setsid nohup $PY -u scripts/hla_popgen/63_sequel2_recall.py --smoke --workers 1 --threads 4 \
    > ~/s05/logs/63_smoke.log 2>&1 < /dev/null & disown
tail -3 ~/s05/logs/63_smoke.log                       # re-run until "done:" appears
column -t -s$'\t' ~/s05/results/63_smoke/run_summary.tsv
column -t -s$'\t' ~/s05/results/63_smoke/trim_summary.tsv
column -t -s$'\t' ~/s05/results/63_smoke/failure_reasons.tsv
# Per-haplotype detail (INTERNAL, look only, do not copy out): sizes, candidates, timings, genes
column -t -s$'\t' ~/s05/internal/63_smoke/hla/immuannot_timing.sequel2.tsv | cut -c1-250
column -t -s$'\t' ~/s05/internal/63_smoke/kir/immuannot_timing.sequel2.kir.tsv | cut -c1-250
cut -f2 ~/s05/internal/63_smoke/hla/immuannot_calls.sequel2.tsv | sort | uniq -c | head -40
grep -h "stage 1\|stage 2\|SKIP" ~/s05/internal/63_smoke/*/pipeline_outputs/*/immuannot_output/.63_last_run.log | cut -c1-200
# GO criteria (all): both people done for both loci; trim_method self_align_confirmed; HLA padded
# ~3-8 MB/hap and KIR ~0.5-3 MB/hap; HLA-A/B/C/DRB1/DQB1 called on both haps; stage-1 candidates
# tens, not thousands. Record the median stage-1/stage-2/immuannot seconds in LOG.md and replace
# the runtime estimate in ledger_parts/63.md. If a hap SKIPs, the reason is in failure_reasons.tsv
# and in the .63_last_run.log above: STOP and report, do not tune thresholds blind.

# --- 4. Full run (991 people; estimate ~500-800 h on 4 vCPU at 1 worker x 4 threads, see ledger) --
#     Only after the smoke passes AND Marc has chosen the box (stay on 4 vCPU, or resize: his call).
#     On 4 vCPU / 31 GB: --workers 1 --threads 4 (CPU-bound; 2 workers x 2 threads is equivalent).
#     On a resized app: --workers ~ min(vCPU/4, RAM_GB/8), --threads 4.
setsid nohup $PY -u scripts/hla_popgen/63_sequel2_recall.py --workers 1 --threads 4 \
    > ~/s05/logs/63.log 2>&1 < /dev/null & disown

# --- 5. Progress (any time; resumable: after a VM stop, re-run the block-4 command unchanged) ----
tail -2 ~/s05/logs/63.log
column -t -s$'\t' ~/s05/results/63/run_summary.tsv      # refreshed every 50 jobs

# --- 6. After it finishes: merge is automatic; re-run on its own if needed ----------------------
$PY scripts/hla_popgen/63_sequel2_recall.py --merge-only 2>&1 | tail -4
ls -la ~/s05/internal/63/merged/
md5sum -c ~/s05/logs/63_prod_md5_pre.txt               # production tables must still be OK
# Retry the failed person-loci once more (after reading failure_reasons.tsv):
# $PY -u scripts/hla_popgen/63_sequel2_recall.py --retry-failed --workers 1 --threads 4 > ~/s05/logs/63_retry.log 2>&1

# --- 7. Pull back PUBLIC only: ~/s05/results/63/{run_summary,trim_summary,failure_reasons}.tsv ---
#     Never ~/s05/internal/63. Counts 1-19 are already "<20".
