# 59 — depth x outcome join + coverage cutoff (VM, Workbench JupyterLab terminal).
# Paste one block at a time; `type` and `Return` are separate actions when driven by computer-use.
# Prerequisite: 58's full run has finished (~/s05/internal/58/per_person_depth.tsv complete, no
# 58.lock present), and the repo on the VM has 59 pulled.

# 0. Session log + sanity (assume the VM restarted; ENVIRONMENT.md quirks #14, #35)
exec > >(tee -a ~/.claude_session.log) 2>&1
cd ~/repos/pilot-validation && git log -1 --oneline
ls scripts/hla_popgen/59_depth_outcome_join.py
timeout 10 ls ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv
ls -la ~/s05/internal/58/per_person_depth.tsv ~/pipeline_outputs/hla_calls_rich.tsv \
       ~/pipeline_outputs/cohort_membership.tsv
test -e ~/s05/internal/58/58.lock && echo "58 STILL RUNNING - wait" || echo "58 done"
ls ~/pipeline_outputs_kir | head -3
mkdir -p ~/s05/logs ~/s05/results ~/s05/internal
# statsmodels is not in the spechla env (Python 3.8); 59 falls back to its own logistic fit.
pixi run -e spechla -- python3 -c "import statsmodels; print('statsmodels', statsmodels.__version__)" \
    || echo "no statsmodels: numpy fallback will be used (engine column says numpy_irls)"

# 1. Smoke test (own output + internal dirs, quirk #22), ~1-2 min
pixi run -e spechla -- python3 -u scripts/hla_popgen/59_depth_outcome_join.py \
    --people-limit 300 --workers 4 \
    --out-dir ~/s05/results/59_smoke --internal-dir ~/s05/internal/59_smoke \
    > ~/s05/logs/59_smoke.log 2>&1; echo "exit $?"
tail -5 ~/s05/logs/59_smoke.log
column -t -s$'\t' ~/s05/results/59_smoke/run_summary.tsv | head -50
column -t -s$'\t' ~/s05/results/59_smoke/proposed_cutoff.tsv | cut -c1-220 | head -30

# 2. Full run (detached; KIR GTF parse of ~12k people is the slow part, a few minutes)
setsid nohup pixi run -e spechla -- python3 -u scripts/hla_popgen/59_depth_outcome_join.py \
    --out-dir ~/s05/results/59 --internal-dir ~/s05/internal/59 --workers 8 \
    > ~/s05/logs/59.log 2>&1 < /dev/null & disown

# 3. Progress / completion (re-run until "done:" appears; a rerun after a VM stop resumes the KIR cache)
tail -3 ~/s05/logs/59.log
ls -la ~/s05/results/59/

# 4. Pull back PUBLIC tables only (never ~/s05/internal/59): results/59/*.tsv
#    outcome_by_depth_bin.tsv hla_gene_absent_by_bin.tsv cutoff_scan.tsv logistic_or.tsv
#    contigs_linear.tsv proposed_cutoff.tsv run_summary.tsv
