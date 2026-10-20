# S06 VM commands (run in order; each step checkpoints and can be re-run)

Assume a fresh VM: remount the bucket first (ENVIRONMENT quirk #14) and check with `ls`. Start the terminal with
`exec > >(tee -a ~/.claude_session.log) 2>&1` (quirk #31).

```bash
cd ~/repos/pilot-validation && git fetch origin && git checkout s06-call-quality && git pull --ff-only
mkdir -p ~/s06/{logs,results,internal}
PY="pixi run -e spechla -- python3"

# 0. the truth-set people, so 70 measures them first (needs 56/57 inputs; ~1 min)
$PY scripts/hla_popgen/71_call_quality_table.py --emit-priority

# 1. exon / CDS features from GENCODE (public coordinates) -- check the printed "not on chr6 primary" list
$PY scripts/hla_popgen/70_call_exon_depth.py --build-features \
    --gencode-gtf ~/ref/gencode.v46.basic.annotation.gtf.gz --out-dir ~/s06/results/70

# 2. smoke test, 20 people, own internal path (quirk #22); look at ~/s06/internal/70_smoke/*.tsv
$PY scripts/hla_popgen/70_call_exon_depth.py --out-dir ~/s06/results/70 \
    --internal-dir ~/s06/internal/70_smoke --limit 20 --workers 8

# 3. full run (checkpointed, resumable; truth-set people first)
setsid nohup $PY -u scripts/hla_popgen/70_call_exon_depth.py --out-dir ~/s06/results/70 \
    --internal-dir ~/s06/internal/70 --priority-file ~/s06/internal/71/priority_ids.txt --workers 8 \
    > ~/s06/logs/70.log 2>&1 < /dev/null & disown
tail -2 ~/s06/logs/70.log

# 4. population depth profile along each gene (300 random people, ~10 min; separate invocation)
$PY scripts/hla_popgen/70_call_exon_depth.py --out-dir ~/s06/results/70 --profile 300 --workers 8

# 5. E7: the extra genes for the trio-truth contig matrix (reruns yak per child; run after 3 or while it runs if CPUs allow)
setsid nohup $PY -u scripts/hla_popgen/57b_trio_gene_origin.py --all-trios --extra-genes DRA,C4A --workers 3 \
    > ~/s06/logs/57b_extra.log 2>&1 < /dev/null & disown

# 6. analysis (partial 70 output is fine; rerun when 70 finishes)
$PY scripts/hla_popgen/71_call_quality_table.py
$PY scripts/hla_popgen/71b_call_quality_stats.py
$PY scripts/hla_popgen/71c_contig_truth_matrix.py
```

Aggregate outputs to pull (small TSVs): `~/s06/results/{70,71,71c}/*.tsv`. Never pull `~/s06/internal`.
