# 60 — TRBV germline copy number vs "never used" (VM, Workbench JupyterLab terminal).
# Paste one block at a time; `type` and `Return` are separate actions when driven by computer-use.
# Prerequisite: 58's full run finished (~/s05/internal/58/trb_cn_per_person.tsv exists) and/or
# 65's run finished (~/s05/internal/65/per_person_gene.tsv); the repo on the VM has 60 pulled.

# 0. Session log + sanity
exec > >(tee -a ~/.claude_session.log) 2>&1
cd ~/repos/pilot-validation && git log -1 --oneline
ls scripts/hla_popgen/60_trb_germline_usage.py
ls -la ~/s05/internal/58/trb_cn_per_person.tsv ~/s05/internal/65/per_person_gene.tsv 2>&1
mkdir -p ~/s05/logs ~/s05/results ~/s05/internal

# 1. Find Aleix's per-person TRBV usage table. The script searches ~/workspace/aleix_rnaseq/
#    (depth <= 4) for *usage*; with 0 or >1 matches it exits and lists the candidates.
#    Expected: the W1 cache vj_usage.pkl ({"V": person x gene clonotype counts}); it may only exist
#    on Aleix's own app instance (~/pipeline_outputs/rnaseq/analysis/), in which case ask him to copy
#    it to the shared bucket.
timeout 20 ls ~/workspace/aleix_rnaseq/ | head
pixi run -e spechla -- python3 -u scripts/hla_popgen/60_trb_germline_usage.py \
    --cn-table ~/s05/internal/58/trb_cn_per_person.tsv --people-limit 300 \
    --out-dir ~/s05/results/60_smoke --internal-dir ~/s05/internal/60_smoke 2>&1 | tail -30
# -> if it ran, the single match was used; otherwise set USAGE to the right file from the list:
USAGE=~/workspace/aleix_rnaseq/REPLACE_WITH_PATH
# For a long table also add: --id-col <id> --gene-col <gene> --value-col <count>

# 2. Smoke test, 58 copy numbers (+ the labelled insertion proxy), own dirs (quirk #22)
pixi run -e spechla -- python3 -u scripts/hla_popgen/60_trb_germline_usage.py \
    --cn-table ~/s05/internal/58/trb_cn_per_person.tsv --usage "$USAGE" --people-limit 300 \
    --out-dir ~/s05/results/60_smoke --internal-dir ~/s05/internal/60_smoke \
    > ~/s05/logs/60_smoke.log 2>&1; echo "exit $?"
tail -8 ~/s05/logs/60_smoke.log
column -t -s$'\t' ~/s05/results/60_smoke/run_summary.tsv
# check: people_in_both_tables > 0 (else the usage id is not the AoU person_id), genes_joined > 0

# 3. Full run, 58 copy numbers (a few minutes; foreground is fine)
pixi run -e spechla -- python3 -u scripts/hla_popgen/60_trb_germline_usage.py \
    --cn-table ~/s05/internal/58/trb_cn_per_person.tsv --usage "$USAGE" \
    --out-dir ~/s05/results/60 --internal-dir ~/s05/internal/60 > ~/s05/logs/60.log 2>&1; echo "exit $?"

# 4. Full run, 65 assembly copy numbers (direct calls for TRBV4-3/3-2/6-3/7-8), when 65 is done.
#    Separate dirs so the two sources never overwrite each other.
pixi run -e spechla -- python3 -u scripts/hla_popgen/60_trb_germline_usage.py \
    --cn-table ~/s05/internal/65/per_person_gene.tsv --cn-source 65 --usage "$USAGE" \
    --haplotype-genes "TRBV4-3,TRBV3-2,TRBV6-3,TRBV7-8,TRBV6-2,TRBV12-4,TRBV6-2/6-3,TRBV6-2/6-3[copy2]" \
    --out-dir ~/s05/results/60_cn65 --internal-dir ~/s05/internal/60_cn65 > ~/s05/logs/60_cn65.log 2>&1; echo "exit $?"

# 5. Pull back PUBLIC tables only (never ~/s05/internal/60*): results/60/*.tsv and results/60_cn65/*.tsv
#    contingency.tsv diagnostic.tsv dosage.tsv haplotype_genes.tsv run_summary.tsv
