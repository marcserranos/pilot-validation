# S05 65 — TRB insertion genotype (reads mode on everyone; assembly mode on a stratified 2,000; 65b).
# Paste ONE block at a time and check its output before the next (ENVIRONMENT.md quirk #2).
# Not meant to run as a whole script. Filesystem probes of the mount are wrapped in `timeout 10`
# (quirk #35). PUBLIC -> ~/s05/results/65*, INTERNAL -> ~/s05/internal/65* (never pulled).
# 58 may still be running: check `pgrep -af 58_regional_depth` and keep total workers within the vCPUs.

# --- 0. Environment, mount, code ----------------------------------------------------------------
PY=~/repos/pilot-validation/.pixi/envs/spechla/bin/python3
ENVBIN=~/repos/pilot-validation/.pixi/envs/spechla/bin          # samtools + minimap2 2.28
REPO=~/repos/pv-65            # the 65 worktree (10-11); ~/repos/pilot-validation-main if 65 was merged there
cd $REPO && git log -1 --format='%h %s'                          # note the commit in LOG.md
timeout 10 ls ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv
mkdir -p ~/s05/results ~/s05/internal ~/s05/logs ~/ref/trb

# --- 1. Aligner ---------------------------------------------------------------------------------
# Option A: mappy for python 3.8 in its own env, put on PYTHONPATH (done 10-11; re-check):
ls ~/envs/mappy38/lib/python3.8/site-packages | grep -i mappy || \
  mamba create -y -p ~/envs/mappy38 -c conda-forge -c bioconda python=3.8 mappy
export PYTHONPATH=~/envs/mappy38/lib/python3.8/site-packages
$PY -c "import mappy; print('mappy', mappy.__version__)"
# Option B (no mappy): add to every 65 command below:  --aligner minimap2 --minimap2 $ENVBIN/minimap2
$ENVBIN/minimap2 --version && $ENVBIN/samtools --version | head -1

# --- 2. Reference: ALT contig + UCSC chr7 (there is no GRCh38 FASTA on the VM) -------------------
curl -L -m 300 -s https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/altLocations.txt.gz | zcat \
  | awk -F'\t' '$2=="chr7" && $3<142900000 && $4>142200000'    # expect chr7_KI270803v1_alt, 142038120-143088503
test -s ~/ref/trb/trb_alt.fa || { curl -L -m 300 -s -o ~/ref/trb/chr7_KI270803v1_alt.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7_KI270803v1_alt.fa.gz && \
  zcat ~/ref/trb/chr7_KI270803v1_alt.fa.gz > ~/ref/trb/trb_alt.fa; }
curl -L -m 300 -s -o ~/ref/trb/chr7.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
ls -la ~/ref/trb/ && zcat ~/ref/trb/chr7.fa.gz | head -c 100

# --- 3. Rebuild probe metadata (same probe set as 10-11) + breakpoint -----------------------------
cp -n ~/s05/results/65/probes.fa ~/ref/trb/probes_imgt_20261011.fa
cd $REPO/scripts/hla_popgen
$PY 65_trb_assembly_genotype.py --build-probes --probe-fasta ~/ref/trb/probes_imgt_20261011.fa \
    --out-dir ~/s05/results/65 --alt-fasta ~/ref/trb/trb_alt.fa \
    --ref-fasta ~/ref/trb/chr7.fa.gz --ref-region chr7:142038120-143088503
cat ~/s05/results/65/reference_check.tsv ~/s05/results/65/probe_meta.tsv
# PASS = TRBV4-3/3-2/7-8 present on alt, absent on primary; 6-2/6-3 n_loci 2 vs 1; a breakpoint line
# with ins_len_estimate ~ 20-22 kb, block_classes TRBV3-2,TRBV4-3 and classes_elsewhere TRBV7-8
# (GRCh38's ALT holds TRBV7-8 on a separate ALT-only segment, so reads mode leaves it out).
# If no breakpoint: stop and report (do not guess one).

# --- 4. Which BAM column has HP tags? (needed for per-haplotype calls / co-presence) ---------------
M=~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv
timeout 10 head -1 $M | tr '\t' '\n' | grep -n -E 'bam'
B=$(timeout 10 awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="grch38_bam")c=i;next}$c!=""{print $c;exit}' $M \
    | sed 's#gs://[^/]*/#'$HOME'/workspace/vwb-aou-datasets-controlled-v9/#')
echo $B && timeout 60 $ENVBIN/samtools view $B chr7:142400000-142400100 | head -3 | grep -o 'HP:i:[0-9]' | head -1
# no HP printed -> add  --bam-column grch38_haplotagged_bam  to the reads commands below (repeat the check on it)

# --- 5. Reads mode smoke (20 people, own paths, quirk #22) ----------------------------------------
cd $REPO/scripts/hla_popgen
$PY -u 65_trb_assembly_genotype.py --mode reads --samtools $ENVBIN/samtools --limit 20 --workers 4 \
    --out-dir ~/s05/results/65_reads_smoke --internal-dir ~/s05/internal/65_reads_smoke 2>&1 | tail -8
cat ~/s05/results/65_reads_smoke/run_summary.tsv ~/s05/results/65_reads_smoke/failures_summary.tsv
# On the VM only (INTERNAL): reads per person and seconds -> sets --workers for the full run
cut -f2,6,8,11-13,15 ~/s05/internal/65_reads_smoke/per_haplotype.tsv | head -12
awk -F'\t' '$2=="insertion_any"' ~/s05/internal/65_reads_smoke/per_person_gene.tsv | cut -f2,7,10,13-15 | head -20

# --- 6. Assembly mode on a stratified 2,000 (validation; I/O bound, ~1-2 days at 4 workers) -------
cd $REPO/scripts/hla_popgen
$PY 65_trb_assembly_genotype.py --mode assembly --sample-stratified 2000 --seed 65 \
    --people-file ~/s05/internal/65_people_strat2000.tsv
setsid nohup $PY -u 65_trb_assembly_genotype.py --mode assembly --people-file ~/s05/internal/65_people_strat2000.tsv \
    --workers 4 > ~/s05/logs/65_assembly.log 2>&1 < /dev/null & disown
tail -2 ~/s05/logs/65_assembly.log

# --- 7. Reads mode on everyone (resumable: rerun the same line after a stop) ----------------------
cd $REPO/scripts/hla_popgen
setsid nohup $PY -u 65_trb_assembly_genotype.py --mode reads --samtools $ENVBIN/samtools \
    --workers 12 > ~/s05/logs/65_reads.log 2>&1 < /dev/null & disown
tail -2 ~/s05/logs/65_reads.log

# --- 8. After both finish: unrelated tables + 65b concordance -------------------------------------
cd $REPO/scripts/hla_popgen
$PY 65_trb_assembly_genotype.py --mode reads --aggregate-only --unrelated --out-dir ~/s05/results/65_reads_unrelated
$PY 65_trb_assembly_genotype.py --mode assembly --aggregate-only --people-file ~/s05/internal/65_people_strat2000.tsv
$PY 65_trb_assembly_genotype.py --compare ~/s05/internal/65_reads/per_person_gene.tsv \
    ~/s05/internal/65_assembly/per_person_gene.tsv --out-dir ~/s05/results/65b
column -t ~/s05/results/65b/65b_method_concordance.tsv | head -30
ls -la ~/s05/results/65_reads ~/s05/results/65_reads_unrelated ~/s05/results/65_assembly ~/s05/results/65b
# Pull ONLY ~/s05/results/65* (PUBLIC). Check every <20 is in place before committing (CLAUDE.md:
# small-cell disclosure is an open item in context/DECISIONS.md). Never pull ~/s05/internal/65*.
