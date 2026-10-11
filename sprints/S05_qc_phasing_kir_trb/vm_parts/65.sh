# S05 65 — TRB insertion/deletion genotype from phased assemblies. VM commands.
# Paste ONE block at a time and check its output before the next (ENVIRONMENT.md quirk #2).
# Not meant to run as a whole script. Every exploratory filesystem probe is wrapped in `timeout 10`
# (quirk #35). PUBLIC -> ~/s05/results/65*, INTERNAL -> ~/s05/internal/65* (never pulled).

# --- 0. Environment, mount, code ----------------------------------------------------------------
PY=~/repos/pilot-validation/.pixi/envs/spechla/bin/python3
ENVBIN=~/repos/pilot-validation/.pixi/envs/spechla/bin
cd ~/repos/pilot-validation-main && git log -1 --format='%h %s'      # note the commit in LOG.md
timeout 10 ls ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv
timeout 10 head -1 ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv \
  | tr '\t' '\n' | grep -n -E 'assembly_hap[12]_(fa|aln2_hg38_paf)|platform|research_id'
mkdir -p ~/s05/results ~/s05/internal ~/s05/logs ~/ref/trb

# --- 1. Aligner: mappy in the env (preferred), else the minimap2 binary --------------------------
$PY -c "import mappy; print('mappy', mappy.__version__)" || $PY -m pip install mappy
$PY -c "import mappy; print('mappy', mappy.__version__)"
ls $ENVBIN/minimap2 && $ENVBIN/minimap2 --version            # pixi.toml lists minimap2 in spechla
# Fallback if neither works (same pattern as ~/envs/yak):
#   mamba create -y -p ~/envs/mm2 -c conda-forge -c bioconda minimap2 samtools
#   then add to every 65 command below:  --aligner minimap2 --minimap2 ~/envs/mm2/bin/minimap2

# --- 2. Reference check inputs: the GRCh38 ALT contig carrying the segment + primary TRB window ---
# Expected ALT name: chr7_KI270803v1_alt. CONFIRM it here (the row whose placement spans the TRB locus):
curl -s https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/altLocations.txt.gz | zcat \
  | awk -F'\t' '$2=="chr7" && $3<142900000 && $4>142200000'
ALT=chr7_KI270803v1_alt     # <- replace with the name printed above if different
curl -s -o ~/ref/trb/$ALT.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$ALT.fa.gz
zcat ~/ref/trb/$ALT.fa.gz > ~/ref/trb/trb_alt.fa && grep -c . ~/ref/trb/trb_alt.fa
# Primary window from the GRCh38 FASTA already on the VM (find it first):
ls ~/ref/
REF=~/ref/REPLACE_WITH_GRCh38.fa   # <- the noalt/primary GRCh38 FASTA listed above (needs its .fai)
$ENVBIN/samtools faidx $REF chr7:142250000-142860000 > ~/ref/trb/trb_primary.fa && head -c 200 ~/ref/trb/trb_primary.fa

# --- 3. Probes (IMGT V-REGION + Ensembl TRBC1 genomic) + reference check -------------------------
# First eyeball the IMGT page: headers must look like >acc|TRBV4-3*01|Homo sapiens|F|V-REGION|...
curl -s 'https://www.imgt.org/genedb/GENElect?query=7.2+TRBV&species=Homo+sapiens' | grep -m3 '>'
cd ~/repos/pilot-validation-main/scripts/hla_popgen
$PY 65_trb_assembly_genotype.py --build-probes --out-dir ~/s05/results/65 \
    --alt-fasta ~/ref/trb/trb_alt.fa --ref-fasta ~/ref/trb/trb_primary.fa
cat ~/s05/results/65/reference_check.tsv
grep -E 'TRBV(4-3|3-2|6-2|6-3|7-8)' ~/s05/results/65/probe_classes.tsv
# PASS = TRBV4-3 and TRBV3-2 present on 'alt', absent on 'primary'; TRBV6-2/TRBV6-3 n_loci 2 on alt, 1 on primary.
# If IMGT fails: download the TRBV V-REGION FASTA by hand, append TRBC1 genomic, then
#   $PY 65_trb_assembly_genotype.py --build-probes --probe-fasta ~/ref/trb/my_probes.fa --out-dir ~/s05/results/65 ...

# --- 4. Smoke test (own paths, quirk #22) --------------------------------------------------------
cd ~/repos/pilot-validation-main/scripts/hla_popgen
$PY -u 65_trb_assembly_genotype.py --probes ~/s05/results/65/probes.fa --limit 20 --workers 4 \
    --out-dir ~/s05/results/65_smoke --internal-dir ~/s05/internal/65_smoke 2>&1 | tail -20
cat ~/s05/results/65_smoke/run_summary.tsv ~/s05/results/65_smoke/failures_summary.tsv
column -t ~/s05/results/65_smoke/haplotype_completeness.tsv
# per-haplotype seconds (INTERNAL, read on the VM only) -> sets --workers/--timeout for the full run:
cut -f1-7,11 ~/s05/internal/65_smoke/per_haplotype.tsv | head -25

# --- 5. Full run (background, survives the terminal; resumable: rerun the same line after a stop) -
cd ~/repos/pilot-validation-main/scripts/hla_popgen
setsid nohup $PY -u 65_trb_assembly_genotype.py --probes ~/s05/results/65/probes.fa \
    --workers 12 > ~/s05/logs/65.log 2>&1 < /dev/null & disown
tail -3 ~/s05/logs/65.log

# --- 6. After it finishes: unrelated-only tables (no re-genotyping) ------------------------------
$PY 65_trb_assembly_genotype.py --probes ~/s05/results/65/probes.fa --aggregate-only --unrelated \
    --out-dir ~/s05/results/65_unrelated
ls -la ~/s05/results/65 ~/s05/results/65_unrelated
# Pull ONLY ~/s05/results/65* (PUBLIC). Check every <20 is in place before committing (CLAUDE.md:
# small-cell disclosure is an open item in context/DECISIONS.md). Never pull ~/s05/internal/65*.
