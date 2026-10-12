# S05 64 — Novel HLA alleles: current IPD-IMGT/HLA release + read-level support. VM commands.
# Paste ONE block at a time and check its output before the next (ENVIRONMENT.md quirk #2).
# Not meant to run as a whole script. Every exploratory filesystem probe is wrapped in `timeout 10`
# (quirk #35). PUBLIC -> ~/s05/results/64*, INTERNAL -> ~/s05/internal/64* (never pulled).

# --- 0. Environment, mount, code ----------------------------------------------------------------
PY=~/repos/pilot-validation/.pixi/envs/spechla/bin/python3
ENVBIN=~/repos/pilot-validation/.pixi/envs/spechla/bin
cd ~/repos/pilot-validation-main && git log -1 --format='%h %s'      # note the commit in LOG.md
ls scripts/hla_popgen/64_novel_allele_validation.py
timeout 10 ls ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv
timeout 10 head -1 ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv \
  | tr '\t' '\n' | grep -n -E 'grch38_bam|platform|research_id'
ls -la ~/pipeline_outputs/hla_calls_rich.tsv ~/pipeline_outputs/cohort_membership.tsv
ls -la ~/s05/results/58/regions.bed      # GRCh38 gene spans from 58 (build with 58 --build-regions if absent)
mkdir -p ~/s05/results ~/s05/internal ~/s05/logs ~/s05/ref

# --- 1. Tools: mappy in the env, samtools in the env ---------------------------------------------
$PY -c "import mappy; print('mappy', mappy.__version__)" || $PY -m pip install mappy
$PY -c "import mappy; a=mappy.Aligner(seq='ACGT'*500, preset='map-hifi'); print('map-hifi ok' if a else 'map-hifi MISSING (script falls back to map-pb)')"
$ENVBIN/samtools --version | head -1

# --- 2. Immuannot per-haplotype files the script reads (one person) ------------------------------
P=$(ls ~/pipeline_outputs/people | head -1); D=~/pipeline_outputs/people/$P/immuannot_output
ls -la $D/hap1.trimmed.fa $D/hap1.gtf.gz $D/hap1/cds.fa.gz
zcat $D/hap1.gtf.gz | awk -F'\t' '$3=="CDS" && /gene_name "HLA-A"/' | head -3     # CDS rows exist
zcat $D/hap1/cds.fa.gz | grep '>' | head -3                                       # >{contig}_{gene}_{i}

# --- 3. IPD-IMGT/HLA current release (GitHub ANHIG/IMGTHLA, ref Latest) --------------------------
cd ~/repos/pilot-validation-main/scripts/hla_popgen
$PY 64_novel_allele_validation.py --fetch-ipd --ipd-release Latest --ipd-dir ~/s05/ref/ipd_latest
ls -la ~/s05/ref/ipd_latest/
head -c 300 ~/s05/ref/ipd_latest/hla_nuc.fasta; echo        # must start '>HLA:HLA00001 A*01:01:01:01 ...'
grep -m2 '^#' ~/s05/ref/ipd_latest/Allelelist.txt            # version line, e.g. '# version: IPD-IMGT/HLA 3.6x.0'
grep -v '^#' ~/s05/ref/ipd_latest/Allelelist_history.txt | head -2 | cut -c1-200   # header must list 3550
# If a file is a Git-LFS pointer the script already re-fetched it from media.githubusercontent.com.
# If GitHub fails entirely: download hla_nuc.fasta / hla_gen.fasta / Allelelist*.txt by hand into
# ~/s05/ref/ipd_latest/ and continue (the script only reads files from --ipd-dir).

# --- 4. Re-match only, classical genes, no BAM access (fast; answers part a) --------------------
$PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest --no-reads \
    --genes A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 --workers 8 \
    --out-dir ~/s05/results/64_ipd --internal-dir ~/s05/internal/64_ipd 2>&1 | tail -8
cat ~/s05/results/64_ipd/run_summary.tsv ~/s05/results/64_ipd/failures_summary.tsv
column -t ~/s05/results/64_ipd/ipd_release.tsv | head -15
grep -E 'classical_8' ~/s05/results/64_ipd/funnel_with_evidence.tsv | column -t
# CHECK: failures_summary call_status 'ok' should dominate; many 'cds_disagree' means the GTF-rebuilt
# CDS and cds.fa.gz differ systematically -> inspect one case before going on (do NOT add
# --ignore-cds-check blindly).

# --- 5. Smoke test with reads (own paths, quirk #22) ---------------------------------------------
cd ~/repos/pilot-validation-main/scripts/hla_popgen
$PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest \
    --genes A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 --limit 30 --workers 4 --samtools $ENVBIN/samtools \
    --out-dir ~/s05/results/64_smoke --internal-dir ~/s05/internal/64_smoke 2>&1 | tail -8
cat ~/s05/results/64_smoke/run_summary.tsv ~/s05/results/64_smoke/failures_summary.tsv
# per-call timing and read counts (INTERNAL, read on the VM only) -> sets --workers/--timeout:
grep -v "^#" ~/s05/internal/64_smoke/per_call_evidence.tsv | cut -f4,12,14,15,26,28-37,48 | column -t | head -25

# --- 6. Optional: short-read check on the smoke set (needs ~/ref/Homo_sapiens_assembly38.fasta) --
timeout 10 head -2 ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/cram/manifest.csv
ls -la ~/ref/Homo_sapiens_assembly38.fasta ~/ref/Homo_sapiens_assembly38.fasta.fai
$PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest --sr \
    --genes A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 --limit 30 --workers 4 --samtools $ENVBIN/samtools \
    --out-dir ~/s05/results/64_smoke_sr --internal-dir ~/s05/internal/64_smoke_sr 2>&1 | tail -5

# --- 7. Full run, all genes, depths 2,3 (background; resumable: rerun the same line after a stop) -
cd ~/repos/pilot-validation-main/scripts/hla_popgen
setsid nohup $PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest --workers 12 \
    --samtools $ENVBIN/samtools > ~/s05/logs/64.log 2>&1 < /dev/null & disown
tail -3 ~/s05/logs/64.log
# Add --sr to the same line for the short-read check (own --internal-dir if run separately).
# Depth 4 (non-coding) is a separate, optional run with its own paths:
#   ... --depths 4 --out-dir ~/s05/results/64_d4 --internal-dir ~/s05/internal/64_d4

# --- 8. After it finishes ------------------------------------------------------------------------
ls -la ~/s05/results/64
grep -E '^(protein|synonymous)\tclassical_8' ~/s05/results/64/funnel_with_evidence.tsv | column -t
# Pull ONLY ~/s05/results/64* (PUBLIC). Before committing: check every <20 is in place, and decide
# with Marc/Aleix whether now_known_alleles.tsv (named alleles, review_flag) is publishable
# (CLAUDE.md: small-cell disclosure is an open item in context/DECISIONS.md). Never pull ~/s05/internal/64*.
