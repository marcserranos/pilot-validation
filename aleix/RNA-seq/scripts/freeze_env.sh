#!/bin/bash
# Record the exact software environment that produced the results, for the methods section.
# Prints to stdout; run_polish.sh writes it to reports/08_methods_audit/environment.txt.
# Nothing participant-level is printed: BAM headers are reduced to the @PG program name and
# version fields (no CL command lines, no @RG sample tags, no paths).
#
#   bash scripts/freeze_env.sh > reports/08_methods_audit/environment.txt
set +e
cd "$(dirname "$0")/.."
COHORT=~/pipeline_outputs/rnaseq/cohort_full.tsv
BUCKET=gs://vwb-aou-datasets-controlled

echo "# Environment freeze -- $(date -u +%Y-%m-%dT%H:%MZ)"
echo "repo_commit: $(git rev-parse HEAD)"
echo "machine_type: $(curl -s -m 3 -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/machine-type | awk -F/ '{print $NF}')"
echo "vcpus: $(nproc)   mem: $(free -g | awk '/Mem:/{print $2" GB"}')"
echo "os: $(. /etc/os-release && echo "$PRETTY_NAME")"
echo "pixi: $(pixi --version 2>/dev/null)"
echo "gcloud: $(gcloud --version 2>/dev/null | head -1)"

echo; echo "## pixi environment (aleix/RNA-seq/pixi.toml)"
pixi list 2>/dev/null | grep -iE '^(Package|python |trust4|samtools|pandas|scipy|numpy|perl) '

echo; echo "## pip packages in the pixi environment"
pixi run pip freeze 2>/dev/null | grep -iE '^(sceptr|torch|olga|umap-learn|pynndescent|scikit-learn|cnsplots|numpy|scipy|pandas|matplotlib|esm|transformers|accelerate|numba|statsmodels)(==| @)'

echo; echo "## TRUST4"
echo "run-trust4 usage header:"; pixi run run-trust4 2>&1 | grep -iE 'v[0-9]+\.[0-9]+|version' | head -3
echo "trust4 binary:"; pixi run trust4 2>&1 | grep -iE 'v[0-9]+\.[0-9]+|version' | head -3
if [[ -d ~/tools/TRUST4/.git ]]; then
  echo "reference source clone (~/tools/TRUST4): $(git -C ~/tools/TRUST4 log -1 --format='%H  %cd' --date=short)"
else
  echo "reference source clone: ~/tools/TRUST4 not present on this VM"
fi
for f in reference/hg38_bcrtcr.fa "reference/human_IMGT+C.fa"; do
  [[ -s "$f" ]] && echo "md5 $(md5sum "$f" | cut -d' ' -f1)  $(basename "$f")  $(wc -c < "$f") bytes"
done

echo; echo "## Aligner recorded in BAM headers (3 BAMs, one each from the first 3 ancestry groups)"
if [[ -s "$COHORT" ]]; then
  awk -F'\t' 'NR>1 && !seen[$2]++ {print $3}' "$COHORT" | head -3 | while read -r rel; do
    gcloud storage cat --billing-project "$GOOGLE_PROJECT" --range=0-8000000 "$BUCKET/$rel" 2>/dev/null \
      | pixi run samtools view -H - 2>/dev/null | grep '^@PG' \
      | awk -F'\t' '{o=""; for(i=1;i<=NF;i++) if($i ~ /^(ID|PN|VN):/) o=o" "$i; print "  @PG"o}' | sort -u
    echo "  --"
  done
else
  echo "  cohort file missing; skipped"
fi
