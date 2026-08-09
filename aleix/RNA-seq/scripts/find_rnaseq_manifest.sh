#!/bin/bash
# Discovery helper -- the RNA-seq manifest (research_id -> RNA BAM path) has not been
# located yet. Every other AoU data type in this project resolves through a
# v9/...manifest.* file (see ../../../context/ENVIRONMENT.md for the srWGS/lrWGS ones) --
# this script browses the bucket the same way those were found. It does NOT assume a path.
#
# [MED] billing project below is this workspace's, reused from the HLA workstream
# (../../../context/ENVIRONMENT.md) -- confirm it's still correct if this VM is in a
# different workspace.
#
# Usage:  bash find_rnaseq_manifest.sh
set -uo pipefail

BILLING_PROJECT="wb-glacial-potato-8710"
BUCKET="vwb-aou-datasets-controlled"

echo "==== top level of v9/ ===="
gcloud storage ls -u "$BILLING_PROJECT" "gs://$BUCKET/v9/" 2>&1 || \
  gsutil -u "$BILLING_PROJECT" ls "gs://$BUCKET/v9/"

echo ""
echo "==== anywhere under v9/ with 'rna' in the path (case-insensitive) ===="
gcloud storage ls -u "$BILLING_PROJECT" -r "gs://$BUCKET/v9/**" 2>&1 | grep -i rna || \
  gsutil -u "$BILLING_PROJECT" ls -r "gs://$BUCKET/v9/**" 2>&1 | grep -i rna

echo ""
echo "Once a manifest.csv/tsv shows up above, inspect its header with:"
echo "  gcloud storage cat -u $BILLING_PROJECT <path> | head -3"
echo "and record the confirmed path + schema in ../README.md's status checklist."
