#!/bin/bash
# One-shot setup for a BRAND-NEW VM: primes everything needed for TRUST4 repertoire
# calling + both embedding models (ESMC + catELMo), so the actual run commands (see
# ../README.md and the chat log) hit zero environment surprises.
#
# WHY THREE SEPARATE ENVIRONMENTS: TRUST4 lives in this folder's pixi env (python 3.12).
# ESMC needs a modern transformers/torch -- also fine in that same pixi env, or plain
# system python3. catELMo needs python 3.6 + allennlp 0.9.0 + torch 1.9.1 + an old
# tensorflow-gpu build (per its own embedders/README.md, verified live 2026-09-07 --
# NOT the TF2.6/Keras2.6 stack an earlier web search summary wrongly suggested; that one
# is for catELMo's separate downstream binding-affinity training pipeline, not embedding).
# Those two dependency sets cannot coexist in one environment -- hence a dedicated conda
# env for catELMo, isolated from everything else, matching this project's own established
# pattern (spechla/specimmune as separate pixi envs for the same reason -- see
# context/ENVIRONMENT.md quirk #17).
#
# Usage:  bash setup_repertoire_and_embedding_envs.sh [--skip-catelmo]
#   --skip-catelmo : set up TRUST4 + ESMC only (catELMo's env is the slowest, riskiest
#                    part -- skip it if time is tight and come back for it separately).
set -uo pipefail

SKIP_CATELMO=0
[[ "${1:-}" == "--skip-catelmo" ]] && SKIP_CATELMO=1

REPO=~/repos/pilot-validation
RNASEQ_DIR="$REPO/aleix/RNA-seq"
CATELMO_DIR=~/tools/catELMo
CATELMO_WEIGHTS_DIR=~/tools/catELMo_weights

echo "==== [1/4] TRUST4 pixi env ===="
cd "$RNASEQ_DIR"
pixi install
bash scripts/setup_trust4_refs.sh || echo "  (refs may already be fetched -- check reference/*.fa)"
echo "  OK -- verify with: pixi run run-trust4 --help"

echo ""
echo "==== [2/4] ESMC + SCEPTR deps (same pixi env, python 3.12) ===="
pixi run pip install --quiet sceptr transformers torch
pixi run python3 -c "import sceptr, transformers, torch; print('  OK -- sceptr', sceptr.__version__ if hasattr(sceptr,'__version__') else '(installed)', '| transformers', transformers.__version__, '| torch', torch.__version__)"

echo ""
echo "==== [3/4] disease-count script -- no install needed, just confirm it's current ===="
cd "$REPO" && git pull --quiet
echo "  deep_immune_breakdown.py has: HLA-linked, autoimmune, immunodeficiency, tumor "
echo "  (lymphoid+solid), plus this week's 7-condition section (Alzheimer's/Parkinson's/"
echo "  atherosclerosis-CAD/Long COVID/RLS/allergic rhinitis/atopic dermatitis)."
echo "  Ready to run once query_overlap_phenotypes.py's pheno cache exists (see README)."

if [[ "$SKIP_CATELMO" -eq 1 ]]; then
  echo ""
  echo "==== [4/4] catELMo -- SKIPPED (--skip-catelmo). Run this script again without the"
  echo "     flag when there's time for it, or run its section manually (see below)."
  echo ""
  echo "==== READY (TRUST4 + ESMC) -- catELMo pending ===="
  exit 0
fi

echo ""
echo "==== [4/4] catELMo -- separate conda env (python 3.6), per embedders/README.md ===="
mkdir -p ~/tools
if [[ ! -d "$CATELMO_DIR" ]]; then
  git clone --quiet https://github.com/Lee-CBG/catELMo.git "$CATELMO_DIR"
else
  echo "  $CATELMO_DIR already exists, skipping clone"
fi

if ! conda env list 2>/dev/null | grep -q "^catELMo "; then
  conda create -y -n catELMo python=3.6 >/dev/null
fi
# per embedders/README.md's exact recipe -- NOT the repo-root catELMo.yml (that's for the
# downstream binding-affinity trainer, different deps, different python version story)
conda run -n catELMo pip install --quiet \
  allennlp==0.9.0 torch==1.9.1 pandas==1.1.5 transformers==4.11.2 \
  tensorflow-gpu==1.14.0 overrides==3.1.0
echo "  catELMo env installed -- verify with: conda run -n catELMo python -c \"import allennlp; print('OK')\""

mkdir -p "$CATELMO_WEIGHTS_DIR"
if [[ ! -s "$CATELMO_WEIGHTS_DIR/weights.hdf5" ]]; then
  echo "  Downloading catELMo pretrained weights (4-layer BiLSTM, 1024-dim, ~few hundred MB)..."
  curl -sL "https://www.dropbox.com/sh/jpw6z71bsn1t7ev/AADRiL7_amT0vQrpep45PcOPa?dl=1" \
    -o "$CATELMO_WEIGHTS_DIR/catELMo.zip"
  (cd "$CATELMO_WEIGHTS_DIR" && unzip -q -o catELMo.zip && rm -f catELMo.zip)
  find "$CATELMO_WEIGHTS_DIR" -maxdepth 2 -iname "*.hdf5" -o -iname "*.json" 2>/dev/null
else
  echo "  weights already present at $CATELMO_WEIGHTS_DIR, skipping download"
fi

if [[ ! -s "$CATELMO_WEIGHTS_DIR/weights.hdf5" ]]; then
  echo "  !! weights.hdf5 not found after unzip -- Dropbox layout may have changed."
  echo "     Check $CATELMO_WEIGHTS_DIR by hand: ls -R $CATELMO_WEIGHTS_DIR"
  echo "     and point --weights-dir at wherever weights.hdf5 + options.json actually landed."
fi

echo ""
echo "==== READY -- TRUST4 + ESMC + catELMo all primed ===="
echo "Run order (see chat log / README for the full sequence):"
echo "  1. build_rnaseq_cohort.py + run_rnaseq_batch.sh          -- call repertoires"
echo "  2. pixi run python3 embed_cdr3s.py <cohort.tsv>          -- SCEPTR + ESMC"
echo "  3. conda run -n catELMo python3 embed_catelmo.py \\"
echo "       ~/pipeline_outputs/rnaseq/embeddings/cdr3_pool_used_for_embedding.tsv \\"
echo "       --weights-dir $CATELMO_WEIGHTS_DIR                  -- catELMo, third arm"
echo "  4. query_overlap_phenotypes.py + deep_immune_breakdown.py -- fresh disease counts"
