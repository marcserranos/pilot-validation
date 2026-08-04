#!/usr/bin/env bash
# One-shot bootstrap for a FRESH Workbench VM, from nothing to ready-to-run-Immuannot.
# Purpose (2026-08-02): before spending real money on the ~13,000-14,500 person production run on
# a large machine, rehearse the exact "spin up a box and get it working" flow on a cheaper test VM
# first -- this script IS that rehearsal, and doubles as the production VM's own setup script later.
# Idempotent -- safe to re-run if it fails partway through or you resize/restart the VM.
#
# Only sets up what Immuannot itself needs (the `specimmune` pixi env has minimap2/samtools/python3,
# same env SpecImmune already uses -- see pixi.toml). Does NOT build the `spechla` env or SpecHLA --
# not needed for this test.
#
# Usage (paste into a fresh Jupyter terminal on the new VM -- nothing needs to pre-exist):
#   curl -fsSL https://raw.githubusercontent.com/marcserranos/pilot-validation/main/scripts/bootstrap_vm.sh | bash
# or, if you'd rather clone first and inspect before running:
#   git clone https://github.com/marcserranos/pilot-validation.git ~/repos/pilot-validation
#   cd ~/repos/pilot-validation && bash scripts/bootstrap_vm.sh
#
# After this completes, run ONE smoke-test person before trusting the machine for anything bigger:
#   cd ~/repos/pilot-validation && pixi shell -e specimmune
#   python3 scripts/run_immuannot_person.py <a_known_good_person_id> --threads 4 --force
set -uo pipefail

REPO_URL="https://github.com/marcserranos/pilot-validation.git"
REPO_DIR="$HOME/repos/pilot-validation"
# NEVER hardcode this -- this project now spans multiple Workbench workspaces (personal +
# institutional/Stanford-pod), each its own GCP project with its own billing-project id
# (confirmed 2026-08-02: hardcoding the personal workspace's id here would have silently pointed
# gcsfuse at the wrong billing project on any other workspace's VM). Always derive it live.
BILLING_PROJECT="$(gcloud config get-value project 2>/dev/null)"
BUCKET="vwb-aou-datasets-controlled"
MOUNT_DIR="$HOME/mnt/aou-controlled"

fail=0
step() { echo ""; echo "== $1 ==" >&2; }
ok()   { echo "  OK: $1" >&2; }
bad()  { echo "  FAILED: $1" >&2; fail=1; }

# --- 0. Sanity: confirm we're actually on the VM, not a Mac (same class of mistake ENVIRONMENT.md
#     flags as the single most expensive lesson of this whole project -- check before anything else). ---
step "0. Confirm this is the VM, not a local machine"
if [[ "$(uname)" != "Linux" ]]; then
  bad "uname reports '$(uname)', not Linux -- this does NOT look like the Workbench VM. STOPPING before doing anything. If you meant to run this locally, you didn't."
  exit 1
fi
ok "Linux confirmed. Hostname: $(hostname). Whoami: $(whoami)."
if [ -z "$BILLING_PROJECT" ]; then
  bad "gcloud config get-value project returned empty -- can't determine which workspace/billing project this VM belongs to. STOPPING (mounting against the wrong/no billing project is exactly the VPC-SC mess from earlier this session)."
  exit 1
fi
ok "Billing project (live-detected, not hardcoded): $BILLING_PROJECT"

# --- 1. pixi ---
step "1. pixi"
if command -v pixi >/dev/null 2>&1; then
  ok "pixi already on PATH ($(pixi --version))."
else
  echo "  Installing pixi..." >&2
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="$HOME/.pixi/bin:$PATH"
  if command -v pixi >/dev/null 2>&1; then
    ok "pixi installed ($(pixi --version)). Add \$HOME/.pixi/bin to PATH in future shells if it's not already in your shell rc."
  else
    bad "pixi install ran but 'pixi' still not found on PATH."
  fi
fi

# --- 2. Clone this repo ---
step "2. Clone pilot-validation repo"
if [ -d "$REPO_DIR/.git" ]; then
  ok "Already cloned at $REPO_DIR -- skipping clone (not pulling automatically; run 'git pull' yourself if you want the latest)."
else
  mkdir -p "$(dirname "$REPO_DIR")"
  if git clone "$REPO_URL" "$REPO_DIR"; then
    ok "Cloned to $REPO_DIR."
  else
    bad "git clone failed."
  fi
fi

# --- 3. pixi env (specimmune only -- Immuannot's actual runtime deps) ---
step "3. pixi install -e specimmune"
if [ -d "$REPO_DIR" ]; then
  ( cd "$REPO_DIR" && pixi install -e specimmune )
  if ( cd "$REPO_DIR" && pixi run -e specimmune -- which minimap2 samtools python3 ) >/dev/null 2>&1; then
    ok "minimap2/samtools/python3 all resolve inside the specimmune env."
  else
    bad "one or more of minimap2/samtools/python3 not found inside the specimmune env after install."
  fi
else
  bad "repo dir missing, cannot pixi install."
fi

# --- 4. Immuannot itself + reference data ---
step "4. Immuannot install (scripts/setup_immuannot.sh)"
if [ -d "$REPO_DIR" ]; then
  ( cd "$REPO_DIR" && bash scripts/setup_immuannot.sh )
  if [ -n "$(find "$HOME/tools/Immuannot" -maxdepth 3 -name immuannot.sh 2>/dev/null)" ]; then
    ok "immuannot.sh found."
  else
    bad "immuannot.sh not found after setup_immuannot.sh ran."
  fi
else
  bad "repo dir missing, cannot run setup_immuannot.sh."
fi

# --- 5. gcsfuse mount (ENVIRONMENT.md quirk #11 -- required flags, --implicit-dirs is NOT optional) ---
step "5. gcsfuse mount"
mkdir -p "$MOUNT_DIR"
if mountpoint -q "$MOUNT_DIR" 2>/dev/null; then
  ok "$MOUNT_DIR already mounted."
else
  gcsfuse --billing-project "$BILLING_PROJECT" --implicit-dirs "$BUCKET" "$MOUNT_DIR"
  sleep 2  # ENVIRONMENT.md quirk #14: give the FUSE layer a beat before trusting it in the same paste
fi
if [ -f "$MOUNT_DIR/v9/wgs/long_read/manifest.tsv" ]; then
  ok "mount verified -- v9 lrWGS manifest is readable."
else
  bad "v9/wgs/long_read/manifest.tsv not found under the mount -- gcsfuse may not actually be up."
fi

# --- 6. hg38 reference (only needed by the Tier 3 self-align fallback) ---
# Added 2026-08-05: run_immuannot_person.py's Tier 3 (--enable-self-align-fallback, used for the
# ~991 sequel2 people who have assembly FASTAs but no aln-to-hg38 BAM/PAF) carves a chr6-only
# slice out of this FASTA to self-align against. Without it, the production orchestrator's phase 2
# skips every one of those people -- quietly, one at a time, ~2 days into an unattended run.
# Public Broad bucket, so plain gs:// works with no requester-pays flags (ENVIRONMENT.md quirk #11).
# Skipped automatically if a chr6 slice was already carved by a previous run.
step "6. hg38 reference (for the Tier 3 self-align fallback)"
mkdir -p "$HOME/ref"
if [ -s "$HOME/ref/chr6.fasta" ]; then
  ok "chr6 slice already carved at ~/ref/chr6.fasta -- full hg38 FASTA not needed."
elif [ -s "$HOME/ref/Homo_sapiens_assembly38.fasta" ]; then
  ok "hg38 reference already present at ~/ref/Homo_sapiens_assembly38.fasta."
else
  echo "  Downloading hg38 reference (~3GB, public bucket)..." >&2
  if gcloud storage cp \
      gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
      gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
      "$HOME/ref/" 2>/dev/null; then
    ok "hg38 reference downloaded to ~/ref/."
  else
    # Deliberately NOT a hard failure: phase 1 (~12,261 people, the bulk of the run) needs none of
    # this, and the orchestrator re-checks and warns at launch anyway.
    echo "  WARNING: hg38 reference download failed. Phase 1 is unaffected; fix before phase 2." >&2
  fi
fi

# --- Summary ---
step "Summary"
if [ "$fail" -eq 0 ]; then
  echo "All checks passed. This VM is ready. Next: run ONE smoke-test person before anything bigger:" >&2
  echo "  cd $REPO_DIR && pixi shell -e specimmune" >&2
  echo "  python3 scripts/run_immuannot_person.py <a_known_good_person_id> --threads 4 --force" >&2
  exit 0
else
  echo "One or more steps FAILED (see 'FAILED:' lines above). Fix those before running anything at scale." >&2
  exit 1
fi
