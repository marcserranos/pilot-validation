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
BILLING_PROJECT="wb-glacial-potato-8710"
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
