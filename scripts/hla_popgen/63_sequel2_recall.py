#!/usr/bin/env python3
"""63 -- HLA and KIR calls for the 991 `sequel2` people (S05 WS-K), into NEW output directories.

Why: the production HLA run (2026-08) and the KIR run (S03 43) were both launched with
`--single-phase --skip-trim-tier self_align_needed`, so the `sequel2` people (assembly FASTA present,
no assembly-to-hg38 PAF/BAM, trim_tier = self_align_needed) were never called. The Tier-3
self-align path that would have handled them failed its only real test (EXPERIMENTS.md
2026-08-05). Diagnosis: sprints/S05_qc_phasing_kir_trb/SEQUEL2_DIAGNOSIS.md. The fix lives in
scripts/run_immuannot_person.py ("Tier 3 v2"); this driver runs that script for the sequel2
people only.

Steps (each idempotent):
  --prepare-refs  builds, under <internal>/ref/ (public reference data, small):
                    hla_stage1.fa   chr6:28,500,001-34,500,000  (HLA window +/- 1 Mb)
                    kir_stage1.fa   chr19:53,600,001-55,920,000 (KIR window +/- 1 Mb)
                    hg38_primary.fa chr1-22, X, Y, M only (no ALT / HLA / decoy contigs, so a true
                                    MHC or KIR contig is not split between chr6 and its ALT copies)
                    hg38_primary.asm5.mmi  minimap2 -x asm5 index of hg38_primary.fa
  (default)       for each selected person, for each locus in --loci (hla, kir): runs
                  run_immuannot_person.py <pid> --enable-self-align-fallback with the refs above,
                  into <internal>/<locus>/pipeline_outputs/<pid>/immuannot_output/.
                  Done = hap1.gtf.gz and hap2.gtf.gz both present with >= 1 real call (the
                  production orchestrator's rule). Resume skips done people; a person-locus that
                  failed --max-attempts times is not retried unless --retry-failed.
  merge           (end of every run, or --merge-only) fragments -> <internal>/<locus>/
                  immuannot_calls.sequel2<sfx>.tsv; then the existing canonical table (READ ONLY)
                  + these rows -> <internal>/merged/immuannot_calls<sfx>.with_sequel2.tsv. The
                  existing ~/pipeline_outputs* trees are never written: every output path is
                  checked at start-up and the script refuses to run if one falls inside them.

Selection: --cohort (default ~/pipeline_outputs/immuannot_cohort_full.tsv, read only), rows with
trim_tier == self_align_needed (--tier). The platform mix of the selection is logged.

Outputs
  INTERNAL (--internal-dir, default ~/s05/internal/63; never pulled; per-person):
    <locus>/pipeline_outputs/<pid>/immuannot_output/{hap1,hap2}.gtf.gz, trimmed FASTAs, PAFs, logs
    <locus>/status.tsv          one row per attempt: person, locus, attempt, rc, seconds, done, reason
    <locus>/immuannot_calls.sequel2<sfx>.tsv, immuannot_timing.sequel2<sfx>.tsv
    merged/immuannot_calls.with_sequel2.tsv, merged/immuannot_calls.kir.with_sequel2.tsv
  PUBLIC (--out-dir, default ~/s05/results/63; aggregate only, counts 1-19 -> "<20"):
    run_summary.tsv             per locus: selected, done, failed, pending
    trim_summary.tsv            per locus x trim_method: haplotypes, median trimmed MB, median
                                stage-1 / stage-2 / immuannot seconds, median candidate contigs
    failure_reasons.tsv         per locus x reason class: people

Robustness: PID lock in <internal>; --limit; --workers (people in parallel); --threads per person
(default: CPUs // workers); --timeout per person-locus; --smoke (2 people, both loci, own
<internal>_smoke / <out>_smoke dirs, quirk #22).

VM (4 vCPU; env ~/repos/pilot-validation/.pixi/envs/specimmune, checkout
~/repos/pilot-validation-main; full list in sprints/S05_qc_phasing_kir_trb/vm_parts/63.sh):
    PY=~/repos/pilot-validation/.pixi/envs/specimmune/bin/python3
    cd ~/repos/pilot-validation-main
    $PY scripts/hla_popgen/63_sequel2_recall.py --prepare-refs
    $PY scripts/hla_popgen/63_sequel2_recall.py --smoke
    setsid nohup $PY -u scripts/hla_popgen/63_sequel2_recall.py --workers 1 \\
        > ~/s05/logs/63.log 2>&1 < /dev/null & disown

Tests: scripts/hla_popgen/tests/test_63_sequel2_recall.py (synthetic plain-gzip assemblies, mocked
minimap2 / immuannot.sh).
"""
import argparse
import concurrent.futures
import glob
import gzip
import os
import re
import subprocess
import sys
import threading
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
try:
    from _s05_common import mask_n  # noqa: E402
except Exception:  # pragma: no cover - flattened deploy
    def mask_n(n, mode="public"):
        n = int(n)
        return "<20" if mode == "public" and 0 < n < 20 else str(n)

HOME = os.path.expanduser("~")
PERSON_SCRIPT = os.path.join(REPO, "scripts", "run_immuannot_person.py")
DEFAULT_COHORT = os.path.join(HOME, "pipeline_outputs", "immuannot_cohort_full.tsv")
DEFAULT_INTERNAL = os.path.join(HOME, "s05", "internal", "63")
DEFAULT_OUT = os.path.join(HOME, "s05", "results", "63")
DEFAULT_MOUNT = os.path.join(HOME, "workspace", "vwb-aou-datasets-controlled-v9")
DEFAULT_HG38 = os.path.join(HOME, "ref", "Homo_sapiens_assembly38.fasta")

FLANK = 1_000_000
LOCI = {
    # region: the exact windows of the production HLA run and of the S03 KIR run (43).
    # sfx: the canonical-file suffix those runs used (immuannot_calls<sfx>.tsv).
    "hla": {"region": "chr6:29500000-33500000", "sfx": "",
            "existing": os.path.join(HOME, "pipeline_outputs", "immuannot_calls.tsv")},
    "kir": {"region": "chr19:54600000-54920000", "sfx": ".kir",
            "existing": os.path.join(HOME, "pipeline_outputs_kir", "immuannot_calls.kir.tsv")},
}
PRIMARY_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
_CALL_RE = re.compile(r'(consensus|allele) "[^"]+"')
_SKIP_RE = re.compile(r"(hap[12]|hap\?)?:? ?SKIP(?::| --) (.*)")


def log(msg):
    print(f"[63 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------------------------
# Paths and safety
# ---------------------------------------------------------------------------------------------
def protected_roots(home=HOME):
    """Every ~/pipeline_outputs* directory: the production trees this script must never write."""
    return [os.path.realpath(p) for p in glob.glob(os.path.join(home, "pipeline_outputs*"))] + [
        os.path.realpath(os.path.join(home, "pipeline_outputs"))]


def assert_safe(path, home=HOME):
    rp = os.path.realpath(path)
    for root in protected_roots(home) + [os.path.realpath(REPO)]:
        if rp == root or rp.startswith(root + os.sep):
            raise SystemExit(f"FATAL: refusing to write under {root} ({path}). 63 writes only to new "
                             f"directories (default ~/s05/internal/63, ~/s05/results/63).")
    base = os.path.basename(rp.rstrip(os.sep))
    if base.startswith("pipeline_outputs") and os.path.dirname(rp) == os.path.realpath(home):
        raise SystemExit(f"FATAL: refusing to write to {path} (a ~/pipeline_outputs* path).")
    return path


def locus_outroot(internal, locus):
    return os.path.join(internal, locus, "pipeline_outputs")


def ref_paths(internal):
    d = os.path.join(internal, "ref")
    return {"dir": d, "hla": os.path.join(d, "hla_stage1.fa"), "kir": os.path.join(d, "kir_stage1.fa"),
            "primary": os.path.join(d, "hg38_primary.fa"),
            "mmi": os.path.join(d, "hg38_primary.asm5.mmi")}


def slice_region(region, flank=FLANK):
    chrom, coords = region.split(":")
    s, e = (int(x) for x in coords.split("-"))
    return f"{chrom}:{max(1, s - flank + 1)}-{e + flank}"


# ---------------------------------------------------------------------------------------------
# Lock
# ---------------------------------------------------------------------------------------------
def acquire_lock(internal_dir):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, "63.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"FATAL: another 63_sequel2_recall.py (pid {pid}) holds {lock}")
        except (ProcessLookupError, ValueError):
            log(f"removing stale lock {lock}")
        except PermissionError:
            raise SystemExit(f"FATAL: {lock} is held by a live process owned by another user")
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


# ---------------------------------------------------------------------------------------------
# Reference preparation (once, single process)
# ---------------------------------------------------------------------------------------------
def _atomic_cmd_to_file(cmd, out_path):
    tmp = f"{out_path}.tmp.{os.getpid()}"
    with open(tmp, "w") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0 or os.path.getsize(tmp) == 0:
        os.remove(tmp)
        raise SystemExit(f"FATAL: {' '.join(cmd)} failed: {proc.stderr[-800:]}")
    os.replace(tmp, out_path)


def prepare_refs(internal, hg38, minimap2="minimap2", samtools="samtools"):
    r = ref_paths(internal)
    os.makedirs(r["dir"], exist_ok=True)
    if not os.path.exists(hg38):
        raise SystemExit(f"FATAL: {hg38} not found (ENVIRONMENT.md: ~/ref/Homo_sapiens_assembly38.fasta)")
    if not os.path.exists(hg38 + ".fai"):
        raise SystemExit(f"FATAL: {hg38}.fai not found; refusing to build an index next to the shared "
                         f"reference. Run `samtools faidx {hg38}` yourself if intended.")
    names = {line.split("\t")[0] for line in open(hg38 + ".fai")}
    for locus in ("hla", "kir"):
        if not (os.path.exists(r[locus]) and os.path.getsize(r[locus]) > 0):
            reg = slice_region(LOCI[locus]["region"])
            log(f"building {r[locus]} = {reg}")
            _atomic_cmd_to_file([samtools, "faidx", hg38, reg], r[locus])
    if not (os.path.exists(r["primary"]) and os.path.getsize(r["primary"]) > 0):
        chroms = [c for c in PRIMARY_CHROMS if c in names]
        missing = sorted(set(PRIMARY_CHROMS) - set(chroms))
        if len(chroms) < 24:
            raise SystemExit(f"FATAL: {hg38}.fai lacks primary chromosomes {missing}")
        log(f"building {r['primary']} ({len(chroms)} primary sequences; dropped "
            f"{len(names) - len(chroms)} ALT/HLA/decoy/unplaced)")
        _atomic_cmd_to_file([samtools, "faidx", hg38] + chroms, r["primary"])
    if not (os.path.exists(r["mmi"]) and os.path.getsize(r["mmi"]) > 0):
        log(f"building {r['mmi']} (minimap2 -x asm5 -d; several minutes, ~5-8 GB RAM)")
        tmp = f"{r['mmi']}.tmp.{os.getpid()}"
        proc = subprocess.run([minimap2, "-x", "asm5", "-d", tmp, r["primary"]],
                              stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0 or not os.path.exists(tmp):
            raise SystemExit(f"FATAL: minimap2 index build failed: {proc.stderr[-800:]}")
        os.replace(tmp, r["mmi"])
    log("refs ready: " + ", ".join(f"{k}={v}" for k, v in r.items() if k != "dir"))
    return r


# ---------------------------------------------------------------------------------------------
# Selection, done-check, status
# ---------------------------------------------------------------------------------------------
def select_people(cohort_path, tier="self_align_needed"):
    df = pd.read_csv(cohort_path, sep="\t", dtype=str)
    for c in ("person_id", "trim_tier"):
        if c not in df.columns:
            raise SystemExit(f"FATAL: {cohort_path} lacks column {c}")
    sel = df[df["trim_tier"] == tier].drop_duplicates("person_id")
    if "platform" in sel.columns:
        log(f"selected {len(sel)} people with trim_tier={tier}; platforms "
            f"{sel['platform'].fillna('NA').value_counts().to_dict()}")
    return list(sel["person_id"].astype(str))


def gtf_has_real_call(path):
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                if not line.startswith("#") and _CALL_RE.search(line):
                    return True
    except (OSError, EOFError):
        return False
    return False


def person_done(outroot, pid):
    d = os.path.join(outroot, str(pid), "immuannot_output")
    return all(gtf_has_real_call(os.path.join(d, f"{h}.gtf.gz")) for h in ("hap1", "hap2"))


STATUS_COLS = ["person_id", "locus", "attempt", "rc", "seconds", "done", "reason"]
_status_lock = threading.Lock()


def read_status(path):
    if not os.path.exists(path):
        return pd.DataFrame(columns=STATUS_COLS)
    return pd.read_csv(path, sep="\t", dtype=str)


def append_status(path, row):
    with _status_lock:
        new = not os.path.exists(path)
        with open(path, "a") as fh:
            if new:
                fh.write("\t".join(STATUS_COLS) + "\n")
            fh.write("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in STATUS_COLS) + "\n")
            fh.flush()
            os.fsync(fh.fileno())


def classify_reason(text):
    """Short, ID-free reason class from the per-person log (used in the PUBLIC failure table)."""
    reasons = [m.group(2) for m in _SKIP_RE.finditer(text or "")]
    if not reasons:
        return "no_skip_line"
    r = reasons[0]
    for key, label in [("needs --self-align-confirm-ref", "no_confirm_ref"),
                       ("not found", "ref_or_file_missing"),
                       ("stage-1 self-align failed", "stage1_minimap2_failed"),
                       ("no contig has a primary block", "stage1_no_candidates"),
                       ("candidates total", "stage1_whole_assembly_pullin"),
                       ("could not extract candidate", "candidate_extract_failed"),
                       ("stage-2 whole-genome alignment failed", "stage2_minimap2_failed"),
                       ("no candidate contig has a primary", "stage2_no_confirmed_contig"),
                       ("confirmed sub-range totals", "size_guard"),
                       ("not chr", "stage1_ref_wrong_chrom"),
                       ("holds", "stage1_ref_wrong_chrom"),
                       ("samtools faidx trim failed", "trim_failed"),
                       ("missing assembly FASTA", "no_assembly_fasta")]:
        if key in r:
            return label
    return "other"


# ---------------------------------------------------------------------------------------------
# One person, one locus
# ---------------------------------------------------------------------------------------------
def person_cmd(pid, locus, args, refs):
    L = LOCI[locus]
    return [args.python, args.person_script, str(pid),
            "--mount", args.mount, "--outroot", locus_outroot(args.internal_dir, locus),
            "--immuannot-dir", args.immuannot_dir, "--refdir", args.refdir,
            "--threads", str(args.threads), "--region", L["region"], "--pad", str(args.pad),
            "--out-suffix", f"{L['sfx']}.{pid}", "--force",
            "--enable-self-align-fallback",
            "--self-align-ref", refs[locus], "--self-align-confirm-ref", refs["mmi"],
            "--hg38-ref", args.hg38_ref, "--minimap2", args.minimap2]


def child_env(python):
    """Puts the env's bin dir (minimap2, samtools, immuannot deps) first on PATH, which is what
    `pixi run -e specimmune` would do, so the env's python can be called directly."""
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(os.path.abspath(python)) + os.pathsep + env.get("PATH", "")
    env["PYTHONUNBUFFERED"] = "1"
    return env


def run_person_locus(pid, locus, args, refs, attempt, runner=subprocess.run):
    outroot = locus_outroot(args.internal_dir, locus)
    pdir = os.path.join(outroot, str(pid), "immuannot_output")
    os.makedirs(pdir, exist_ok=True)
    cmd = person_cmd(pid, locus, args, refs)
    t0 = time.time()
    try:
        proc = runner(cmd, capture_output=True, text=True, timeout=args.timeout,
                      env=child_env(args.python))
        rc, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as e:
        rc, out, err = "timeout", str(e.stdout or ""), str(e.stderr or "") + "\nhap?: SKIP -- timeout"
    secs = round(time.time() - t0, 1)
    with open(os.path.join(pdir, ".63_last_run.log"), "w") as fh:
        fh.write(f"# cmd: {' '.join(cmd)}\n# rc: {rc}\n# seconds: {secs}\n--- stdout ---\n{out}\n"
                 f"--- stderr ---\n{err}\n")
    done = person_done(outroot, pid)
    reason = "" if done else ("timeout" if rc == "timeout" else classify_reason(err))
    append_status(os.path.join(args.internal_dir, locus, "status.tsv"),
                  {"person_id": pid, "locus": locus, "attempt": attempt, "rc": rc,
                   "seconds": secs, "done": int(done), "reason": reason})
    return done, reason, secs


# ---------------------------------------------------------------------------------------------
# Merge + public summaries
# ---------------------------------------------------------------------------------------------
def _read_tsvs(paths):
    frames = []
    for p in paths:
        try:
            frames.append(pd.read_csv(p, sep="\t", dtype=str))
        except pd.errors.EmptyDataError:
            continue
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _fragments(outroot, base, sfx):
    pat = re.compile(rf"^{re.escape(base)}{re.escape(sfx)}\.[^.]+\.tsv$")
    return sorted(p for p in glob.glob(os.path.join(outroot, f"{base}{sfx}.*.tsv"))
                  if pat.match(os.path.basename(p)))


def _write_atomic(df, path):
    assert_safe(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def merge_locus(internal, locus, existing_path=None):
    """Fragments -> <internal>/<locus>/immuannot_{calls,timing}.sequel2<sfx>.tsv (fragments kept);
    existing canonical (read only) + sequel2 calls -> <internal>/merged/...with_sequel2.tsv.
    Returns a dict of counts for the log."""
    L = LOCI[locus]
    outroot = locus_outroot(internal, locus)
    res = {"locus": locus}
    for base, key in (("immuannot_calls", ["person_id", "gene"]), ("immuannot_timing", ["person_id", "hap"])):
        df = _read_tsvs(_fragments(outroot, base, L["sfx"]))
        if df.empty:
            res[base] = 0
            continue
        df = df.drop_duplicates(subset=key, keep="last")
        _write_atomic(df, os.path.join(internal, locus, f"{base}.sequel2{L['sfx']}.tsv"))
        res[base] = len(df)
    calls_path = os.path.join(internal, locus, f"immuannot_calls.sequel2{L['sfx']}.tsv")
    existing_path = existing_path if existing_path is not None else L["existing"]
    if os.path.exists(calls_path) and existing_path and os.path.exists(existing_path):
        new = pd.read_csv(calls_path, sep="\t", dtype=str)
        old = pd.read_csv(existing_path, sep="\t", dtype=str)
        overlap = set(old["person_id"]) & set(new["person_id"])
        if overlap:
            log(f"WARNING [{locus}]: {len(overlap)} sequel2 people already in {existing_path}; their "
                f"existing rows are kept and the new rows dropped from the combined table.")
            new = new[~new["person_id"].isin(overlap)]
        both = pd.concat([old, new], ignore_index=True)
        out = os.path.join(internal, "merged", f"immuannot_calls{L['sfx']}.with_sequel2.tsv")
        _write_atomic(both, out)
        res.update(existing_people=old["person_id"].nunique(), added_people=new["person_id"].nunique(),
                   combined_people=both["person_id"].nunique(), overlap=len(overlap), combined_path=out)
    return res


def _median(s):
    if s is None:
        return ""
    s = pd.to_numeric(s, errors="coerce").dropna()
    return f"{s.median():.2f}" if len(s) else ""


def write_public(internal, out_dir, people, loci):
    assert_safe(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    run_rows, trim_rows, fail_rows = [], [], []
    for locus in loci:
        outroot = locus_outroot(internal, locus)
        st = read_status(os.path.join(internal, locus, "status.tsv"))
        done = {p for p in people if person_done(outroot, p)}
        attempted = set(st["person_id"]) & set(people) if len(st) else set()
        failed = attempted - done
        run_rows.append({"locus": locus, "region": LOCI[locus]["region"],
                         "people_selected": mask_n(len(people)), "people_done": mask_n(len(done)),
                         "people_failed_so_far": mask_n(len(failed)),
                         "people_not_attempted": mask_n(len(set(people) - attempted - done))})
        if len(st):
            last = st[st["person_id"].isin(failed)].groupby("person_id").tail(1)
            for reason, n in last["reason"].fillna("").value_counts().items():
                fail_rows.append({"locus": locus, "reason": reason or "unknown", "people": mask_n(n)})
        tim_path = os.path.join(internal, locus, f"immuannot_timing.sequel2{LOCI[locus]['sfx']}.tsv")
        if os.path.exists(tim_path):
            t = pd.read_csv(tim_path, sep="\t", dtype=str)
            t["trim_method"] = t.get("trim_method", pd.Series(dtype=str)).fillna("none")
            for m, g in t.groupby("trim_method"):
                trim_rows.append({"locus": locus, "trim_method": m, "haplotypes": mask_n(len(g)),
                                  "median_trimmed_mb": _median(g.get("trimmed_mb")),
                                  "median_candidates": _median(g.get("self_align_candidates")),
                                  "median_stage1_s": _median(g.get("self_align_seconds")),
                                  "median_stage2_s": _median(g.get("confirm_seconds")),
                                  "median_immuannot_s": _median(g.get("immuannot_seconds")),
                                  "median_hap_total_s": _median(g.get("hap_total_seconds"))})
    pd.DataFrame(run_rows).to_csv(os.path.join(out_dir, "run_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(trim_rows, columns=["locus", "trim_method", "haplotypes", "median_trimmed_mb",
                                     "median_candidates", "median_stage1_s", "median_stage2_s",
                                     "median_immuannot_s", "median_hap_total_s"]).to_csv(
        os.path.join(out_dir, "trim_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(fail_rows, columns=["locus", "reason", "people"]).to_csv(
        os.path.join(out_dir, "failure_reasons.tsv"), sep="\t", index=False)


# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default=DEFAULT_COHORT)
    ap.add_argument("--tier", default="self_align_needed")
    ap.add_argument("--people", help="Optional file of person_ids (one per line) to restrict to.")
    ap.add_argument("--loci", default="hla,kir")
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL)
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    ap.add_argument("--mount", default=DEFAULT_MOUNT)
    ap.add_argument("--hg38-ref", default=DEFAULT_HG38)
    ap.add_argument("--immuannot-dir", default=os.path.join(HOME, "tools", "Immuannot"))
    ap.add_argument("--refdir", default=os.path.join(HOME, "tools", "Immuannot_refdata"))
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--person-script", default=PERSON_SCRIPT)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--pad", type=int, default=100_000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--threads", type=int, default=None, help="per person (default CPUs // workers)")
    ap.add_argument("--timeout", type=int, default=3 * 3600, help="seconds per person-locus")
    ap.add_argument("--max-attempts", type=int, default=2)
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--smoke", action="store_true",
                    help="2 people (or --limit), both loci, <internal>_smoke and <out>_smoke dirs")
    ap.add_argument("--prepare-refs", action="store_true")
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--no-merge", action="store_true")
    return ap


def main(argv=None, runner=subprocess.run):
    args = build_parser().parse_args(argv)
    if args.smoke:
        args.internal_dir = args.internal_dir.rstrip("/") + "_smoke"
        args.out_dir = args.out_dir.rstrip("/") + "_smoke"
        args.limit = args.limit or 2
    for p in (args.internal_dir, args.out_dir):
        assert_safe(p)
    loci = [x.strip() for x in args.loci.split(",") if x.strip()]
    for x in loci:
        if x not in LOCI:
            raise SystemExit(f"FATAL: unknown locus {x}; choose from {sorted(LOCI)}")
    args.threads = args.threads or max(1, (os.cpu_count() or 4) // max(1, args.workers))

    lock = acquire_lock(args.internal_dir)
    try:
        refs = ref_paths(args.internal_dir)
        if args.smoke:  # reuse the full run's refs if they exist (public data, no need to rebuild)
            full = ref_paths(args.internal_dir[:-len("_smoke")])
            if all(os.path.exists(full[k]) for k in ("hla", "kir", "mmi")):
                refs = full
        if args.prepare_refs:
            prepare_refs(args.internal_dir, args.hg38_ref, args.minimap2, args.samtools)
            return 0
        people = select_people(args.cohort, args.tier)
        if args.people:
            keep = {line.strip() for line in open(args.people) if line.strip()}
            people = [p for p in people if p in keep]
        if args.limit:
            people = people[:args.limit]
        if args.merge_only:
            for locus in loci:
                log(f"merge {merge_locus(args.internal_dir, locus)}")
            write_public(args.internal_dir, args.out_dir, people, loci)
            return 0
        missing = [k for k in ("hla", "kir", "mmi") if not os.path.exists(refs[k])]
        if missing:
            raise SystemExit(f"FATAL: refs {missing} not built; run with --prepare-refs first")
        if not os.path.isfile(os.path.join(args.mount, "v9", "wgs", "long_read", "manifest.tsv")):
            raise SystemExit(f"FATAL: {args.mount}/v9/wgs/long_read/manifest.tsv not found -- "
                             f"controlled-data automount down? (ENVIRONMENT.md quirk #35)")

        todo = []
        for locus in loci:
            outroot = locus_outroot(args.internal_dir, locus)
            st = read_status(os.path.join(args.internal_dir, locus, "status.tsv"))
            tries = st.groupby("person_id").size().to_dict() if len(st) else {}
            n_done = n_given_up = 0
            for p in people:
                if person_done(outroot, p):
                    n_done += 1
                elif tries.get(p, 0) >= args.max_attempts and not args.retry_failed:
                    n_given_up += 1
                else:
                    todo.append((p, locus, tries.get(p, 0) + 1))
            log(f"[{locus}] {len(people)} selected, {n_done} done (resume), {n_given_up} given up, "
                f"{sum(1 for t in todo if t[1] == locus)} to run")
        # person-major order: HLA then KIR for the same person, so a partial run covers whole people
        order = {p: i for i, p in enumerate(people)}
        todo.sort(key=lambda t: (order[t[0]], loci.index(t[1])))
        log(f"{len(todo)} person-locus jobs, {args.workers} worker(s) x {args.threads} thread(s)")

        t0, n_ok, n_fail = time.time(), 0, 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(run_person_locus, p, l, args, refs, a, runner): (p, l)
                    for p, l, a in todo}
            for i, fut in enumerate(concurrent.futures.as_completed(futs), 1):
                done, reason, secs = fut.result()
                n_ok += int(done)
                n_fail += int(not done)
                el = time.time() - t0
                log(f"{i}/{len(todo)} jobs: ok {n_ok}, failed {n_fail} (last {futs[fut][1]} "
                    f"{secs:.0f}s{', ' + reason if reason else ''}); {el/3600:.2f} h elapsed, "
                    f"ETA {el / i * (len(todo) - i) / 3600:.1f} h")
                if not args.no_merge and i % 50 == 0:
                    for locus in loci:
                        merge_locus(args.internal_dir, locus)
                    write_public(args.internal_dir, args.out_dir, people, loci)
        if not args.no_merge:
            for locus in loci:
                log(f"merge {merge_locus(args.internal_dir, locus)}")
        write_public(args.internal_dir, args.out_dir, people, loci)
        log(f"done: {n_ok} ok, {n_fail} failed this invocation, {(time.time() - t0)/3600:.2f} h")
        return 0
    finally:
        try:
            os.remove(lock)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
