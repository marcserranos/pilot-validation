#!/usr/bin/env python3
"""WS5 (Sprint S03, call #8 summary section 7) -- KIR scoping + pilot.

Cole (call #8, 2026-09-22): the production Immuannot run only ever trimmed assemblies to the
chr6 HLA window (`run_immuannot_person.py`'s `DEFAULT_REGION`) -- KIR genes (chr19, the leukocyte
receptor complex / LRC) were **never called**, even though Immuannot's own reference data set
(IPD-KIR) and its own detection logic support them natively (reference/README_Immuannot.md's
"Gene coverage" section lists all 17 KIR genes; no extra flag selects the family -- Immuannot just
annotates whatever gene families its refdata alleles map onto, HLA or KIR, wherever the trimmed
contig happens to align). So the only change needed to call KIR is trimming to the RIGHT region --
`run_immuannot_person.py --region` already exists and is fully backward-compatible (default
unchanged, chr6:29,500,000-33,500,000) -- this script does not modify that file at all, only
supplies a different --region value and a person-selection wrapper.

KIR window, GRCh38 (justification -- see README.md's "Coordinates" section for the primary-source
citations, since a wrong region here silently ships an empty KIR panel, not a loud error):
  KIR3DL3 (centromeric-most framework gene): chr19:54,724,442-54,736,632
  KIR3DL2 (telomeric-most framework gene):   chr19:54,850,443-54,867,207
  -> core KIR tandem array spans ~143 kb, consistent with the literature's "~150 kb" figure.
  Anchors used to set the outer trim boundary, since the KIR array itself is NOT a safe boundary
  (large structural CNVs -- insertion/deletion haplotypes reported up to and beyond the array's own
  span -- can push a real gene's alignment past the array's textbook edges):
    LILRB1 (immediately centromeric of the array): chr19:54,616,309-54,638,022
    FCAR    (immediately telomeric of the array):   chr19:54,874,235-54,891,420
  KIR_REGION_DEFAULT below = chr19:54,600,000-54,920,000 (320 kb) -- ~110 kb margin on the
  centromeric side (past LILRB1's start) and ~30 kb margin on the telomeric side (past FCAR's end).
  This region is only the OUTER candidate-contig-selection boundary passed to
  run_immuannot_person.py's existing --region/--pad machinery: the .paf-based Tier-1 trim
  (regions_from_paf) then adds its own further +/-100 kb pad around whatever actually aligns here
  (DEFAULT_PAD, unchanged, already proven safe for the analogous HLA case) -- so the real
  structural margin is wider than 320 kb once a person's own contig alignment is factored in.

Tier 3 (self-align) chr6-only-reference bug (run_immuannot_person.py's `ensure_chr6_ref()` calls
`samtools faidx <hg38 ref> chr6` with the chromosome name LITERALLY hardcoded -- not parameterized
by the target --region's chromosome). Not fixed in that file (out of this script's backward-compat
scope). Worked around here WITHOUT touching run_immuannot_person.py: `ensure_chr6_ref()` is a
no-op the moment its cache path already exists and is non-empty, so `ensure_kir_ref_cache()` below
pre-extracts chr19 (not chr6) into a *separate* cache file and this script passes that path via
run_immuannot_person.py's own `--chr6-ref-cache` flag -- the misleading flag NAME still means "the
Tier-3 self-align reference cache path", it is not chr6-specific in what it accepts, only in its
default. Confirmed by reading ensure_chr6_ref()'s body (scripts/run_immuannot_person.py) before
relying on this, not assumed.

Usage (on the VM, from ~/s03, --threads capped at 2 per CLAUDE.md's pilot-only compute budget):
  PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 41_kir_pilot.py --n-people 20 --threads 2 \\
      --immuannot-script ~/repos/pilot-validation/scripts/run_immuannot_person.py \\
      --outroot ~/s03/results/41/pipeline_outputs --out-suffix .kir41

Disclosure: person_ids are selected and used ONLY inside this VM process and its own log/tsv
files under --outroot (never printed to a channel meant to be read back through the JupyterLab
REST/websocket bridge -- VM_CHANNEL.md's hard rule). Progress lines printed to stdout are
aggregate-only (counts, indices, timings) by construction -- grep this file for `print(` to audit.
"""
import argparse
import gzip
import os
import random
import re
import subprocess
import sys
import time
from collections import defaultdict

import pandas as pd

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
    "samples_relatedness.tsv")
DEFAULT_IMMUANNOT_SCRIPT = os.path.expanduser("~/repos/pilot-validation/scripts/run_immuannot_person.py")
DEFAULT_HG38_REF = os.path.expanduser("~/ref/Homo_sapiens_assembly38.fasta")
DEFAULT_KIR_REF_CACHE = os.path.expanduser("~/ref/chr19_kir.fasta")

KIR_REGION_DEFAULT = "chr19:54600000-54920000"
KIR_PAD_DEFAULT = 100_000

# Same kinship threshold as scripts/hla_popgen/24_novelty_by_field.py / 29_hla_ld_by_ancestry.py
# (KIN_MIN = 0.0442, third-degree-or-closer per KING) -- copied, not re-derived, so "unrelated"
# means exactly what it means in the rest of this sub-project (VM_CHANNEL.md's instruction).
KIN_MIN = 0.0442
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ASSEMBLY_PLATFORMS = {"revio", "sequel2e", "sequel2"}  # the only platforms with hap-assembly FASTA

# 17 KIR genes Immuannot's refdata covers (reference/README_Immuannot.md "Gene coverage").
KIR_GENES = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DP1",
             "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2",
             "KIR3DL3", "KIR3DP1", "KIR3DS1"]
FRAMEWORK_GENES = ["KIR3DL3", "KIR2DL4", "KIR3DP1", "KIR3DL2"]  # present on ~all haplotypes
# "B-content" genes: any one of these on a haplotype marks it cB (variable/activating-rich)
# rather than cA (fixed, mostly-inhibitory) content, standard KIR nomenclature (Uhrberg et al 1997,
# Hsu et al 2002 haplotype-A/B framework already used informally in SCHEMA.md's kir gene_class row).
B_CONTENT_GENES = {"KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5",
                    "KIR3DS1"}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Step 1: pilot cohort selection (unrelated + ancestry/platform-stratified).
# Kinship logic copied verbatim in spirit from scripts/hla_popgen/24_novelty_by_field.py's
# greedy_unrelated()/load_relatedness_pairs() so "unrelated" carries the same meaning across the
# whole hla_popgen sub-project, without importing that module's heavier plotting dependencies
# (matplotlib/scipy/umap) into a scoping/timing script that doesn't need them.
# ---------------------------------------------------------------------------
def load_relatedness_pairs(path):
    if not os.path.exists(path):
        log(f"WARNING: relatedness table not found at {path!r} -- proceeding WITHOUT a kinship "
            f"filter (pilot only; a full run must not skip this).")
        return []
    df = pd.read_csv(path, sep="\t", dtype={"i.s": str, "j.s": str})
    missing = {"i.s", "j.s", "kin"} - set(df.columns)
    if missing:
        log(f"WARNING: relatedness table missing {missing}; columns {list(df.columns)} -- "
            f"proceeding WITHOUT a kinship filter.")
        return []
    return list(zip(df["i.s"], df["j.s"], df["kin"]))


def greedy_unrelated(person_ids, pairs, kin_min=KIN_MIN):
    ids = set(map(str, person_ids))
    adj = defaultdict(set)
    for i, j, k in pairs:
        i, j = str(i), str(j)
        try:
            k = float(k)
        except (TypeError, ValueError):
            continue
        if i == j or k < kin_min or i not in ids or j not in ids:
            continue
        adj[i].add(j)
        adj[j].add(i)
    removed = []
    while True:
        live = [(len(v), n) for n, v in adj.items() if v]
        if not live:
            break
        maxdeg = max(d for d, _ in live)
        victim = min(n for d, n in live if d == maxdeg)
        for nb in adj[victim]:
            adj[nb].discard(victim)
        adj[victim] = set()
        removed.append(victim)
    return ids - set(removed), removed


def pick_pilot_cohort(cohort_path, relatedness_path, n_people, seed):
    """Returns (person_ids list, selection_summary dict). selection_summary is aggregate-only
    (counts per ancestry/platform bucket actually picked) -- safe to print/return."""
    cohort = pd.read_csv(cohort_path, sep="\t", dtype=str, low_memory=False)
    cohort = cohort[cohort["in_lr"].astype(str) == "True"]
    cohort = cohort[cohort["platform"].isin(ASSEMBLY_PLATFORMS)]
    cohort["ancestry_pred"] = cohort["ancestry_pred"].fillna("UNASSIGNED").str.upper()

    pairs = load_relatedness_pairs(relatedness_path)
    kept, removed = greedy_unrelated(cohort["person_id"], pairs)
    cohort = cohort[cohort["person_id"].isin(kept)]
    log(f"[41] candidate pool after unrelated+assembly-platform filter: {len(cohort)} "
        f"({len(removed)} relatives dropped)")

    rng = random.Random(seed)
    buckets = defaultdict(list)
    for _, row in cohort.iterrows():
        buckets[(row["ancestry_pred"], row["platform"])].append(row["person_id"])
    for k in buckets:
        rng.shuffle(buckets[k])

    # Round-robin across (ancestry, platform) buckets so both dimensions get mixed, prioritizing
    # ancestries present in ANCESTRY_ORDER first (skip UNASSIGNED unless we run short), and always
    # trying to include >=1 sequel2 person (no aln-to-hg38 BAM -- exercises the Tier-3 self-align
    # path with the chr19 reference-cache workaround, per this script's own module docstring).
    ordered_keys = sorted(buckets.keys(), key=lambda k: (k[0] not in ANCESTRY_ORDER, k[0], k[1]))
    picked, summary = [], defaultdict(int)
    idx = 0
    while len(picked) < n_people and any(buckets[k] for k in ordered_keys):
        k = ordered_keys[idx % len(ordered_keys)]
        if buckets[k]:
            pid = buckets[k].pop()
            picked.append(pid)
            summary[k] += 1
        idx += 1
    return picked, {f"{a}/{p}": n for (a, p), n in sorted(summary.items())}


# ---------------------------------------------------------------------------
# Step 2: Tier-3 chr19 self-align reference cache (workaround for the hardcoded-chr6 bug -- see
# module docstring). Mirrors run_immuannot_person.py's own ensure_chr6_ref() logic exactly
# (idempotent: a second call with an existing non-empty cache is a no-op) but extracts chr19.
# ---------------------------------------------------------------------------
def ensure_kir_ref_cache(hg38_ref, cache_path):
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        return cache_path
    if not os.path.exists(hg38_ref):
        log(f"WARNING: {hg38_ref} not found -- Tier-3 self-align unavailable for sequel2 people "
            f"in this pilot (they will SKIP with a clear reason, not silently produce 0 calls).")
        return None
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    tmp = cache_path + ".tmp"
    with open(tmp, "w") as out:
        proc = subprocess.run(["samtools", "faidx", hg38_ref, "chr19"], stdout=out,
                               stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0 or os.path.getsize(tmp) == 0:
        log(f"WARNING: could not extract chr19 from {hg38_ref}: {proc.stderr}")
        os.remove(tmp)
        return None
    os.replace(tmp, cache_path)
    return cache_path


# ---------------------------------------------------------------------------
# Step 3: run run_immuannot_person.py per person (sequential -- respects the <=2-core pilot
# budget: hap1/hap2 run sequentially by default inside that script, each getting --threads).
# ---------------------------------------------------------------------------
def run_one_person(pid, args, kir_ref_cache, platform_of):
    cmd = [sys.executable, args.immuannot_script, str(pid),
           "--region", args.region, "--pad", str(args.pad),
           "--outroot", args.outroot, "--threads", str(args.threads),
           "--out-suffix", args.out_suffix, "--time-budget-min", str(args.time_budget_min)]
    if args.refdir:
        cmd += ["--refdir", args.refdir]
    if args.mount:
        cmd += ["--mount", args.mount]
    if platform_of.get(str(pid)) == "sequel2" and kir_ref_cache:
        cmd += ["--enable-self-align-fallback", "--chr6-ref-cache", kir_ref_cache,
                "--hg38-ref", args.hg38_ref]

    t0 = time.perf_counter()
    cpu0 = os.times()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    cpu1 = os.times()
    cpu = (cpu1.children_user - cpu0.children_user) + (cpu1.children_system - cpu0.children_system)
    ok = proc.returncode == 0
    log_path = os.path.join(args.outroot, f"41_person_stdout_stderr.{pid}.log")
    try:
        with open(log_path, "w") as f:
            f.write(f"# cmd: {' '.join(cmd)}\n# exit_code: {proc.returncode}\n"
                     f"# wall_seconds: {wall:.1f}\n# cpu_seconds: {cpu:.1f}\n"
                     f"# ---- STDOUT ----\n{proc.stdout}\n# ---- STDERR ----\n{proc.stderr}\n")
    except OSError as e:
        log(f"  WARNING: could not write per-person log {log_path}: {e}")
    return {"ok": ok, "wall_seconds": wall, "cpu_seconds": cpu, "returncode": proc.returncode,
            "log_path": log_path}


# ---------------------------------------------------------------------------
# Step 4: gene-content / novelty / sanity checks, re-parsing the raw per-hap gtf.gz directly
# (immuannot_calls.tsv keeps only ONE row per (person,gene) -- collapses multi-copy genes and
# drops the novelty ("...:new") tag entirely; SCHEMA.md Table 1 is the finer-grained per-copy
# analogue but this script doesn't write into that shared production table).
# ---------------------------------------------------------------------------
def parse_hap_gtf(gtf_path):
    """Returns list of dicts: {gene, consensus, is_novel, is_undetermined}, one per transcript
    row (one per gene-copy call -- multi-copy genes get >1 row, matching Table 1's grain)."""
    rows = []
    if not gtf_path or not os.path.exists(gtf_path):
        return rows
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            gene_m = re.search(r'gene_name "([^"]+)"', attrs)
            cons_m = re.search(r'consensus "([^"]+)"', attrs)
            if not gene_m or not cons_m:
                continue
            consensus = cons_m.group(1)
            rows.append({
                "gene": gene_m.group(1),
                "consensus": consensus,
                # Novelty suffix: confirmed on real pilot output (2026-09-22/23) that KIR
                # consensus strings append "new" DIRECTLY to the allele digits with NO colon
                # (e.g. "KIR3DL2*00201new"), unlike IMMUANNOT_GTF_SPEC.md's HLA-style illustrative
                # example ("HLA-A*01:01:new") -- HLA alleles are colon-field-delimited, KIR's are
                # not, so there is no colon for "new" to follow. Checking endswith("new") alone
                # (not ":new") catches both conventions; a first version of this script checked
                # ":new" only and silently measured 0% novelty on 383 real KIR calls -- caught by
                # spot-checking one real hap1.gtf.gz's raw consensus strings, not by inspection of
                # this code alone. "undetermined" does not end in "new" so it can't collide.
                "is_novel": consensus.rstrip('"').endswith("new") and consensus != "undetermined",
                "is_undetermined": consensus == "undetermined",
            })
    return rows


def aggregate_quality(person_ids, args):
    """Reads back each successfully-run person's hap1/hap2 gtf.gz (still on the VM, under
    --outroot) and computes AGGREGATE-ONLY statistics. No person_id, no allele string, leaves
    this function -- only counts/rates/fractions in the returned dict."""
    n_hap_attempted = 0
    n_hap_with_calls = 0
    genes_per_hap = []
    framework_present_count = defaultdict(int)
    ac_hap_content = defaultdict(int)  # "cA" / "cB" counts
    n_2dl2_and_2dl3 = 0
    n_3dl1_and_3ds1 = 0
    n_haps_with_2dl2_or_2dl3 = 0
    n_haps_with_3dl1_or_3ds1 = 0
    total_calls, total_novel, total_undetermined = 0, 0, 0
    kir_gene_hap_presence = defaultdict(int)

    for pid in person_ids:
        for hap in ("hap1", "hap2"):
            gtf_path = os.path.join(args.outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
            n_hap_attempted += 1
            rows = parse_hap_gtf(gtf_path)
            kir_rows = [r for r in rows if r["gene"] in KIR_GENES]
            if not kir_rows:
                continue
            n_hap_with_calls += 1
            genes_called = sorted({r["gene"] for r in kir_rows})
            genes_per_hap.append(len(genes_called))
            for g in genes_called:
                kir_gene_hap_presence[g] += 1
            for g in FRAMEWORK_GENES:
                if g in genes_called:
                    framework_present_count[g] += 1
            content = "cB" if B_CONTENT_GENES.intersection(genes_called) else "cA"
            ac_hap_content[content] += 1
            has_2dl2, has_2dl3 = "KIR2DL2" in genes_called, "KIR2DL3" in genes_called
            has_3dl1, has_3ds1 = "KIR3DL1" in genes_called, "KIR3DS1" in genes_called
            if has_2dl2 or has_2dl3:
                n_haps_with_2dl2_or_2dl3 += 1
            if has_2dl2 and has_2dl3:
                n_2dl2_and_2dl3 += 1
            if has_3dl1 or has_3ds1:
                n_haps_with_3dl1_or_3ds1 += 1
            if has_3dl1 and has_3ds1:
                n_3dl1_and_3ds1 += 1
            for r in kir_rows:
                total_calls += 1
                total_novel += int(r["is_novel"])
                total_undetermined += int(r["is_undetermined"])

    def pct(n, d):
        return round(100.0 * n / d, 1) if d else None

    return {
        "n_people": len(person_ids),
        "n_hap_attempted": n_hap_attempted,
        "n_hap_with_kir_calls": n_hap_with_calls,
        "hap_success_rate_pct": pct(n_hap_with_calls, n_hap_attempted),
        "mean_kir_genes_per_hap": (round(sum(genes_per_hap) / len(genes_per_hap), 2)
                                    if genes_per_hap else None),
        "framework_gene_presence_pct": {g: pct(framework_present_count[g], n_hap_with_calls)
                                         for g in FRAMEWORK_GENES},
        "hap_content_A_vs_B": dict(ac_hap_content),
        "kir_gene_hap_presence_pct": {g: pct(kir_gene_hap_presence.get(g, 0), n_hap_with_calls)
                                       for g in KIR_GENES},
        "pct_haps_2dl2_and_2dl3_cooccur": pct(n_2dl2_and_2dl3, n_haps_with_2dl2_or_2dl3),
        "pct_haps_3dl1_and_3ds1_cooccur": pct(n_3dl1_and_3ds1, n_haps_with_3dl1_or_3ds1),
        "total_kir_calls": total_calls,
        "pct_novel": pct(total_novel, total_calls),
        "pct_undetermined": pct(total_undetermined, total_calls),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-people", type=int, default=20)
    ap.add_argument("--seed", type=int, default=41)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--immuannot-script", default=DEFAULT_IMMUANNOT_SCRIPT)
    ap.add_argument("--refdir", default=None)
    ap.add_argument("--mount", default=None)
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--out-suffix", default=".kir41")
    ap.add_argument("--region", default=KIR_REGION_DEFAULT)
    ap.add_argument("--pad", type=int, default=KIR_PAD_DEFAULT)
    ap.add_argument("--threads", type=int, default=2,
                    help="Capped at 2 for this pilot per CLAUDE.md ('up to 2 cores for the "
                         "pilot; do not start the full-cohort run yourself').")
    ap.add_argument("--time-budget-min", type=float, default=30)
    ap.add_argument("--hg38-ref", default=DEFAULT_HG38_REF)
    ap.add_argument("--kir-ref-cache", default=DEFAULT_KIR_REF_CACHE)
    ap.add_argument("--dry-run", action="store_true", help="Only select+print cohort summary.")
    ap.add_argument("--force", action="store_true",
                    help="Re-run a person even if both hap1/hap2.gtf.gz already exist under "
                         "--outroot (default: skip -- makes a killed-and-relaunched pilot resumable).")
    args = ap.parse_args()

    if args.threads > 2:
        log(f"WARNING: --threads {args.threads} > 2, clamping to 2 (pilot compute budget).")
        args.threads = 2

    os.makedirs(args.outroot, exist_ok=True)

    log(f"[41] region={args.region} pad={args.pad}")
    person_ids, selection_summary = pick_pilot_cohort(
        args.cohort_membership, args.relatedness_table, args.n_people, args.seed)
    log(f"[41] selected {len(person_ids)} people. ancestry/platform mix: {selection_summary}")

    if args.dry_run:
        log("[41] --dry-run: stopping after selection.")
        return

    cohort_df = pd.read_csv(args.cohort_membership, sep="\t", dtype=str, low_memory=False)
    platform_of = dict(zip(cohort_df["person_id"].astype(str), cohort_df["platform"]))
    n_sequel2 = sum(1 for p in person_ids if platform_of.get(str(p)) == "sequel2")
    kir_ref_cache = ensure_kir_ref_cache(args.hg38_ref, args.kir_ref_cache) if n_sequel2 else None
    if n_sequel2:
        log(f"[41] {n_sequel2} sequel2 people in pilot -- Tier-3 self-align enabled with "
            f"chr19 ref cache: {kir_ref_cache}")

    results = []
    t_run0 = time.perf_counter()
    for i, pid in enumerate(person_ids, 1):
        h1 = os.path.join(args.outroot, str(pid), "immuannot_output", "hap1.gtf.gz")
        h2 = os.path.join(args.outroot, str(pid), "immuannot_output", "hap2.gtf.gz")
        if not args.force and os.path.exists(h1) and os.path.exists(h2):
            log(f"[41] person {i}/{len(person_ids)}: SKIP -- already has both hap gtf.gz "
                f"outputs (resumable run; pass --force to redo).")
            results.append({"ok": True, "wall_seconds": 0.0, "cpu_seconds": 0.0,
                            "returncode": 0, "log_path": None,
                            "platform": platform_of.get(str(pid)), "skipped": True})
            continue
        log(f"[41] person {i}/{len(person_ids)} ...")
        r = run_one_person(pid, args, kir_ref_cache, platform_of)
        r["platform"] = platform_of.get(str(pid))
        results.append(r)
        log(f"[41] person {i}/{len(person_ids)}: {'OK' if r['ok'] else 'FAILED'} "
            f"wall={r['wall_seconds']:.1f}s cpu={r['cpu_seconds']:.1f}s")
    run_wall = time.perf_counter() - t_run0

    n_ok = sum(1 for r in results if r["ok"])
    mean_wall = sum(r["wall_seconds"] for r in results) / len(results) if results else 0
    mean_cpu = sum(r["cpu_seconds"] for r in results) / len(results) if results else 0
    log(f"[41] === Pilot run summary: {n_ok}/{len(results)} people ok, total wall "
        f"{run_wall/60:.1f} min, mean/person wall {mean_wall:.1f}s cpu {mean_cpu:.1f}s ===")

    quality = aggregate_quality(person_ids, args)
    log(f"[41] === Quality (aggregate-only): {quality} ===")

    summary_path = os.path.join(args.outroot, f"41_kir_pilot_summary{args.out_suffix}.tsv")
    pd.DataFrame([{
        "n_people": len(person_ids), "n_people_ok": n_ok,
        "person_success_rate_pct": round(100.0 * n_ok / len(person_ids), 1) if person_ids else None,
        "total_wall_min": round(run_wall / 60, 2),
        "mean_wall_sec_per_person": round(mean_wall, 1),
        "mean_cpu_sec_per_person": round(mean_cpu, 1),
        "selection_summary": str(selection_summary),
        **{f"q_{k}": str(v) for k, v in quality.items()},
    }]).to_csv(summary_path, sep="\t", index=False)
    log(f"[41] Aggregate summary written: {summary_path}")


if __name__ == "__main__":
    main()
