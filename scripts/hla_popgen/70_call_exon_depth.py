#!/usr/bin/env python3
"""S06 70: exon-, CDS- and per-base long-read depth of every HLA gene, per person (VM).

Why: 58/59/67 measured depth over whole gene spans or the 4 Mb MHC window. Two-field typing is decided by the
exons, so the quantity that should predict a typing error is depth over the exons / coding sequence of the called
gene, and how evenly that sequence is covered. This script measures it, per person, in one pass over the
person's GRCh38 long-read BAM.

Features (built from a GENCODE GTF by `--build-features`; GRCh38, 0-based half-open BED):
  gene          the gene row (what 58 measured)
  exons         union of the exons of the gene's protein-coding transcripts
  cds           union of the CDS rows of those transcripts (no UTR, no introns)
  ex<k>         each exon of the reference transcript (MANE Select, else Ensembl canonical, else the transcript
                with the longest CDS), k = GENCODE exon_number; `cds_ord` counts only exons that carry CDS
  ARS flag      the exons that encode the peptide-binding groove: class I (A, B, C, E, F, G) CDS exons 2 and 3,
                class II (DRB*, DQB1, DPB1, DRA, DQA*, DPA1) CDS exon 2. These decide 2-field typing in practice.

Per person (INTERNAL, `--internal-dir`):
  feature_depth.tsv      person x feature x MAPQ floor (0, 20): sum depth, reads, mean depth (samtools bedcov -c)
  cds_perbase.tsv        person x gene x MAPQ floor: over the CDS union and over its ARS part, per-base minimum, 5th and
                         10th percentile, median, mean, share of bases at 0x and under 5/10/20x, MAPQ-agnostic
                         (samtools depth -a over cds.bed)
  failures.tsv           person + reason. A missing or slow BAM is counted, never dropped silently.
Relative depth is NOT computed here: 71 divides by 58's per-person control-window depth (same MAPQ floor), so the
two scripts agree on the denominator.

Public (aggregate reference data only, safe to commit): features.tsv / features.bed (GENCODE coordinates),
profile_by_gene.tsv (`--profile`: depth along each gene in 10 bp bins, median / quartiles / 10th-90th percentile
over a random sample of people, with the exon structure, to draw the methodology figure).

VM (samtools lives in the spechla pixi env, as for 58):
    cd ~/repos/pilot-validation
    pixi run -e spechla -- python3 scripts/hla_popgen/70_call_exon_depth.py --build-features \\
        --gencode-gtf ~/ref/gencode.v46.basic.annotation.gtf.gz --out-dir ~/s06/results/70
    setsid nohup pixi run -e spechla -- python3 -u scripts/hla_popgen/70_call_exon_depth.py \\
        --features-bed ~/s06/results/70/features.bed --out-dir ~/s06/results/70 \\
        --internal-dir ~/s06/internal/70 --priority-file ~/s06/internal/70/priority_ids.txt --workers 8 \\
        > ~/s06/logs/70.log 2>&1 < /dev/null & disown
    Smoke test first: add --limit 20 --internal-dir ~/s06/internal/70_smoke (own path, quirk #22).
    Population profile (sample of 300): add --profile 300 (a separate invocation; no per-person tables are written).

Tests: scripts/hla_popgen/tests/test_70_call_exon_depth.py (synthetic, no samtools).
"""
import argparse
import gzip
import importlib.util
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import _disclosure as _disc  # noqa: E402


def _load58():
    spec = importlib.util.spec_from_file_location("rd58", os.path.join(HERE, "58_regional_depth.py"))
    m = importlib.util.module_from_spec(spec)
    sys.modules["rd58"] = m
    spec.loader.exec_module(m)
    return m


RD = _load58()   # load_people, uri_to_local, parse_bedcov helpers; one definition of "a person's BAM"

MAPQS = (0, 20)
DEFAULT_OUT_DIR = os.path.expanduser("~/s06/results/70")
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s06/internal/70")

# Genes whose typing we score (56's ALLG) plus DRA / C4A for the contig matrix. GENCODE spells HLA genes with the
# HLA- prefix; MICA, TAP1, TAP2, C4A, C4B are bare. DRB3/4/5 and others that GRCh38 carries only on ALT contigs
# are not on a noalt BAM and are reported as absent, not as zero depth.
GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
         "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DQA2", "HLA-DRA", "MICA", "TAP1", "TAP2", "C4A", "C4B"]
CLASS_I = {"HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"}
CLASS_II = {"HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DQB1", "HLA-DPB1",
            "HLA-DRA", "HLA-DQA1", "HLA-DQA2", "HLA-DPA1"}
FEATURE_COLS = ["chrom", "start", "end", "name", "gene", "kind", "exon_number", "cds_ord", "ars",
                "length", "strand", "transcript_id"]


def log(msg):
    print(f"[70 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def bare(gene):
    return gene[4:] if gene.startswith("HLA-") else gene


# ---------------------------------------------------------------------------
# Features from GENCODE
# ---------------------------------------------------------------------------
def _attr(attrs, key):
    m = re.search(key + r' "([^"]+)"', attrs)
    return m.group(1) if m else None


def _tags(attrs):
    return set(re.findall(r'tag "([^"]+)"', attrs))


def merge_intervals(iv):
    """Union of (start, end) half-open intervals."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def parse_gtf(gtf_path, genes=GENES, chroms=("chr6",)):
    """-> {gene: dict(chrom, strand, gene=(s,e), tx={tid: dict(type, tags, exons={n:(s,e)}, cds=[(s,e)])})}.
    Coordinates converted to 0-based half-open."""
    want = set(genes)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    out = {}
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in chroms or f[2] not in ("gene", "transcript", "exon", "CDS"):
                continue
            name = _attr(f[8], "gene_name")
            if name not in want:
                continue
            g = out.setdefault(name, dict(chrom=f[0], strand=f[6], gene=None, tx={}))
            s0, e = int(f[3]) - 1, int(f[4])
            if f[2] == "gene":
                g["gene"] = (s0, e)
                continue
            tid = _attr(f[8], "transcript_id")
            t = g["tx"].setdefault(tid, dict(type=None, tags=set(), exons={}, cds=[]))
            if f[2] == "transcript":
                t["type"] = _attr(f[8], "transcript_type") or _attr(f[8], "transcript_biotype")
                t["tags"] = _tags(f[8])
            elif f[2] == "exon":
                t["exons"][int(_attr(f[8], "exon_number") or len(t["exons"]) + 1)] = (s0, e)
            else:
                t["cds"].append((s0, e))
    return out


def reference_transcript(g):
    """MANE Select, else Ensembl canonical, else the protein-coding transcript with the longest CDS."""
    coding = {tid: t for tid, t in g["tx"].items() if t["cds"]}
    if not coding:
        return None
    for tag in ("MANE_Select", "Ensembl_canonical"):
        hit = [tid for tid, t in coding.items() if tag in t["tags"]]
        if hit:
            return sorted(hit)[0]
    return max(coding, key=lambda tid: (sum(e - s for s, e in coding[tid]["cds"]), tid))


def is_ars(gene, cds_ord):
    if gene in CLASS_I:
        return cds_ord in (2, 3)
    if gene in CLASS_II:
        return cds_ord == 2
    return False


def build_features(gtf_path, genes=GENES):
    """DataFrame FEATURE_COLS, plus the list of requested genes the GTF does not have on chr6 (ALT-only)."""
    parsed = parse_gtf(gtf_path, genes)
    rows, missing = [], []
    for gene in genes:
        g = parsed.get(gene)
        if not g or g["gene"] is None or not any(t["cds"] for t in g["tx"].values()):
            missing.append(gene)
            continue
        coding = {tid: t for tid, t in g["tx"].items() if t["cds"]}
        ref = reference_transcript(g)
        base = dict(chrom=g["chrom"], gene=bare(gene), strand=g["strand"])
        ex_union = merge_intervals([iv for t in coding.values() for iv in t["exons"].values()])
        cds_union = merge_intervals([iv for t in coding.values() for iv in t["cds"]])
        s, e = g["gene"]
        rows.append(dict(base, start=s, end=e, kind="gene", exon_number=0, cds_ord=0, ars=False, transcript_id=""))
        for k, (s, e) in enumerate(ex_union, 1):
            rows.append(dict(base, start=s, end=e, kind="exons", exon_number=k, cds_ord=0, ars=False,
                             transcript_id="union"))
        for k, (s, e) in enumerate(cds_union, 1):
            rows.append(dict(base, start=s, end=e, kind="cds", exon_number=k, cds_ord=0, ars=False,
                             transcript_id="union"))
        # per-exon intervals of the reference transcript; cds_ord counts only exons overlapping a CDS row
        t = g["tx"][ref]
        order = sorted(t["exons"], key=lambda n: n)
        ord_ = 0
        for n in order:
            xs, xe = t["exons"][n]
            cds_here = [(max(xs, a), min(xe, b)) for a, b in t["cds"] if min(xe, b) > max(xs, a)]
            if not cds_here:
                continue          # pure UTR exon: not a coding exon, not measured
            ord_ += 1
            rows.append(dict(base, start=min(a for a, _ in cds_here), end=max(b for _, b in cds_here), kind="ex",
                             exon_number=n, cds_ord=ord_, ars=is_ars(gene, ord_), transcript_id=ref))
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=FEATURE_COLS), missing
    df["length"] = df["end"] - df["start"]
    df["name"] = [f"{g}|{k}" + (f"{n}" if k in ("ex", "exons", "cds") else "")
                  for g, k, n in zip(df["gene"], df["kind"], df["exon_number"])]
    # "ex" rows are clipped to their CDS part, so an ex interval measures coding bases only (UTR bases of the
    # first / last coding exon excluded); the name says which exon_number of the reference transcript it is.
    return df[FEATURE_COLS].reset_index(drop=True), missing


def write_bed(df, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    df[["chrom", "start", "end", "name"]].to_csv(path, sep="\t", header=False, index=False)


def cds_bed(df):
    return df[df["kind"] == "cds"][["chrom", "start", "end", "gene"]]


def ars_intervals(df):
    return df[(df["kind"] == "ex") & df["ars"]][["chrom", "start", "end", "gene"]]


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
def run_bedcov(bam, bed_path, mapq, timeout, samtools="samtools"):
    out = subprocess.run([samtools, "bedcov", "-c", "-Q", str(mapq), bed_path, bam],
                         capture_output=True, text=True, timeout=timeout, check=True).stdout
    return RD.parse_bedcov(out)


def run_depth(bam, bed_path, mapq, timeout, samtools="samtools"):
    """samtools depth -a over a BED -> DataFrame chrom, pos (1-based), depth."""
    out = subprocess.run([samtools, "depth", "-a", "-Q", str(mapq), "-b", bed_path, bam],
                         capture_output=True, text=True, timeout=timeout, check=True).stdout
    return parse_depth(out)


def parse_depth(text):
    rows = [ln.split("\t") for ln in text.splitlines() if ln.strip()]
    if not rows:
        return pd.DataFrame(columns=["chrom", "pos", "depth"])
    return pd.DataFrame({"chrom": [r[0] for r in rows], "pos": np.array([int(r[1]) for r in rows]),
                         "depth": np.array([int(r[2]) for r in rows])})


def perbase_summary(depth_df, cds, ars):
    """depth_df over the CDS union; cds / ars: DataFrames chrom,start,end,gene (0-based half-open).
    -> {gene: stats dict}, in the person's absolute depth units."""
    out = {}
    if depth_df.empty:
        return out
    pos0 = depth_df["pos"].to_numpy() - 1
    dep = depth_df["depth"].to_numpy().astype(float)
    chrom = depth_df["chrom"].to_numpy()
    for gene, sub in cds.groupby("gene"):
        m = np.zeros(len(pos0), bool)
        for s, e in zip(sub["start"], sub["end"]):
            m |= (pos0 >= s) & (pos0 < e) & (chrom == sub["chrom"].iloc[0])
        a = np.zeros(len(pos0), bool)
        for s, e in zip(ars.loc[ars["gene"] == gene, "start"], ars.loc[ars["gene"] == gene, "end"]):
            a |= (pos0 >= s) & (pos0 < e) & (chrom == sub["chrom"].iloc[0])
        d = dep[m]
        if not len(d):
            continue
        row = dict(cds_bases=int(len(d)), cds_min=float(d.min()), cds_p05=float(np.percentile(d, 5)),
                   cds_p10=float(np.percentile(d, 10)), cds_median=float(np.median(d)), cds_mean=float(d.mean()),
                   cds_share_0=float((d == 0).mean()), cds_share_lt5=float((d < 5).mean()),
                   cds_share_lt10=float((d < 10).mean()), cds_share_lt20=float((d < 20).mean()))
        da = dep[a]
        if len(da):
            row.update(ars_bases=int(len(da)), ars_min=float(da.min()), ars_mean=float(da.mean()),
                       ars_share_lt10=float((da < 10).mean()))
        out[gene] = row
    return out


FEATURE_DEPTH_COLS = ["person_id", "feature", "mapq", "sum_depth", "n_reads", "mean_depth"]
PERBASE_COLS = ["person_id", "gene", "mapq", "cds_bases", "cds_min", "cds_p05", "cds_p10", "cds_median", "cds_mean",
                "cds_share_0", "cds_share_lt5", "cds_share_lt10", "cds_share_lt20",
                "ars_bases", "ars_min", "ars_mean", "ars_share_lt10"]


def measure_person(person, feat_bed, cds_path, cds, ars, timeout, samtools="samtools",
                   bedcov=None, depth=None):
    """-> (feature_rows, perbase_rows, None) or (None, None, reason)."""
    bedcov, depth = bedcov or run_bedcov, depth or run_depth
    bam = person["bam"]
    if not bam or not os.path.exists(bam):
        return None, None, "bam_missing"
    frows, prows = [], []
    for q in MAPQS:
        try:
            df = bedcov(bam, feat_bed, q, timeout, samtools)
            dp = depth(bam, cds_path, q, timeout, samtools)
        except subprocess.TimeoutExpired:
            return None, None, "timeout"
        except subprocess.CalledProcessError:
            return None, None, "samtools_error"
        df["mean_depth"] = df["sum_depth"] / (df["end"] - df["start"]).clip(lower=1)
        for r in df.itertuples():
            frows.append(dict(person_id=person["person_id"], feature=r.name, mapq=q, sum_depth=r.sum_depth,
                              n_reads=r.n_reads, mean_depth=round(r.mean_depth, 4)))
        for gene, st in perbase_summary(dp, cds, ars).items():
            prows.append(dict(person_id=person["person_id"], gene=gene, mapq=q,
                              **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in st.items()}))
    return frows, prows, None


# ---------------------------------------------------------------------------
# Checkpoint, lock, resume
# ---------------------------------------------------------------------------
def acquire_lock(internal_dir):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, "70.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"FATAL: another 70_call_exon_depth.py (pid {pid}) holds {lock}")
        except (ValueError, ProcessLookupError):
            log(f"removing stale lock {lock}")
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


def append_rows(path, rows, cols):
    new = not os.path.exists(path)
    if new:
        _disc.assert_internal_path_allowed(path)
    with open(path, "a") as fh:
        if new:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
            fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in cols) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def done_people(internal_dir):
    done = set()
    for name in ("feature_depth.tsv", "failures.tsv"):
        p = os.path.join(internal_dir, name)
        if os.path.exists(p):
            done |= set(pd.read_csv(p, sep="\t", dtype={"person_id": str}, comment="#")["person_id"])
    return done


def order_people(people, priority_ids):
    """Truth-set people (priority) first, so a partial run already supports the error analysis."""
    pri = set(priority_ids)
    key = people["person_id"].map(lambda p: 0 if p in pri else 1)
    return people.assign(_k=key).sort_values(["_k", "person_id"], kind="stable").drop(columns="_k").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Population depth profile along each gene (public aggregate)
# ---------------------------------------------------------------------------
BIN = 10
FLANK = 500


def profile_regions(features):
    g = features[features["kind"] == "gene"]
    return pd.DataFrame({"chrom": g["chrom"], "start": (g["start"] - FLANK).clip(lower=0), "end": g["end"] + FLANK,
                         "gene": g["gene"]})


def bin_profile(dp, regions):
    """Mean depth per BIN-bp bin of each region -> {gene: array}."""
    out = {}
    pos0 = dp["pos"].to_numpy() - 1
    dep = dp["depth"].to_numpy().astype(float)
    for r in regions.itertuples():
        n = int(np.ceil((r.end - r.start) / BIN))
        m = (pos0 >= r.start) & (pos0 < r.end)
        idx = (pos0[m] - r.start) // BIN
        s = np.bincount(idx, weights=dep[m], minlength=n)[:n]
        c = np.bincount(idx, minlength=n)[:n]
        out[r.gene] = np.divide(s, c, out=np.full(n, np.nan), where=c > 0)
    return out


def summarise_profile(stack, regions, features, ctrl_depth):
    """stack: {gene: 2D array people x bins}; depth divided by each person's control depth (relative)."""
    rows = []
    for r in regions.itertuples():
        a = np.vstack(stack[r.gene]) / np.asarray(ctrl_depth)[:, None]
        for b in range(a.shape[1]):
            col = a[:, b]
            col = col[np.isfinite(col)]
            if len(col) < 20:
                continue
            q = np.percentile(col, [10, 25, 50, 75, 90])
            rows.append(dict(gene=r.gene, bin_start=int(r.start + b * BIN), bin_end=int(r.start + (b + 1) * BIN),
                             n_people=len(col), rel_p10=round(q[0], 4), rel_p25=round(q[1], 4),
                             rel_p50=round(q[2], 4), rel_p75=round(q[3], 4), rel_p90=round(q[4], 4)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-features", action="store_true")
    ap.add_argument("--gencode-gtf")
    ap.add_argument("--features-bed")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    ap.add_argument("--cohort", default=RD.DEFAULT_COHORT)
    ap.add_argument("--manifest", default=os.path.join(RD.AOU_ROOT, RD.LR_MANIFEST_REL))
    ap.add_argument("--control58", default=os.path.expanduser("~/s05/internal/58/per_person_depth.tsv"),
                    help="58's per-person table; only used by --profile for the control-depth denominator")
    ap.add_argument("--priority-file", help="person ids (one per line) to measure first")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--profile", type=int, help="measure the depth profile of this many random people (no per-person tables)")
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)

    feat_path = os.path.join(a.out_dir, "features.tsv")
    if a.build_features:
        if not a.gencode_gtf:
            raise SystemExit("--build-features needs --gencode-gtf")
        feats, missing = build_features(a.gencode_gtf)
        feats.to_csv(feat_path, sep="\t", index=False)
        write_bed(feats, os.path.join(a.out_dir, "features.bed"))
        write_bed(feats[feats["kind"] == "cds"].assign(name=lambda d: d["gene"]), os.path.join(a.out_dir, "cds.bed"))
        log(f"features: {len(feats)} rows, genes {sorted(feats['gene'].unique())}; not on chr6 primary: {missing}")
        return 0

    feats = pd.read_csv(feat_path, sep="\t")
    feat_bed = a.features_bed or os.path.join(a.out_dir, "features.bed")
    cds_path = os.path.join(a.out_dir, "cds.bed")
    cds, ars = cds_bed(feats), ars_intervals(feats)

    people = RD.load_people(a.cohort, a.manifest)
    if a.priority_file and os.path.exists(a.priority_file):
        people = order_people(people, [x.strip() for x in open(a.priority_file) if x.strip()])

    if a.profile:
        return run_profile(a, feats, people)

    lock = acquire_lock(a.internal_dir)
    try:
        done = done_people(a.internal_dir)
        todo = people[~people["person_id"].isin(done)]
        if a.limit:
            todo = todo.head(a.limit)
        log(f"people {len(people)}, done {len(done)}, to run {len(todo)}")
        fpath = os.path.join(a.internal_dir, "feature_depth.tsv")
        ppath = os.path.join(a.internal_dir, "cds_perbase.tsv")
        xpath = os.path.join(a.internal_dir, "failures.tsv")
        t0, nfail = time.time(), 0
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(measure_person, p, feat_bed, cds_path, cds, ars, a.timeout, a.samtools): p["person_id"]
                    for p in todo.to_dict("records")}
            for k, fu in enumerate(as_completed(futs), 1):
                pid = futs[fu]
                fr, pr, why = fu.result()
                if why:
                    nfail += 1
                    append_rows(xpath, [dict(person_id=pid, reason=why)], ["person_id", "reason"])
                else:
                    append_rows(fpath, fr, FEATURE_DEPTH_COLS)
                    append_rows(ppath, pr, PERBASE_COLS)
                if k % 100 == 0 or k == len(futs):
                    log(f"{k}/{len(futs)} done, {nfail} failed, {(time.time() - t0) / k:.1f} s/person")
        log("DONE 70")
    finally:
        try:
            os.remove(lock)
        except OSError:
            pass
    return 0


def run_profile(a, feats, people):
    regions = profile_regions(feats)
    rng = np.random.default_rng(70)
    pick = people[people["bam"].map(lambda b: bool(b) and os.path.exists(b))]
    pick = pick.iloc[rng.permutation(len(pick))[: a.profile]]
    reg_path = os.path.join(a.out_dir, "profile_regions.bed")
    regions[["chrom", "start", "end", "gene"]].to_csv(reg_path, sep="\t", header=False, index=False)
    c58 = pd.read_csv(a.control58, sep="\t", dtype={"person_id": str}, comment="#")
    c58 = c58[c58["mapq"] == 0].drop_duplicates("person_id").set_index("person_id")["control_depth"]
    stack = {g: [] for g in regions["gene"]}
    ctrl = []

    def one(p):
        dp = run_depth(p["bam"], reg_path, 0, a.timeout, a.samtools)
        return p["person_id"], bin_profile(dp, regions)

    with ThreadPoolExecutor(a.workers) as ex:
        for k, (pid, prof) in enumerate(ex.map(one, pick.to_dict("records")), 1):
            if pid in c58.index:
                ctrl.append(float(c58[pid]))
                for g, arr in prof.items():
                    stack[g].append(arr)
            if k % 50 == 0:
                log(f"profile {k}/{len(pick)}")
    out = summarise_profile(stack, regions, feats, ctrl)
    out.to_csv(os.path.join(a.out_dir, "profile_by_gene.tsv"), sep="\t", index=False)
    log(f"profile_by_gene.tsv: {len(out)} bins from {len(ctrl)} people")
    return 0


if __name__ == "__main__":
    sys.exit(main())
