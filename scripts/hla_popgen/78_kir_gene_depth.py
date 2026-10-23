#!/usr/bin/env python3
"""S06 78: exon / CDS depth for the KIR genes that exist on primary chr19 of GRCh38 (VM; the KIR counterpart of 70).

Same measurement as 70: ONE index-backed `samtools view -u BAM chr19:54500001-55100000 | samtools depth -a -Q q -b gene_regions.bed -`
per MAPQ floor (0 and 20), every feature derived from the per-base arrays (70's own helpers are imported, not copied).
Genes: KIR2DL1, KIR2DL3, KIR2DL4, KIR2DP1, KIR2DS4, KIR3DL1, KIR3DL2, KIR3DL3, KIR3DP1, KIR3DX1 (GENCODE v46 reference transcript, CDS rows;
KIR genes have no peptide-binding-exon convention, so there is no ARS metric).

WHAT THIS CAN AND CANNOT SAY
  - KIR2DL2, KIR2DL5A/B, KIR2DS1, KIR2DS2, KIR2DS3, KIR2DS5 and KIR3DS1 are not on primary chr19 of GRCh38 (alternative-haplotype contigs only,
    which the alignments do not carry). They get NO depth metric (NA in 78b, never 0). Their reads, when present, map onto their
    paralogues, so a depth value read off a paralogue would be a statement about the wrong gene.
  - The KIR genes are >90% identical in sequence and sit in segmental duplications. Mapping to one reference is reference-biased and multi-mapping:
    reads of a gene absent from the reference pile onto a paralogue (depth above expectation), reads of a divergent allele soak into MAPQ 0.
    The quantitative indicator is the MAPQ share: mapq_frac = mean CDS depth at MAPQ>=20 / mean CDS depth at MAPQ>=0 (1 = every read maps uniquely),
    and the MAPQ>=0 depth relative to the person's own control depth (a value far from 1 means collapsed or lost copies). Both are written.
  - Depth is relative to the person's own control depth: the mean depth in two 100 kb control windows on chr19 away from the LRC
    (default chr19:51,000,001-51,100,000 and chr19:56,000,001-56,100,000; --control-regions to change).

Outputs
  --build-features --gencode-gtf GTF     results/78/{features.tsv, features.bed, gene_regions.bed, control_<i>.bed}
  run (default)                          INTERNAL ~/s06/internal/78/{feature_depth.tsv, cds_perbase.tsv, control_depth.tsv, failures.tsv}   (checkpointed, --resume by default)
  --summarise                            INTERNAL ~/s06/internal/78/gene_depth_metrics.INTERNAL.tsv  (person x gene) + aggregate distributions results/78/

People: by default only those in 77c's INTERNAL records (trio members + pair members), because that is what the error analysis needs and
it is a few thousand BAM reads rather than the cohort; --all measures every long-read person.

    cd ~/repos/pilot-validation-main && P="pixi run -e spechla -- python3 ~/repos/pv-s06/scripts/hla_popgen/78_kir_gene_depth.py"
    $P --build-features --gencode-gtf ~/ref/gencode.v46.basic.annotation.gtf.gz
    $P --limit 20 --internal-dir ~/s06/internal/78_smoke --workers 4        # smoke
    setsid nohup $P --workers 4 > ~/s06/logs/78.log 2>&1 < /dev/null & disown
    $P --summarise
"""
import argparse
import importlib.util
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402
import _disclosure as _disc  # noqa: E402


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


M70 = sys.modules.get("m70") or _load("m70", "70_call_exon_depth.py")
RD = M70.RD

KIR_ON_CHR19 = ["KIR2DL1", "KIR2DL3", "KIR2DL4", "KIR2DP1", "KIR2DS4", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DX1"]
NOT_ON_GRCH38 = ["KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR3DS1"]
KIR_REGION = "chr19:54500001-55100000"
DEFAULT_CONTROL = "chr19:51000001-51100000,chr19:56000001-56100000"
OUT = os.path.expanduser("~/s06/results/78/")
INT = os.path.expanduser("~/s06/internal/78/")
I77C = os.path.expanduser("~/s06/internal/77c/")
I77 = os.path.expanduser("~/s06/internal/77/")
MAPQS = M70.MAPQS


def log(msg):
    C.log("[78] " + msg)


# ---------------------------------------------------------------------------
# features: 70's builder pointed at chr19
# ---------------------------------------------------------------------------
def build_kir_features(gtf_path, genes=KIR_ON_CHR19):
    """70.build_features with the chromosome restricted to chr19 (70's builder looks up the module-level parse_gtf)."""
    orig = M70.parse_gtf

    def on_chr19(path, genes=genes, chroms=("chr19",)):
        return orig(path, genes, chroms=("chr19",))

    M70.parse_gtf = on_chr19
    try:
        return M70.build_features(gtf_path, genes)
    finally:
        M70.parse_gtf = orig


def parse_regions(spec):
    out = []
    for r in re.split(r"[;\s]+|,(?=chr)", spec):
        r = r.strip()
        if not r:
            continue
        chrom, rng = r.split(":")
        s, e = rng.replace(",", "").split("-")
        out.append((chrom, int(s), int(e)))
    return out


def write_control_beds(out_dir, regions):
    """One BED per control window (0-based half-open), so each `samtools view` region is paired with its own -b file."""
    paths = []
    for i, (chrom, s, e) in enumerate(regions):
        p = os.path.join(out_dir, "control_%d.bed" % i)
        with open(p, "w") as fh:
            fh.write("%s\t%d\t%d\tcontrol%d\n" % (chrom, s - 1, e, i))
        paths.append((p, "%s:%d-%d" % (chrom, s, e)))
    return paths


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------
def control_mean(bam, mapq, controls, timeout, samtools, depth):
    """Mean per-base depth over the control windows (positions with no read count as 0 because `samtools depth -a`)."""
    tot, n = 0.0, 0
    for bed, region in controls:
        dp = depth(bam, bed, mapq, timeout, samtools, region)
        tot += float(dp["depth"].sum())
        n += len(dp)
    return tot / n if n else float("nan")


def measure_person(person, regions_bed, regions, features, controls, timeout, samtools="samtools", depth=None):
    """-> (feature_rows, perbase_rows, control_rows, None) or (None, None, None, reason). `depth(bam, bed, mapq, timeout, samtools, region)` is injectable."""
    depth = depth or M70.run_depth
    bam = person["bam"]
    if not bam or not os.path.exists(bam):
        return None, None, None, "bam_missing"
    import subprocess
    frows, prows, crows = [], [], []
    for q in MAPQS:
        try:
            dp = depth(bam, regions_bed, q, timeout, samtools, KIR_REGION)
            cm = control_mean(bam, q, controls, timeout, samtools, depth)
        except subprocess.TimeoutExpired:
            return None, None, None, "timeout"
        except subprocess.CalledProcessError:
            return None, None, None, "samtools_error"
        arrays = M70.gene_arrays(dp, regions)
        for name, sm, mean in M70.feature_sums(arrays, features):
            frows.append(dict(person_id=person["person_id"], feature=name, mapq=q, sum_depth=round(sm, 2), mean_depth=round(mean, 4)))
        for gene, st in M70.perbase_summary(arrays, features).items():
            prows.append(dict(person_id=person["person_id"], gene=gene, mapq=q,
                              **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in st.items()}))
        crows.append(dict(person_id=person["person_id"], mapq=q, control_depth=round(cm, 4)))
    return frows, prows, crows, None


CONTROL_COLS = ["person_id", "mapq", "control_depth"]


def done_people(internal_dir):
    done = set()
    for name in ("control_depth.tsv", "failures.tsv"):
        p = os.path.join(internal_dir, name)
        if os.path.exists(p):
            done |= set(pd.read_csv(p, sep="\t", dtype={"person_id": str}, comment="#")["person_id"])
    return done


def acquire_lock(internal_dir):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, "78.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit("FATAL: another 78_kir_gene_depth.py (pid %d) holds %s" % (pid, lock))
        except (ValueError, ProcessLookupError):
            log("removing stale lock " + lock)
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


def people_of_interest(int77c=I77C):
    """Trio members and pair members of 77c's person-level records."""
    ids = set()
    for fn, cols in (("records_trio_kir.INTERNAL.tsv", ("child", "tx", "other")), ("records_pair_kir.INTERNAL.tsv", ("person", "other"))):
        p = os.path.join(int77c, fn)
        if not os.path.exists(p):
            continue
        d = M70_read_int(p, usecols=list(cols), dtype=str)
        for c in cols:
            ids |= set(d[c].dropna())
    return ids


def M70_read_int(path, **kw):
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0, **kw)


# ---------------------------------------------------------------------------
# summarise: person x gene metrics (relative to the person's own control depth)
# ---------------------------------------------------------------------------
def gene_metrics(fd, perbase, ctrl, features):
    """fd: feature_depth rows (person_id, feature, mapq, sum_depth, mean_depth); perbase: cds_perbase; ctrl: control_depth.
    -> person_id, gene, depth_defined=True plus rel_gene, rel_cds (MAPQ>=20), cds_p10_rel, cds_min_rel, cds_share_lt10 (absolute, MAPQ>=20),
    mapq_frac (CDS depth at MAPQ>=20 over MAPQ>=0), rel_cds_q0 (MAPQ>=0 CDS depth over control), dev_cds_q0 = |log2(rel_cds_q0)|."""
    f = features.set_index("name")
    fd = fd.copy()
    fd["gene"] = fd["feature"].map(f["gene"])
    fd["kind"] = fd["feature"].map(f["kind"])
    c = ctrl.set_index(["person_id", "mapq"])["control_depth"]
    rows = {}

    def put(p, g, key, val):
        rows.setdefault((p, g), dict(person_id=p, gene=g))[key] = val

    for r in fd[fd["kind"].isin(["gene", "cds"])].itertuples():
        cd = c.get((r.person_id, int(r.mapq)), np.nan)
        rel = r.mean_depth / cd if cd and np.isfinite(cd) and cd > 0 else np.nan
        tag = "gene" if r.kind == "gene" else "cds"
        put(r.person_id, r.gene, "rel_%s_q%d" % (tag, int(r.mapq)), rel)
        put(r.person_id, r.gene, "abs_%s_q%d" % (tag, int(r.mapq)), r.mean_depth)
    for r in perbase.itertuples():
        cd = c.get((r.person_id, int(r.mapq)), np.nan)
        if int(r.mapq) == 20 and cd and np.isfinite(cd) and cd > 0:
            put(r.person_id, r.gene, "cds_p10_rel", r.cds_p10 / cd)
            put(r.person_id, r.gene, "cds_min_rel", r.cds_min / cd)
            put(r.person_id, r.gene, "cds_share_lt10", r.cds_share_lt10)
    out = pd.DataFrame(list(rows.values()))
    if out.empty:
        return out
    out["rel_gene"] = out.get("rel_gene_q20")
    out["rel_cds"] = out.get("rel_cds_q20")
    a0, a20 = out.get("abs_cds_q0"), out.get("abs_cds_q20")
    out["mapq_frac"] = np.where(a0 > 0, a20 / a0, np.nan)
    out["rel_cds_q0"] = out.get("rel_cds_q0")
    out["dev_cds_q0"] = np.abs(np.log2(out["rel_cds_q0"].where(out["rel_cds_q0"] > 0)))
    out["depth_defined"] = True
    keep = ["person_id", "gene", "depth_defined", "rel_gene", "rel_cds", "cds_p10_rel", "cds_min_rel", "cds_share_lt10", "mapq_frac",
            "rel_cds_q0", "dev_cds_q0"]
    return out[[k for k in keep if k in out.columns]]


def distribution_table(gm):
    """Aggregate quantiles per gene (public; people counts masked)."""
    rows = []
    for g, d in gm.groupby("gene"):
        for col in ("rel_gene", "rel_cds", "cds_p10_rel", "cds_min_rel", "cds_share_lt10", "mapq_frac", "rel_cds_q0"):
            v = d[col].to_numpy(float)
            v = v[np.isfinite(v)]
            if len(v) < 20:
                continue
            q = np.percentile(v, [5, 25, 50, 75, 95])
            rows.append(dict(gene=g, metric=col, n_people=C.mask_n(len(v)), p05=round(q[0], 4), p25=round(q[1], 4), p50=round(q[2], 4),
                             p75=round(q[3], 4), p95=round(q[4], 4)))
    return pd.DataFrame(rows)


def summarise(a):
    feats = pd.read_csv(os.path.join(a.out_dir, "features.tsv"), sep="\t")
    rd = lambda n: M70_read_int(os.path.join(a.internal_dir, n), dtype={"person_id": str})  # noqa: E731
    gm = gene_metrics(rd("feature_depth.tsv"), rd("cds_perbase.tsv"), rd("control_depth.tsv"), feats)
    path = os.path.join(a.internal_dir, "gene_depth_metrics.INTERNAL.tsv")
    _disc.assert_internal_path_allowed(path)
    with open(path, "w") as fh:
        fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
        gm.to_csv(fh, sep="\t", index=False)
    log("wrote %s (%d person-gene rows, %d people)" % (path, len(gm), gm["person_id"].nunique() if len(gm) else 0))
    C.write_tsv(distribution_table(gm), os.path.join(a.out_dir, "depth_distribution_by_gene.tsv"))
    cov = pd.DataFrame([dict(gene=g, depth_available="yes" if g in KIR_ON_CHR19 else "no (absent from GRCh38 primary chr19)")
                        for g in KIR_ON_CHR19 + NOT_ON_GRCH38])
    C.write_tsv(cov, os.path.join(a.out_dir, "depth_coverage_by_gene.tsv"))
    return 0


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-features", action="store_true")
    ap.add_argument("--gencode-gtf")
    ap.add_argument("--summarise", action="store_true")
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--int77c", default=I77C)
    ap.add_argument("--cohort", default=RD.DEFAULT_COHORT)
    ap.add_argument("--manifest", default=os.path.join(RD.AOU_ROOT, RD.LR_MANIFEST_REL))
    ap.add_argument("--control-regions", default=DEFAULT_CONTROL)
    ap.add_argument("--all", action="store_true", help="every long-read person, not only 77c's trio/pair members")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    feat_path = os.path.join(a.out_dir, "features.tsv")
    if a.build_features:
        if not a.gencode_gtf:
            raise SystemExit("--build-features needs --gencode-gtf")
        feats, missing = build_kir_features(a.gencode_gtf)
        feats.to_csv(feat_path, sep="\t", index=False)
        M70.write_bed(feats, os.path.join(a.out_dir, "features.bed"))
        M70.gene_regions(feats).to_csv(os.path.join(a.out_dir, "gene_regions.bed"), sep="\t", header=False, index=False)
        write_control_beds(a.out_dir, parse_regions(a.control_regions))
        log("features: %d rows, genes %s; not found on chr19 primary: %s; never available (not on GRCh38): %s"
            % (len(feats), sorted(feats["gene"].unique()), missing, NOT_ON_GRCH38))
        return 0
    if a.summarise:
        return summarise(a)

    feats = pd.read_csv(feat_path, sep="\t")
    regions = M70.gene_regions(feats)
    regions_bed = os.path.join(a.out_dir, "gene_regions.bed")
    controls = [(os.path.join(a.out_dir, "control_%d.bed" % i), "%s:%d-%d" % r) for i, r in enumerate(parse_regions(a.control_regions))]
    people = RD.load_people(a.cohort, a.manifest)
    if not a.all:
        want = people_of_interest(a.int77c)
        if not want:
            raise SystemExit("no 77c person-level records under %s (run 77c first, or pass --all)" % a.int77c)
        people = people[people["person_id"].isin(want)].reset_index(drop=True)
    lock = acquire_lock(a.internal_dir)
    try:
        done = done_people(a.internal_dir)
        todo = people[~people["person_id"].isin(done)]
        if a.limit:
            todo = todo.head(a.limit)
        log("people %d, done %d, to run %d" % (len(people), len(done), len(todo)))
        fpath, ppath = os.path.join(a.internal_dir, "feature_depth.tsv"), os.path.join(a.internal_dir, "cds_perbase.tsv")
        cpath, xpath = os.path.join(a.internal_dir, "control_depth.tsv"), os.path.join(a.internal_dir, "failures.tsv")
        t0, nfail = time.time(), 0
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(measure_person, p, regions_bed, regions, feats, controls, a.timeout, a.samtools): p["person_id"]
                    for p in todo.to_dict("records")}
            for k, fu in enumerate(as_completed(futs), 1):
                fr, pr, cr, why = fu.result()
                if why:
                    nfail += 1
                    M70.append_rows(xpath, [dict(person_id=futs[fu], reason=why)], ["person_id", "reason"])
                else:
                    M70.append_rows(fpath, fr, M70.FEATURE_DEPTH_COLS)
                    M70.append_rows(ppath, pr, M70.PERBASE_COLS)
                    M70.append_rows(cpath, cr, CONTROL_COLS)
                if k % 100 == 0 or k == len(futs):
                    log("%d/%d done, %d failed, %.1f s/person" % (k, len(futs), nfail, (time.time() - t0) / k))
        log("DONE 78")
    finally:
        try:
            os.remove(lock)
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
