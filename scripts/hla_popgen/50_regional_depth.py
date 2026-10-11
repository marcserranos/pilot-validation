#!/usr/bin/env python3
"""Regional long-read depth: MHC, KIR and TRB, relative to each person's own genome-wide depth.

Sprint S05, WS-A and WS-B. Supervisor ask (2026-10-06): flag people whose long-read depth over the
MHC (or KIR) is far below their genome-wide depth, e.g. ~8-10x over the MHC at ~17-20x overall,
viewed per region rather than per person. This script
measures it, and in the same pass reads depth over every TRB V/J/C gene, because a germline deletion
of a TRBV gene (the lead from Aleix's W1: TRBV4-3 never used by ~27% of people) shows up in the
GRCh38-aligned long reads as roughly half depth (one copy) or no depth (both copies) over that gene.

What it measures, per person (one pass over that person's GRCh38 long-read BAM):

  - `samtools bedcov -c` over a BED of regions, at two mapping-quality floors:
      MAPQ >= 0  (every aligned read; the honest total, inflates in near-identical paralogs),
      MAPQ >= 20 (uniquely placed reads; deflates in near-identical paralogs).
    bedcov returns the summed per-base depth and the read count for each region.
    mean_depth = sum_depth / region_length.
  - Regions (built by `--build-regions` from a GENCODE GTF, or read from `--regions-bed`):
      control windows   five 2 Mb autosomal windows away from known large structural
                        polymorphism; their median mean depth is the person's genome-wide proxy;
      MHC window        chr6:29,500,000-33,500,000, the exact window the HLA pipeline trimmed to;
      HLA genes         the 65 genes of SCHEMA.md where GENCODE has them on the primary assembly;
      KIR window        first to last KIR gene on chr19 in GENCODE, padded 10 kb;
      KIR genes         every KIR gene GENCODE places on chr19 of the primary assembly (the
                        reference haplotype only; KIR2DS*/3DS1 etc. that GRCh38 carries only on
                        ALT contigs are absent from a noalt BAM and so not measurable here);
      TRB window        first to last TRB gene on chr7, padded 10 kb;
      TRB genes         every TRBV/TRBJ/TRBC/TRBD gene GENCODE places on chr7.
  - rel_depth = mean_depth(region) / median(mean_depth(control windows)), computed at the same
    MAPQ floor. 1.0 means the region is covered like the rest of the genome. For a diploid gene,
    0.5 means one copy, 0 means none.

TRB copy-number call (per person, per TRB gene, MAPQ >= 0 by default):
    copies = 0 if rel < 0.25, 1 if 0.25 <= rel < 0.75, 2 if rel >= 0.75 (`--cn-cuts`).
    A gene is "assessable" only if its median MAPQ>=20 / MAPQ>=0 depth ratio across the cohort is
    >= 0.8 (`--min-unique-frac`). Below that, reads are not uniquely placed (a near-identical
    paralog such as TRBV6-2/6-3 or TRBV12-3/12-4), depth is shared between copies, and a
    per-gene copy number would be wrong; such genes are reported but flagged, never called.
    Thresholds are reported next to the full rel_depth histogram, so the reader can see whether
    the three copy-number modes really separate where the cuts are.

Outputs (`--out-dir`, PUBLIC, aggregate only, safe to commit):
    regions.bed                      the exact regions used (public GENCODE coordinates)
    run_summary.tsv                  people attempted / measured / failed, platform mix, timing
    depth_summary_by_region.tsv      per region class x MAPQ x ancestry x platform: quantiles of
                                     rel_depth and mean_depth, and how many people fall below
                                     0.5 / 0.6 / 0.7 rel_depth (counts masked <20 in PUBLIC)
    trb_gene_assessability.tsv       per TRB gene: median unique fraction, assessable flag
    trb_cn_by_ancestry.tsv           per assessable TRB gene x ancestry: genotype counts 0/1/2,
                                     deletion-allele frequency, HWE exact p
    trb_rel_depth_histogram.tsv      per TRB gene, 0.05-wide rel_depth bins (counts masked <20)
INTERNAL (`--internal-dir`, never in the repo, never pulled):
    per_person_depth.tsv             one row per person x region x MAPQ (person_id included)
    trb_cn_per_person.tsv            one row per person x TRB gene: rel_depth and copies
These two INTERNAL tables are what 50b (outcome join) and 51 (TRB germline vs usage) read.

Robustness (ENVIRONMENT.md quirks #22, #23, #26, #35):
    - hard check that the controlled-data automount and the manifest exist before anything else;
    - PID-file lock in the internal dir, so a second launch fails loudly instead of corrupting;
    - per-person results are appended + fsynced to the INTERNAL table as they finish, and a rerun
      skips people already present (resume);
    - every samtools call has a timeout; a person whose BAM is missing or times out is counted
      as failed with a reason, never silently dropped.

Usage (VM, from the repo; samtools lives in the spechla pixi env):
    cd ~/repos/pilot-validation
    pixi run -e spechla -- python3 scripts/hla_popgen/50_regional_depth.py --build-regions \\
        --gencode-gtf ~/ref/gencode.v46.basic.annotation.gtf.gz --out-dir ~/s05/results/50
    setsid nohup pixi run -e spechla -- python3 -u scripts/hla_popgen/50_regional_depth.py \\
        --regions-bed ~/s05/results/50/regions.bed --out-dir ~/s05/results/50 \\
        --internal-dir ~/s05/internal/50 --workers 24 > ~/s05/logs/50.log 2>&1 < /dev/null & disown
    Smoke test first: add --limit 20 --internal-dir ~/s05/internal/50_smoke (own path, quirk #22).

Tests: scripts/hla_popgen/tests/test_50_regional_depth.py (synthetic, no samtools needed).
"""
import argparse
import gzip
import json
import math
import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _disclosure as _disc  # noqa: E402

AOU_ROOT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
LR_MANIFEST_REL = "v9/wgs/long_read/manifest.tsv"
DEFAULT_COHORT = os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/s05/results/50")
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s05/internal/50")

ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
MAPQS = (0, 20)

# GRCh38, 0-based half-open BED. 2 Mb gene-dense-but-ordinary autosomal windows chosen away from
# the MHC, KIR, TRB, IG loci, centromeres, and the well-known large inversion polymorphisms
# (8p23.1, 17q21.31, 16p11.2). Their median is robust to one window being odd in one person.
CONTROL_WINDOWS = [
    ("chr2", 100_000_000, 102_000_000, "control_chr2"),
    ("chr3", 150_000_000, 152_000_000, "control_chr3"),
    ("chr5", 100_000_000, 102_000_000, "control_chr5"),
    ("chr10", 60_000_000, 62_000_000, "control_chr10"),
    ("chr12", 80_000_000, 82_000_000, "control_chr12"),
]
MHC_WINDOW = ("chr6", 29_500_000, 33_500_000, "window_MHC")
LOCUS_PAD = 10_000

# SCHEMA.md gene list (bare names as GENCODE spells them).
HLA_GENE_NAMES = {
    "HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1",
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DRA",
    "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8",
    "HLA-DRB9", "HLA-DPA2", "HLA-DPB2", "HLA-DQA2", "HLA-DQB2",
    "HLA-E", "HLA-F", "HLA-G",
    "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S", "HLA-T", "HLA-U", "HLA-V",
    "HLA-W", "HLA-Y", "MICA", "MICB", "TAP1", "TAP2", "C4A", "C4B", "HFE",
}
TRB_RE = re.compile(r"^TRB[VJCD]")
KIR_RE = re.compile(r"^KIR[23]D")


def log(msg):
    print(f"[50 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Regions
# ---------------------------------------------------------------------------
def _gtf_attr(attrs, key):
    m = re.search(key + r' "([^"]+)"', attrs)
    return m.group(1) if m else None


def build_regions_from_gtf(gtf_path):
    """Returns a DataFrame chrom,start,end,name,region_class (0-based half-open), from the gene
    rows of a GENCODE GTF. Only primary-assembly chromosomes chr6/chr7/chr19 are scanned."""
    opener = gzip.open if gtf_path.endswith(".gz") else open
    genes = []
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[0] not in ("chr6", "chr7", "chr19"):
                continue
            name = _gtf_attr(f[8], "gene_name")
            if not name:
                continue
            start0, end = int(f[3]) - 1, int(f[4])
            if f[0] == "chr6" and name in HLA_GENE_NAMES:
                genes.append((f[0], start0, end, name, "gene_HLA"))
            elif f[0] == "chr7" and TRB_RE.match(name):
                genes.append((f[0], start0, end, name, "gene_TRB"))
            elif f[0] == "chr19" and KIR_RE.match(name):
                genes.append((f[0], start0, end, name, "gene_KIR"))
    return assemble_regions(genes)


def assemble_regions(gene_rows):
    """gene_rows: iterable of (chrom, start0, end, name, region_class). Adds controls, the MHC
    window, and padded KIR/TRB locus windows spanning the genes found. Duplicate gene names
    (rare, e.g. a gene split in two GENCODE records) are merged to their outer span."""
    df = pd.DataFrame(list(gene_rows), columns=["chrom", "start", "end", "name", "region_class"])
    if not df.empty:
        df = (df.groupby(["chrom", "name", "region_class"], as_index=False)
                .agg(start=("start", "min"), end=("end", "max")))
    rows = [dict(chrom=c, start=s, end=e, name=n, region_class="control")
            for c, s, e, n in CONTROL_WINDOWS]
    c, s, e, n = MHC_WINDOW
    rows.append(dict(chrom=c, start=s, end=e, name=n, region_class="window_MHC"))
    for cls, label in (("gene_KIR", "window_KIR"), ("gene_TRB", "window_TRB")):
        sub = df[df["region_class"] == cls] if not df.empty else df
        if len(sub):
            rows.append(dict(chrom=sub["chrom"].iloc[0], start=max(0, int(sub["start"].min()) - LOCUS_PAD),
                             end=int(sub["end"].max()) + LOCUS_PAD, name=label, region_class=label))
    out = pd.concat([pd.DataFrame(rows), df[["chrom", "start", "end", "name", "region_class"]]],
                    ignore_index=True)
    if out["name"].duplicated().any():
        raise ValueError(f"duplicate region names: {sorted(out.loc[out['name'].duplicated(), 'name'])}")
    return out[["chrom", "start", "end", "name", "region_class"]]


def read_regions_bed(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     names=["chrom", "start", "end", "name", "region_class"])
    return df


def write_regions_bed(df, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    df.to_csv(path, sep="\t", header=False, index=False)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def uri_to_local(uri, root=AOU_ROOT):
    if not isinstance(uri, str) or not uri.startswith(BUCKET_PREFIX):
        return None
    return os.path.join(root, uri[len(BUCKET_PREFIX):])


def load_people(cohort_path, manifest_path, root=AOU_ROOT):
    """One row per long-read person with called HLA: person_id, ancestry, platform, bam path.
    When the manifest has several rows for one person, the row whose platform matches the
    cohort's platform wins, then the first row with an existing BAM."""
    coh = pd.read_csv(cohort_path, sep="\t", dtype={"person_id": str})
    if "in_lr" in coh.columns:
        coh = coh[coh["in_lr"].astype(str).str.lower().isin(["true", "1"])]
    coh["ancestry"] = coh.get("ancestry_pred", pd.Series(index=coh.index, dtype=str)).astype(str).str.upper()
    man = pd.read_csv(manifest_path, sep="\t", dtype=str)
    man = man.rename(columns={"research_id": "person_id"})
    man["bam"] = man["grch38_bam"].map(lambda u: uri_to_local(u, root))
    rows = []
    by_pid = {pid: g for pid, g in man.groupby("person_id")}
    for _, r in coh.iterrows():
        g = by_pid.get(r["person_id"])
        bam, platform = None, r.get("platform")
        if g is not None:
            pref = g[g["platform"].astype(str) == str(platform)] if "platform" in g else g
            for cand in list(pref["bam"]) + list(g["bam"]):
                if cand:
                    bam = cand
                    break
            if (platform is None or (isinstance(platform, float) and math.isnan(platform))) and len(g):
                platform = g["platform"].iloc[0]
        rows.append(dict(person_id=r["person_id"], ancestry=r["ancestry"],
                         platform=str(platform), bam=bam))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
def run_bedcov(bam, bed_path, mapq, timeout, samtools="samtools"):
    """Returns DataFrame name, sum_depth, n_reads. Raises on failure."""
    cmd = [samtools, "bedcov", "-c", "-Q", str(mapq), bed_path, bam]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=True).stdout
    return parse_bedcov(out)


def parse_bedcov(text):
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        # chrom start end name region_class sum_depth n_reads
        rows.append(dict(name=f[3], start=int(f[1]), end=int(f[2]),
                         sum_depth=float(f[-2]), n_reads=int(f[-1])))
    return pd.DataFrame(rows)


def measure_person(person, bed_path, timeout, samtools="samtools", runner=None):
    """Returns (rows, None) or (None, reason). rows: list of dicts with person_id, region,
    mapq, mean_depth, n_reads, rel_depth."""
    runner = runner or run_bedcov
    bam = person["bam"]
    if not bam or not os.path.exists(bam):
        return None, "bam_missing"
    out = []
    for q in MAPQS:
        try:
            df = runner(bam, bed_path, q, timeout, samtools)
        except subprocess.TimeoutExpired:
            return None, "timeout"
        except subprocess.CalledProcessError:
            return None, "samtools_error"
        df["mean_depth"] = df["sum_depth"] / (df["end"] - df["start"]).clip(lower=1)
        ctrl = df.loc[df["name"].str.startswith("control_"), "mean_depth"]
        base = float(np.median(ctrl)) if len(ctrl) else float("nan")
        if not base or not np.isfinite(base) or base <= 0:
            return None, "zero_control_depth"
        for r in df.itertuples():
            out.append(dict(person_id=person["person_id"], region=r.name, mapq=q,
                            mean_depth=round(r.mean_depth, 4), n_reads=r.n_reads,
                            control_depth=round(base, 4),
                            rel_depth=round(r.mean_depth / base, 4)))
    return out, None


# ---------------------------------------------------------------------------
# Checkpoint + lock
# ---------------------------------------------------------------------------
PER_PERSON_COLS = ["person_id", "region", "mapq", "mean_depth", "n_reads", "control_depth",
                   "rel_depth"]
FAIL_COLS = ["person_id", "reason"]


def acquire_lock(internal_dir):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, "50.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"FATAL: another 50_regional_depth.py (pid {pid}) holds {lock}")
        except (ValueError, ProcessLookupError):
            log(f"removing stale lock {lock}")
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


def done_people(internal_dir):
    done = set()
    for name in ("per_person_depth.tsv", "failures.tsv"):
        p = os.path.join(internal_dir, name)
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t", dtype={"person_id": str}, comment="#")
            done |= set(d["person_id"])
    return done


def append_rows(path, rows, cols):
    new = not os.path.exists(path)
    if new:
        _disc.assert_internal_path_allowed(path)
    with open(path, "a") as fh:
        if new:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
            fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


# ---------------------------------------------------------------------------
# Aggregation (PUBLIC)
# ---------------------------------------------------------------------------
def region_class_map(regions):
    return dict(zip(regions["name"], regions["region_class"]))


def depth_summary(per, people, regions, mode):
    """Quantiles of rel_depth per region (windows + every gene) x mapq x ancestry, plus counts
    below 0.5/0.6/0.7. Also an 'ALL' ancestry row."""
    cls = region_class_map(regions)
    d = per.merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    d["region_class"] = d["region"].map(cls)
    out = []
    groups = [("ancestry", a) for a in ANCESTRY_ORDER] + [("ancestry", "ALL")]
    groups += [("platform", p) for p in sorted(d["platform"].dropna().unique())]
    for (region, q), g in d.groupby(["region", "mapq"]):
        for kind, val in groups:
            sub = g if val == "ALL" else g[g[kind] == val]
            n = len(sub)
            if n == 0:
                continue
            rel = sub["rel_depth"].to_numpy()
            row = dict(region=region, region_class=cls.get(region), mapq=q, stratum_type=kind,
                       stratum=val, n_people=_disc.mask(n, _disc.PARTICIPANT, mode).display)
            if mode == _disc.PUBLIC and n < _disc.SUPPRESS_BELOW:
                out.append(row)
                continue
            for p in (1, 5, 25, 50, 75, 95):
                row[f"rel_p{p:02d}"] = round(float(np.percentile(rel, p)), 3)
            row["mean_depth_p50"] = round(float(np.median(sub["mean_depth"])), 2)
            for cut in (0.5, 0.6, 0.7):
                k = int((rel < cut).sum())
                row[f"n_rel_lt_{cut}"] = _disc.mask(k, _disc.PARTICIPANT, mode).display
            out.append(row)
    return pd.DataFrame(out)


def trb_assessability(per, regions, min_unique_frac):
    cls = region_class_map(regions)
    trb = [r for r, c in cls.items() if c == "gene_TRB"]
    p = per[per["region"].isin(trb)].pivot_table(index=["person_id", "region"], columns="mapq",
                                                 values="mean_depth").reset_index()
    if p.empty or 0 not in p.columns or 20 not in p.columns:
        return pd.DataFrame(columns=["gene", "median_unique_frac", "assessable"])
    p["frac"] = np.where(p[0] > 0, p[20] / p[0], np.nan)
    a = p.groupby("region")["frac"].median().reset_index()
    a.columns = ["gene", "median_unique_frac"]
    a["median_unique_frac"] = a["median_unique_frac"].round(3)
    a["assessable"] = a["median_unique_frac"] >= min_unique_frac
    return a.sort_values("gene").reset_index(drop=True)


def call_copies(rel, cuts=(0.25, 0.75)):
    lo, hi = cuts
    return np.where(rel < lo, 0, np.where(rel < hi, 1, 2))


def hwe_exact_p(n_aa, n_ab, n_bb):
    """Wigginton et al. 2005 exact HWE test (two-sided). n_ab heterozygotes."""
    n = n_aa + n_ab + n_bb
    if n == 0:
        return float("nan")
    n_a = 2 * n_aa + n_ab
    n_b = 2 * n - n_a
    if n_a > n_b:
        n_a, n_b = n_b, n_a
    probs = {}
    mid = int(n_a * n_b / (2 * n))
    if (n_a - mid) % 2:
        mid += 1
    probs[mid] = 1.0
    het = mid
    while het > 1:
        hom_r = (n_a - het) // 2
        hom_c = n - het - hom_r
        probs[het - 2] = probs[het] * het * (het - 1) / (4.0 * (hom_r + 1) * (hom_c + 1))
        het -= 2
    het = mid
    while het <= n_a - 2:
        hom_r = (n_a - het) // 2
        hom_c = n - het - hom_r
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2) * (het + 1))
        het += 2
    total = sum(probs.values())
    obs = probs.get(n_ab, 0.0)
    return min(1.0, sum(v for v in probs.values() if v <= obs * (1 + 1e-9)) / total)


def trb_tables(per, people, regions, assess, cuts, mapq, mode):
    """Returns (cn_by_ancestry, histogram, cn_per_person[INTERNAL])."""
    cls = region_class_map(regions)
    trb = [r for r, c in cls.items() if c == "gene_TRB"]
    d = per[(per["region"].isin(trb)) & (per["mapq"] == mapq)].copy()
    d = d.merge(people[["person_id", "ancestry"]], on="person_id", how="left")
    d["copies"] = call_copies(d["rel_depth"].to_numpy(), cuts)
    ok = set(assess.loc[assess["assessable"], "gene"])
    d["assessable"] = d["region"].isin(ok)

    cn_rows = []
    for gene, g in d[d["assessable"]].groupby("region"):
        for anc in ANCESTRY_ORDER + ["ALL"]:
            sub = g if anc == "ALL" else g[g["ancestry"] == anc]
            n = len(sub)
            if n == 0:
                continue
            c0, c1, c2 = (int((sub["copies"] == k).sum()) for k in (0, 1, 2))
            del_af = (2 * c0 + c1) / (2 * n)
            row = dict(gene=gene, ancestry=anc,
                       n_people=_disc.mask(n, _disc.PARTICIPANT, mode).display,
                       n_0copies=_disc.mask(c0, _disc.PARTICIPANT, mode).display,
                       n_1copy=_disc.mask(c1, _disc.PARTICIPANT, mode).display,
                       n_2copies=_disc.mask(c2, _disc.PARTICIPANT, mode).display)
            small = mode == _disc.PUBLIC and any(0 < k < _disc.SUPPRESS_BELOW for k in (n, c0, c1, c2))
            row["pct_0copies"] = "" if small else round(100 * c0 / n, 2)
            row["deletion_allele_freq"] = "" if small else round(del_af, 4)
            row["hwe_p"] = "" if small else float(f"{hwe_exact_p(c0, c1, c2):.3g}")
            cn_rows.append(row)

    hist_rows = []
    bins = np.round(np.arange(0, 1.55, 0.05), 2)
    for gene, g in d.groupby("region"):
        counts, edges = np.histogram(g["rel_depth"].clip(upper=1.5 - 1e-9), bins=bins)
        for k, lo in zip(counts, edges[:-1]):
            hist_rows.append(dict(gene=gene, bin_lo=float(lo), bin_hi=round(float(lo) + 0.05, 2),
                                  n_people=_disc.mask(int(k), _disc.PARTICIPANT, mode).display,
                                  assessable=gene in ok))
    per_person = d[["person_id", "region", "rel_depth", "copies", "assessable"]].rename(
        columns={"region": "gene"})
    return pd.DataFrame(cn_rows), pd.DataFrame(hist_rows), per_person


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-regions", action="store_true",
                    help="build regions.bed from --gencode-gtf into --out-dir, then exit")
    ap.add_argument("--gencode-gtf")
    ap.add_argument("--regions-bed")
    ap.add_argument("--cohort", default=DEFAULT_COHORT)
    ap.add_argument("--aou-root", default=AOU_ROOT)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--timeout", type=int, default=600, help="seconds per samtools call")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--cn-cuts", default="0.25,0.75")
    ap.add_argument("--cn-mapq", type=int, default=0, choices=list(MAPQS))
    ap.add_argument("--min-unique-frac", type=float, default=0.8)
    ap.add_argument("--aggregate-only", action="store_true",
                    help="skip measurement; rebuild PUBLIC tables from the INTERNAL table")
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args(argv)

    if args.build_regions:
        if not args.gencode_gtf:
            ap.error("--build-regions needs --gencode-gtf")
        regions = build_regions_from_gtf(args.gencode_gtf)
        path = os.path.join(args.out_dir, "regions.bed")
        write_regions_bed(regions, path)
        log(f"wrote {path}: " + json.dumps(regions["region_class"].value_counts().to_dict()))
        return 0

    if not args.regions_bed:
        ap.error("--regions-bed is required (build it first with --build-regions)")
    regions = read_regions_bed(args.regions_bed)
    cuts = tuple(float(x) for x in args.cn_cuts.split(","))

    manifest = os.path.join(args.aou_root, LR_MANIFEST_REL)
    if not args.aggregate_only and not os.path.exists(manifest):
        raise SystemExit(f"FATAL: {manifest} not found -- is the controlled-data automount up? "
                         f"(ENVIRONMENT.md quirks #26/#35). Refusing to start.")
    people = load_people(args.cohort, manifest, args.aou_root)
    if args.limit:
        people = people.head(args.limit)
    log(f"{len(people)} long-read people; platforms {people['platform'].value_counts().to_dict()}")

    per_path = os.path.join(args.internal_dir, "per_person_depth.tsv")
    fail_path = os.path.join(args.internal_dir, "failures.tsv")
    t0 = time.time()
    if not args.aggregate_only:
        lock = acquire_lock(args.internal_dir)
        try:
            done = done_people(args.internal_dir)
            todo = [r for r in people.to_dict("records") if r["person_id"] not in done]
            log(f"{len(done)} already done (resume), {len(todo)} to measure, {args.workers} workers")
            n_ok = n_fail = 0
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                futs = {ex.submit(measure_person, p, args.regions_bed, args.timeout, args.samtools): p
                        for p in todo}
                for i, fut in enumerate(as_completed(futs), 1):
                    p = futs[fut]
                    rows, reason = fut.result()
                    if rows:
                        append_rows(per_path, rows, PER_PERSON_COLS)
                        n_ok += 1
                    else:
                        append_rows(fail_path, [dict(person_id=p["person_id"], reason=reason)], FAIL_COLS)
                        n_fail += 1
                    if i % 200 == 0 or i == len(todo):
                        rate = i / max(1e-9, time.time() - t0)
                        log(f"{i}/{len(todo)} done ({n_ok} ok, {n_fail} failed), {rate:.2f} people/s")
        finally:
            os.remove(lock)

    per = pd.read_csv(per_path, sep="\t", dtype={"person_id": str}, comment="#")
    fails = (pd.read_csv(fail_path, sep="\t", dtype={"person_id": str}, comment="#")
             if os.path.exists(fail_path) else pd.DataFrame(columns=FAIL_COLS))
    per = per[per["person_id"].isin(people["person_id"])]
    mode = args.disclosure
    os.makedirs(args.out_dir, exist_ok=True)
    write_regions_bed(regions, os.path.join(args.out_dir, "regions.bed"))

    summ = depth_summary(per, people, regions, mode)
    _disc.write_tsv(summ, os.path.join(args.out_dir, "depth_summary_by_region.tsv"), mode)
    assess = trb_assessability(per, regions, args.min_unique_frac)
    _disc.write_tsv(assess, os.path.join(args.out_dir, "trb_gene_assessability.tsv"), mode)
    cn, hist, cn_pp = trb_tables(per, people, regions, assess, cuts, args.cn_mapq, mode)
    _disc.write_tsv(cn, os.path.join(args.out_dir, "trb_cn_by_ancestry.tsv"), mode)
    _disc.write_tsv(hist, os.path.join(args.out_dir, "trb_rel_depth_histogram.tsv"), mode)
    _disc.write_tsv(cn_pp, os.path.join(args.internal_dir, "trb_cn_per_person.tsv"), _disc.INTERNAL)

    measured = per["person_id"].nunique()
    fail_reasons = fails[fails["person_id"].isin(people["person_id"])]["reason"].value_counts().to_dict()
    run = [dict(metric="people_in_cohort", value=len(people)),
           dict(metric="people_measured", value=measured),
           dict(metric="people_failed", value=int(sum(fail_reasons.values())))]
    run += [dict(metric=f"failed_{k}", value=_disc.mask(v, _disc.QC_TALLY, mode).display)
            for k, v in sorted(fail_reasons.items())]
    run += [dict(metric="regions", value=len(regions)),
            dict(metric="trb_genes_assessable", value=int(assess["assessable"].sum()) if len(assess) else 0),
            dict(metric="cn_cuts", value=args.cn_cuts), dict(metric="cn_mapq", value=args.cn_mapq),
            dict(metric="min_unique_frac", value=args.min_unique_frac),
            dict(metric="wall_seconds_this_invocation", value=int(time.time() - t0))]
    _disc.write_tsv(pd.DataFrame(run), os.path.join(args.out_dir, "run_summary.tsv"), mode)
    log(f"done: {measured} people measured; tables in {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
