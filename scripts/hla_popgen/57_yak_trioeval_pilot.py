#!/usr/bin/env python3
"""S05 57: trio k-mer phasing benchmark (yak trioeval / yak qv) on AoU long-read trios (VM; run from ~/s05).

Why: 56/56b measure typing and phasing error indirectly (relatives share haplotypes). A trio gives the
direct truth used by the HPRC/hifiasm papers: parental k-mers label every stretch of the child's
haplotype contigs as paternal or maternal -> switch error (adjacent hap-specific k-mers flip parent)
and hamming error (k-mers on the minority parent of their contig).

Trio = long-read child with two first-degree relatives (2^-2.5 <= kin <= 2^-1.5) that are absent from
the relatedness table as a pair (unrelated to each other). Any other configuration (sib + parent,
child + parent, two half-sibs) would put the pair in the table, so these are parent-parent-child.
Which parent is the father does not matter (switch/hamming are symmetric).

Region design (keeps memory/time small on a 4-vCPU app): parental and child k-mer databases are built
from the reads aligned to a padded DB window only; the child's hap1/hap2 contigs are trimmed to the
target window with the hifiasm->hg38 PAF and evaluated against those databases.
  MHC (chr6:29.5-33.5 Mb), class I A..B span, class II DR-DQ span, KIR (chr19 LRC), chr20 control.
Caveats: k-mers whose reads map outside the DB window are missing from a parent (paralogs) -> a small
upward bias, identical for all regions' platforms; region QV uses child reads from the DB window only.

Outputs (aggregate PUBLIC): ~/s05/results/57/trioeval_summary.tsv, trio_selection.tsv, run.log
INTERNAL (never leaves the VM): ~/s05/internal/57/per_trio_region.INTERNAL.tsv (checkpoint/resume),
raw yak output under ~/s05/work57/ (deleted per trio unless --keep).
"""
import argparse
import csv
import os
import random
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

import _s05_common as C

MOUNT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/")
LR_MANIFEST = MOUNT + "v9/wgs/long_read/manifest.tsv"
YAK = os.path.expanduser("~/envs/yak/bin/yak")
SAM = os.path.expanduser("~/envs/yak/bin/samtools")
OUT = os.path.expanduser("~/s05/results/57/")
INT = os.path.expanduser("~/s05/internal/57/")
WORK = os.path.expanduser("~/s05/work57/")
HIFI = ("revio", "sequel2e", "sequel2")
ASM = ("revio", "sequel2e")
# DB windows (reads -> k-mer DBs) and evaluation targets (child contig pieces) on GRCh38
DB = {"MHC": "chr6:28500000-34500000", "KIR": "chr19:54500000-55100000", "chr20": "chr20:1-64444167"}
TARGETS = [  # name, chrom, start, end, DB key
    ("MHC", "chr6", 29_500_000, 33_500_000, "MHC"),
    ("MHC_classI_A_to_B", "chr6", 29_900_000, 31_400_000, "MHC"),
    ("MHC_classII_DR_DQ", "chr6", 32_400_000, 32_700_000, "MHC"),
    ("KIR", "chr19", 54_700_000, 54_900_000, "KIR"),
    ("chr20_control", "chr20", 5_000_000, 60_000_000, "chr20"),
]
QV_CAP = 60.0
FIELDS = ["child", "child_platform", "target", "hap", "n_pieces", "bp",
          "sw_err", "sw_tot", "hm_err", "hm_tot", "qv"]


def sh(cmd):
    r = subprocess.run(cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[:160]} -> rc {r.returncode}: {r.stderr[-400:]}")
    return r.stdout


def mpath(uri):
    """gs://bucket/pooled/... -> mount path; '' if missing on the mount."""
    if not uri or uri == "NA":
        return ""
    p = MOUNT + uri.split("/", 3)[-1] if uri.startswith("gs://") else MOUNT + uri
    return p if os.path.exists(p) else ""


def manifest_index(people):
    """person -> {'bam','plat'} (any resolving HiFi row) and {'h1','h2','p1','p2','plat'} (assembly row).
    Only trio members are existence-checked (mount stats are slow)."""
    bam, asm = {}, {}
    with open(LR_MANIFEST) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pid, pl = r["research_id"], r["platform"].lower()
            if pl not in HIFI or pid not in people:
                continue
            if pid not in bam:
                b = mpath(r.get("grch38_bam", ""))
                if b and os.path.exists(b + ".bai"):
                    bam[pid] = {"bam": b, "plat": pl}
            if pid not in asm and pl in ASM:
                d = {k: mpath(r.get(c, "")) for k, c in (
                    ("h1", "assembly_hap1_fa"), ("h2", "assembly_hap2_fa"),
                    ("p1", "assembly_hap1_aln2_hg38_paf"), ("p2", "assembly_hap2_aln2_hg38_paf"))}
                if all(d.values()) and os.path.exists(d["h1"] + ".gzi"):
                    d["plat"] = pl
                    asm[pid] = d
    return bam, asm


def find_trios(lr_people):
    rel = pd.read_csv(C.DEFAULT_RELATEDNESS, sep="\t", dtype={"i.s": str, "j.s": str})
    first = rel[(rel["kin"] >= 2 ** -2.5) & (rel["kin"] <= 2 ** -1.5)]
    any_pair = set(zip(rel["i.s"], rel["j.s"])) | set(zip(rel["j.s"], rel["i.s"]))
    nb = defaultdict(set)
    for i, j in zip(first["i.s"], first["j.s"]):
        nb[i].add(j)
        nb[j].add(i)
    out = []
    for c in sorted(lr_people):
        ps = sorted(nb.get(c, ()))
        for a in range(len(ps)):
            for b in range(a + 1, len(ps)):
                if (ps[a], ps[b]) not in any_pair:
                    out.append((c, ps[a], ps[b]))
    return out


def region_pieces(paf, chrom, s, e, minlen=1000):
    """Child contig intervals aligned to chrom:s-e (primary PAF lines, tp:A:P), trimmed linearly to the
    target window, merged per contig. Returns [(contig, qs, qe)] 1-based inclusive for faidx."""
    iv = defaultdict(list)
    with open(paf) as f:
        for ln in f:
            x = ln.split("\t")
            if x[5] != chrom or "tp:A:P" not in ln:
                continue
            ql, qs, qe, strand, ts, te = int(x[1]), int(x[2]), int(x[3]), x[4], int(x[7]), int(x[8])
            a, b = max(ts, s), min(te, e)
            if b - a < minlen:
                continue
            f0, f1 = (a - ts) / (te - ts), (b - ts) / (te - ts)
            if strand == "+":
                q0, q1 = qs + f0 * (qe - qs), qs + f1 * (qe - qs)
            else:
                q0, q1 = qe - f1 * (qe - qs), qe - f0 * (qe - qs)
            iv[x[0]].append((max(0, int(q0)), min(ql, int(q1))))
    out = []
    for ctg, lst in iv.items():
        lst.sort()
        cur = list(lst[0])
        for a, b in lst[1:]:
            if a <= cur[1]:
                cur[1] = max(cur[1], b)
            else:
                out.append((ctg, cur[0] + 1, cur[1]))
                cur = [a, b]
        out.append((ctg, cur[0] + 1, cur[1]))
    return out


def yak_db(bam, region, path, threads):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return
    tmp, fq = path + ".tmp", path + ".fq"  # yak count with -b reads its input twice -> no pipe
    sh(f"set -o pipefail; {SAM} view -u -F 0x900 '{bam}' {region} | {SAM} fastq - 2>/dev/null > {fq} "
       f"&& {YAK} count -k31 -b30 -t{threads} -o {tmp} {fq} 2>/dev/null; rm -f {fq}")
    if os.path.getsize(tmp) < 100_000:
        raise RuntimeError(f"near-empty k-mer DB for {region}")
    os.replace(tmp, path)


def parse_trioeval(txt):
    sw = hm = (0, 0)
    for ln in txt.splitlines():
        x = ln.split("\t")
        if x[0] == "W" and len(x) >= 4:
            sw = (int(x[1]), int(x[2]))
        elif x[0] == "H" and len(x) >= 4:
            hm = (int(x[1]), int(x[2]))
    return sw, hm


def parse_qv(txt):
    q = np.nan
    for ln in txt.splitlines():
        x = ln.split("\t")
        if x[0] == "QV" and len(x) >= 3:
            q = float(x[-1])
    return q


def run_trio(trio, bam, asm, done, threads, keep):
    c, p, m = trio
    wd = WORK + c + "/"
    os.makedirs(wd, exist_ok=True)
    rows = []
    try:
        for h in ("1", "2"):  # local copy so samtools faidx can build its .fai/.gzi (mount copies are read-only)
            fa = wd + f"h{h}.fa.gz"
            if not os.path.exists(fa + ".fai"):
                sh(f"cp '{asm[c]['h' + h]}' {fa}.part && mv {fa}.part {fa} && {SAM} faidx {fa}")
        for key, region in DB.items():
            todo = [t for t in TARGETS if t[4] == key and (c, t[0]) not in done]
            if not todo:
                continue
            for who, tag in ((p, "pat"), (m, "mat"), (c, "child")):
                yak_db(bam[who]["bam"], region, wd + f"{tag}.{key}.yak", threads)
            for name, chrom, s, e, _ in todo:
                for h in ("1", "2"):
                    pcs = region_pieces(asm[c]["p" + h], chrom, s, e)
                    seq = wd + f"h{h}.{name}.fa"
                    bp = sum(b - a + 1 for _, a, b in pcs)
                    if not pcs:
                        rows.append([c, asm[c]["plat"], name, h, 0, 0, 0, 0, 0, 0, np.nan])
                        continue
                    with open(seq + ".reg", "w") as f:
                        f.write("".join(f"{ctg}:{a}-{b}\n" for ctg, a, b in pcs))
                    sh(f"{SAM} faidx {wd}h{h}.fa.gz -r {seq}.reg > {seq}")
                    te = sh(f"{YAK} trioeval -t{threads} {wd}pat.{key}.yak {wd}mat.{key}.yak {seq} 2>/dev/null")
                    qv = sh(f"{YAK} qv -t{threads} {wd}child.{key}.yak {seq} 2>/dev/null")
                    if keep:
                        open(seq + ".trioeval.txt", "w").write(te)
                    (sw, hm), q = parse_trioeval(te), parse_qv(qv)
                    rows.append([c, asm[c]["plat"], name, h, len(pcs), bp, sw[0], sw[1], hm[0], hm[1], q])
            if not keep:
                sh(f"rm -f {wd}*.{key}.yak")
        return rows, None
    except Exception as ex:  # one bad trio must not stop the run
        return rows, f"{c}: {ex}"
    finally:
        if not keep:
            sh(f"rm -rf {wd}")


def summarise(df):
    out = []
    df = df.copy()
    df["qv"] = pd.to_numeric(df["qv"], errors="coerce").clip(upper=QV_CAP)  # yak reports inf when no k-mer is missing
    for (tgt, plat), g in [((t, "all"), g) for t, g in df.groupby("target")] + list(df.groupby(["target", "child_platform"])):
        per = g.groupby("child").agg(sw_err=("sw_err", "sum"), sw_tot=("sw_tot", "sum"),
                                     hm_err=("hm_err", "sum"), hm_tot=("hm_tot", "sum"),
                                     qv=("qv", "median"), n_pieces=("n_pieces", "sum"), bp=("bp", "sum"))
        n = len(per)
        small = 0 < n < 20
        row = {"target": tgt, "child_platform": plat, "n_trios": C.mask_n(n)}
        for k in ("sw", "hm"):
            e, t = per[k + "_err"].sum(), per[k + "_tot"].sum()
            row[k + "_pooled_pct"] = round(100 * e / t, 3) if t and not small else ""
            r = (100 * per[k + "_err"] / per[k + "_tot"].replace(0, np.nan)).dropna()
            for qn, qq in (("p25", .25), ("median", .5), ("p75", .75)):
                row[f"{k}_{qn}_pct"] = "" if small or r.empty else round(float(r.quantile(qq)), 3)
            # no "% trios with zero error" column: with 25-50 trios it encodes counts of 1-19 people
        row["qv_median"] = "" if small else round(float(per["qv"].median()), 1)
        row["pieces_per_trio_median"] = "" if small else float(per["n_pieces"].median())
        row["mb_per_trio_median"] = "" if small else round(float(per["bp"].median()) / 1e6, 2)
        out.append(row)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-platform", type=int, default=20)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--summarise-only", action="store_true")
    a = ap.parse_args()
    if not os.path.isdir(MOUNT + "v9"):
        sys.exit("FATAL: AoU bucket mount missing")
    for d in (OUT, INT, WORK):
        os.makedirs(d, exist_ok=True)
    ck = INT + "per_trio_region.INTERNAL.tsv"
    if not a.summarise_only:
        cohort = C.load_cohort(C.DEFAULT_COHORT)
        lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
        trios = find_trios(lr)
        bam, asm = manifest_index({x for t in trios for x in t})
        ok = [t for t in trios if t[0] in asm and t[1] in bam and t[2] in bam]
        seen, uniq = set(), []
        for t in ok:  # one trio per child
            if t[0] not in seen:
                seen.add(t[0])
                uniq.append(t)
        rng = random.Random(57)
        pick = []
        for pl in ASM:
            cand = [t for t in uniq if asm[t[0]]["plat"] == pl]
            rng.shuffle(cand)
            pick += cand[:a.per_platform]
        sel = [{"stage": "candidate trios (kinship)", "n": C.mask_n(len(trios))},
               {"stage": "child assembly + both parents HiFi BAM on mount", "n": C.mask_n(len(ok))},
               {"stage": "distinct children", "n": C.mask_n(len(uniq))}]
        sel += [{"stage": f"selected, child {pl}", "n": C.mask_n(sum(asm[t[0]]["plat"] == pl for t in pick))} for pl in ASM]
        sel += [{"stage": f"available, child {pl}", "n": C.mask_n(sum(asm[t[0]]["plat"] == pl for t in uniq))} for pl in ASM]
        C.write_tsv(pd.DataFrame(sel), OUT + "trio_selection.tsv")
        done = set()
        if os.path.exists(ck):
            d0 = pd.read_csv(ck, sep="\t", dtype={"child": str})
            done = set(zip(d0["child"], d0["target"]))
        else:
            with open(ck, "w") as f:
                f.write("\t".join(FIELDS) + "\n")
        todo = [t for t in pick if any((t[0], x[0]) not in done for x in TARGETS)]
        C.log(f"trios selected {len(pick)}, to run {len(todo)}")
        with ThreadPoolExecutor(a.workers) as ex:
            futs = [ex.submit(run_trio, t, bam, asm, done, a.threads, a.keep) for t in todo]
            for k, fu in enumerate(futs, 1):
                rows, err = fu.result()
                with open(ck, "a") as f:
                    for r in rows:
                        f.write("\t".join(map(str, r)) + "\n")
                    f.flush()
                    os.fsync(f.fileno())
                C.log(f"trio {k}/{len(todo)} rows {len(rows)}" + (f" ERROR {err}" if err else ""))
    df = pd.read_csv(ck, sep="\t", dtype={"child": str}).drop_duplicates(["child", "target", "hap"], keep="last")
    C.write_tsv(summarise(df), OUT + "trioeval_summary.tsv")
    C.log("DONE 57")


if __name__ == "__main__":
    main()
