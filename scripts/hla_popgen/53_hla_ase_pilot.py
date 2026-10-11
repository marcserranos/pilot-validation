#!/usr/bin/env python3
"""S05 pilot: haplotype-resolved allele-specific expression (ASE) of HLA-A/B/C from each person's
own blood RNA-seq, aligned to a PERSONAL two-haplotype transcript reference built from their own
phased assembly. Question: do protein-identical alleles that differ in mRNA sequence (UTR or
synonymous) show allelic imbalance, i.e. is there expression-relevant variation invisible to
protein-level HLA typing?

Per person (unrelated, LR assembly + RNA-seq):
  1. mRNA of HLA-A/B/C on each haplotype = exons from hap{1,2}.gtf.gz (copy 1) spliced out of
     hap{1,2}.trimmed.fa (reverse-complemented on '-' strand). All other annotated MHC genes of both
     haps are added as decoys so cross-mapping reads have somewhere else to go.
  2. RNA reads over chr6:29.6-31.6 Mb pulled from the AoU STAR BAM (gcloud cp -> samtools view ->
     fastq; the local BAM is deleted immediately).
  3. minimap2 -ax sr against the personal reference; a read is INFORMATIVE for a gene if its
     primary alignment is to that gene's hap1 or hap2 mRNA with MAPQ >= 10 (the two copies differ,
     so the aligner can tell them apart).
  4. Mapping-bias control: 4,000 reads simulated equally from both haplotypes' mRNAs (150 bp,
     0.2% substitutions) are aligned the same way; their log2(h1/h2) is the person's bias.
Pair classes (per gene): mrna_identical (cannot be measured; counted only), protein_identical_
mrna_diff (same 2-field, different mRNA), protein_diff.

Outputs: INTERNAL per-person table stays on the VM (~/s05/internal/53/); the PUBLIC summary is
aggregate per gene x class with "<20" masking.

    python3 53_hla_ase_pilot.py --n-per-class 30 --workers 4 --out-dir ~/s05/results/53
"""
import argparse
import gzip
import os
import random
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import binomtest

import _s05_common as C

GENES = ["A", "B", "C"]
REGION = "chr6:29600000-31600000"
PEOPLE_ROOT = os.path.expanduser("~/pipeline_outputs/people")
DEFAULT_MINIMAP2 = os.path.expanduser("~/repos/pilot-validation/.pixi/envs/spechla/bin/minimap2")
DEFAULT_RNA_MANIFEST = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq/manifest.tsv")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def gtf_mrnas(gtf_path, contigs):
    """{gene_id: (gene_name, mRNA)} from exon rows; gene_id keeps the '.N' copy suffix."""
    ex = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] != "exon" or x[0] not in contigs:
                continue
            gid = re.search(r'gene_id "([^"]+)"', x[8])
            gname = re.search(r'gene_name "([^"]+)"', x[8])
            if not gid or not gname:
                continue
            ex.setdefault(gid.group(1), [C.norm_gene(gname.group(1)), x[6], []])[2].append(
                (x[0], int(x[3]), int(x[4])))
    out = {}
    for gid, (gname, strand, exons) in ex.items():
        exons.sort(key=lambda e: e[1])
        s = "".join(contigs[c][a - 1:b] for c, a, b in exons)
        if strand == "-":
            s = s.translate(COMP)[::-1]
        if len(s) >= 300:
            out[gid] = (gname, s.upper())
    return out


def is_copy1(gid):
    return not re.search(r"\.(\d+)$", gid) or gid.endswith(".1")


def simulate(seqs, n, rng, rl=150, err=0.002):
    out = []
    for name, s in seqs:
        for i in range(n):
            if len(s) <= rl:
                continue
            p = rng.randrange(0, len(s) - rl)
            r = list(s[p:p + rl])
            for j in range(rl):
                if rng.random() < err:
                    r[j] = rng.choice("ACGT")
            out.append((f"sim_{name}_{i}", "".join(r)))
    return out


def count_informative(sam_lines, targets):
    cnt = {t: 0 for t in targets}
    for line in sam_lines:
        if line.startswith("@"):
            continue
        x = line.split("\t", 6)
        flag, ref, mapq = int(x[1]), x[2], int(x[4])
        if flag & 0x904 or ref not in cnt or mapq < 10:
            continue
        cnt[ref] += 1
    return cnt


def run_person(task):
    pid, bam_gs, bai_gs, classes, mm2, threads, seed = task
    bp = ["--billing-project", os.environ.get("GOOGLE_CLOUD_PROJECT", "")]
    rng = random.Random(seed)
    base = os.path.join(PEOPLE_ROOT, pid, "immuannot_output")
    tmp = tempfile.mkdtemp(prefix=f"ase_{pid}_", dir=os.path.expanduser("~/s05/tmp"))
    try:
        refs, focal = [], {}
        for hap in ("hap1", "hap2"):
            contigs = read_fasta(os.path.join(base, f"{hap}.trimmed.fa"))
            for gid, (gname, s) in gtf_mrnas(os.path.join(base, f"{hap}.gtf.gz"), contigs).items():
                tag = f"{hap}|{gid}"
                refs.append((tag, s))
                if gname in GENES and is_copy1(gid):
                    focal.setdefault(gname, {})[hap] = (tag, s)
        with open(f"{tmp}/ref.fa", "w") as f:
            for t, s in refs:
                f.write(f">{t}\n{s}\n")
        bam = f"{tmp}/rna.bam"
        subprocess.run(["gcloud", "storage", "cp", "-q"] + bp + [bam_gs, bam], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["gcloud", "storage", "cp", "-q"] + bp + [bai_gs, bam + ".bai"], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if not os.path.exists(bam + ".bai"):
            subprocess.run(["samtools", "index", bam], check=True)
        reg = subprocess.run(["samtools", "view", "-b", "-F", "0x904", bam, REGION], check=True,
                             capture_output=True).stdout
        os.remove(bam)
        fq = subprocess.run(["samtools", "fastq", "-F", "0x980", "-"], input=reg, check=True,
                            capture_output=True).stdout
        with open(f"{tmp}/r.fq", "wb") as f:
            f.write(fq)
        sam = subprocess.run([mm2, "-ax", "sr", "-t", str(threads), "--secondary=no", f"{tmp}/ref.fa",
                              f"{tmp}/r.fq"], check=True, capture_output=True, text=True).stdout.splitlines()
        out = []
        for gene, d in focal.items():
            if "hap1" not in d or "hap2" not in d:
                continue
            (t1, s1), (t2, s2) = d["hap1"], d["hap2"]
            cls = classes.get(gene)
            if cls is None:
                continue
            if s1 == s2:
                out.append({"gene": gene, "pair_class": "mrna_identical", "n1": 0, "n2": 0,
                            "sim1": 0, "sim2": 0})
                continue
            cnt = count_informative(sam, [t1, t2])
            sim = simulate([(t1, s1), (t2, s2)], 2000, rng)
            with open(f"{tmp}/sim.fq", "w") as f:
                for nm, r in sim:
                    f.write(f"@{nm}\n{r}\n+\n{'I' * len(r)}\n")
            ssam = subprocess.run([mm2, "-ax", "sr", "-t", str(threads), "--secondary=no", f"{tmp}/ref.fa",
                                   f"{tmp}/sim.fq"], check=True, capture_output=True, text=True).stdout.splitlines()
            sc = count_informative(ssam, [t1, t2])
            out.append({"gene": gene, "pair_class": cls, "n1": cnt[t1], "n2": cnt[t2],
                        "sim1": sc[t1], "sim2": sc[t2], "mrna_len_diff": abs(len(s1) - len(s2))})
        for r in out:
            r["person_id"] = pid
        return out, None
    except Exception as e:  # noqa: BLE001 -- one failed person must not kill the pilot
        return [], f"{type(e).__name__}"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def pair_classes(t1, unrelated):
    """{pid: {gene: class}} for genes with exactly 2 copy-1 calls."""
    sub = t1[(t1["copy_index"] == 1) & t1["gene"].isin(GENES) & t1["person_id"].isin(unrelated)]
    res = {}
    for (pid, gene), g in sub.groupby(["person_id", "gene"]):
        if len(g) != 2:
            continue
        c1, c2 = list(g["consensus"])
        d1, d2 = list(g["novelty_depth"])
        p1, p2 = C.two_field(c1), C.two_field(c2)
        if p1 is None or p2 is None or d1 == 2 or d2 == 2:
            continue
        cls = "protein_diff" if p1 != p2 else (
            "protein_identical_mrna_diff" if (c1 != c2 or d1 >= 3 or d2 >= 3) else "identical_call")
        res.setdefault(pid, {})[gene] = cls
    return res


def summarize(df, mode):
    rows = []
    for (gene, cls), g in df.groupby(["gene", "pair_class"]):
        row = {"gene": gene, "pair_class": cls, "n_people": C.mask_n(len(g), mode)}
        if cls != "mrna_identical" and len(g):
            n = g["n1"] + g["n2"]
            ok = g[n >= 30]
            lo = np.log2((ok["n1"] + .5) / (ok["n2"] + .5))
            ls = np.log2((ok["sim1"] + .5) / (ok["sim2"] + .5))
            pv = [binomtest(int(a), int(a + b), 0.5).pvalue for a, b in zip(ok["n1"], ok["n2"])]
            sig = (np.array(pv) < 0.01) & (np.abs(lo - ls) > 0.5) if len(ok) else np.array([])
            row.update({"n_people_ge30_informative": C.mask_n(len(ok), mode),
                        "median_informative_reads": int(n.median()),
                        "median_abs_log2_obs": round(float(np.abs(lo).median()), 2) if len(ok) else "",
                        "median_abs_log2_sim_bias": round(float(np.abs(ls).median()), 2) if len(ok) else "",
                        "median_abs_log2_obs_minus_bias": round(float(np.abs(lo - ls).median()), 2) if len(ok) else "",
                        "pct_imbalanced_p01_gt_1.41x": C.rate_cells(int(sig.sum()), len(ok), mode, 1)["pct"] if len(ok) else ""})
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--rna-manifest", default=DEFAULT_RNA_MANIFEST)
    ap.add_argument("--minimap2", default=DEFAULT_MINIMAP2)
    ap.add_argument("--n-per-class", type=int, default=30)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s05/results/53"))
    ap.add_argument("--internal-dir", default=os.path.expanduser("~/s05/internal/53"))
    ap.add_argument("--disclosure", default="public", choices=["public", "internal"])
    ap.add_argument("--seed", type=int, default=20261006)
    args = ap.parse_args()
    os.makedirs(os.path.expanduser("~/s05/tmp"), exist_ok=True)
    t1 = C.load_table1(args.table1)
    unrelated = C.greedy_unrelated(set(t1["person_id"]), C.load_relatedness_pairs(args.relatedness_table))
    man = pd.read_csv(args.rna_manifest, sep="\t", dtype=str).drop_duplicates("research_id")
    bam_of = dict(zip(man["research_id"], man["markduplicates_bam_file_path"]))
    bai_of = dict(zip(man["research_id"], man["markduplicates_bam_index_path"]))
    pcls = {p: v for p, v in pair_classes(t1, unrelated).items() if p in bam_of}
    rng = random.Random(args.seed)
    focal = sorted(p for p, v in pcls.items() if "protein_identical_mrna_diff" in v.values())
    other = sorted(p for p, v in pcls.items() if p not in set(focal))
    pick = rng.sample(focal, min(args.n_per_class, len(focal))) + rng.sample(other, min(args.n_per_class, len(other)))
    C.log(f"eligible focal={len(focal)} other={len(other)}; running {len(pick)} people")
    tasks = [(p, bam_of[p], bai_of[p], pcls[p], args.minimap2, 1, args.seed + i) for i, p in enumerate(pick)]
    rows, errs = [], {}
    with ProcessPoolExecutor(args.workers) as ex:
        for i, (out, err) in enumerate(ex.map(run_person, tasks)):
            rows.extend(out)
            if err:
                errs[err] = errs.get(err, 0) + 1
            if i % 10 == 0:
                C.log(f"  done {i + 1}/{len(tasks)} errors={errs}")
    df = pd.DataFrame(rows)
    os.makedirs(args.internal_dir, exist_ok=True)
    df.to_csv(os.path.join(args.internal_dir, "ase_per_person.INTERNAL.tsv"), sep="\t", index=False)
    C.write_tsv(summarize(df, args.disclosure), os.path.join(os.path.expanduser(args.out_dir), "ase_by_gene_class.tsv"))
    C.log(f"errors: {errs}")
    C.log("DONE 53")


if __name__ == "__main__":
    main()
