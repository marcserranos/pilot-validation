"""96 RNA-seq allele-specific expression of recurrent HLA-F signatures (VM side; think-tank experiment 'independent technology').

For carriers of a signature who also have blood RNA-seq, reuse 53's machinery: a PERSONAL two-haplotype transcript reference built from the
carrier's own phased assembly (so a variant read is not penalised against GRCh38/IPD), RNA reads over chr6:29.6-31.6 Mb from the AoU STAR BAM,
minimap2 -ax sr, and a read is informative for HLA-F if it maps to F of hap1 or hap2 with MAPQ >= 10. A 4,000-read simulated mapping-bias control
is aligned the same way. Variant-bearing haplotype = the haplotype 91 recorded for the person.
Readouts per signature (aggregate): carriers analysed (>= --min-reads informative reads), pooled variant-haplotype read fraction, per-person
bias-corrected log2(variant/other) median, and how many are imbalanced (binomial p < 0.01 and |log2 obs - log2 sim| > 0.5).
Expectation if the allele is real and normally expressed: variant fraction ~ 0.5 after bias correction. Variant reads absent (fraction ~ 0) while
the other haplotype is well covered = the variant is not in the transcript (artefact or silenced).

Private per-person results are appended to <out>/ase96.private.jsonl (mode 600, resumable); the public table masks counts < 20.
"""
import argparse
import csv
import importlib.util
import json
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def log(m):
    print(m, file=sys.stderr, flush=True)


def load53():
    spec = importlib.util.spec_from_file_location("a53", os.path.join(HERE, "53_hla_ase_pilot.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.GENES = ["F"]
    return m


def summarise(rows, supp, show_reads):
    """rows: dicts with sig, var_n, oth_n, var_sim, oth_sim, min_reads -> text."""
    from math import log2
    from scipy.stats import binomtest
    by = defaultdict(list)
    for r in rows:
        by[r["sig"]].append(r)
    out = []
    for sig, rs in sorted(by.items()):
        ok = [r for r in rs if r["var_n"] + r["oth_n"] >= r["min_reads"]]
        out.append("== %s: carriers run %s, with >= min informative reads %s" % (sig, supp(len(rs)), supp(len(ok))))
        if not ok:
            out.append("")
            continue
        tv, to = sum(r["var_n"] for r in ok), sum(r["oth_n"] for r in ok)
        out.append("pooled variant-haplotype read fraction %.3f%s" % (tv / max(1, tv + to), " (%d variant + %d other reads)" % (tv, to) if show_reads else ""))
        lo = [log2((r["var_n"] + .5) / (r["oth_n"] + .5)) - log2((r["var_sim"] + .5) / (r["oth_sim"] + .5)) for r in ok]
        med = sorted(lo)[len(lo) // 2]
        out.append("bias-corrected median log2(variant/other) %.2f" % med)
        sig_n = sum(1 for r, l in zip(ok, lo) if binomtest(r["var_n"], r["var_n"] + r["oth_n"], 0.5).pvalue < 0.01 and abs(l) > 0.5)
        zero = sum(1 for r in ok if r["var_n"] == 0)
        out.append("imbalanced carriers %s; carriers with zero variant reads %s" % (supp(sig_n), supp(zero)))
        out.append("")
    return "\n".join(out) + "\n"


def run(a):
    a53 = load53()
    spec = importlib.util.spec_from_file_location("r91", os.path.join(HERE, "91_novel_recurrence.py"))
    r91 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(r91)
    os.environ.setdefault("GOOGLE_CLOUD_PROJECT", a.project)
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    os.makedirs(os.path.expanduser("~/s05/tmp"), exist_ok=True)
    src = os.path.expanduser(a.tables)
    raw = list(csv.DictReader(open(os.path.join(src, "recurrence_raw.private.tsv")), delimiter="\t"))
    carriers = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(src, "carriers.private.tsv")), delimiter="\t"):
        carriers[(r["gene"], r["ref_allele"], r["cds_events"])][r["person_id"]] = r["haps"].split(",")
    man = {}
    with open(os.path.expanduser(a.rna_manifest)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            man.setdefault(r["research_id"], (r["markduplicates_bam_file_path"], r["markduplicates_bam_index_path"]))
    jp = os.path.join(out, "ase96.private.jsonl")
    done = set()
    if os.path.exists(jp):
        done = {(d["sig"], d["person_id"]) for d in map(json.loads, open(jp))}
    rnd = random.Random(11)
    jobs = []
    for prot in a.signatures.split(","):
        cand = [r for r in raw if r["gene"] == "F" and r["protein"] == prot]
        if not cand:
            continue
        r = max(cand, key=lambda x: int(x["persons"]))
        byp = carriers[(r["gene"], r["ref_allele"], r["cds_events"])]
        ps = sorted(p for p, h in byp.items() if p in man and len(h) == 1)   # heterozygous carriers with RNA-seq
        if len(ps) > a.max_per_signature:
            ps = sorted(rnd.sample(ps, a.max_per_signature))
        log("%s: heterozygous carriers with RNA-seq selected: %d" % (prot, len(ps)))
        jobs += [(prot, p, byp[p][0]) for p in ps if (prot, p) not in done]
    for i, (prot, p, var_hap) in enumerate(jobs):
        bam, bai = man[p]
        rows, err = a53.run_person((p, bam, bai, {"F": "x"}, os.path.expanduser(a.minimap2), 1, 7 + i))
        rec = {"sig": prot, "person_id": p, "min_reads": a.min_reads, "err": err}
        for r in rows:
            if r["pair_class"] == "mrna_identical":
                rec["err"] = "identical_mrna"
                continue
            v, o = ("1", "2") if var_hap == "hap1" else ("2", "1")
            rec.update({"var_n": r["n" + v], "oth_n": r["n" + o], "var_sim": r["sim" + v], "oth_sim": r["sim" + o]})
        with open(jp, "a") as f:
            f.write(json.dumps(rec) + "\n")
        os.chmod(jp, 0o600)
        log("%d/%d %s done err=%s" % (i + 1, len(jobs), prot, rec["err"]))
    res = [d for d in map(json.loads, open(jp)) if "var_n" in d]
    with open(os.path.join(out, "ase96_pub.txt"), "w") as f:
        f.write(summarise(res, r91.supp, False))
    rawp = os.path.join(out, "ase96_raw.private.txt")
    with open(rawp, "w") as f:
        f.write(summarise(res, str, True))
    os.chmod(rawp, 0o600)
    log("wrote ase96 summaries to %s (%d analysable carrier-runs)" % (out, len(res)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--tables", default="~/s08/r366")
    p.add_argument("--out-dir", default="~/s08/r366")
    p.add_argument("--rna-manifest", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq/manifest.tsv")
    p.add_argument("--minimap2", default="~/repos/pilot-validation-main/.pixi/envs/spechla/bin/minimap2")
    p.add_argument("--project", default="wb-cordial-leechee-9743", help="requester-pays billing project if GOOGLE_CLOUD_PROJECT is unset")
    p.add_argument("--signatures", default="T237I,P205T,Y106C")
    p.add_argument("--max-per-signature", type=int, default=13)
    p.add_argument("--min-reads", type=int, default=20)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
