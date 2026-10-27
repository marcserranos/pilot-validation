"""92 novel burden by sequencing platform and ancestry (VM side, aggregate output).

Skeptic test for 'novel' calls: a real new allele follows population structure, an assembly/consensus error follows the sequencing
platform or read depth. One pass over hla_calls_rich.tsv; per person counts of novel calls, split into indel-bearing and
substitution-only (replayed from the public IPD CDS + Immuannot cs, as 91), per classical-gene calls, joined to the lrWGS manifest
platform and to ancestry_pred. Output: per (platform, ancestry) group: persons, mean novel calls per person, share of persons with
>= 1 substitution-only / indel-bearing novel call; counts < 20 persons suppressed. No person ids in the output.
"""
import argparse
import csv
import gzip
import importlib.util
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CLASSICAL = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
NONCLASSICAL = ["E", "F", "G"]


def log(m):
    print(m, file=sys.stderr, flush=True)


def load_91():
    spec = importlib.util.spec_from_file_location("r91", os.path.join(HERE, "91_novel_recurrence.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read_platforms(path):
    plat = defaultdict(set)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            plat[r["research_id"]].add(r.get("platform") or "NA")
    return {p: (next(iter(v)) if len(v) == 1 else "multi") for p, v in plat.items()}


def edit_kind(core, ref, cds_mut):
    first = cds_mut.split(",")[0].split("|")
    if ref is None or len(first) < 2:
        return "unknown"
    rows = core.replay_from_ref(ref, first[1])
    if rows is None:
        return "unknown"
    st = core.analyse(rows[0], rows[1])["stats"]
    ind = st["ins_events"] + st["del_events"]
    return "indel" if ind else ("snv" if st["snv"] else "none")


def run(a):
    r91 = load_91()
    core = r91.load_core()
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    plat = read_platforms(os.path.expanduser(a.manifest))
    refs = {}
    per = defaultdict(lambda: defaultdict(int))
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            p = r["person_id"]
            if p not in lr or g not in CLASSICAL + NONCLASSICAL:
                continue
            grp = "classical" if g in CLASSICAL else "nonclassical"
            per[p]["calls_" + grp] += 1
            nov = core.novelty_class(r["consensus"])
            if nov == "known":
                continue
            per[p]["novel_" + grp] += 1
            if nov == "novel-protein" and r["cds_mut"]:
                allele = r["cds_mut"].split(",")[0].split("|")[0]
                gname = allele.split("*")[0]
                if gname not in refs:
                    fp = os.path.join(os.path.expanduser(a.refdata), "CDSseq", gname + ".fa.gz")
                    refs[gname] = r91.read_fasta_gz(fp) if os.path.exists(fp) else {}
                k = edit_kind(core, refs[gname].get(allele), r["cds_mut"])
                per[p]["%s_%s_%s" % (k, "protein", grp)] += 1
    groups = defaultdict(list)
    for p in lr:
        groups[(plat.get(p, "NA"), anc.get(p, "NA"))].append(per.get(p, {}))
        groups[(plat.get(p, "NA"), "ALL")].append(per.get(p, {}))
        groups[("ALL", anc.get(p, "NA"))].append(per.get(p, {}))
    cols = ["calls_classical", "novel_classical", "snv_protein_classical", "indel_protein_classical",
            "calls_nonclassical", "novel_nonclassical", "snv_protein_nonclassical", "indel_protein_nonclassical"]
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "novel_burden_pub.tsv"), "w") as f:
        f.write("platform\tancestry\tpersons\t" + "\t".join("mean_" + c for c in cols) +
                "\tshare_persons_snv_classical\tshare_persons_indel_classical\n")
        for (pl, an), rows in sorted(groups.items()):
            n = len(rows)
            if n < 20:
                f.write("%s\t%s\t<20\n" % (pl, an))
                continue
            means = ["%.3f" % (sum(x.get(c, 0) for x in rows) / n) for c in cols]
            s1 = sum(1 for x in rows if x.get("snv_protein_classical", 0)) / n
            s2 = sum(1 for x in rows if x.get("indel_protein_classical", 0)) / n
            f.write("%s\t%s\t%d\t%s\t%.4f\t%.4f\n" % (pl, an, n, "\t".join(means), s1, s2))
    log("done: %d persons" % len(lr))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--manifest", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv")
    p.add_argument("--refdata", default="~/tools/Immuannot_refdata")
    p.add_argument("--out-dir", default="~/s08")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
