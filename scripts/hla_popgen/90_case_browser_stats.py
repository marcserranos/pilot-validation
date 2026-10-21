"""90 stats (VM side, aggregate output): how are Immuannot's novel calls made up, over EVERY call in the cohort?

One streaming pass over hla_calls_rich.tsv. For each call with a CDS alignment (cds_mut) it rebuilds the alignment from
the public IPD reference allele + Immuannot's own cs string (no BAM, no per-person file), classifies it with the same code
the viewer uses, and counts. Output is counts per gene x novelty class x edit-distance bin: no ids, no sequences.

  cohort_stats.json   rows: gene, novelty, bin, n, and how many are indel_only / snv_only / mixed, how many carry a
                      protein-altering substitution, are synonymous-only, have a frameshift, are homopolymer-or-tandem
                      indels only; plus n_partial (clipped alignments that cannot be rebuilt from cs alone).
  novel_snv_index.tsv (only with --index) row ids of novel calls with no indel, for the exhaustive case extraction.
"""
import csv
import gzip
import json
import os
import sys
from collections import defaultdict

BINS = [(1, 1, "1"), (2, 2, "2"), (3, 5, "3-5"), (6, 10, "6-10"), (11, 30, "11-30"), (31, 100, "31-100"), (101, 10 ** 9, ">100")]
HLA_GENES = ["A", "B", "C", "E", "F", "G", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "DRB3", "DRB4", "DRB5"]
FIELDS = ["n", "indel_only", "snv_only", "mixed", "protein_snv", "syn_snv_only", "frameshift", "hp_only_indel", "n_partial"]


def bin_of(nm):
    for lo, hi, label in BINS:
        if lo <= nm <= hi:
            return label
    return None


def log(m):
    print(m, file=sys.stderr, flush=True)


def read_fasta_gz(path):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def classify(core, ref, cds_mut):
    """Stats dict for one call, or None when it cannot be rebuilt from the cs alone."""
    first = cds_mut.split(",")[0].split("|")
    if len(first) < 2:
        return None
    rows = core.replay_from_ref(ref, first[1])
    if rows is None:
        return None
    return core.analyse(rows[0], rows[1])


def run(a, core):
    only_novel = not a.all_calls
    refs, cells, idx = {}, defaultdict(lambda: dict.fromkeys(FIELDS, 0)), []
    total = skipped = 0
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if g not in HLA_GENES or not r["cds_mut"]:
                continue
            nov = core.novelty_class(r["consensus"])
            if only_novel and nov == "known":
                continue
            allele = r["cds_mut"].split(",")[0].split("|")[0]
            gname = allele.split("*")[0]
            if gname not in refs:
                p = os.path.join(os.path.expanduser(a.refdata), "CDSseq", gname + ".fa.gz")
                refs[gname] = read_fasta_gz(p) if os.path.exists(p) else {}
            ref = refs[gname].get(allele)
            total += 1
            res = classify(core, ref, r["cds_mut"]) if ref else None
            try:
                nm = int(float(r["cds_distance"]))
            except ValueError:
                nm = None
            cell = cells[(g, nov, bin_of(nm) if nm else "?")]
            cell["n"] += 1
            if res is None:
                cell["n_partial"] += 1
                skipped += 1
                continue
            st, ev = res["stats"], res["events"]
            k = core.edit_class(st)
            if k == "indel-only":
                cell["indel_only"] += 1
            elif k == "snv-only":
                cell["snv_only"] += 1
            elif k == "mixed":
                cell["mixed"] += 1
            cell["protein_snv"] += (st["missense"] + st["nonsense"]) > 0
            cell["syn_snv_only"] += (k == "snv-only" and st["missense"] + st["nonsense"] == 0)
            cell["frameshift"] += st["frameshifts"] > 0
            ind = [e for e in ev if e["type"] != "snv"]
            cell["hp_only_indel"] += bool(ind) and k == "indel-only" and all(e.get("homopolymer") or e.get("tandem") for e in ind)
            if a.index and k in ("snv-only",):
                idx.append((r["person_id"], r["hap"], g, r["copy_index"], nov, nm))
    rows = [dict(zip(("gene", "novelty", "bin"), key), **v) for key, v in sorted(cells.items(), key=lambda x: str(x[0]))]
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    meta = {"scope": "all calls with a CDS alignment" if a.all_calls else "every novel call (Immuannot 'new') with a CDS alignment",
            "unit": "calls (haplotype x gene), not people", "n_calls": total, "n_partial": skipped,
            "note": "aggregate counts, internal: apply the <20 rule before any publication"}
    with open(os.path.join(out, "cohort_stats.json"), "w") as f:
        json.dump({"meta": meta, "rows": rows}, f, separators=(",", ":"))
    if a.index:
        p = os.path.join(out, "novel_snv_index.private.tsv")
        with open(p, "w") as f:
            f.write("person\thap\tgene\tcopy\tnovelty\tnm\n")
            f.writelines("\t".join(map(str, x)) + "\n" for x in idx)
        os.chmod(p, 0o600)
    log("calls considered: %d, not rebuildable from cs alone: %d, cells: %d" % (total, skipped, len(rows)))
    return rows


def register(sub):
    p = sub.add_parser("stats", help="VM: cohort-wide composition of novel calls -> <out-dir>/cohort_stats.json (aggregate)")
    p.add_argument("--out-dir", default="~/case_browser")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--refdata", default="~/tools/Immuannot_refdata")
    p.add_argument("--all-calls", action="store_true", help="include known-allele calls too (slower)")
    p.add_argument("--index", action="store_true", help="also write the private index of SNV-only novel calls")
    p.set_defaults(func=lambda a: run(a, _CORE))


_CORE = None  # set by 90_case_browser.main()
