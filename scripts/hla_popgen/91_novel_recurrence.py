"""91 novel-call recurrence (VM side, aggregate output): which 'novel' substitution-only HLA calls recur, in whom, and are they inherited?

One streaming pass over hla_calls_rich.tsv. Every novel call (Immuannot 'new') whose CDS alignment is rebuildable from the public IPD
reference allele + Immuannot's cs string, and that has NO indel, is keyed by (gene, reference allele, set of CDS substitutions). Per key:
  calls / persons, persons left after greedy removal of relatives (>= 2nd degree, kin >= 2^-3.5), ancestry of the unrelated carriers
  (lowercase ancestry_pred upper-cased, ancestry denominators included), homozygous persons (both haplotypes carry it),
  background consensus alleles, and transmission: carriers with two LR-cohort parents (first-degree pair with no kinship between
  them, as 90 extract's load_trios) split into 'a parent also carries it' / 'neither parent carries it'.
A real germline allele is transmitted; a calling artefact appears in the child alone. Parents' calls are only 'seen' when they have
the same key among novel calls: a parent whose call was typed as the known allele counts as non-carrier (conservative toward 'not transmitted').

Outputs (--out-dir, default ~/s08):
  recurrence_pub.tsv           counts <20 written '<20' (safe for the report); signatures with < --min-persons persons dropped
  recurrence_raw.private.tsv   same with exact counts (mode 600, stays on the VM)
  carriers.private.tsv         key -> person ids (mode 600, VM only; for later phenotype/EHR work on the VM)
Checkpointing is not needed: one pass, ~minutes, no partial output.
"""
import argparse
import csv
import gzip
import importlib.util
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_GENES = ["A", "B", "C", "E", "F", "G", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "DRB3", "DRB4", "DRB5"]
FIRST_LO, FIRST_HI = 2 ** -2.5, 2 ** -1.5
SECOND_MIN = 2 ** -3.5


def log(m):
    print(m, file=sys.stderr, flush=True)


def supp(n, floor=20):
    return "<%d" % floor if 0 < n < floor else str(n)


def load_core():
    spec = importlib.util.spec_from_file_location("cb90", os.path.join(HERE, "90_case_browser.py"))
    m = importlib.util.module_from_spec(spec)
    sys.modules["cb90"] = m
    spec.loader.exec_module(m)
    return m


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


def greedy_unrelated(persons, kin_pairs):
    """Largest-ish subset of `persons` with no pair in kin_pairs (dict person -> set of relatives). Deterministic."""
    keep, banned = [], set()
    for p in sorted(persons, key=lambda x: (len(kin_pairs.get(x, ())), x)):
        if p in banned:
            continue
        keep.append(p)
        banned.add(p)
        banned.update(kin_pairs.get(p, ()))
    return keep


def read_membership(path):
    anc, lr = {}, set()
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("in_lr", "").lower() in ("true", "1"):
                lr.add(r["person_id"])
                anc[r["person_id"]] = (r.get("ancestry_pred") or "NA").upper()
    return anc, lr


def read_relatives(path, universe):
    """(first_degree, any_kin >= 2nd degree) adjacency, restricted to `universe` persons."""
    first, anyk = defaultdict(set), defaultdict(set)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            i, j = r["i.s"], r["j.s"]
            if i not in universe or j not in universe:
                continue
            kin = float(r["kin"])
            if kin >= SECOND_MIN:
                anyk[i].add(j)
                anyk[j].add(i)
            if FIRST_LO <= kin <= FIRST_HI:
                first[i].add(j)
                first[j].add(i)
    return first, anyk


def parents_of(c, first, anyk):
    ps = sorted(first.get(c, ()))
    for a in range(len(ps)):
        for b in range(a + 1, len(ps)):
            if ps[b] not in anyk.get(ps[a], ()):
                return (ps[a], ps[b])
    return None


def signature(core, ref, cds_mut):
    """(ref_allele, tuple of 'refpos:ref>obs'), or None when not rebuildable or any indel is present."""
    first = cds_mut.split(",")[0].split("|")
    if len(first) < 2:
        return None
    rows = core.replay_from_ref(ref, first[1])
    if rows is None:
        return None
    res = core.analyse(rows[0], rows[1])
    if res["stats"]["ins_events"] or res["stats"]["del_events"] or not res["stats"]["snv"]:
        return None
    prot = []
    for e in res["events"]:
        if e["type"] == "snv" and e.get("effect") in ("missense", "nonsense"):
            prot.append("%s%d%s" % (e["aa_ref"], e["codon"], e["aa_obs"]))
    ev = tuple(sorted("%d:%s>%s" % (e["refpos"], e["ref"], e["obs"]) for e in res["events"] if e["type"] == "snv"))
    return first[0], ev, ";".join(dict.fromkeys(prot)) or "synonymous"


def run(a):
    core = load_core()
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    anc, lr = read_membership(os.path.expanduser(a.membership))
    log("LR persons: %d" % len(lr))
    refs = {}
    calls = defaultdict(lambda: defaultdict(list))   # key -> person -> [(hap, consensus)]
    prot_of = {}
    n_rows = n_sig = 0
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if g not in HLA_GENES or not r["cds_mut"] or core.novelty_class(r["consensus"]) == "known":
                continue
            n_rows += 1
            allele = r["cds_mut"].split(",")[0].split("|")[0]
            gname = allele.split("*")[0]
            if gname not in refs:
                p = os.path.join(os.path.expanduser(a.refdata), "CDSseq", gname + ".fa.gz")
                refs[gname] = read_fasta_gz(p) if os.path.exists(p) else {}
            ref = refs[gname].get(allele)
            sig = signature(core, ref, r["cds_mut"]) if ref else None
            if sig is None:
                continue
            key = (g, sig[0], ";".join(sig[1]))
            prot_of[key] = sig[2]
            calls[key][r["person_id"]].append((r["hap"], r["consensus"]))
            n_sig += 1
    log("novel rows %d, substitution-only rebuilt %d, distinct keys %d" % (n_rows, n_sig, len(calls)))

    cand = {k: v for k, v in calls.items() if len(v) >= a.min_persons}
    carriers_all = set().union(*[set(v) for v in cand.values()]) if cand else set()
    first, anyk = read_relatives(os.path.expanduser(a.relatedness), lr)
    log("relatives loaded for %d persons" % len(anyk))
    anc_den = Counter(anc.values())
    unr_all = greedy_unrelated(lr, anyk)
    anc_den_unrel = Counter(anc[p] for p in unr_all)
    trio = {c: parents_of(c, first, anyk) for c in carriers_all}

    rows = []
    for key, byp in cand.items():
        persons = sorted(byp)
        unrel = greedy_unrelated(persons, anyk)
        ac = Counter(anc.get(p, "NA") for p in unrel)
        hom = sum(1 for p in persons if len({h for h, _ in byp[p]}) > 1)
        bg = Counter(c for p in persons for _, c in byp[p])
        tr_with = tr_without = tr_denovo_only = 0
        for p in persons:
            pr = trio.get(p)
            if not pr:
                continue
            if any(q in byp for q in pr):
                tr_with += 1
            else:
                tr_without += 1
        rows.append({"gene": key[0], "ref_allele": key[1], "protein": prot_of[key], "cds_events": key[2],
                     "calls": sum(len(v) for v in byp.values()), "persons": len(persons), "unrelated": len(unrel),
                     "homozygous_persons": hom, "n_background_alleles": len(bg),
                     "top_background": ";".join("%s=%d" % kv for kv in bg.most_common(3)),
                     "anc_unrelated": ";".join("%s=%d" % kv for kv in sorted(ac.items())),
                     "carrier_children_with_2_parents": tr_with + tr_without,
                     "child_parent_also_carries": tr_with, "child_neither_parent_carries": tr_without})
    rows.sort(key=lambda r: (-r["unrelated"], r["gene"]))
    cols = list(rows[0]) if rows else []
    rawp = os.path.join(out, "recurrence_raw.private.tsv")
    with open(rawp, "w") as f:
        f.write("\t".join(cols) + "\n")
        f.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    os.chmod(rawp, 0o600)
    numeric = {"calls", "persons", "unrelated", "homozygous_persons", "carrier_children_with_2_parents",
               "child_parent_also_carries", "child_neither_parent_carries"}

    def pub_cell(c, v):
        if c in numeric:
            return supp(v)
        if c in ("anc_unrelated", "top_background"):
            return ";".join(p.split("=")[0] + "=" + supp(int(p.split("=")[1])) for p in v.split(";") if "=" in p)
        return str(v)
    with open(os.path.join(out, "recurrence_pub.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        f.writelines("\t".join(pub_cell(c, r[c]) for c in cols) + "\n" for r in rows)
    cp = os.path.join(out, "carriers.private.tsv")
    with open(cp, "w") as f:
        f.write("gene\tref_allele\tcds_events\tperson_id\thaps\n")
        for key, byp in cand.items():
            for p, v in byp.items():
                f.write("%s\t%s\t%s\t%s\t%s\n" % (key[0], key[1], key[2], p, ",".join(sorted({h for h, _ in v}))))
    os.chmod(cp, 0o600)
    with open(os.path.join(out, "ancestry_denominators.tsv"), "w") as f:
        f.write("ancestry\tlr_persons\tlr_unrelated\n")
        for k in sorted(anc_den):
            f.write("%s\t%s\t%s\n" % (k, supp(anc_den[k]), supp(anc_den_unrel.get(k, 0))))
    log("wrote %d signature rows to %s" % (len(rows), out))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
    p.add_argument("--refdata", default="~/tools/Immuannot_refdata")
    p.add_argument("--out-dir", default="~/s08")
    p.add_argument("--min-persons", type=int, default=2)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
