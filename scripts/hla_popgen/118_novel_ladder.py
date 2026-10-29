"""118 novel-allele confidence ladder (VM side): the full distribution of "minimally novel" HLA alleles and every confidence field behind it.

Pool (the "minimally novel" set, decided by Marc 2026-10-29): Immuannot-novel calls whose CDS is rebuildable from the public IPD reference allele + cs string,
carry NO indel, and whose CDS is not an exact named allele in the CURRENT IPD release (--ipd-nuc). One row per distinct (gene, observed CDS) = one candidate allele
(script 91 keys on reference allele + substitution set instead, so the same CDS reached from two templates is two signatures there and one allele here).
Per allele: carriers (all / unrelated), ancestries, platforms (Sequel II vs Revio) among unrelated carriers, homozygous carriers, trio transmission, nearest-allele
distance, effect (synonymous / missense / nonsense), reference-template flags (suffix N/L/S/C/A/Q, truncated), IPD match status and a SHORT-READ confidence:
for every checked site of every carrier the srWGS reads (90 sr: realigned to the person's own called haplotype) give k reads with the called base out of n, and
  LOD = k*log10(p1/e) + (n-k)*log10((1-p1)/(1-e)),  p1 = 0.5 (het site) or 0.98 (hom), e = 0.02 background (error + mis-assignment + reference bias)
so LOD >= 3 means the read counts are 1000x likelier under "the variant is on this haplotype" than under "absent". A carrier's score is its WEAKEST site; an allele's
level comes from its best carrier and its independent-support count from carriers with score >= 3 (relatives collapsed with the same greedy rule as 91).

Inputs (VM): hla_calls_rich.tsv, cohort_membership.tsv, relatedness, Immuannot_refdata, IPD hla_nuc.fasta, lrWGS manifest (platform), and from the case browser
(~/case_browser) case_map.private.tsv + sr.jsonl produced by `90_case_browser_extract.py --exhaustive-snv <big cap>` then `90_case_browser_sr.py sr`.
Outputs (--out-dir, default ~/s10):
  ladder_data.private.json   EXACT counts (mode 600; stays in the Workbench; feeds the exact-count report built by 119 on the VM)
  ladder_data.pub.json       counts of 1-19 persons replaced by classes (1 / 2-19 / >=20; exact only from 20 up), funnel and sr counts <20 -> null; cloud-safe
  novel_ladder.private.tsv   one row per allele, exact (mode 600)
No person ids are written anywhere.
"""
import argparse
import csv
import importlib.util
import json
import math
import os
import random
import re
import sys
from bisect import bisect_right
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
E_BG, P_HET, P_HOM = 0.02, 0.5, 0.98
LOD_CAP = 20.0
NEAR_MAX = 15            # substitutions; beyond this the reference template is probably the wrong one (S08 log 10-23)
SUFFIX = re.compile(r"\d[NLSCAQ]$")
STOP = set("*X.")
CLASSICAL = {"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"}
SR_LEVELS = ["not evaluated", "contradicted (LOD<0)", "weak (0-3)", "moderate (3-6)", "strong (>=6)"]
BUCKETS = ["1 ready to submit", "2 submittable", "3 uncertain", "4 unlikely"]


def log(m):
    print("[118] " + m, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ pure scoring
def site_lod(obs, alt, oth, het):
    """LOD of 'variant on this haplotype' vs 'absent' for one site; None without reads. het None (other haplotype uncalled) is treated as hom."""
    n = obs + alt + oth
    if n <= 0:
        return None
    p1 = P_HET if het else P_HOM
    return obs * math.log10(p1 / E_BG) + (n - obs) * math.log10((1 - p1) / (1 - E_BG))


def case_lod(sites):
    """Weakest site wins. Sites with an extraction error or no reads are ignored; None when nothing is usable."""
    v = []
    for s in sites or []:
        if str(s.get("verdict", "")).startswith("error"):
            continue
        x = site_lod(s.get("obs", 0), s.get("alt", 0), s.get("oth", 0), s.get("het"))
        if x is not None:
            v.append(x)
    return None if not v else max(-LOD_CAP, min(LOD_CAP, min(v)))


def sr_level(best_lod):
    if best_lod is None:
        return 0
    return 1 if best_lod < 0 else 2 if best_lod < 3 else 3 if best_lod < 6 else 4


def effect_of(prot):
    if prot == "synonymous":
        return 0
    return 2 if any(p and p[-1] in STOP for p in prot.split(";")) else 1


def ref_flags(ref_allele, ref_len, gene_median_len):
    """(clean, why): clean = no S/N/Q/L/C/A suffix and not a truncated template."""
    if SUFFIX.search(ref_allele or ""):
        return False, "suffix"
    if gene_median_len and ref_len and ref_len < 0.97 * gene_median_len:
        return False, "truncated"
    return True, ""


def bucket(r):
    """Marc-adjustable default ladder (same rule is re-implemented in 119's JS and tested for parity).
    r: rc (clean template), near (<=15 edits), srl (0-4), srn (independent carriers with LOD>=3), xu (unrelated carriers), tr (transmitted), mp (>=2 platforms)."""
    if not (r["rc"] and r["near"]) or r["srl"] == 1:
        return 3
    if r["srl"] >= 3 and r["srn"] >= 2 and r["xu"] >= 2 and (r["tr"] or r["mp"]):
        return 0
    if r["srl"] >= 3:
        return 1
    return 2


def xclass(n):
    return 0 if n <= 1 else 1 if n < 20 else 2


def build_cube(alleles, masked):
    """alleles -> list of cells {flags..., w}. masked: exact carrier counts only from 20 up (xu/xa null below), class always present."""
    cells = Counter()
    for a in alleles:
        xu, xa = a["xu"], a["xa"]
        key = (a["gc"], a["eff"], a["ipd"], int(a["rc"]), int(a["near"]), a["srl"], min(a["srn"], 2), int(a["tr"]), int(a["mp"]), int(a["hm"]),
               xclass(xu), xclass(xa), (xu if (not masked or xu >= 20) else None), (xa if (not masked or xa >= 20) else None))
        cells[key] += 1
    names = ["gc", "eff", "ipd", "rc", "near", "srl", "srn", "tr", "mp", "hm", "xuc", "xac", "xu", "xa"]
    out = []
    for k, w in sorted(cells.items(), key=lambda kv: str(kv[0])):
        d = dict(zip(names, k))
        d["w"] = w
        out.append(d)
    return out


# ------------------------------------------------------------------ saturation
def saturation(carriers, persons, perms, points, seed=1):
    """carriers: {allele: set(person)} (only persons in `persons` count). Random order of `persons`; for every permutation the number of alleles
    seen at each point (distinct alleles with >=1 carrier among the first n persons). Returns {"n": points, "med": [], "lo": [], "hi": []} (5th/95th)."""
    rng = random.Random(seed)
    plist = sorted(persons)
    keep = {a: [p for p in c if p in set(plist)] for a, c in carriers.items()}
    keep = {a: c for a, c in keep.items() if c}
    runs = []
    for _ in range(perms):
        rng.shuffle(plist)
        pos = {p: i + 1 for i, p in enumerate(plist)}
        first = sorted(min(pos[p] for p in c) for c in keep.values())
        runs.append([bisect_right(first, n) for n in points])
    med, lo, hi = [], [], []
    for j in range(len(points)):
        col = sorted(r[j] for r in runs)
        med.append(col[len(col) // 2])
        lo.append(col[int(0.05 * (len(col) - 1))])
        hi.append(col[int(math.ceil(0.95 * (len(col) - 1)))])
    return {"n": list(points), "med": med, "lo": lo, "hi": hi, "n_alleles": len(keep), "n_alleles_no_carrier_in_set": len(carriers) - len(keep)}


def sat_points(n, k=40):
    if n <= 0:
        return []
    pts = sorted({max(1, int(round(n ** (i / (k - 1))))) for i in range(k)} | {n})
    return pts


# ------------------------------------------------------------------ sr join
def read_sr(case_map_path, sr_path):
    """(gene, person, hap) -> case LOD, from the case browser files. hap is the last digit as in 90's hap_no."""
    pm = {}
    with open(case_map_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pm[r["case"]] = (r["gene"], r["person"], str(r["hap"])[-1])
    out = {}
    with open(sr_path) as f:
        for ln in f:
            if not ln.strip():
                continue
            b = json.loads(ln)
            k = pm.get(b.get("id"))
            if not k:
                continue
            lod = case_lod(b.get("sites"))
            if k not in out or (lod is not None and (out[k] is None or lod > out[k])):
                out[k] = lod
    return out


def load(name):
    spec = importlib.util.spec_from_file_location("m" + name[:2], os.path.join(HERE, name))
    m = importlib.util.module_from_spec(spec)
    sys.modules["m" + name[:2]] = m
    spec.loader.exec_module(m)
    return m


def supp(n, floor=20):
    return None if 0 < n < floor else n


# ------------------------------------------------------------------ main pass
def classify_row(core, ref, cds_mut):
    first = cds_mut.split(",")[0].split("|")
    if len(first) < 2:
        return "no_cs", None
    rows = core.replay_from_ref(ref, first[1])
    if rows is None:
        return "not_rebuildable", None
    res = core.analyse(rows[0], rows[1])
    st = res["stats"]
    if st["ins_events"] or st["del_events"]:
        return ("indel_plus_subs" if st["snv"] else "indel_only"), None
    if not st["snv"]:
        return "cds_identical", None
    prot = []
    for e in res["events"]:
        if e["type"] == "snv" and e.get("effect") in ("missense", "nonsense"):
            prot.append("%s%d%s" % (e["aa_ref"], e["codon"], e["aa_obs"]))
    return "ok", (first[0], int(st["snv"]) if not isinstance(st["snv"], bool) else 1, ";".join(dict.fromkeys(prot)) or "synonymous",
                  rows[1].replace("-", "").upper())


def run(a):
    r91, r92 = load("91_novel_recurrence.py"), load("92_novel_burden.py")
    core = r91.load_core()
    nuc = r91.read_ipd_nuc(os.path.expanduser(a.ipd_nuc))
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    plat = r92.read_platforms(os.path.expanduser(a.manifest))
    sr = read_sr(os.path.join(os.path.expanduser(a.sr_dir), "case_map.private.tsv"), os.path.join(os.path.expanduser(a.sr_dir), "sr.jsonl"))
    log("persons %d, short-read case scores %d" % (len(lr), len(sr)))
    refs, medlen = {}, {}
    funnel = Counter()
    al = {}                                   # (gene, cds) -> dict
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if g not in r91.HLA_GENES or r["person_id"] not in lr:
                continue
            funnel["calls_in_lr_cohort"] += 1
            nov = core.novelty_class(r["consensus"])
            if nov == "known":
                continue
            funnel["novel:" + nov] += 1
            if not r["cds_mut"]:
                funnel["drop:no_cds_mut"] += 1
                continue
            allele = r["cds_mut"].split(",")[0].split("|")[0]
            gname = allele.split("*")[0]
            if gname not in refs:
                p = os.path.join(os.path.expanduser(a.refdata), "CDSseq", gname + ".fa.gz")
                refs[gname] = r91.read_fasta_gz(p) if os.path.exists(p) else {}
                L = sorted(len(s) for s in refs[gname].values())
                medlen[gname] = L[len(L) // 2] if L else 0
            ref = refs[gname].get(allele)
            if not ref:
                funnel["drop:ref_missing"] += 1
                continue
            st, sig = classify_row(core, ref, r["cds_mut"])
            if sig is None:
                funnel["drop:" + st] += 1
                continue
            funnel["pool_calls"] += 1
            k = (g, sig[3])
            d = al.setdefault(k, {"gene": g, "cds": sig[3], "refs": Counter(), "nedit": sig[1], "prot": sig[2], "by": defaultdict(list)})
            d["refs"][allele] += 1
            d["by"][r["person_id"]].append(str(r["hap"])[-1])
    log("distinct pool alleles before IPD filter: %d" % len(al))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unr_all = set(r91.greedy_unrelated(lr, anyk))
    rows, rec = [], []
    for k, d in al.items():
        g = d["gene"]
        m, names = r91.ipd_match(nuc, g, d["cds"])
        funnel["alleles:ipd_" + m] += 1
        if m == "exact":
            continue
        persons = sorted(d["by"])
        unrel = r91.greedy_unrelated(persons, anyk)
        ref_allele = d["refs"].most_common(1)[0][0]
        rc, why = ref_flags(ref_allele, len(refs[ref_allele.split("*")[0]].get(ref_allele, "")), medlen[ref_allele.split("*")[0]])
        sc = {p: max((sr[(g, p, h)] for h in set(d["by"][p]) if (g, p, h) in sr and sr[(g, p, h)] is not None), default=None) for p in persons}
        evald = [p for p in persons if any((g, p, h) in sr for h in set(d["by"][p]))]
        good = [p for p in persons if sc[p] is not None and sc[p] >= 3]
        srn = len(r91.greedy_unrelated(good, anyk))
        best = max((v for v in sc.values() if v is not None), default=None)
        trio = {p: r91.parents_of(p, first, anyk) for p in persons}
        tr = sum(1 for p in persons if trio[p] and any(q in d["by"] for q in trio[p]))
        trn = sum(1 for p in persons if trio[p])
        pl = {plat.get(p, "NA") for p in unrel} - {"NA"}
        hom = sum(1 for p in persons if len(set(d["by"][p])) > 1)
        x = {"gc": 0 if g in CLASSICAL else 1, "eff": effect_of(d["prot"]), "ipd": 0 if m == "none" else 1, "rc": rc, "near": d["nedit"] <= NEAR_MAX,
             "srl": sr_level(best), "srn": srn, "tr": tr > 0, "mp": len(pl) >= 2, "hm": hom > 0, "xu": len(unrel), "xa": len(persons)}
        x["bucket"] = bucket(x)
        rec.append(x)
        rows.append(dict(gene=g, protein=d["prot"], n_edits=d["nedit"], ref=ref_allele, ref_flag=why or "clean", ipd=m, ipd_hit=";".join(names),
                         carriers=len(persons), unrelated=len(unrel), hom=hom, platforms=len(pl), trio_children=trn, trio_transmitted=tr,
                         sr_evaluated=len(evald), sr_best_lod="" if best is None else "%.1f" % best, sr_indep_support=srn, bucket=BUCKETS[x["bucket"]],
                         ancestries=";".join("%s=%d" % kv for kv in sorted(Counter(anc.get(p, "NA") for p in unrel).items()))))
        d["_unrel"] = set(unrel)
        d["_rec"] = x
    rows.sort(key=lambda r: (-r["unrelated"], r["gene"]))
    tp = os.path.join(out, "novel_ladder.private.tsv")
    with open(tp, "w") as f:
        f.write("\t".join(rows[0]) + "\n" if rows else "")
        f.writelines("\t".join(str(r[c]) for c in r) + "\n" for r in rows)
    os.chmod(tp, 0o600)

    pool = [d for d in al.values() if "_rec" in d]
    sat = {}
    pts_u, pts_a = sat_points(len(unr_all)), sat_points(len(lr))
    sets = {"all pool": lambda x: True, "buckets 1-2": lambda x: x["bucket"] <= 1, "bucket 1": lambda x: x["bucket"] == 0,
            "classical genes, buckets 1-2": lambda x: x["bucket"] <= 1 and x["gc"] == 0, "non-classical (E/F/G), buckets 1-2": lambda x: x["bucket"] <= 1 and x["gc"] == 1}
    for name, fn in sets.items():
        sel = [d for d in pool if fn(d["_rec"])]
        sat[name] = {"unrelated": saturation({id(d): d["_unrel"] for d in sel}, unr_all, a.perms, pts_u),
                     "all": saturation({id(d): set(d["by"]) for d in sel}, set(lr), a.perms, pts_a)}
        log("saturation %s: %d alleles" % (name, len(sel)))

    lods = Counter()
    for v in sr.values():
        if v is not None:
            lods[int(math.floor(max(-5, min(15, v))))] += 1
    meta = {"ipd_release": a.ipd_release, "n_lr_persons": len(lr), "n_unrelated_persons": len(unr_all), "e_bg": E_BG, "p_het": P_HET, "p_hom": P_HOM,
            "near_max": NEAR_MAX, "sr_levels": SR_LEVELS, "buckets": BUCKETS, "perms": a.perms}
    for masked, name in ((False, "private"), (True, "pub")):
        data = {"meta": dict(meta, masked=masked), "cube": build_cube(rec, masked),
                "funnel": {k: (v if not masked else supp(v)) for k, v in sorted(funnel.items())},
                "sr_lod_hist": {str(k): (v if not masked else supp(v)) for k, v in sorted(lods.items())}, "saturation": sat}
        p = os.path.join(out, "ladder_data.%s.json" % name)
        with open(p, "w") as f:
            json.dump(data, f, separators=(",", ":"))
        if not masked:
            os.chmod(p, 0o600)
    log("wrote %d alleles to %s" % (len(rows), out))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
    p.add_argument("--manifest", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv")
    p.add_argument("--refdata", default="~/tools/Immuannot_refdata")
    p.add_argument("--ipd-nuc", required=True, help="hla_nuc.fasta of the NEWEST IPD-IMGT/HLA release")
    p.add_argument("--ipd-release", required=True, help="release label recorded in the output, e.g. 3.66.0")
    p.add_argument("--sr-dir", default="~/case_browser", help="holds case_map.private.tsv and sr.jsonl")
    p.add_argument("--out-dir", default="~/s10")
    p.add_argument("--perms", type=int, default=100)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
