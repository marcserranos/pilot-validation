"""80 extract (VM side): pick example calls and pull their sequences into one anonymised cases.json.

Read-only on shared paths (pipeline outputs, Immuannot reference, 57b's existing origin table); writes only
under --out-dir. Never runs yak, minimap2 or any BAM access: every alignment shown is replayed from files
Immuannot already wrote (mm2.ipd.cds.paf.gz) and the parent-vs-child comparison is a small direct alignment.

Selection
  ladder   per gene x distance bucket (Immuannot cds_distance 1, 2, ~5, ~10, ~50, ~100, >=150), a few random
           calls each, mixed across novelty classes. HLA from hla_calls_rich.tsv, KIR from a sample of the
           per-person KIR GTFs.
  mendel   57b origin table + the called alleles: child haplotypes whose 2-field allele is absent from the
           parent they were assigned to (plus a few concordant controls).
  switch   57b origin table: contigs whose neighbouring genes were assigned to different parents (plus controls).
"""
import csv
import gzip
import json
import math
import os
import random
import re
import sys
import traceback
from collections import defaultdict

TIERS = {1: (1, 1), 2: (2, 2), 5: (4, 6), 10: (8, 12), 50: (40, 60), 100: (80, 120), 200: (150, 10 ** 6)}
HLA_GENES = ["A", "B", "C", "E", "F", "G", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "DRB3", "DRB4", "DRB5"]
MENDEL_GENES = {"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"}
ATTR_RE = re.compile(r'(\w+)\s+("([^"]*)"|[^;\s]+)')


def where(e):
    """Aggregate-safe description of an exception: type, short message, file:line of the failing frame."""
    fr = traceback.extract_tb(e.__traceback__)[-1]
    return "%s:%s@%s:%d" % (type(e).__name__, str(e)[:40], os.path.basename(fr.filename), fr.lineno)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def bare(g):
    return g[4:] if g.startswith("HLA-") else g


def tier_of(nm):
    for t, (lo, hi) in TIERS.items():
        if lo <= nm <= hi:
            return t
    return None


# ------------------------------------------------------------------ Immuannot file readers
def parse_attrs(s):
    return {m.group(1): (m.group(3) if m.group(3) is not None else m.group(2)) for m in ATTR_RE.finditer(s)}


def read_gtf(path):
    """Immuannot GTF -> {(gene_name, copy): call dict}. Copy index comes from the '.N' suffix of the gene_id."""
    genes, trans, cds, exons = {}, {}, defaultdict(list), defaultdict(list)
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            x = ln.rstrip("\n").split("\t")
            if len(x) < 9:
                continue
            a = parse_attrs(x[8])
            ft = x[2]
            if ft == "gene":
                genes[a["gene_id"]] = dict(a, contig=x[0], start=int(x[3]), end=int(x[4]), strand=x[6])
            elif ft == "transcript":
                trans[a["transcript_id"]] = dict(a, contig=x[0])
            elif ft == "CDS":
                cds[a["transcript_id"]].append((int(x[3]), int(x[4])))
            elif ft == "exon":
                m = re.match(r"\s*(\d+)", a.get("exon_number", ""))
                exons[a["transcript_id"]].append((int(m.group(1)) if m else 0, int(x[3]), int(x[4])))
    out = {}
    for tid, t in trans.items():
        g = genes.get(t.get("gene_id"))
        if g and exons.get(tid) and not all(n for n, _, _ in exons[tid]):      # blank exon_number: number 5'->3'
            rev = g["strand"] == "-"
            exons[tid] = [(k, s_, e_) for k, (_, s_, e_) in enumerate(sorted(exons[tid], key=lambda z: z[1], reverse=rev), 1)]
        if not g:
            continue
        m = re.search(r"\.(\d+)$", t["gene_id"])
        copy = int(m.group(1)) if m else 1
        out[(g["gene_name"], copy)] = {
            "gene_name": g["gene_name"], "copy": copy, "contig": g["contig"], "start": g["start"], "end": g["end"],
            "strand": g["strand"], "template_allele": g.get("template_allele", ""),
            "template_distance": g.get("template_distance", ""), "consensus": t.get("consensus", ""),
            "alleles": t.get("alleles", ""), "cds_distance": t.get("cds_distance", ""),
            "cds_mut": t.get("cds_mut", ""), "warn": t.get("template_warning", ""),
            "cds": sorted(cds.get(tid, [])), "exons": exons.get(tid, [])}
    return out


def exon_lengths(call, obs_len):
    """[(exon_number, bases)] tiling the observed CDS 5'->3' (CDS rows exclude the stop codon: it is added to
    the last exon so the list sums to the CDS record's length)."""
    segs = sorted(call["cds"], reverse=(call["strand"] == "-"))
    out = []
    for s, e in segs:
        n = next((k for k, es, ee in call["exons"] if es <= s and e <= ee), None)
        out.append([n, e - s + 1])
    tot = sum(l for _, l in out)
    if out and 0 < obs_len - tot <= 3:
        out[-1][1] += obs_len - tot
    return [tuple(x) for x in out]


def read_fasta_gz(path, want=None):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None and (want is None or name in want):
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
    if name is not None and (want is None or name in want):
        seqs[name] = "".join(buf)
    return seqs


class Store:
    """Caches per-person Immuannot files and the public IPD reference CDS per gene."""

    def __init__(self, hla_root, kir_root, refdata):
        self.roots = {"HLA": hla_root, "KIR": kir_root}
        self.refdata = refdata
        self._gtf, self._ref = {}, {}

    def odir(self, species, person, hap):
        return os.path.join(self.roots[species], person, "immuannot_output")

    def calls(self, species, person, hap):
        k = (species, person, hap)
        if k not in self._gtf:
            self._gtf[k] = read_gtf(os.path.join(self.odir(species, person, hap), "hap%s.gtf.gz" % hap))
            if len(self._gtf) > 400:
                self._gtf.pop(next(iter(self._gtf)))
        return self._gtf[k]

    def ref_cds(self, gene_name):
        if gene_name not in self._ref:
            p = os.path.join(self.refdata, "CDSseq", gene_name + ".fa.gz")
            self._ref[gene_name] = read_fasta_gz(p) if os.path.exists(p) else {}
        return self._ref[gene_name]

    def observed(self, species, person, hap, call):
        """Observed CDS record for one gene copy, and the PAF lines Immuannot aligned it with."""
        d = os.path.join(self.odir(species, person, hap), "hap%s" % hap)
        suffix = "_%s_%d" % (call["gene_name"], call["copy"])
        rec_id, buf, done = None, [], False
        with gzip.open(os.path.join(d, "cds.fa.gz"), "rt") as f:
            for ln in f:
                if ln.startswith(">"):
                    if rec_id is not None:
                        done = True
                        break
                    name = ln[1:].split()[0]
                    if name.endswith(suffix) and name.startswith(call["contig"]):
                        rec_id = name
                elif rec_id is not None:
                    buf.append(ln.strip())
        seq = "".join(buf) if rec_id is not None else None
        if seq is None:
            raise KeyError("no CDS record")
        paf = []
        pp = os.path.join(d, "mm2.ipd.cds.paf.gz")
        if os.path.exists(pp):
            with gzip.open(pp, "rt") as f:
                for ln in f:
                    x = ln.rstrip("\n").split("\t")
                    if len(x) > 12 and x[5] == rec_id:
                        tags = {t[:2]: t[5:] for t in x[12:]}
                        paf.append({"q": x[0], "qlen": int(x[1]), "qs": int(x[2]), "qe": int(x[3]), "strand": x[4],
                                    "tlen": int(x[6]), "ts": int(x[7]), "te": int(x[8]), "AS": int(tags.get("AS", 0)),
                                    "NM": int(tags.get("NM", 0)), "cs": tags.get("cs", "")})
        return seq, paf


# ------------------------------------------------------------------ alignment bundles
def immuannot_alignment(core, store, species, person, hap, call, ars):
    """Alignment object for one call exactly as Immuannot saw it (PAF replay), or an identity block when
    Immuannot found a perfect CDS match and therefore ran no CDS search."""
    obs, paf = store.observed(species, person, hap, call)
    exons = exon_lengths(call, len(obs))
    first = (call["cds_mut"].split(",")[0].split("|")[0]) if call["cds_mut"] else ""
    ref_all = store.ref_cds(call["gene_name"])
    lines = [p for p in paf if p["q"] == first] or sorted(paf, key=lambda p: -p["AS"])[:1]
    if not lines:
        allele = (call["template_allele"] or first)
        ref = ref_all.get(allele, obs)
        note = "no CDS search was needed: identical to %s" % allele if ref == obs else \
            "no PAF line found for this call; showing observed against itself"
        a = core.make_alignment(allele, "observed", obs, obs, "immuannot-paf", exons=exons, ars_exons=ars,
                                note=note, problems=0 if ref == obs else 1)
        return a, allele
    p = max(lines, key=lambda q: q["AS"])
    ref = ref_all.get(p["q"])
    if ref is None:
        raise KeyError("reference allele CDS missing: %s" % p["q"])
    aref, aobs, c5, c3, bad = core.alignment_from_paf(ref, obs, p["qs"], p["qe"], p["ts"], p["te"], p["cs"], p["strand"])
    bad += int(p["tlen"] != len(obs))
    note = "" if p["q"] == first or not first else "PAF line for %s (cds_mut named %s)" % (p["q"], first)
    return core.make_alignment(p["q"], "observed", aref, aobs, "immuannot-paf", c5, c3, exons, ars, note, problems=bad), p["q"]


def pick_diverse(core, rows, k, rng):
    """Up to k rows, one person at most once, round-robin over novelty classes."""
    rng.shuffle(rows)
    by = defaultdict(list)
    for r in rows:
        by[core.novelty_class(r["consensus"])].append(r)
    out, used = [], set()
    while len(out) < k and any(by.values()):
        for cls in sorted(by):
            while by[cls]:
                r = by[cls].pop()
                if r["person"] not in used:
                    used.add(r["person"])
                    out.append(r)
                    break
            if len(out) >= k:
                break
    return out


# ------------------------------------------------------------------ selection from tables
def load_rich(path, genes):
    keep = set(genes)
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = bare(r["gene"])
            if g in keep:
                rows.append({"person": r["person_id"], "hap": r["hap"], "contig": r["contig"], "gene": g,
                             "copy": int(float(r["copy_index"] or 1)), "consensus": r["consensus"],
                             "cds_distance": r["cds_distance"], "cds_mut": r["cds_mut"],
                             "template_distance": r["template_distance"], "start": r["gene_start"], "end": r["gene_end"],
                             "n_tied": r.get("n_tied", "")})
    return rows


def hap_no(h):
    return str(h)[-1]


def build_ladder(core, store, rows, per_cell, rng, species, report):
    cells = defaultdict(list)
    for r in rows:
        try:
            nm = int(float(r["cds_distance"]))
        except (TypeError, ValueError):
            continue
        t = tier_of(nm)
        if t and r["cds_mut"]:
            cells[(r["gene"], t)].append(r)
    cases = []
    for (gene, t), rs in sorted(cells.items()):
        for r in pick_diverse(core, rs, per_cell, rng):
            try:
                person, hap = r["person"], hap_no(r["hap"])
                call = store.calls(species, person, hap).get((("HLA-" if species == "HLA" else "") + gene, r["copy"]))
                if call is None:
                    report["no_gtf_call"] += 1
                    continue
                a, allele = immuannot_alignment(core, store, species, person, hap, call, core.ars_for(gene) if species == "HLA" else ())
                a["title"] = "Immuannot alignment"
                cases.append(({"tab": "ladder", "species": species, "gene": gene, "tier": t,
                               "novelty": core.novelty_class(r["consensus"]), "consensus": r["consensus"],
                               "template_distance": r["template_distance"], "cds_distance": r["cds_distance"],
                               "warn": call["warn"], "alignments": [a]},
                              {"person": person, "hap": hap, "contig": r["contig"], "gene": gene, "copy": r["copy"]}))
            except Exception as e:  # keep going: report the reason class, never the person
                report["errors:" + where(e)] += 1
    return cases


def kir_rows(store, root, n_people, rng, report):
    people = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d, "immuannot_output")))
    rng.shuffle(people)
    rows = []
    for pid in people[:n_people]:
        for h in ("1", "2"):
            try:
                calls = store.calls("KIR", pid, h)
            except Exception as e:
                report["errors:kir_gtf:" + where(e)] += 1
                continue
            for (g, copy), c in calls.items():
                if g.startswith("KIR") and c["cds_mut"] and c["cds_distance"]:
                    rows.append({"person": pid, "hap": h, "contig": c["contig"], "gene": g, "copy": copy,
                                 "consensus": c["consensus"], "cds_distance": c["cds_distance"], "cds_mut": c["cds_mut"],
                                 "template_distance": c["template_distance"]})
    return rows


# ------------------------------------------------------------------ trios (57's convention: sorted parents, first = 'pat')
def load_trios(relatedness, children):
    first, anyp = defaultdict(set), set()
    lo, hi = 2 ** -2.5, 2 ** -1.5
    with open(relatedness) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            i, j, kin = r["i.s"], r["j.s"], float(r["kin"])
            anyp.add((i, j))
            anyp.add((j, i))
            if lo <= kin <= hi:
                first[i].add(j)
                first[j].add(i)
    out = {}
    for c in children:
        ps = sorted(first.get(c, ()))
        for a in range(len(ps)):
            for b in range(a + 1, len(ps)):
                if (ps[a], ps[b]) not in anyp and c not in out:
                    out[c] = (ps[a], ps[b])
    return out


def read_origin(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["gene"] = bare(r["gene"])
            r["hap"] = hap_no(r["hap"])
            r["pat_k"], r["mat_k"] = int(float(r["pat_k"])), int(float(r["mat_k"]))
            rows.append(r)
    return rows


def build_mendel(core, store, rich_by, origin, trio_of, n_bad, n_ok, rng, report):
    tf = core_two_field
    pool_bad, pool_ok = [], []
    for o in origin:
        if o["origin"] not in ("pat", "mat") or o["gene"] not in MENDEL_GENES or o["child"] not in trio_of:
            continue
        child_rows = [r for r in rich_by.get((o["child"], o["gene"]), []) if hap_no(r["hap"]) == o["hap"] and r["copy"] == 1]
        if not child_rows:
            continue
        pid = trio_of[o["child"]][0 if o["origin"] == "pat" else 1]
        prows = [r for r in rich_by.get((pid, o["gene"]), []) if r["copy"] == 1]
        ctf = tf(child_rows[0]["consensus"])
        if not prows or ctf is None:
            continue
        item = (o, child_rows[0], pid, prows)
        (pool_ok if ctf in {tf(r["consensus"]) for r in prows} else pool_bad).append(item)
    rng.shuffle(pool_bad)
    rng.shuffle(pool_ok)
    cases = []
    for item, tag in [(i, "bad") for i in pool_bad[:n_bad]] + [(i, "ok") for i in pool_ok[:n_ok]]:
        o, cr, pid, prows = item
        try:
            ars = core.ars_for(o["gene"])
            gname = "HLA-" + o["gene"]
            child = o["child"]
            ccall = store.calls("HLA", child, o["hap"])[(gname, 1)]
            aln, cobs = [], store.observed("HLA", child, o["hap"], ccall)[0]
            a, _ = immuannot_alignment(core, store, "HLA", child, o["hap"], ccall, ars)
            a["title"] = "child hap vs its reference (Immuannot)"
            aln.append(a)
            pobs = []
            for r in sorted(prows, key=lambda r: hap_no(r["hap"])):
                h = hap_no(r["hap"])
                pcall = store.calls("HLA", pid, h)[(gname, 1)]
                pa, _ = immuannot_alignment(core, store, "HLA", pid, h, pcall, ars)
                pa["title"] = "parent hap%s vs its reference (Immuannot)" % h
                aln.append(pa)
                pobs.append((h, r, store.observed("HLA", pid, h, pcall)[0], pcall))
            for h, r, ps, pcall in pobs:
                ra, rb = core.banded_align(ps, cobs)
                d = core.make_alignment("parent hap%s" % h, "child", ra, rb, "direct", exons=exon_lengths(ccall, len(cobs)), ars_exons=ars,
                                        note="direct alignment of the two called CDS (not Immuannot's); parent call %s" % r["consensus"])
                d["title"] = "child vs parent hap%s (direct)" % h
                aln.append(d)
            who = "Parent A" if o["origin"] == "pat" else "Parent B"
            best = min((x["stats"]["nm"] for x in aln if x["source"] == "direct"), default=None)
            cases.append(({"tab": "mendel", "species": "HLA", "gene": o["gene"], "tier": "-", "novelty": core.novelty_class(cr["consensus"]),
                           "consensus": cr["consensus"], "verdict": ("child allele absent from assigned parent" if tag == "bad" else "control: allele present in assigned parent"),
                           "facts": [["child call", cr["consensus"]], ["assigned parent", "%s (k-mers %d vs %d)" % (who, o["pat_k"], o["mat_k"]) if o["origin"] == "pat" else "%s (k-mers %d vs %d)" % (who, o["mat_k"], o["pat_k"])],
                                     ["parent calls", " / ".join(r["consensus"] for r in prows)], ["closest parent haplotype", "%s edits" % best]],
                           "note": "Assignment is 57b's k-mer parent of origin. 'Closest parent haplotype' is the edit count of the direct CDS alignment.",
                           "alignments": aln},
                          {"person": child, "hap": o["hap"], "contig": o["contig"], "gene": o["gene"], "copy": 1, "parent": pid}))
        except Exception as e:
            report["errors:mendel:" + where(e)] += 1
    return cases


def core_two_field(consensus):
    if not consensus or "*" not in consensus:
        return None
    g, rest = bare(consensus).split("*", 1)
    f = rest.split(":")
    return None if len(f) < 2 or "new" in f[:2] or not f[0] or not f[1] else "%s*%s:%s" % (g, f[0], f[1])


def build_switch(origin, rich_by, n_sw, n_ok, rng):
    by = defaultdict(list)
    for o in origin:
        by[(o["child"], o["hap"], o["contig"])].append(o)
    sw, ok = [], []
    for key, rs in by.items():
        trk = []
        for o in rs:
            rr = [r for r in rich_by.get((o["child"], o["gene"]), []) if hap_no(r["hap"]) == o["hap"] and r["copy"] == 1 and r["contig"] == o["contig"]]
            if rr:
                trk.append((int(float(rr[0]["start"])), o, rr[0]["consensus"]))
        trk.sort(key=lambda x: x[0])
        if len(trk) < 4:
            continue
        res = [i for i, (_, o, _) in enumerate(trk) if o["origin"] in ("pat", "mat")]
        flips = []
        for a, b in zip(res, res[1:]):
            if trk[a][1]["origin"] != trk[b][1]["origin"]:
                flips.append({"i": a, "gap_kb": round((trk[b][0] - trk[a][0]) / 1000.0, 1), "to": b})
        (sw if flips else ok).append((key, trk, flips))
    rng.shuffle(sw)
    rng.shuffle(ok)
    cases = []
    for key, trk, flips in sw[:n_sw] + ok[:n_ok]:
        x0 = trk[0][0]
        for f in flips:       # marker halfway between the two flanking resolved genes (they may not be adjacent)
            f["i"] = (f["i"] + f["to"] - 1) / 2.0
        cases.append(({"tab": "switch", "species": "HLA", "gene": None, "tier": "-", "novelty": "known", "n_switch": len(flips),
                       "switches": [{"i": f["i"], "gap_kb": f["gap_kb"]} for f in flips],
                       "track": [{"gene": o["gene"], "pos_kb": round((s - x0) / 1000.0, 1), "origin": o["origin"], "pat_k": o["pat_k"], "mat_k": o["mat_k"], "allele": cons}
                                 for s, o, cons in trk]},
                      {"person": key[0], "hap": key[1], "contig": key[2], "gene": "", "copy": ""}))
    return cases


# ------------------------------------------------------------------ driver
def run(a, core):
    rng = random.Random(a.seed)
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    store = Store(os.path.expanduser(a.hla_root), os.path.expanduser(a.kir_root), os.path.expanduser(a.refdata))
    report = defaultdict(int)
    allc = []
    rich = load_rich(os.path.expanduser(a.rich), HLA_GENES)
    log("rich rows: %d" % len(rich))
    allc += build_ladder(core, store, [r for r in rich], a.per_cell, rng, "HLA", report)
    if a.kir_people:
        krows = kir_rows(store, os.path.expanduser(a.kir_root), a.kir_people, rng, report)
        log("KIR rows with CDS differences: %d" % len(krows))
        allc += build_ladder(core, store, krows, a.kir_per_cell, rng, "KIR", report)
    origin_path = os.path.expanduser(a.origin)
    if os.path.exists(origin_path) and (a.mendel or a.switch):
        origin = read_origin(origin_path)
        kids = {o["child"] for o in origin}
        by = defaultdict(list)
        for r in rich:
            by[(r["person"], r["gene"])].append(r)
        if a.mendel:
            rel = os.path.expanduser(a.relatedness)
            if os.path.exists(rel):
                trio_of = load_trios(rel, kids)
                log("origin rows %d, children %d, trios %d" % (len(origin), len(kids), len(trio_of)))
                allc += build_mendel(core, store, by, origin, trio_of, a.mendel, a.mendel_controls, rng, report)
            else:
                log("relatedness file not found (bucket unmounted? ENVIRONMENT quirks #14/#35): skipping Mendelian cases")
        if a.switch:
            allc += build_switch(origin, by, a.switch, a.switch_controls, rng)
    cases, private = [], []
    for k, (c, priv) in enumerate(allc, 1):
        c["id"] = "c%03d" % k
        cases.append(c)
        private.append(dict(priv, case=c["id"]))
    with open(os.path.join(out, "cases.json"), "w") as f:
        json.dump({"meta": {"title": "Case browser — real data, keep inside the Workbench", "seed": a.seed}, "cases": cases}, f, separators=(",", ":"))
    pp = os.path.join(out, "case_map.private.tsv")
    with open(pp, "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=["case", "person", "hap", "contig", "gene", "copy", "parent"], restval="")
        w.writeheader()
        w.writerows(private)
    os.chmod(pp, 0o600)
    # aggregate-only summary to the console (no ids)
    by_tab = defaultdict(int)
    bad = 0
    for c in cases:
        by_tab[c["tab"] + "/" + c["species"]] += 1
        bad += sum(1 for x in c.get("alignments", []) if x.get("problems"))
    log("cases by tab: %s" % dict(by_tab))
    log("alignments whose cs replay disagreed with the sequences: %d" % bad)
    log("skipped/errors: %s" % dict(report))
    return cases


def register(sub):
    p = sub.add_parser("extract", help="VM: pull cases into <out-dir>/cases.json (+ private id map)")
    p.add_argument("--out-dir", default="~/case_browser")
    p.add_argument("--hla-root", default="~/pipeline_outputs/people")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--refdata", default="~/tools/Immuannot_refdata")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--origin", default="~/s05/internal/57b/gene_origin.INTERNAL.tsv")
    p.add_argument("--relatedness", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
    p.add_argument("--per-cell", type=int, default=3)
    p.add_argument("--kir-people", type=int, default=300)
    p.add_argument("--kir-per-cell", type=int, default=2)
    p.add_argument("--mendel", type=int, default=30)
    p.add_argument("--mendel-controls", type=int, default=6)
    p.add_argument("--switch", type=int, default=25)
    p.add_argument("--switch-controls", type=int, default=6)
    p.add_argument("--seed", type=int, default=1)
    p.set_defaults(func=lambda a: run(a, _CORE))


_CORE = None  # set by 80_case_browser.main()
