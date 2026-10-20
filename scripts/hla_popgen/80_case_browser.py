#!/usr/bin/env python3
"""S08 80: base-level case browser for HLA/KIR call errors (Immuannot alignment, Mendelian, switches).

A deliberately small, dependency-free (stdlib, Python >= 3.8) pair of steps:

  extract  (VM only)   read already-computed outputs, write ONE anonymised JSON of cases
                       (observed CDS + Immuannot's own minimap2 cs string, trio haplotypes, origin tracks).
                       Read-only on shared paths; writes only under --out-dir. Never runs yak.
  build    (anywhere)  turn that JSON into one self-contained index.html (80_case_browser_viewer.html).
  demo     (laptop)    synthetic cases -> index.html, to try the viewer without any data.

The alignment shown for the "Immuannot" tabs is NOT recomputed: it is replayed from the `cds_mut` cs
string Immuannot wrote (reference allele = minimap2 target, observed CDS = query), so what you see is
exactly what the caller saw. Only the parent-vs-child comparisons use a direct banded alignment
(labelled as such in the viewer).

Disclosure: the JSON/HTML hold participant-derived sequence. They stay in the Workbench (open the HTML
from JupyterLab); case ids are `c001`-style, the id->person map is a separate *.private.tsv. Nothing
here is committed (see .gitignore) and the repo holds code only.

    python3 80_case_browser.py demo --out-dir /tmp/cb_demo
    python3 80_case_browser.py build --cases ~/case_browser/cases.json --out-dir ~/case_browser
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VIEWER = os.path.join(HERE, "80_case_browser_viewer.html")

# ---------------------------------------------------------------------------------------- genetics
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(
    (a, b, c) for a in _BASES for b in _BASES for c in _BASES)}


def translate(seq):
    seq = seq.upper()
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))


# ---------------------------------------------------------------------------------------- Immuannot alignment
CS_RE = re.compile(r"(:\d+|=[acgtnACGTN]+|\*[acgtnACGTN]{2}|\+[acgtnACGTN]+|-[acgtnACGTN]+|~[acgtn]{2}\d+[acgtn]{2})")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def parse_cs(cs):
    """minimap2 cs string -> [(op, payload)]; ops '=' (int run or str), '*' (target_base, query_base), '+', '-'."""
    ops, pos, cs = [], 0, cs.strip()
    while pos < len(cs):
        m = CS_RE.match(cs, pos)
        if not m:
            raise ValueError("bad cs near %r" % cs[pos:pos + 12])
        t, pos = m.group(1), m.end()
        if t[0] == ":":
            ops.append(("=", int(t[1:])))
        elif t[0] == "=":
            ops.append(("=", t[1:].upper()))
        elif t[0] == "*":
            ops.append(("*", (t[1].upper(), t[2].upper())))
        elif t[0] in "+-":
            ops.append((t[0], t[1:].upper()))
        else:
            raise ValueError("intron token in a CDS cs string: %r" % t)
    return ops


def alignment_from_paf(ref, obs, qs, qe, ts, te, cs, strand="+"):
    """Replay Immuannot's own CDS alignment (mm2.ipd.cds.paf: query = IPD reference allele, target = the
    observed CDS) as gapped rows. cs semantics are target->query, so '*xy' = observed x / reference y,
    '+s' = bases only in the reference, '-s' = bases only in the observed CDS. Unaligned ends are kept as
    'clip' columns. Returns (aref, aobs, clip5, clip3, problems); problems counts any disagreement between
    the cs string and the two sequences (0 means the replay is exactly what minimap2 aligned)."""
    ref = ref.upper()
    obs = obs.upper()
    if strand == "-":
        ref = revcomp(ref)
        qs, qe = len(ref) - qe, len(ref) - qs
    ti, qi, bad, R, O = ts, qs, 0, [], []
    for op, p in parse_cs(cs):
        if op == "=":
            n = p if isinstance(p, int) else len(p)
            if obs[ti:ti + n] != ref[qi:qi + n] or (not isinstance(p, int) and obs[ti:ti + n] != p):
                bad += 1
            R.append(ref[qi:qi + n])
            O.append(obs[ti:ti + n])
            ti, qi = ti + n, qi + n
        elif op == "*":
            bad += (obs[ti:ti + 1] != p[0]) + (ref[qi:qi + 1] != p[1])
            R.append(p[1])
            O.append(p[0])
            ti, qi = ti + 1, qi + 1
        elif op == "+":
            bad += ref[qi:qi + len(p)] != p
            R.append(p)
            O.append("-" * len(p))
            qi += len(p)
        else:
            bad += obs[ti:ti + len(p)] != p
            R.append("-" * len(p))
            O.append(p)
            ti += len(p)
    bad += (ti != te) + (qi != qe)

    def flank(r, o):
        w = max(len(r), len(o))
        return r.ljust(w, "-"), o.ljust(w, "-")

    lr, lo = flank(ref[:qs], obs[:ts])
    rr, ro = flank(ref[qe:], obs[te:])
    return lr + "".join(R) + rr, lo + "".join(O) + ro, len(lr), len(rr), bad


# ---------------------------------------------------------------------------------------- direct aligner
def banded_align(a, b, match=2, mis=-3, gap_open=-5, gap_ext=-1, band=40):
    """Global affine-gap (Gotoh) alignment of `a` (reference row) vs `b` (observed row), restricted to a band
    around the main diagonal (widened by the length difference). Used only where no Immuannot alignment
    exists (parent-vs-child haplotype CDS). Returns gapped strings (ra, rb)."""
    if a == b:
        return a, b
    n, m = len(a), len(b)
    dlo, dhi = min(0, m - n) - band, max(0, m - n) + band
    W = dhi - dlo + 1
    NEG = -10 ** 9
    M = [[NEG] * W for _ in range(n + 1)]
    X = [[NEG] * W for _ in range(n + 1)]   # a[i-1] against a gap in b
    Y = [[NEG] * W for _ in range(n + 1)]   # b[j-1] against a gap in a
    PM = [[0] * W for _ in range(n + 1)]
    PX = [[0] * W for _ in range(n + 1)]
    PY = [[0] * W for _ in range(n + 1)]
    M[0][-dlo] = 0
    for i in range(n + 1):
        for k in range(W):
            j = i + k + dlo
            if j < 0 or j > m or (i == 0 and j == 0):
                continue
            if i > 0 and j > 0:
                s = match if a[i - 1] == b[j - 1] else mis
                c = (M[i - 1][k], X[i - 1][k], Y[i - 1][k])
                st = c.index(max(c))
                if c[st] > NEG // 2:
                    M[i][k], PM[i][k] = c[st] + s, st
            if i > 0 and k + 1 < W:
                o, e = M[i - 1][k + 1] + gap_open, X[i - 1][k + 1] + gap_ext
                if max(o, e) > NEG // 2:
                    X[i][k], PX[i][k] = (e, 1) if e > o else (o, 0)
            if j > 0 and k > 0:
                o, e = M[i][k - 1] + gap_open, Y[i][k - 1] + gap_ext
                if max(o, e) > NEG // 2:
                    Y[i][k], PY[i][k] = (e, 1) if e > o else (o, 0)
    i, j = n, m
    k = j - i - dlo
    c = (M[i][k], X[i][k], Y[i][k])
    st = c.index(max(c))
    ra, rb = [], []
    while i > 0 or j > 0:
        k = j - i - dlo
        if st == 0:
            ra.append(a[i - 1])
            rb.append(b[j - 1])
            st, i, j = PM[i][k], i - 1, j - 1
        elif st == 1:
            ra.append(a[i - 1])
            rb.append("-")
            st, i = (1 if PX[i][k] else 0), i - 1
        else:
            ra.append("-")
            rb.append(b[j - 1])
            st, j = (2 if PY[i][k] else 0), j - 1
    return "".join(reversed(ra)), "".join(reversed(rb))


# ---------------------------------------------------------------------------------------- events
def _run_len(seq, pos, base):
    """Length of the run of `base` in `seq` that covers/abuts index `pos`."""
    s = e = pos
    while s > 0 and seq[s - 1] == base:
        s -= 1
    while e < len(seq) and seq[e] == base:
        e += 1
    return e - s


def analyse(aref, aobs, clip5=0, clip3=0, exons=None, ars_exons=()):
    """Classify every difference in a gapped pair. `exons` = [(exon_number, n_query_bases)] in 5'->3' order
    of the OBSERVED CDS (they tile it); `ars_exons` = exon numbers forming the peptide-binding groove.
    Columns in the 5'/3' clip flanks are ignored. Returns dict(events, stats, aa_ref, aa_obs, exon_cols)."""
    L = len(aref)
    core0, core1 = clip5, L - clip3
    # per-column ungapped indexes
    ridx, qidx = [-1] * L, [-1] * L
    r = q = 0
    for c in range(L):
        if aref[c] != "-":
            ridx[c] = r
            r += 1
        if aobs[c] != "-":
            qidx[c] = q
            q += 1
    ref_u = aref.replace("-", "")
    obs_u = aobs.replace("-", "")
    # exon of every column (via observed index; deleted columns inherit the previous base's exon)
    bounds, acc = [], 0
    for n, ln in (exons or []):
        bounds.append((n, acc, acc + ln))
        acc += ln

    def exon_of_q(qi):
        for n, s, e in bounds:
            if s <= qi < e:
                return n
        return None

    col_exon, last = [], None
    for c in range(L):
        ex = exon_of_q(qidx[c]) if qidx[c] >= 0 else last
        col_exon.append(ex)
        last = ex if ex is not None else last
    # events
    events, c = [], core0
    shift_mod = 0   # cumulative (ref - obs) length difference mod 3 before this column
    cum = 0
    while c < core1:
        a, b = aref[c], aobs[c]
        if a == b:
            c += 1
            continue
        if a != "-" and b != "-":
            ev = {"type": "snv", "col": c, "len": 1, "ref": a, "obs": b, "refpos": ridx[c] + 1,
                  "qpos": qidx[c] + 1, "frame": cum % 3}
            c += 1
        else:
            kind = "ins" if a == "-" else "del"
            s = c
            while c < core1 and ((aref[c] == "-") if kind == "ins" else (aobs[c] == "-")) and aref[c] != aobs[c]:
                c += 1
            seq = (aobs if kind == "ins" else aref)[s:c]
            ev = {"type": kind, "col": s, "len": c - s, "ref": "" if kind == "ins" else seq,
                  "obs": seq if kind == "ins" else "", "refpos": (ridx[s] if ridx[s] >= 0 else
                                                                   next((ridx[x] for x in range(s, L) if ridx[x] >= 0), 0)) + 1,
                  "qpos": (qidx[s] if qidx[s] >= 0 else
                           next((qidx[x] for x in range(s, L) if qidx[x] >= 0), 0)) + 1,
                  "frame": cum % 3}
            # homopolymer / tandem slippage, judged in the sequence that carries the extra bases
            long_row = (aobs if kind == "ins" else aref).replace("-", "")
            lp = (qidx[s] if kind == "ins" else ridx[s])
            ev["homopolymer"] = False
            ev["tandem"] = False
            if len(set(seq)) == 1:
                rl = _run_len(long_row, lp, seq[0])
                ev["run"] = rl
                ev["homopolymer"] = rl >= 3
            k = len(seq)
            if not ev["homopolymer"] and k >= 2 and lp >= 0:
                if long_row[max(0, lp - k):lp] == seq or long_row[lp + k:lp + 2 * k] == seq:
                    ev["tandem"] = True
            cum += (c - s) * (1 if kind == "del" else -1)
            ev["inframe"] = (ev["len"] % 3 == 0)
        ev["exon"] = col_exon[ev["col"]]
        ev["ars"] = ev["exon"] in ars_exons if ev["exon"] is not None else False
        events.append(ev)
    # amino-acid consequence of substitutions (only while ref/obs codons are in register)
    for ev in events:
        if ev["type"] != "snv":
            continue
        ri, qi = ev["refpos"] - 1, ev["qpos"] - 1
        if (ri - qi) % 3 == 0 and ri // 3 == qi // 3:
            cr, co = ref_u[3 * (ri // 3):3 * (ri // 3) + 3], obs_u[3 * (qi // 3):3 * (qi // 3) + 3]
            if len(cr) == 3 and len(co) == 3:
                ev["codon"] = ri // 3 + 1
                ev["cod_ref"], ev["cod_obs"] = cr, co
                ev["aa_ref"], ev["aa_obs"] = CODON.get(cr, "X"), CODON.get(co, "X")
                ev["effect"] = ("synonymous" if ev["aa_ref"] == ev["aa_obs"] else
                                "nonsense" if ev["aa_obs"] == "*" else "missense")
                continue
        ev["effect"] = "frameshifted"
    for ev in events:
        if ev["type"] != "snv":
            ev["effect"] = "inframe indel" if ev.get("inframe") else "frameshift"
    # per-column amino-acid tracks (letter at the middle base of each codon)
    aa_ref, aa_obs = [" "] * L, [" "] * L
    for c in range(L):
        if ridx[c] >= 0 and ridx[c] % 3 == 1 and 3 * (ridx[c] // 3) + 3 <= len(ref_u):
            aa_ref[c] = CODON.get(ref_u[ridx[c] - 1:ridx[c] + 2], "X")
        if qidx[c] >= 0 and qidx[c] % 3 == 1 and 3 * (qidx[c] // 3) + 3 <= len(obs_u):
            aa_obs[c] = CODON.get(obs_u[qidx[c] - 1:qidx[c] + 2], "X")
    # exon column spans
    exon_cols = []
    for n, s, e in bounds:
        cols = [c for c in range(L) if qidx[c] >= s and qidx[c] < e]
        if cols:
            exon_cols.append({"n": n, "c0": cols[0], "c1": cols[-1] + 1, "ars": n in ars_exons})
    n_snv = sum(1 for e in events if e["type"] == "snv")
    ins_b = sum(e["len"] for e in events if e["type"] == "ins")
    del_b = sum(e["len"] for e in events if e["type"] == "del")
    nm = n_snv + ins_b + del_b
    cds_len = len(ref_u) if len(ref_u) else 1
    indels = [e for e in events if e["type"] != "snv"]
    stats = {
        "nm": nm, "snv": n_snv,
        "ins_events": sum(1 for e in events if e["type"] == "ins"), "ins_bases": ins_b,
        "del_events": sum(1 for e in events if e["type"] == "del"), "del_bases": del_b,
        "homopolymer_indels": sum(1 for e in indels if e.get("homopolymer")),
        "tandem_indels": sum(1 for e in indels if e.get("tandem")),
        "frameshifts": sum(1 for e in indels if not e.get("inframe")),
        "missense": sum(1 for e in events if e.get("effect") == "missense"),
        "synonymous": sum(1 for e in events if e.get("effect") == "synonymous"),
        "nonsense": sum(1 for e in events if e.get("effect") == "nonsense"),
        "in_ars": sum(1 for e in events if e.get("ars")),
        "cds_len": len(ref_u), "obs_len": len(obs_u),
        "pct_cds": round(100.0 * nm / cds_len, 2),
        "clip5": clip5, "clip3": clip3,
    }
    return {"events": events, "stats": stats, "aa_ref": "".join(aa_ref), "aa_obs": "".join(aa_obs),
            "exon_cols": exon_cols}


def make_alignment(label_ref, label_obs, aref, aobs, source, clip5=0, clip3=0, exons=None, ars_exons=(), note="",
                   problems=0):
    """One viewer alignment. `problems` = disagreements found while replaying Immuannot's cs (0 = exact)."""
    a = analyse(aref, aobs, clip5, clip3, exons, ars_exons)
    return {"source": source, "ref_label": label_ref, "obs_label": label_obs, "aref": aref, "aobs": aobs,
            "clip5": clip5, "clip3": clip3, "note": note, "problems": problems, **a}


def novelty_class(consensus):
    """Immuannot 'new' depth: 2nd field -> protein-altering, 3rd -> synonymous, 4th -> beyond-CDS."""
    if not consensus or "new" not in consensus:
        return "known"
    fields = consensus.split("*", 1)[-1].split(":")
    d = fields.index("new") + 1 if "new" in fields else 0
    return {1: "novel-protein", 2: "novel-protein", 3: "novel-synonymous", 4: "novel-beyond-CDS"}.get(d, "novel-other")


ARS = {"classI": (2, 3), "classII": (2,)}
CLASS_I = {"A", "B", "C", "E", "F", "G"}


def ars_for(gene):
    g = gene.replace("HLA-", "")
    if g in CLASS_I:
        return ARS["classI"]
    if g.startswith(("DR", "DQ", "DP")):
        return ARS["classII"]
    return ()


def describe_alignment(a):
    """One-line summary of an alignment: edit count and what kind (for the Mendelian case header)."""
    st, ev = a["stats"], a["events"]
    if not st["nm"]:
        return "identical CDS"
    parts = []
    if st["snv"]:
        parts.append("%d substitution%s (%d missense, %d synonymous)" % (
            st["snv"], "" if st["snv"] == 1 else "s", st["missense"], st["synonymous"]))
    ind = [e for e in ev if e["type"] != "snv"]
    if ind:
        parts.append("%d indel%s (%d in homopolymers, %d frameshifting)" % (
            len(ind), "" if len(ind) == 1 else "s", sum(1 for e in ind if e.get("homopolymer")),
            sum(1 for e in ind if not e.get("inframe"))))
    return "%d edit%s: %s" % (st["nm"], "" if st["nm"] == 1 else "s", "; ".join(parts))


def normalise_cases(cases):
    """Idempotent touch-ups so a cases.json written by an older extract renders the same as a new one:
    fill unmapped exon numbers from their neighbours, and give Mendelian cases the descriptive summary."""
    for c in cases:
        for a in c.get("alignments", []):
            ex = a.get("exon_cols") or []
            for i, e in enumerate(ex):
                if e.get("n") is None:
                    e["n"] = ex[i - 1]["n"] + 1 if i and ex[i - 1]["n"] else i + 1
        if c.get("tab") == "mendel":
            direct = [a for a in c["alignments"] if a["source"] == "direct"]
            if direct:
                best = describe_alignment(min(direct, key=lambda a: a["stats"]["nm"]))
                c["facts"] = [[k, best if k == "closest parent haplotype" else v] for k, v in c.get("facts", [])]
    return cases


# ---------------------------------------------------------------------------------------- build / demo
def build_html(cases, out_dir, meta=None):
    cases = normalise_cases(cases)
    os.makedirs(out_dir, exist_ok=True)
    with open(VIEWER) as f:
        tpl = f.read()
    payload = json.dumps({"meta": meta or {}, "cases": cases}, separators=(",", ":")).replace("</", "<\\/")
    html = tpl.replace("/*__DATA__*/null", payload)
    out = os.path.join(out_dir, "index.html")
    with open(out, "w") as f:
        f.write(html)
    return out


def cmd_build(a):
    with open(a.cases) as f:
        d = json.load(f)
    out = build_html(d["cases"], a.out_dir, d.get("meta"))
    print("wrote %s (%d cases, %.1f KB)" % (out, len(d["cases"]), os.path.getsize(out) / 1024.0))


def cmd_demo(a):
    from random import Random
    rng = Random(7)
    cases = []
    gene_len = {"DEMO-A": 1098, "DEMO-B": 1089}
    for gene, ln in gene_len.items():
        ref = "ATG" + "".join(rng.choice("ACGT") for _ in range(ln - 6)) + "TGA"
        # make some homopolymers so indels have somewhere to slip
        ref = ref[:300] + "AAAAAA" + ref[306:600] + "GGGGG" + ref[605:]
        for tier in (1, 2, 5, 10, 50, 100):
            for rep in range(2):
                obs, ops = list(ref), []
                for _ in range(tier):
                    p = rng.randrange(3, ln - 6)
                    t = rng.random()
                    if t < 0.7:
                        obs[p] = rng.choice([x for x in "ACGT" if x != obs[p]])
                    elif t < 0.85:
                        obs.insert(p, obs[p] if rng.random() < .6 else rng.choice("ACGT"))
                    else:
                        del obs[p]
                aref, aobs = banded_align(ref, "".join(obs))
                novel = "known" if tier <= 2 and rep == 0 else "novel-protein" if rep else "novel-synonymous"
                cases.append({
                    "id": "c%03d" % (len(cases) + 1), "tab": "ladder", "species": "DEMO", "gene": gene, "tier": tier,
                    "novelty": novel, "consensus": "%s*01:01" % gene, "template_distance": tier * 3,
                    "alignments": [make_alignment(gene + " reference (synthetic)", "observed (synthetic)", aref, aobs,
                                                  "synthetic", exons=[(1, 300), (2, 270), (3, ln - 570)], ars_exons=(2, 3))]})
    ref = "ATG" + "".join(rng.choice("ACGT") for _ in range(ln - 6)) + "TGA"
    child = ref[:400] + ("G" if ref[400] != "G" else "A") + ref[401:]
    other = ref[:700] + ref[701:]
    mk = lambda lbl, x, y: make_alignment(lbl[0], lbl[1], *banded_align(x, y), "direct", exons=[(1, 300), (2, 270), (3, ln - 570)], ars_exons=(2, 3))
    cases.append({
        "id": "c%03d" % (len(cases) + 1), "tab": "mendel", "species": "DEMO", "gene": "DEMO-B", "tier": "-", "novelty": "known",
        "verdict": "child allele absent from assigned parent", "consensus": "DEMO-B*07:02",
        "facts": [["child call", "B*07:02"], ["assigned parent", "Parent B (k-mers 41 vs 2)"], ["parent calls", "B*07:03 / B*15:01"]],
        "note": "SYNTHETIC: one substitution against parent hap 1, one 1-bp deletion against hap 2.",
        "alignments": [dict(mk(("parent hap1", "child"), ref, child), title="child vs parent hap 1"),
                       dict(mk(("parent hap2", "child"), other, child), title="child vs parent hap 2")]})
    genes = ["A", "C", "B", "MICA", "DRB3", "DRB1", "DQA1", "DQB1", "DQA2", "TAP2", "TAP1", "DPA1", "DPB1"]
    ori = ["pat"] * 10 + ["mat", "pat", "amb"]
    cases.append({
        "id": "c%03d" % (len(cases) + 1), "tab": "switch", "species": "DEMO", "gene": None, "tier": "-", "novelty": "known",
        "n_switch": 2, "switches": [{"i": 9.5, "gap_kb": 16.2}, {"i": 10.5, "gap_kb": 219.2}],
        "track": [{"gene": g, "pos_kb": 30000 + 400 * i, "origin": o, "pat_k": 3 if o != "pat" else 90, "mat_k": 80 if o == "mat" else 2, "allele": g + "*01:01"}
                  for i, (g, o) in enumerate(zip(genes, ori))]})
    out = build_html(cases, a.out_dir, {"title": "SYNTHETIC DEMO — no real data", "synthetic": True})
    print("wrote", out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--cases", required=True)
    b.add_argument("--out-dir", required=True)
    d = sub.add_parser("demo")
    d.add_argument("--out-dir", required=True)
    import importlib.util as iu
    spec = iu.spec_from_file_location("cb_extract", os.path.join(HERE, "80_case_browser_extract.py"))
    ext = iu.module_from_spec(spec)
    spec.loader.exec_module(ext)
    ext._CORE = sys.modules[__name__]
    ext.register(sub)
    a = ap.parse_args(argv)
    if a.cmd == "build":
        cmd_build(a)
    elif a.cmd == "demo":
        cmd_demo(a)
    else:
        a.func(a)


if __name__ == "__main__":
    main()
