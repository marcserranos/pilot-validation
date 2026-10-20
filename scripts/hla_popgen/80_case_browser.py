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


# ---------------------------------------------------------------------------------------- cs replay
CS_RE = re.compile(r"(:\d+|=[acgtnACGTN]+|\*[acgtnACGTN]{2}|\+[acgtnACGTN]+|-[acgtnACGTN]+|~[acgtn]{2}\d+[acgtn]{2})")


def parse_cs(cs):
    """cs string -> [(op, payload)] ; ops: '=' (n or seq), '*' (ref, obs), '+' seq, '-' seq."""
    ops, pos = [], 0
    cs = cs.strip()
    while pos < len(cs):
        m = CS_RE.match(cs, pos)
        if not m:
            raise ValueError("bad cs near %r" % cs[pos:pos + 12])
        t = m.group(1)
        pos = m.end()
        c = t[0]
        if c == ":":
            ops.append(("=", int(t[1:])))
        elif c == "=":
            ops.append(("=", t[1:].upper()))
        elif c == "*":
            ops.append(("*", (t[1].upper(), t[2].upper())))
        elif c in "+-":
            ops.append((c, t[1:].upper()))
        else:
            raise ValueError("intron token in a CDS cs string: %r" % t)
    return ops


def cs_extent(ops):
    """(target bases consumed, query bases consumed)."""
    tl = ql = 0
    for op, p in ops:
        if op == "=":
            n = p if isinstance(p, int) else len(p)
            tl += n
            ql += n
        elif op == "*":
            tl += 1
            ql += 1
        elif op == "+":
            ql += len(p)
        else:
            tl += len(p)
    return tl, ql


def replay_cs(query, cs, qoff=0):
    """Rebuild the gapped (target, query) rows that minimap2 saw. Matches are read from `query`.
    Returns (aref, aobs, problems) where problems counts query-vs-cs inconsistencies."""
    ops = parse_cs(cs)
    qi, ref, obs, bad = qoff, [], [], 0
    for op, p in ops:
        if op == "=":
            seg = query[qi:qi + p] if isinstance(p, int) else p
            if not isinstance(p, int) and query[qi:qi + len(p)] != p:
                bad += 1
            ref.append(seg)
            obs.append(seg)
            qi += len(seg)
        elif op == "*":
            if query[qi:qi + 1] != p[1]:
                bad += 1
            ref.append(p[0])
            obs.append(p[1])
            qi += 1
        elif op == "+":
            if query[qi:qi + len(p)] != p:
                bad += 1
            ref.append("-" * len(p))
            obs.append(p)
            qi += len(p)
        else:
            ref.append(p)
            obs.append("-" * len(p))
    return "".join(ref), "".join(obs), bad


def resolve_alignment(query, cs, target=None):
    """Replay `cs` against the observed CDS (`query`) and, when the reference allele CDS (`target`) is known,
    place the aligned block inside it (minimap2 may clip ends). Returns dict(aref, aobs, tstart, qstart, ...)
    with unaligned flanks appended as 'clip' columns (shown but never counted as differences)."""
    ops = parse_cs(cs)
    tl, ql = cs_extent(ops)
    best = None
    for qoff in range(0, max(0, len(query) - ql) + 1):
        aref, aobs, bad = replay_cs(query, cs, qoff)
        tstart, tbad = 0, 0
        if target:
            ru = aref.replace("-", "")
            tstart = target.find(ru)
            if tstart < 0:
                tstart, tbad = 0, 1 + sum(1 for x, y in zip(ru, target) if x != y)
        score = bad * 1000 + tbad
        if best is None or score < best[0]:
            best = (score, qoff, tstart, aref, aobs, bad)
        if score == 0:
            break
    _, qoff, tstart, aref, aobs, bad = best
    t_end, q_end = tstart + tl, qoff + ql
    lt = target[:tstart] if target else ""
    rt = target[t_end:] if target else ""
    lq, rq = query[:qoff], query[q_end:]
    pad = lambda a, b: (a.ljust(max(len(a), len(b)), "-"), b.ljust(max(len(a), len(b)), "-"))
    lref, lobs = pad(lt, lq)
    rref, robs = pad(rt, rq)
    return {"aref": lref + aref + rref, "aobs": lobs + aobs + robs,
            "clip5": len(lref), "clip3": len(rref), "tstart": tstart - 0, "qstart": qoff,
            "problems": bad, "ref_known": bool(target)}


# ---------------------------------------------------------------------------------------- direct aligner
def banded_align(a, b, match=2, mis=-3, gap_open=-5, gap_ext=-1, band=None):
    """Global affine-gap (Gotoh) alignment of `a` (reference row) vs `b` (observed row), banded around the
    diagonal. Ties prefer diagonal moves, then gaps placed 3' (right-most), the usual indel convention.
    Returns gapped strings (ra, rb)."""
    n, m = len(a), len(b)
    if a == b:
        return a, b
    band = band or (abs(n - m) + 40)
    NEG = -10 ** 9
    W = 2 * band + 1
    # rows indexed 0..n, columns j = i + (k - band) for k in 0..W-1 (shifted for length difference)
    shift = (m - n) // 2

    def idx(i, j):
        k = j - i - shift + band
        return k if 0 <= k < W and 0 <= j <= m else -1

    M = [[NEG] * W for _ in range(n + 1)]
    X = [[NEG] * W for _ in range(n + 1)]   # gap in b (deletion: a base consumed)
    Y = [[NEG] * W for _ in range(n + 1)]   # gap in a (insertion: b base consumed)
    tb = [[(0, 0, 0)] * W for _ in range(n + 1)]
    k0 = idx(0, 0)
    if k0 < 0:
        return banded_align(a, b, band=band * 2)
    M[0][k0] = 0
    for i in range(n + 1):
        for j in range(max(0, i + shift - band), min(m, i + shift + band) + 1):
            k = idx(i, j)
            if k < 0 or (i == 0 and j == 0):
                continue
            best_m = (NEG, 0)
            if i > 0 and j > 0:
                kp = idx(i - 1, j - 1)
                if kp >= 0:
                    s = match if a[i - 1] == b[j - 1] else mis
                    cand = ((M[i - 1][kp], 0), (X[i - 1][kp], 1), (Y[i - 1][kp], 2))
                    v, f = max(cand, key=lambda t: (t[0], -t[1]))
                    if v > NEG // 2:
                        M[i][k] = v + s
                        tb[i][k] = (0, f, tb[i][k][2])
            if i > 0:
                kp = idx(i - 1, j)
                if kp >= 0:
                    o, e = M[i - 1][kp] + gap_open, X[i - 1][kp] + gap_ext
                    if max(o, e) > NEG // 2:
                        X[i][k] = max(o, e)
                        tb[i][k] = (tb[i][k][0], tb[i][k][1], 0) if False else tb[i][k]
            if j > 0:
                kp = idx(i, j - 1)
                if kp >= 0:
                    o, e = M[i][kp] + gap_open, Y[i][kp] + gap_ext
                    if max(o, e) > NEG // 2:
                        Y[i][k] = max(o, e)
    # traceback by score recomputation (keeps the code short; matrices are small)
    i, j = n, m
    k = idx(i, j)
    if k < 0:
        return banded_align(a, b, band=band * 2)
    state = max((M[i][k], 0), (X[i][k], 1), (Y[i][k], 2), key=lambda t: (t[0], -t[1]))[1]
    ra, rb = [], []
    while i > 0 or j > 0:
        k = idx(i, j)
        if state == 0:
            ra.append(a[i - 1])
            rb.append(b[j - 1])
            s = match if a[i - 1] == b[j - 1] else mis
            kp = idx(i - 1, j - 1)
            tgt = M[i][k] - s
            state = 0 if M[i - 1][kp] == tgt else (1 if X[i - 1][kp] == tgt else 2)
            i, j = i - 1, j - 1
        elif state == 1:
            ra.append(a[i - 1])
            rb.append("-")
            kp = idx(i - 1, j)
            state = 1 if X[i - 1][kp] + gap_ext == X[i][k] else 0
            i -= 1
        else:
            ra.append("-")
            rb.append(b[j - 1])
            kp = idx(i, j - 1)
            state = 2 if Y[i][kp] + gap_ext == Y[i][k] else 0
            j -= 1
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


def make_alignment(label_ref, label_obs, aref, aobs, source, clip5=0, clip3=0, exons=None, ars_exons=(), note=""):
    a = analyse(aref, aobs, clip5, clip3, exons, ars_exons)
    return {"source": source, "ref_label": label_ref, "obs_label": label_obs, "aref": aref, "aobs": aobs,
            "clip5": clip5, "clip3": clip3, "note": note, **a}


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


# ---------------------------------------------------------------------------------------- build / demo
def build_html(cases, out_dir, meta=None):
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
    try:
        from importlib import import_module  # noqa: F401
        import importlib.util as iu
        p = os.path.join(HERE, "80_case_browser_extract.py")
        if os.path.exists(p):
            spec = iu.spec_from_file_location("cb_extract", p)
            mod = iu.module_from_spec(spec)
            spec.loader.exec_module(mod)
            mod.register(sub)
    except ImportError:
        pass
    a = ap.parse_args(argv)
    if a.cmd == "build":
        cmd_build(a)
    elif a.cmd == "demo":
        cmd_demo(a)
    else:
        a.func(a)


if __name__ == "__main__":
    main()
