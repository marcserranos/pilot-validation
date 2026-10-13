#!/usr/bin/env python3
"""S05 66: complementary suppression ACROSS the committed PUBLIC tables (laptop; no participant data).

Why: `_disclosure.suppress_equations` only sees the equations inside one table. A reader can also
subtract across tables (65 genotype counts = sums of allele-fraction bins; 58 copy-number counts =
sums of relative-depth bins), along a chain (a waterfall step minus the next), or against a total
published elsewhere (64 call outcomes sum to the calls tested). The 2026-10-12 S05 fact-check found
1-19 people recoverable that way. This tool reads the PUBLIC tables as published (display strings
only, no raw values), writes down every sum a reader can rebuild, and masks more cells until no
count of 1-19 people can be recovered.

Model. Every relation is a linear equation `sum(coef * cell) = const` over published count cells
and "residuals" (unpublished, non-negative parts a reader can still compute: the 23 people with no
listed ancestry, people measured but not called, people lost at a waterfall step, haplotypes not
split, ...). A cell is shown if its text is a plain non-negative integer; anything else (`<20`,
`<20*`, blank) is unknown. Unknowns are solved for per connected component (SVD): an unknown is
recoverable when the equations fix its value. It leaks when that value is 1-19. A second check
looks at one equation at a time with the bounds `<20` = 1..19, residual >= 0 (e.g. two `<20` cells
that sum to 2).

Fix. For a leaking unknown, take its equation with the fewest unknowns and mask the smallest shown
non-zero, non-protected count in it as `<20*` (the complementary mark), blanking the rate columns
that would give it back. Repeat until nothing leaks. Protected cells are values published in many
places (cohort sizes, run-summary totals); masking them in one table would not hide them. Then
percentages are checked on their own: a shown pct with a shown denominator that pins the numerator
(or its complement) to 1-19 is blanked.

Tables (edited in place, untouched cells byte-identical):
    reports/hla_popgen/65_trb_assembly_genotype/reads_{unrelated,all}/  genotype_by_stratum x allele_fraction_histogram
    reports/hla_popgen/58_regional_depth/          trb_cn_by_ancestry x trb_rel_depth_histogram x depth_summary_by_region
    reports/hla_popgen/50_qc_foundation/           waterfall, dq_ascertainment
    reports/hla_popgen/64_novel_allele_validation/tables/  failures_summary, summary_by_stratum (x run_summary)
    reports/hla_popgen/51_kir_structure_qc/        allelic_pair_events
    reports/hla_popgen/53_hla_ase_pilot/           ase_by_gene_class (80 people per gene, README)
    reports/hla_popgen/55_structure_expression/    structure_expression
    reports/hla_popgen/56_phasing_typing_benchmark/ refined_sweep

    python3 scripts/hla_popgen/66_cross_table_suppression.py --check   # report leaks, exit 1 if any
    python3 scripts/hla_popgen/66_cross_table_suppression.py           # mask and write back
    python3 scripts/hla_popgen/66_cross_table_suppression.py --dir65 ~/s05/results/65_reads   # one 65 PUBLIC dir
    python3 scripts/hla_popgen/66_cross_table_suppression.py --dir58 ~/s05/results/58         # one 58 PUBLIC dir

Not modelled (see the S05 report's residual list): differences between two cohorts that are both
published (65 reads_all minus reads_unrelated = related people; 50 waterfall MID 495 -> 487), and
the cumulative score-threshold sweeps (56 refined_sweep min_score t vs t+1, 58 n_rel_lt_0.5/0.6/0.7).
"""
import argparse
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np

SUPPRESS_BELOW = 20
PRIMARY = f"<{SUPPRESS_BELOW}"
SECONDARY = f"<{SUPPRESS_BELOW}*"
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
_INT = re.compile(r"^\d+$")
_NUM = re.compile(r"^-?\d+(\.\d+)?$")

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.path.join(REPO, "reports", "hla_popgen")


# ---------------------------------------------------------------------------
# Tables as published: rows of strings, written back byte-identical where untouched.
# ---------------------------------------------------------------------------
class Table:
    def __init__(self, name, path):
        self.name, self.path = name, path
        with open(path) as f:
            text = f.read()
        self.trailing_nl = text.endswith("\n")
        lines = text.split("\n")
        if self.trailing_nl:
            lines = lines[:-1]
        self.header = lines[0].split("\t")
        self.col = {c: i for i, c in enumerate(self.header)}
        self.rows = [ln.split("\t") for ln in lines[1:]]
        self.masked = 0       # cells masked <20* by this tool
        self.blanked = 0      # rate cells blanked by this tool
        self.dirty = False

    def __len__(self):
        return len(self.rows)

    def get(self, r, c):
        i = self.col[c]
        row = self.rows[r]
        return row[i] if i < len(row) else ""

    def set(self, r, c, v):
        i = self.col[c]
        row = self.rows[r]
        while len(row) <= i:
            row.append("")
        if row[i] != v:
            row[i] = v
            self.dirty = True

    def where(self, **eq):
        """Row indices whose columns equal the given strings (a list value = any of)."""
        out = []
        for r in range(len(self.rows)):
            ok = True
            for c, v in eq.items():
                x = self.get(r, c)
                if (x not in v) if isinstance(v, (list, tuple, set)) else (x != v):
                    ok = False
                    break
            if ok:
                out.append(r)
        return out

    def text(self):
        body = "\n".join("\t".join(r) for r in [self.header] + self.rows)
        return body + ("\n" if self.trailing_nl else "")

    def save(self):
        with open(self.path, "w") as f:
            f.write(self.text())


def cell_value(s):
    s = str(s).strip()
    return int(s) if _INT.match(s) else None


def cell_bounds(s):
    s = str(s).strip()
    if s == PRIMARY:
        return 1, SUPPRESS_BELOW - 1
    if s.startswith("<"):
        return 1, math.inf
    return 0, math.inf


# ---------------------------------------------------------------------------
# The problem: equations over cells ("c", table, row, col) and residuals ("r", name).
# ---------------------------------------------------------------------------
class Problem:
    def __init__(self):
        self.tables = {}
        self.eqs = []                       # (terms [(key, coef)], const, label)
        self.rate_links = defaultdict(dict)  # table -> {col: [rate cols]}
        self.protected = set()
        self.pct_checks = []                # (table, pct_col, denom_col, decimals, also_blank)
        self.unfixable = {}
        self.log = []

    # -- building ---------------------------------------------------------
    def load(self, name, path, required=True):
        if not os.path.exists(path):
            if required:
                raise FileNotFoundError(path)
            return None
        t = Table(name, path)
        self.tables[name] = t
        return t

    @staticmethod
    def k(t, r, c):
        return ("c", t.name, r, c)

    def eq(self, terms, const=0, label=""):
        terms = [(k, c) for k, c in terms if k is not None]
        if len(terms) >= 2 or (terms and const):
            self.eqs.append((terms, const, label))

    def sum_eq(self, total, parts, resid=None, label=""):
        """total = sum(parts) [+ resid]. total may be an int constant."""
        terms = [(p, -1) for p in parts]
        if resid is not None:
            terms.append((("r", resid), -1))
        if isinstance(total, int):
            self.eq(terms, -total, label)
        else:
            self.eq([(total, 1)] + terms, 0, label)

    def protect(self, t, rows, cols):
        for r in rows:
            for c in cols:
                if c in t.col:
                    self.protected.add(self.k(t, r, c))

    # -- reading ------------------------------------------------------------
    def disp(self, key):
        return self.tables[key[1]].get(key[2], key[3])

    def val(self, key):
        return None if key[0] == "r" else cell_value(self.disp(key))

    def bounds(self, key):
        return (0, math.inf) if key[0] == "r" else cell_bounds(self.disp(key))

    # -- analysis -----------------------------------------------------------
    def analyse(self, keep=False):
        """-> (leaks, inconsistent). leak = (var, value, how). keep=True also stores, per component,
        what solve() needs to pick a cell that really un-pins a leak (self._comp)."""
        inconsistent, reduced = [], []
        self._known_eqs = set()
        for ei, (terms, const, label) in enumerate(self.eqs):
            b, u = const, defaultdict(int)
            for key, c in terms:
                v = self.val(key)
                if v is None:
                    u[key] += c
                else:
                    b -= c * v
            u = {k: c for k, c in u.items() if c}
            if not u:
                self._known_eqs.add(ei)
                if b != 0:
                    inconsistent.append((label, b))
                continue
            reduced.append((u, b, label, ei))
        # connected components over unknowns
        parent = {}

        def find(x):
            while parent.setdefault(x, x) != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for u, _b, _l, _e in reduced:
            ks = list(u)
            for k in ks[1:]:
                parent[find(k)] = find(ks[0])
            find(ks[0])
        comps = defaultdict(list)
        for i, (u, _b, _l, _e) in enumerate(reduced):
            comps[find(next(iter(u)))].append(i)
        leaks = []
        self._comp = {}
        for root, idx in comps.items():
            vars_ = sorted({k for i in idx for k in reduced[i][0]}, key=str)
            vi = {k: j for j, k in enumerate(vars_)}
            A = np.zeros((len(idx), len(vars_)))
            b = np.zeros(len(idx))
            for row, i in enumerate(idx):
                for k, c in reduced[i][0].items():
                    A[row, vi[k]] = c
                b[row] = reduced[i][1]
            # eigen-decomposition of A'A (n x n, n = unknowns): null space and least-squares solution
            lam, V = np.linalg.eigh(A.T @ A)
            pos = lam > 1e-9 * max(1.0, lam.max(initial=0))
            Vp = V[:, pos]
            x = Vp @ ((Vp.T @ (A.T @ b)) / lam[pos])
            if np.abs(A @ x - b).max(initial=0) > 1e-6:
                inconsistent.append((f"component of {len(vars_)} unknowns ({reduced[idx[0]][2]})", None))
            null = V[:, ~pos]
            found = False
            for j, k in enumerate(vars_):
                if null.shape[1] == 0 or np.abs(null[j]).max() < 1e-8:
                    v = x[j]
                    rv = round(v)
                    # a masked cell must never be recoverable (even a <20* holding >= 20); an
                    # unpublished residual only matters when it is 1-19 people
                    if abs(v - rv) < 1e-6 and (1 <= rv < SUPPRESS_BELOW or (k[0] == "c" and rv > 0)):
                        leaks.append((k, int(rv), "solved"))
                        found = True
            if keep and found:
                rows = [reduced[i][3] for i in idx]
                comp = dict(vi=vi, A=A, rows={e: r for r, e in enumerate(rows)}, Vp=Vp, lam=lam[pos])
                for k in vars_:
                    self._comp[k] = comp
        # one equation at a time, with bounds
        for u, b, label, ei in reduced:
            if len(u) < 2:
                continue
            bd = {k: self.bounds(k) for k in u}
            for k, c in u.items():
                lo_r = hi_r = 0.0
                for k2, c2 in u.items():
                    if k2 == k:
                        continue
                    lo2, hi2 = bd[k2]
                    lo_r += c2 * lo2 if c2 > 0 else c2 * hi2
                    hi_r += c2 * hi2 if c2 > 0 else c2 * lo2
                if math.isinf(lo_r) or math.isinf(hi_r):
                    continue
                a, z = sorted(((b - hi_r) / c, (b - lo_r) / c))
                a, z = max(a, bd[k][0]), min(z, bd[k][1])
                if a == z and 1 <= a < SUPPRESS_BELOW:
                    leaks.append((k, int(a), "bounds"))
        return leaks, inconsistent

    def _unpins(self, var, cell):
        """Would masking `cell` leave `var` undetermined? Exact within var's component: the new
        unknown's column must lie in the component's column space with a non-zero `var` part. A cell
        that also sits in a fully shown equation would be pinned again by it -> False."""
        comp = self._comp.get(var)
        if comp is None:
            return False
        eq_of = self._index()
        a = np.zeros(comp["A"].shape[0])
        for ei in eq_of[cell]:
            if ei in self._known_eqs:
                return False
            r = comp["rows"].get(ei)
            if r is not None:
                a[r] += sum(c for k, c in self.eqs[ei][0] if k == cell)
        nz = np.flatnonzero(a)
        if not len(nz):
            return False
        A, Vp = comp["A"], comp["Vp"]
        xp = Vp @ ((Vp.T @ (A[nz].T @ a[nz])) / comp["lam"])
        if np.abs(A @ xp - a).max() > 1e-6:
            return False
        return abs(xp[comp["vi"][var]]) > 1e-6

    def _maskable(self, k):
        return k[0] == "c" and k not in self.protected and (self.val(k) or 0) > 0

    def _ranked_candidates(self, var, max_depth=6):
        """Shown, maskable cells near `var`: (depth, value, key), breadth-first over equations."""
        eq_of = self._index()
        seen_eq, frontier, out, seen = set(), {var}, [], set()
        for depth in range(1, max_depth + 1):
            nxt = set()
            for u in frontier:
                for ei in eq_of[u]:
                    if ei in seen_eq:
                        continue
                    seen_eq.add(ei)
                    for k, _c in self.eqs[ei][0]:
                        if self._maskable(k) and k not in seen:
                            seen.add(k)
                            out.append((depth, self.val(k), str(k), k))
                        elif self.val(k) is None and k != var:
                            nxt.add(k)
            frontier = nxt
            if not frontier:
                break
        out.sort()
        return out

    def mask(self, key):
        t = self.tables[key[1]]
        t.set(key[2], key[3], SECONDARY)
        t.masked += 1
        for rc in self.rate_links[t.name].get(key[3], ()):
            if rc in t.col and t.get(key[2], rc) != "":
                t.set(key[2], rc, "")
                t.blanked += 1

    def _index(self):
        if getattr(self, "_eq_of", None) is None or self._eq_n != len(self.eqs):
            self._eq_of = defaultdict(list)
            for i, (terms, _c, _l) in enumerate(self.eqs):
                for k in {k for k, _ in terms}:
                    self._eq_of[k].append(i)
            self._eq_n = len(self.eqs)
        return self._eq_of

    def solve(self, max_passes=300, max_tests=40):
        """Greedy. Per pass, for each leak whose equations were not touched yet in this pass: mask the
        nearest (then smallest) shown cell that provably un-pins it; if none does on its own, the
        smallest shown cell of its own equations. Then re-analyse. -> True when nothing leaks except
        leaks with no maskable cell left (self.unfixable)."""
        eq_of = self._index()
        for _ in range(max_passes):
            leaks, _inc = self.analyse(keep=True)
            leaks = [lk for lk in leaks if lk[0] not in self.unfixable]
            if not leaks:
                return not self.unfixable
            touched = set()
            for var, v, how in leaks:
                if any(i in touched for i in eq_of[var]):
                    continue
                cands = self._ranked_candidates(var)
                pick = None
                if how == "solved":
                    pick = next((k for _d, _v, _s, k in cands[:max_tests] if self._unpins(var, k)), None)
                if pick is None:
                    near = [c for c in cands if c[0] == 1] or cands
                    pick = near[0][3] if near else None
                if pick is None:
                    self.unfixable[var] = (v, how)
                    continue
                self.log.append(f"mask {fmt_key(self, pick)}={self.val(pick)} (protects {fmt_key(self, var)}, {how})")
                self.mask(pick)
                touched.update(eq_of[var])
                touched.update(eq_of[pick])
        return False

    def check_pcts(self, apply=True):
        """A shown pct with a shown denominator pins the numerator; blank it when the numerator or its
        complement is 1-19. -> list of (table, row, pct_col)."""
        hits = []
        for tname, pc, dc, dec, also in self.pct_checks:
            t = self.tables.get(tname)
            if t is None or pc not in t.col or dc not in t.col:
                continue
            half = 0.5 * 10 ** -dec
            for r in range(len(t)):
                p, n = t.get(r, pc).strip(), cell_value(t.get(r, dc))
                if not _NUM.match(p) or not n:
                    continue
                p = float(p)
                lo, hi = math.ceil(n * (p - half) / 100 - 1e-9), math.floor(n * (p + half) / 100 + 1e-9)
                ks = [k for k in range(max(lo, 0), min(hi, n) + 1)]
                if ks and (all(1 <= k < SUPPRESS_BELOW for k in ks) or all(1 <= n - k < SUPPRESS_BELOW for k in ks)):
                    hits.append((tname, r, pc))
                    if apply:
                        for c in [pc] + list(also):
                            if c in t.col and t.get(r, c) != "":
                                t.set(r, c, "")
                                t.blanked += 1
        return hits


def fmt_key(P, key):
    if key[0] == "r":
        return "residual " + "/".join(map(str, key[1]))
    t = P.tables[key[1]]
    lab = "/".join(t.get(key[2], c) for c in t.header[:3])
    return f"{key[1]}[{lab}].{key[3]}"


def _anc_resid_needed(P, t, all_r, parts, col):
    """Is there an unpublished 'other ancestry' part? True unless ALL equals the ancestry sum."""
    tot = P.val(P.k(t, all_r, col))
    vals = [P.val(P.k(t, r, col)) for r in parts]
    if tot is None or any(v is None for v in vals):
        return True
    return tot != sum(vals)


# ---------------------------------------------------------------------------
# Equation sets per script.
# ---------------------------------------------------------------------------
GENO_COLS = ["n_people_measured", "n_called", "n_0copies", "n_1copy", "n_2copies"]


def spec_65(P, d, tag="65"):
    """65 reads-mode PUBLIC dir: genotype_by_stratum x allele_fraction_histogram (x run_summary)."""
    G = P.load(f"{tag}/genotype_by_stratum", os.path.join(d, "genotype_by_stratum.tsv"))
    H = P.load(f"{tag}/allele_fraction_histogram", os.path.join(d, "allele_fraction_histogram.tsv"), required=False)
    RS = P.load(f"{tag}/run_summary", os.path.join(d, "run_summary.tsv"), required=False)
    meta = {RS.get(r, "metric"): RS.get(r, "value") for r in range(len(RS))} if RS else {}
    cuts = [float(x) for x in meta.get("af_cuts", "0.2,0.8").split(",")]
    cuts_rule = meta.get("genotyper", "cuts") == "cuts"
    est = ["insertion_allele_freq", "pct_0copies", "hwe_p"]
    P.rate_links[G.name] = {"n_people_measured": ["pct_called"], "n_called": ["pct_called"] + est,
                            "n_0copies": est, "n_1copy": est, "n_2copies": est}
    P.protect(G, range(len(G)), ["n_people_measured"])
    P.pct_checks.append((G.name, "pct_called", "n_people_measured", 1, []))
    k = P.k
    feats = sorted({G.get(r, "feature") for r in range(len(G))})
    oth = {}
    for f in feats:
        all_r = G.where(feature=f, stratum_type="ancestry", stratum="ALL")
        if not all_r:
            continue
        a = all_r[0]
        anc = G.where(feature=f, stratum_type="ancestry", stratum=ANCESTRY_ORDER)
        plat = G.where(feature=f, stratum_type="platform")
        oth[f] = _anc_resid_needed(P, G, a, anc, "n_people_measured")
        for c in GENO_COLS:
            P.sum_eq(k(G, a, c), [k(G, r, c) for r in anc], (G.name, f, "oth", c) if oth[f] else None,
                     f"{G.name} {f} {c}: ancestries -> ALL")
            if plat:
                P.sum_eq(k(G, a, c), [k(G, r, c) for r in plat], None, f"{G.name} {f} {c}: platforms -> ALL")
    for r in range(len(G)):
        lab = f"{G.get(r, 'feature')}/{G.get(r, 'stratum')}"
        P.sum_eq(k(G, r, "n_called"), [k(G, r, c) for c in ("n_0copies", "n_1copy", "n_2copies")], None,
                 f"{G.name} {lab}: genotypes -> called")
        P.sum_eq(k(G, r, "n_people_measured"), [k(G, r, "n_called")], (G.name, lab, "not_called"),
                 f"{G.name} {lab}: called + not called -> measured")
    if RS:
        rs = {RS.get(r, "metric"): r for r in range(len(RS))}
        P.protect(RS, range(len(RS)), ["value"])
        # people_called_any_feature = insertion_any n_called only under the cuts rule: the mixture
        # genotyper can leave a person ambiguous on one feature and call them on another.
        for metric, feat, col in (("people_called_any_feature", "insertion_any", "n_called"),
                                  ("people_genotyped", None, "n_people_measured")):
            if metric not in rs or (metric == "people_called_any_feature" and not cuts_rule):
                continue
            for f in ([feat] if feat else feats):
                a = G.where(feature=f, stratum_type="ancestry", stratum="ALL")
                if a:
                    P.sum_eq(k(RS, rs[metric], "value"), [k(G, a[0], col)], None, f"{RS.name} {metric} = {G.name}")
    if H is None:
        return
    groups = defaultdict(list)
    for r in range(len(H)):
        groups[(H.get(r, "feature"), H.get(r, "stratum_type"), H.get(r, "stratum"))].append(r)
    for (f, st, s), rows in groups.items():
        g = G.where(feature=f, stratum_type=st, stratum=s)
        if not g:
            continue
        bins = [k(H, r, "n_people") for r in rows]
        if cuts_rule:   # mixture genotyper: the histogram also holds 'ambiguous' people
            P.sum_eq(k(G, g[0], "n_called"), bins, None, f"{H.name} {f}/{s}: bins -> {G.name} n_called")
            lo_edges = [float(H.get(r, "af_lo")) for r in rows]
            for col, lo, hi in (("n_0copies", -1, cuts[0]), ("n_1copy", cuts[0], cuts[1]), ("n_2copies", cuts[1], 9)):
                part = [k(H, r, "n_people") for r, e in zip(rows, lo_edges) if lo - 1e-9 <= e < hi - 1e-9]
                P.sum_eq(k(G, g[0], col), part, None, f"{H.name} {f}/{s}: bins [{max(lo, 0)},{min(hi, 1)}) -> {G.name} {col}")
    by_bin = defaultdict(dict)
    for (f, st, s), rows in groups.items():
        for r in rows:
            by_bin[(f, H.get(r, "af_lo"))][(st, s)] = r
    for (f, lo), cells in by_bin.items():
        if ("ancestry", "ALL") not in cells:
            continue
        a = cells[("ancestry", "ALL")]
        anc = [cells[("ancestry", x)] for x in ANCESTRY_ORDER if ("ancestry", x) in cells]
        plat = [r for (st, _s), r in cells.items() if st == "platform"]
        P.sum_eq(k(H, a, "n_people"), [k(H, r, "n_people") for r in anc],
                 (H.name, f, lo, "oth") if oth.get(f, True) else None, f"{H.name} {f} bin {lo}: ancestries -> ALL")
        if plat:
            P.sum_eq(k(H, a, "n_people"), [k(H, r, "n_people") for r in plat], None,
                     f"{H.name} {f} bin {lo}: platforms -> ALL")


def spec_58(P, d, tag="58"):
    """58: trb_cn_by_ancestry x trb_rel_depth_histogram x depth_summary_by_region."""
    CN = P.load(f"{tag}/trb_cn_by_ancestry", os.path.join(d, "trb_cn_by_ancestry.tsv"))
    HB = P.load(f"{tag}/trb_rel_depth_histogram", os.path.join(d, "trb_rel_depth_histogram.tsv"))
    DS = P.load(f"{tag}/depth_summary_by_region", os.path.join(d, "depth_summary_by_region.tsv"), required=False)
    RS = P.load(f"{tag}/run_summary", os.path.join(d, "run_summary.tsv"), required=False)
    meta = {RS.get(r, "metric"): RS.get(r, "value") for r in range(len(RS))} if RS else {}
    cuts = [float(x) for x in meta.get("cn_cuts", "0.25,0.75").split(",")]
    cn_mapq = meta.get("cn_mapq", "0")
    k = P.k
    cols = ["n_people", "n_0copies", "n_1copy", "n_2copies"]
    P.rate_links[CN.name] = {c: ["pct_0copies", "deletion_allele_freq", "hwe_p"] for c in cols}
    P.protect(CN, range(len(CN)), ["n_people"])
    genes = sorted({CN.get(r, "gene") for r in range(len(CN))})
    cn_all = {}
    for gname in genes:
        a = CN.where(gene=gname, ancestry="ALL")
        anc = CN.where(gene=gname, ancestry=ANCESTRY_ORDER)
        if not a:
            continue
        cn_all[gname] = a[0]
        resid = _anc_resid_needed(P, CN, a[0], anc, "n_people")
        for c in cols:
            P.sum_eq(k(CN, a[0], c), [k(CN, r, c) for r in anc], (CN.name, gname, "oth", c) if resid else None,
                     f"{CN.name} {gname} {c}: ancestries -> ALL")
    for r in range(len(CN)):
        P.sum_eq(k(CN, r, "n_people"), [k(CN, r, c) for c in cols[1:]], None,
                 f"{CN.name} {CN.get(r, 'gene')}/{CN.get(r, 'ancestry')}: copies -> people")
    ds_all = {}
    if DS is not None:
        dcols = [c for c in DS.header if c.startswith("n_rel_lt_")] + ["n_people"]
        P.protect(DS, range(len(DS)), ["n_people"])
        regions = defaultdict(list)
        for r in range(len(DS)):
            regions[(DS.get(r, "region"), DS.get(r, "mapq"))].append(r)
        for (reg, mq), rows in regions.items():
            a = [r for r in rows if DS.get(r, "stratum_type") == "ancestry" and DS.get(r, "stratum") == "ALL"]
            if not a:
                continue
            if mq == cn_mapq:
                ds_all[reg] = a[0]
            anc = [r for r in rows if DS.get(r, "stratum_type") == "ancestry" and DS.get(r, "stratum") in ANCESTRY_ORDER]
            plat = [r for r in rows if DS.get(r, "stratum_type") == "platform"]
            resid = _anc_resid_needed(P, DS, a[0], anc, "n_people")
            for c in dcols:
                P.sum_eq(k(DS, a[0], c), [k(DS, r, c) for r in anc], (DS.name, reg, mq, "oth", c) if resid else None,
                         f"{DS.name} {reg} mapq{mq} {c}: ancestries -> ALL")
                if plat:
                    P.sum_eq(k(DS, a[0], c), [k(DS, r, c) for r in plat], None,
                             f"{DS.name} {reg} mapq{mq} {c}: platforms -> ALL")
    hg = defaultdict(list)
    for r in range(len(HB)):
        hg[HB.get(r, "gene")].append(r)
    for gname, rows in hg.items():
        bins = [(float(HB.get(r, "bin_lo")), float(HB.get(r, "bin_hi")), k(HB, r, "n_people")) for r in rows]
        if gname in cn_all:
            tot = k(CN, cn_all[gname], "n_people")
        elif gname in ds_all:
            tot = k(DS, ds_all[gname], "n_people")
        else:
            tot = None
        if tot is not None:
            P.sum_eq(tot, [b[2] for b in bins], None, f"{HB.name} {gname}: bins -> people")
        if gname in cn_all:
            for col, lo, hi in (("n_0copies", -1, cuts[0]), ("n_1copy", cuts[0], cuts[1]), ("n_2copies", cuts[1], 99)):
                part = [b[2] for b in bins if lo - 1e-9 <= b[0] < hi - 1e-9]
                P.sum_eq(k(CN, cn_all[gname], col), part, None, f"{HB.name} {gname}: bins [{max(lo, 0)},{hi}) -> {CN.name} {col}")
        if gname in ds_all:
            for c in DS.header:
                if c.startswith("n_rel_lt_"):
                    x = float(c[len("n_rel_lt_"):])
                    part = [b[2] for b in bins if b[1] <= x + 1e-9]
                    P.sum_eq(k(DS, ds_all[gname], c), part, None, f"{HB.name} {gname}: bins < {x} -> {DS.name} {c}")


WATERFALL_CHAIN = ["all_cohort_membership_rows", "in_lr", "assembly_platform", "has_hla_calls",
                   "unrelated_kin_lt_0.0442", "unrelated_and_ancestry_p_ge_0.90",
                   "unrelated_and_ancestry_p_ge_0.95", "unrelated_and_ancestry_p_ge_0.98"]
# Cohort denominators published in many other tables (and in prose): masking them here hides nothing.
WATERFALL_PROTECTED_ROWS = ["in_lr", "assembly_platform", "has_hla_calls", "unrelated_kin_lt_0.0442"]


def spec_50(P, d):
    k = P.k
    W = P.load("50/waterfall", os.path.join(d, "waterfall.tsv"))
    anc = [c for c in ANCESTRY_ORDER if c in W.col]
    plat = [c for c in W.header if c.startswith("platform_")]
    rows = {W.get(r, "step"): r for r in range(len(W))}
    chain = [rows[s] for s in WATERFALL_CHAIN if s in rows]
    P.protect(W, [rows[s] for s in WATERFALL_PROTECTED_ROWS if s in rows], ["total"] + anc + plat)
    P.protect(W, range(len(W)), ["total", "platform_sequel2"])
    for s in WATERFALL_CHAIN:
        if s not in rows:
            continue
        r = rows[s]
        # people with no listed ancestry drop out at the ancestry-probability filters
        resid = None if s.startswith("unrelated_and_ancestry_p") else ("50/waterfall", s, "oth")
        P.sum_eq(k(W, r, "total"), [k(W, r, c) for c in anc], resid, f"waterfall {s}: ancestries -> total")
        P.sum_eq(k(W, r, "total"), [k(W, r, c) for c in plat], None, f"waterfall {s}: platforms -> total")
    lost0 = all(W.get(r, "total") == "0" for s, r in rows.items() if s.startswith("lost_no_hla_calls"))
    for a, b in zip(chain, chain[1:]):
        sa, sb = W.get(a, "step"), W.get(b, "step")
        exact = lost0 and sa in ("in_lr", "assembly_platform") and sb in ("assembly_platform", "has_hla_calls")
        for c in ["total"] + anc + plat:
            P.sum_eq(k(W, a, c), [k(W, b, c)], None if exact else ("50/waterfall", sa, sb, c),
                     f"waterfall {c}: {sa} -> {sb}")

    Q = P.load("50/dq_ascertainment", os.path.join(d, "dq_ascertainment.tsv"))
    P.rate_links[Q.name] = {"n_haps": ["pct_split", "lo", "hi"], "n_split": ["pct_split", "lo", "hi"]}
    P.pct_checks.append((Q.name, "pct_split", "n_haps", 2, ["lo", "hi"]))
    cells = {}
    for r in range(len(Q)):
        cells[(Q.get(r, "scope"), Q.get(r, "DQA1_group"), Q.get(r, "DQB1_group"))] = r
    scopes = sorted({s for s, _a, _b in cells if s == "ALL" or s.startswith(("anc_", "plat_"))})
    pairs = [("G1", "G1"), ("G1", "G2"), ("G2", "G1"), ("G2", "G2")]
    for s in scopes:
        for ga in ("G1", "G2"):
            if (s, ga, "any") not in cells:
                continue
            parts = [cells[(s, ga, gb)] for gb in ("G1", "G2") if (s, ga, gb) in cells]
            for c in ("n_haps", "n_split"):
                P.sum_eq(k(Q, cells[(s, ga, "any")], c), [k(Q, r, c) for r in parts], None, f"dq {s} {ga}: groups -> any")
        for ga, gb in pairs + [("G1", "any"), ("G2", "any")]:
            if (s, ga, gb) in cells:
                r = cells[(s, ga, gb)]
                P.sum_eq(k(Q, r, "n_haps"), [k(Q, r, "n_split")], ("dq", s, ga, gb, "not_split"),
                         f"dq {s} {ga}/{gb}: split + one contig -> haps")
    for ga, gb in pairs + [("G1", "any"), ("G2", "any")]:
        if ("ALL", ga, gb) not in cells:
            continue
        a = cells[("ALL", ga, gb)]
        for pre in ("anc_", "plat_"):
            parts = [r for (s, x, y), r in cells.items() if s.startswith(pre) and (x, y) == (ga, gb)]
            for c in ("n_haps", "n_split"):
                P.sum_eq(k(Q, a, c), [k(Q, r, c) for r in parts], None, f"dq {ga}/{gb} {c}: {pre}* -> ALL")
    sp = Q.where(scope="split_haps_cross_group")
    if sp and ("ALL", "G1", "any") in cells:
        g1a, g2a = cells[("ALL", "G1", "any")], cells[("ALL", "G2", "any")]
        P.sum_eq(k(Q, sp[0], "n_haps"), [k(Q, g1a, "n_split"), k(Q, g2a, "n_split")], None, "dq split haps total")
        cross = [cells[("ALL", "G1", "G2")], cells[("ALL", "G2", "G1")]]
        P.sum_eq(k(Q, sp[0], "n_split"), [k(Q, r, "n_split") for r in cross], None, "dq split cross-group")
    sc = Q.where(scope="same_contig_haps_cross_group")
    if sc:
        P.eq([(k(Q, sc[0], "n_split"), 1), (("r", ("dq", "ALL", "G1", "G2", "not_split")), -1),
              (("r", ("dq", "ALL", "G2", "G1", "not_split")), -1)], 0, "dq same-contig cross-group = one-contig G1/G2 + G2/G1")
        if ("ALL", "G1", "any") in cells:
            P.eq([(k(Q, sc[0], "n_haps"), 1), (("r", ("dq", "ALL", "G1", "any", "not_split")), -1),
                  (("r", ("dq", "ALL", "G2", "any", "not_split")), -1)], 0, "dq same-contig haps total")


def spec_64(P, d, with_strata=False):
    """failures_summary x run_summary. with_strata adds summary_by_stratum: its ~2,000-unknown
    component is too slow for the greedy solver as written (hours), so it is off by default and the
    stratum table is NOT yet re-masked (open item)."""
    k = P.k
    F = P.load("64/failures_summary", os.path.join(d, "failures_summary.tsv"))
    RS = P.load("64/run_summary", os.path.join(d, "run_summary.tsv"))
    rs = {RS.get(r, "metric"): r for r in range(len(RS))}
    P.protect(RS, range(len(RS)), ["value"])
    st = F.where(stage="call_status")
    lr = F.where(stage="long_read")
    ok = F.where(stage="call_status", reason="ok")
    P.protect(F, ok, ["n_calls"])  # = calls tested by reads, quoted as the long-read base
    if "calls_selected" in rs:
        P.sum_eq(k(RS, rs["calls_selected"], "value"), [k(F, r, "n_calls") for r in st], None, "64 call_status -> calls_selected")
    if ok:
        P.sum_eq(k(F, ok[0], "n_calls"), [k(F, r, "n_calls") for r in lr], None, "64 long_read reasons -> ok")
    if not with_strata:
        return
    S = P.load("64/summary_by_stratum", os.path.join(d, "summary_by_stratum.tsv"))
    verdicts = ["supported", "contradicted", "ambiguous"]
    cnt = ["n_calls", "n_calls_tested"] + ["n_" + v for v in verdicts]
    P.rate_links[S.name] = {"n_" + v: ["pct_" + v] for v in verdicts}
    P.rate_links[S.name]["n_calls_tested"] = ["pct_" + v for v in verdicts]
    idx = {}
    for r in range(len(S)):
        idx[tuple(S.get(r, c) for c in ("novelty_class", "gene", "platform", "ancestry", "hp_split"))] = r
    for r in range(len(S)):
        lab = "/".join(S.get(r, c) for c in ("novelty_class", "gene", "platform", "ancestry", "hp_split"))
        P.sum_eq(k(S, r, "n_calls_tested"), [k(S, r, "n_" + v) for v in verdicts], None, f"64 {lab}: verdicts -> tested")
        P.sum_eq(k(S, r, "n_calls"), [k(S, r, "n_calls_tested")], ("64", lab, "untested"), f"64 {lab}: tested + untested -> calls")
    keys = list(idx)
    ncls = sorted({x[0] for x in keys})
    genes = sorted({x[1] for x in keys})
    single = [g for g in genes if g not in ("classical_8", "ALL")]
    plats = sorted({x[2] for x in keys} - {"ALL"})

    for (nc, g, pl, an, sp) in keys:
        if sp == "all":
            parts = [idx[(nc, g, pl, an, s2)] for s2 in ("hp_only", "other") if (nc, g, pl, an, s2) in idx]
            for c in cnt:
                P.sum_eq(k(S, idx[(nc, g, pl, an, sp)], c), [k(S, r, c) for r in parts], None, f"64 {nc}/{g}/{pl}/{an}: hp split")
        if pl == "ALL":
            parts = [idx[(nc, g, p2, an, sp)] for p2 in plats if (nc, g, p2, an, sp) in idx]
            for c in cnt:
                P.sum_eq(k(S, idx[(nc, g, pl, an, sp)], c), [k(S, r, c) for r in parts], None, f"64 {nc}/{g}/{an}/{sp}: platforms")
        if an == "ALL":
            parts = [idx[(nc, g, pl, a2, sp)] for a2 in ANCESTRY_ORDER if (nc, g, pl, a2, sp) in idx]
            for c in cnt:
                P.sum_eq(k(S, idx[(nc, g, pl, an, sp)], c), [k(S, r, c) for r in parts], ("64", nc, g, pl, sp, "oth", c),
                         f"64 {nc}/{g}/{pl}/{sp}: ancestries")
        if g == "classical_8":
            parts = [idx[(nc, g2, pl, an, sp)] for g2 in single if (nc, g2, pl, an, sp) in idx]
            for c in cnt:
                P.sum_eq(k(S, idx[(nc, g, pl, an, sp)], c), [k(S, r, c) for r in parts], None, f"64 {nc}/{pl}/{an}/{sp}: genes")
        if g == "ALL" and (nc, "classical_8", pl, an, sp) in idx:
            for c in cnt:
                P.sum_eq(k(S, idx[(nc, g, pl, an, sp)], c), [k(S, idx[(nc, "classical_8", pl, an, sp)], c)], None,
                         f"64 {nc}/{pl}/{an}/{sp}: ALL = classical_8 (8 genes only)")
    tops = [idx[(nc, "ALL", "ALL", "ALL", "all")] for nc in ncls if (nc, "ALL", "ALL", "ALL", "all") in idx]
    if "novel_calls_unrelated" in rs:
        P.sum_eq(k(RS, rs["novel_calls_unrelated"], "value"), [k(S, r, "n_calls") for r in tops], None, "64 classes -> novel calls")
    if ok:
        P.sum_eq(k(F, ok[0], "n_calls"), [k(S, r, "n_calls_tested") for r in tops], None, "64 classes tested -> ok calls")


def spec_51(P, d):
    k = P.k
    T = P.load("51/allelic_pair_events", os.path.join(d, "allelic_pair_events.tsv"))
    P.rate_links[T.name] = {"co_called": ["co_called_pct", "overlap_share_pct", "same_contig_share_pct", "diff_contig_share_pct"],
                            "haps_with_either": ["co_called_pct"],
                            "overlap": ["overlap_share_pct"], "same_contig": ["same_contig_share_pct"],
                            "diff_contig": ["diff_contig_share_pct"]}
    P.protect(T, range(len(T)), ["haps_with_either"])
    parts = ["overlap", "same_contig", "diff_contig"]
    for r in range(len(T)):
        P.sum_eq(k(T, r, "co_called"), [k(T, r, c) for c in parts], None, f"51 {T.get(r, 'pair')}/{T.get(r, 'scope')}: event types")
    for pair in sorted({T.get(r, "pair") for r in range(len(T))}):
        a, c_ = T.where(pair=pair, scope="all_haps"), T.where(pair=pair, scope="complete_haps")
        if a and c_:
            for c in ["haps_with_either", "co_called"] + parts:
                P.sum_eq(k(T, a[0], c), [k(T, c_[0], c)], ("51", pair, "incomplete", c), f"51 {pair} {c}: complete + incomplete")
    for c in parts:
        P.pct_checks.append((T.name, f"{c}_share_pct", "co_called", 1, []))


def spec_53(P, d, n_people=80):
    k = P.k
    T = P.load("53/ase_by_gene_class", os.path.join(d, "ase_by_gene_class.tsv"))
    stats = [c for c in T.header if c not in ("gene", "pair_class", "n_people")]
    P.rate_links[T.name] = {"n_people": [], "n_people_ge30_informative": ["pct_imbalanced_p01_gt_1.41x"]}
    for gname in sorted({T.get(r, "gene") for r in range(len(T))}):
        rows = T.where(gene=gname)
        # every person of the pilot has exactly one pair class per gene (README: n = 80)
        P.sum_eq(n_people, [k(T, r, "n_people") for r in rows], None, f"53 {gname}: classes -> {n_people} people")
    for r in range(len(T)):
        if T.get(r, "pair_class") != "mrna_identical":
            P.sum_eq(k(T, r, "n_people"), [k(T, r, "n_people_ge30_informative")], ("53", r, "lt30"), "53 >=30 informative")
    P.pct_checks.append((T.name, "pct_imbalanced_p01_gt_1.41x", "n_people_ge30_informative", 1, []))
    _ = stats


def spec_55(P, d):
    k = P.k
    T = P.load("55/structure_expression", os.path.join(d, "structure_expression.tsv"))
    P.rate_links[T.name] = {"n_people": ["pct_tpm_lt1"]}
    P.pct_checks.append((T.name, "pct_tpm_lt1", "n_people", 1, []))
    sets = defaultdict(dict)
    for r in range(len(T)):
        sets[(T.get(r, "gene"), T.get(r, "set"))][T.get(r, "copies")] = r
    for (gname, s), cp in sets.items():
        if "rho" in cp:
            P.sum_eq(k(T, cp["rho"], "n_people"), [k(T, cp[c], "n_people") for c in ("0", "1", "2") if c in cp], None,
                     f"55 {gname}/{s}: copies -> people")
        if s == "kir_all" and (gname, "kir_both_haps_complete") in sets:
            sub = sets[(gname, "kir_both_haps_complete")]
            for c, r in cp.items():
                if c in sub:
                    P.sum_eq(k(T, r, "n_people"), [k(T, sub[c], "n_people")], ("55", gname, c, "incomplete"),
                             f"55 {gname} copies {c}: complete + incomplete -> all")
    rho = [r for (_g, _s), cp in sets.items() for c, r in cp.items() if c == "rho"]
    P.protect(T, rho, ["n_people"])


def spec_56(P, d):
    k = P.k
    T = P.load("56/refined_sweep", os.path.join(d, "refined_sweep.tsv"))
    P.rate_links[T.name] = {"n": ["pct", "lo", "hi"], "k": ["pct", "lo", "hi"]}
    for chk in sorted({T.get(r, "check") for r in range(len(T))}):
        a0 = T.where(check=chk, item="ALL", min_score="0")
        pl = [r for r in T.where(check=chk, item="ALL") if T.get(r, "min_score").startswith("platform:")]
        if a0 and pl:
            for c in ("n", "k"):
                P.sum_eq(k(T, a0[0], c), [k(T, r, c) for r in pl], None, f"56 {chk}: platform pairs -> ALL")
        scores = sorted({T.get(r, "min_score") for r in T.where(check=chk)} - {T.get(r, "min_score") for r in pl})
        for t in scores:
            a = T.where(check=chk, item="ALL", min_score=t)
            parts = [r for r in T.where(check=chk, min_score=t) if T.get(r, "item") != "ALL"]
            if a and parts:
                for c in ("n", "k"):
                    P.sum_eq(k(T, a[0], c), [k(T, r, c) for r in parts], None, f"56 {chk} score>={t}: items -> ALL")


# ---------------------------------------------------------------------------
def build(which=None, dir65=None, dir58=None):
    P = Problem()
    if dir65 or dir58:
        if dir65:
            spec_65(P, dir65, "65")
        if dir58:
            spec_58(P, dir58, "58")
        return P
    d65 = os.path.join(R, "65_trb_assembly_genotype")
    specs = {
        "65_unrelated": lambda: spec_65(P, os.path.join(d65, "reads_unrelated"), "65u"),
        "65_all": lambda: spec_65(P, os.path.join(d65, "reads_all"), "65a"),
        "58": lambda: spec_58(P, os.path.join(R, "58_regional_depth")),
        "50": lambda: spec_50(P, os.path.join(R, "50_qc_foundation")),
        "64": lambda: spec_64(P, os.path.join(R, "64_novel_allele_validation", "tables")),
        "51": lambda: spec_51(P, os.path.join(R, "51_kir_structure_qc")),
        "53": lambda: spec_53(P, os.path.join(R, "53_hla_ase_pilot")),
        "55": lambda: spec_55(P, os.path.join(R, "55_structure_expression")),
        "56": lambda: spec_56(P, os.path.join(R, "56_phasing_typing_benchmark")),
    }
    for name, f in specs.items():
        if which is None or name in which:
            f()
    return P


def report(P, leaks, inconsistent, pct_hits):
    lines = []
    for label, b in inconsistent:
        lines.append(f"INCONSISTENT equation (model wrong?): {label} {'' if b is None else f'(off by {b})'}")
    for var, v, how in leaks:
        lines.append(f"LEAK {fmt_key(P, var)} = {v} ({how})")
    for t, r, c in pct_hits:
        lines.append(f"LEAK via pct {t}[{r}].{c}")
    return lines


def run(P, apply=True, verbose=False):
    """Solve, check pcts, write back. -> dict summary."""
    _leaks0, inconsistent = P.analyse()
    ok = P.solve() if apply else None
    pct_hits = P.check_pcts(apply=apply)
    leaks, inconsistent2 = P.analyse()
    if apply:
        for t in P.tables.values():
            if t.dirty:
                t.save()
    out = dict(leaks=leaks, inconsistent=inconsistent2, pct_hits=pct_hits, solved=ok,
               unfixable=dict(P.unfixable),
               masked={n: t.masked for n, t in P.tables.items() if t.masked or t.blanked},
               blanked={n: t.blanked for n, t in P.tables.items() if t.blanked})
    if verbose:
        for ln in P.log:
            print(ln)
    return out


def apply_dir(dir65=None, dir58=None):
    """For 58/65 PUBLIC aggregation on the VM: run the cross-table pass over one output dir."""
    return run(build(dir65=dir65, dir58=dir58), apply=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report leaks only; exit 1 if any")
    ap.add_argument("--only", nargs="*", help="subset of: 65_unrelated 65_all 58 50 64 51 53 55 56")
    ap.add_argument("--dir65", help="a 65 reads-mode PUBLIC output dir (instead of the committed tables)")
    ap.add_argument("--dir58", help="a 58 PUBLIC output dir (instead of the committed tables)")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    P = build(a.only, a.dir65, a.dir58)
    if a.check:
        leaks, inconsistent = P.analyse()
        pct_hits = P.check_pcts(apply=False)
        for ln in report(P, leaks, inconsistent, pct_hits):
            print(ln)
        print(f"{len(P.eqs)} equations over {len(P.tables)} tables: {len(leaks)} leaks, "
              f"{len(pct_hits)} pct leaks, {len(inconsistent)} inconsistent")
        return 1 if (leaks or pct_hits or inconsistent) else 0
    res = run(P, apply=True, verbose=a.verbose)
    for n in sorted(set(res["masked"]) | set(res["blanked"])):
        t = P.tables[n]
        print(f"{n}: {t.masked} cells masked <20*, {t.blanked} rate cells blanked")
    for var, (v, how) in res["unfixable"].items():
        print(f"UNFIXABLE (only protected cells left): {fmt_key(P, var)} = {v} ({how})")
    for ln in report(P, res["leaks"], res["inconsistent"], []):
        print(ln)
    return 1 if (res["leaks"] or res["inconsistent"]) else 0


if __name__ == "__main__":
    sys.exit(main())
