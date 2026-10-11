#!/usr/bin/env python3
"""61 -- KIR gene content and allele frequencies by ancestry, as one Figure 1 panel (S05 WS-F).

Four views, all per predicted ancestry (AFR/AMR/EAS/EUR/MID/SAS), all on COMPLETE KIR haplotypes
only (51's rule: every KIR call of the haplotype on ONE contig, and both ends, KIR3DL3 and
KIR3DL2, present -- `51_kir_structure_qc.classify_hap`, reused, not re-derived):
  a. gene frequency      % of complete haplotypes carrying each of the 17 KIR genes
                         (+ person-level carrier frequency in the table, among people whose two
                         haplotypes are both complete)
  b. top alleles         top N alleles per gene per ancestry (named only when >= 20 haplotypes)
  c. cA / cB share       cB = any B-content gene on the haplotype (41/43/51 definition)
  d. genes per haplotype distribution of the number of distinct KIR genes per complete haplotype

Three modes:
  --mode vm     Reads the per-person KIR GTFs (~/pipeline_outputs_kir/<pid>/immuannot_output/
                hap{1,2}.gtf.gz) with 44's contig-aware parser, keeps unrelated people (greedy
                kin >= 0.0442, _s05_common) minus the people flagged exclude=True in 59's
                INTERNAL exclude_people.tsv (--exclude-people; the coverage filter),
                classifies every haplotype with 51's rule, and writes PUBLIC aggregates to
                --out-dir (default ~/s05/results/61). Draws the figure too when matplotlib is
                importable; otherwise run --mode figure after pulling the tables.
  --mode local  Builds a DRAFT of the same tables and figure from committed PUBLIC aggregates only:
                  51_kir_structure_qc/gene_content_haplotypes.tsv   (content strings with >= 20
                      complete haplotypes, % per ancestry)  -> a, d (estimated, with bounds)
                  51_kir_structure_qc/cB_share_by_completeness.tsv   -> c (exact)
                  43_kir_full_cohort/kir_allele_freq.tsv             -> b (pooled, known calls,
                      all haplotypes; no per-ancestry allele table is committed)
                  43_kir_full_cohort/kir_gene_by_ancestry.tsv        -> reference column only
                Every draft table carries a `basis` column and the figure is titled DRAFT.
  --mode figure Re-draws the figure from a tables directory written by --mode vm.

Allele resolution (vm mode; both written, figure uses --panel-d-resolution):
  protein     first three digits of the KIR allele name (the protein field). Calls of tier known,
              novel_genomic_known_cds or novel_cds_synonymous keep the protein field of their
              nearest named allele (their CDS/protein is unchanged); novel_protein calls are pooled
              as "novel protein"; undetermined calls are dropped.
  full_known  full allele name, known (exact IPD-KIR match) calls only -- 43's convention.

Disclosure (PUBLIC, the only mode of this script): participant/haplotype counts 1-19 -> "<20";
a rate is blanked when its numerator or denominator is 1-19 (_s05_common.rate_cells); an allele is
named only when >= 20 haplotypes carry it, the rest pooled as "other (<20 each)", the pooled row
itself masked if 1-19. The figure hatches every masked cell. No person_id is written or printed.

Usage
  local draft (laptop, from the repo root):
    python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode local
  VM (from the repo clone; see sprints/S05_qc_phasing_kir_trb/vm_parts/61.sh):
    python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode vm --out-dir ~/s05/results/61 \\
        --exclude-people ~/s05/internal/59/exclude_people.tsv --workers 8
  figure from pulled VM tables:
    python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode figure \\
        --tables-dir reports/hla_popgen/61_kir_by_ancestry_panel

Tests: scripts/hla_popgen/tests/test_61_kir_by_ancestry_panel.py (synthetic GTFs, no AoU data).
"""
import argparse
import importlib.util
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _s05_common as C  # noqa: E402

_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
REPORT_ROOT = os.path.join(_REPO_ROOT, "reports", "hla_popgen")
DEFAULT_LOCAL_OUT = os.path.join(REPORT_ROOT, "61_kir_by_ancestry_panel")
DEFAULT_VM_OUT = os.path.expanduser("~/s05/results/61")
SUPPRESS = 20

# cen -> tel display order (51's KIR_GENES order)
KIR_GENES = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1",
             "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DS3", "KIR2DS5", "KIR2DS1",
             "KIR2DS4", "KIR2DL5B", "KIR3DL2"]
B_CONTENT = {"KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5",
             "KIR3DS1"}
FRAMEWORK = {"KIR3DL3", "KIR2DL4", "KIR3DP1", "KIR3DL2"}
DEFAULT_PANEL_D_GENES = ["KIR2DL1", "KIR2DL4", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR2DS4"]
ANC = C.ANCESTRY_ORDER
NOVEL_PROTEIN_LABEL = "novel protein"


# ---------------------------------------------------------------------------
# Module reuse (41: novelty tiers; 44: contig-aware GTF parse; 51: completeness rule)
# ---------------------------------------------------------------------------
_mods = {}


def _load(filename, modname):
    if modname not in _mods:
        spec = importlib.util.spec_from_file_location(modname, os.path.join(_THIS_DIR, filename))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[modname] = mod
        spec.loader.exec_module(mod)
        _mods[modname] = mod
    return _mods[modname]


def m41():
    return _load("41_kir_pilot.py", "kir_pilot_for_61")


def m44():
    return _load("44_kir_recurrence_saturation.py", "kir_recurrence_for_61")


def m51():
    return _load("51_kir_structure_qc.py", "kir_structure_qc_for_61")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def protein_field(consensus):
    """'KIR2DL1*0030204' -> 'KIR2DL1*003'; 'KIR3DL2*00201new' -> 'KIR3DL2*002'. None if no digits."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    gene, rest = consensus.split("*", 1)
    digits = re.match(r"\d+", rest.replace("new", ""))
    if not digits:
        return None
    return f"{gene}*{digits.group(0)[:3]}"


def allele_label(row, resolution):
    """Allele label of one KIR call at `resolution`, or None when the call is not counted."""
    tier = row["novelty_tier"]
    if resolution == "full_known":
        return row["consensus"] if tier == "known" else None
    if tier in ("known", "novel_genomic_known_cds", "novel_cds_synonymous"):
        return protein_field(row["consensus"])
    if tier == "novel_protein":
        return NOVEL_PROTEIN_LABEL
    return None


def read_exclude_people(path):
    """Reads 59's INTERNAL exclude_people.tsv directly (`# INTERNAL` first line; one row per person
    per region; a person is excluded when its `exclude` column is True). Also accepts any TSV with a
    `person_id` column (+ optional `exclude`; without it every listed person is excluded), or a plain
    list with one person_id per line. Returns a set of str."""
    if not path:
        return set()
    if not os.path.exists(path):
        raise SystemExit(f"FATAL: --exclude-people {path!r} does not exist")
    first = ""
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                first = line.rstrip("\n")
                break
    if "person_id" in first.split("\t"):
        df = pd.read_csv(path, sep="\t", dtype={"person_id": str}, comment="#")
        if "exclude" in df.columns:
            df = df[df["exclude"].astype(str).str.lower().isin(["true", "1", "yes"])]
        return set(df["person_id"].astype(str))
    ids = set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.add(s.split("\t")[0])
    return ids


def rate(k, n, digits=1):
    return C.rate_cells(k, n, "public", digits=digits)


def suppress_partition(counts, floor=SUPPRESS):
    """Complementary suppression for categories that partition a published total: primary cells
    are 1-19; while fewer than two cells are masked or the masked cells sum to < floor (so the
    total minus the shown cells would recover a small count), also mask the smallest shown
    non-zero cell. Returns the set of masked categories (primary + secondary)."""
    masked = {c for c, v in counts.items() if 0 < v < floor}
    if not masked:
        return masked
    while len(masked) < 2 or sum(counts[c] for c in masked) < floor:
        rest = [(v, str(c), c) for c, v in counts.items() if c not in masked and v > 0]
        if not rest:
            break
        masked.add(min(rest)[2])
    return masked


def partition_rates(counts, total, digits=1):
    """rate() per category with complementary suppression; secondary cells show '<20*'."""
    sup = suppress_partition(counts)
    out = {}
    for c, k in counts.items():
        r = rate(k, total, digits)
        if c in sup and not (0 < k < SUPPRESS):
            r = dict(r, n="<20*", pct="", lo="", hi="")
        out[c] = r
    return out


# ---------------------------------------------------------------------------
# VM mode: parse + aggregate
# ---------------------------------------------------------------------------
def parse_person(args):
    """(pid, outroot) -> list of (pid, hap, calls); calls = list of dicts with contig, gene,
    consensus, novelty_tier (KIR genes only). A missing GTF is skipped (not a haplotype)."""
    pid, outroot = args
    k41, k44 = m41(), m44()
    out = []
    for hap in ("hap1", "hap2"):
        path = os.path.join(outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
        if not os.path.exists(path):
            continue
        try:
            rows = k44.parse_hap_gtf_full(path, k41)
        except (OSError, EOFError):
            continue
        calls = [{"contig": r["contig"], "gene": r["gene"], "consensus": r["consensus"],
                  "novelty_tier": r["novelty_tier"]} for r in rows if r["gene"] in KIR_GENES]
        out.append((str(pid), hap, calls))
    return out


def classify(calls):
    """51's completeness rule on our call dicts."""
    return m51().classify_hap([(c["contig"], c["gene"], 0, 0) for c in calls])


def aggregate(haps, ancestry_of, keep_people, top_n=5):
    """haps: list of (pid, hap, calls). keep_people: unrelated, not excluded. Returns dict of
    PUBLIC DataFrames. Only complete haplotypes of kept people enter a-d."""
    recs = []
    for pid, hap, calls in haps:
        if pid not in keep_people:
            continue
        anc = ancestry_of.get(pid, "NA")
        recs.append({"pid": pid, "hap": hap, "anc": anc, "cls": classify(calls), "calls": calls,
                     "genes": frozenset(c["gene"] for c in calls)})
    strata = [(a, [r for r in recs if r["anc"] == a]) for a in ANC] + [("ALL", recs)]
    tabs = {}

    # completeness by ancestry (is the rule ancestry-neutral?)
    rows = []
    for a, rr in strata:
        n = len(rr)
        row = {"ancestry": a, "n_haps": C.mask_n(n)}
        classes = ("complete", "one_contig_open_end", "multi_contig", "empty")
        pr = partition_rates({k: sum(x["cls"] == k for x in rr) for k in classes}, n)
        for k in classes:
            row[f"{k}_pct"] = pr[k]["pct"]
        rows.append(row)
    tabs["completeness_by_ancestry"] = pd.DataFrame(rows)

    comp = {a: [r for r in rr if r["cls"] == "complete"] for a, rr in strata}

    # a. gene frequency (haplotype) + person-level carrier frequency (both haps complete)
    by_person = defaultdict(list)
    for r in recs:
        by_person[r["pid"]].append(r)
    persons_both = {a: [] for a, _ in strata}
    for pid, rr in by_person.items():
        if len(rr) == 2 and all(x["cls"] == "complete" for x in rr):
            g = rr[0]["genes"] | rr[1]["genes"]
            if rr[0]["anc"] in persons_both:
                persons_both[rr[0]["anc"]].append(g)
            persons_both["ALL"].append(g)
    rows = []
    for gene in KIR_GENES:
        for a, _ in strata:
            cc = comp[a]
            k = sum(gene in x["genes"] for x in cc)
            rh = partition_rates({"yes": k, "no": len(cc) - k}, len(cc))["yes"]
            pp = persons_both[a]
            k = sum(gene in g for g in pp)
            rp = partition_rates({"yes": k, "no": len(pp) - k}, len(pp))["yes"]
            rows.append({"gene": gene, "ancestry": a, "n_complete_haps": rh["d"],
                         "n_haps_carrying": rh["n"], "hap_freq_pct": rh["pct"],
                         "hap_freq_lo": rh["lo"], "hap_freq_hi": rh["hi"],
                         "n_people_both_complete": rp["d"], "n_people_carrying": rp["n"],
                         "carrier_freq_pct": rp["pct"], "carrier_freq_lo": rp["lo"],
                         "carrier_freq_hi": rp["hi"], "basis": "vm_complete_haplotypes"})
    tabs["gene_freq_by_ancestry"] = pd.DataFrame(rows)

    # b. cA/cB
    rows = []
    for a, _ in strata:
        cc = comp[a]
        nb = sum(bool(x["genes"] & B_CONTENT) for x in cc)
        pr = partition_rates({"cB": nb, "cA": len(cc) - nb}, len(cc))
        r = pr["cB"]
        rows.append({"ancestry": a, "n_complete_haps": r["d"], "n_cB": r["n"],
                     "n_cA": pr["cA"]["n"], "cB_pct": r["pct"], "cB_lo": r["lo"],
                     "cB_hi": r["hi"], "basis": "vm_complete_haplotypes"})
    tabs["content_by_ancestry"] = pd.DataFrame(rows)

    # c. genes per haplotype
    rows = []
    for a, _ in strata:
        cc = comp[a]
        k = Counter(len(x["genes"]) for x in cc)
        n = len(cc)
        bins = {g: k.get(g, 0) for g in range(min(k, default=0), max(k, default=-1) + 1)}
        pr = partition_rates(bins, n)
        for g in bins:
            r = pr[g]
            rows.append({"ancestry": a, "n_genes": g, "n_haps": r["n"], "n_complete_haps": r["d"],
                         "pct": r["pct"], "basis": "vm_complete_haplotypes"})
        if n >= SUPPRESS:
            v = np.array([len(x["genes"]) for x in cc])
            rows.append({"ancestry": a, "n_genes": "mean", "n_haps": "", "n_complete_haps": str(n),
                         "pct": f"{v.mean():.2f}", "basis": "vm_complete_haplotypes"})
            rows.append({"ancestry": a, "n_genes": "median", "n_haps": "", "n_complete_haps": str(n),
                         "pct": f"{np.median(v):.1f}", "basis": "vm_complete_haplotypes"})
    tabs["genes_per_hap_by_ancestry"] = pd.DataFrame(rows)

    # d. top alleles per gene per ancestry
    rows = []
    for resolution in ("protein", "full_known"):
        for gene in KIR_GENES:
            for a, _ in strata:
                cnt = Counter()
                for x in comp[a]:
                    for c in x["calls"]:
                        if c["gene"] != gene:
                            continue
                        lab = allele_label(c, resolution)
                        if lab:
                            cnt[lab] += 1
                rows += allele_rows(gene, a, cnt, resolution, top_n)
    tabs["top_alleles_by_ancestry"] = pd.DataFrame(rows)

    # run summary
    n_people = len({r["pid"] for r in recs})
    rows = [{"metric": "people_kept_unrelated_not_excluded", "value": C.mask_n(n_people)},
            {"metric": "haplotypes_parsed_kept_people", "value": C.mask_n(len(recs))},
            {"metric": "haplotypes_complete", "value": C.mask_n(len(comp["ALL"]))}]
    for a in ANC:
        rows.append({"metric": f"complete_haplotypes_{a}", "value": C.mask_n(len(comp[a]))})
    tabs["run_summary"] = pd.DataFrame(rows)
    return tabs


def allele_rows(gene, ancestry, cnt, resolution, top_n, basis="vm_complete_haplotypes"):
    """Up to top_n alleles named, each only when >= 20 haplotypes carry it; the novel-protein pool
    is a row of its own when >= 20. Everything else is pooled as "other". If the pooled remainder
    would be 1-19 (recoverable as total minus the shown rows), the smallest shown row is folded
    into it until it is 0 or >= 20."""
    cnt = Counter(cnt)
    total = sum(cnt.values())
    out = []
    if total == 0:
        return out
    named = [(lab, n) for lab, n in cnt.most_common()
             if n >= SUPPRESS and lab not in (NOVEL_PROTEIN_LABEL, "other")]
    shown = named[:top_n]
    if cnt.get(NOVEL_PROTEIN_LABEL, 0) >= SUPPRESS:
        shown.append((NOVEL_PROTEIN_LABEL, cnt[NOVEL_PROTEIN_LABEL]))
    rest = total - sum(n for _, n in shown)
    while 0 < rest < SUPPRESS and shown:
        smallest = min(range(len(shown)), key=lambda i: shown[i][1])
        rest += shown.pop(smallest)[1]
    rank = 0
    for lab, n in shown:
        if lab != NOVEL_PROTEIN_LABEL:
            rank += 1
        r = rate(n, total)
        out.append({"gene": gene, "ancestry": ancestry, "resolution": resolution,
                    "rank": rank if lab != NOVEL_PROTEIN_LABEL else "", "allele": lab,
                    "n_haps": r["n"], "n_calls_gene": r["d"], "freq_pct": r["pct"],
                    "freq_lo": r["lo"], "freq_hi": r["hi"], "basis": basis})
    if rest > 0:
        r = rate(rest, total)
        out.append({"gene": gene, "ancestry": ancestry, "resolution": resolution, "rank": "",
                    "allele": "other", "n_haps": r["n"], "n_calls_gene": r["d"],
                    "freq_pct": r["pct"], "freq_lo": r["lo"], "freq_hi": r["hi"], "basis": basis})
    return out


def run_vm(args):
    root = os.path.expanduser(args.kir_outroot)
    if not os.path.isdir(root):
        raise SystemExit(f"FATAL: KIR outroot {root!r} not found")
    pids = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    if args.limit:
        pids = pids[:args.limit]
    C.log(f"[61] KIR person dirs: {len(pids)}")
    haps = []
    if args.workers > 1:
        with ProcessPoolExecutor(args.workers) as ex:
            for i, recs in enumerate(ex.map(parse_person, [(p, root) for p in pids], chunksize=64)):
                haps.extend(recs)
                if i % 2000 == 0:
                    C.log(f"[61]   parsed {i}/{len(pids)}")
    else:
        for p in pids:
            haps.extend(parse_person((p, root)))
    cohort = C.load_cohort(args.cohort_membership)
    ancestry_of = dict(zip(cohort["person_id"].astype(str), cohort["ancestry_pred"]))
    people = {p for p, _, _ in haps}
    unrelated = C.greedy_unrelated(people, C.load_relatedness_pairs(args.relatedness_table))
    excl = read_exclude_people(args.exclude_people)
    keep = (unrelated - excl) if args.unrelated_only else (people - excl)
    C.log(f"[61] people with KIR haps {len(people)}; unrelated {len(unrelated)}; "
          f"excluded by list {len(people & excl)}; kept {len(keep)}")
    tabs = aggregate(haps, ancestry_of, keep, top_n=args.top_n)
    tabs["run_summary"] = pd.concat([tabs["run_summary"], pd.DataFrame([
        {"metric": "people_with_kir_haplotypes", "value": C.mask_n(len(people))},
        {"metric": "people_excluded_by_list", "value": C.mask_n(len(people & excl))},
        {"metric": "exclude_list_used", "value": bool(args.exclude_people)},
        {"metric": "unrelated_only", "value": bool(args.unrelated_only)},
    ])], ignore_index=True)
    out = os.path.expanduser(args.out_dir)
    for name, df in tabs.items():
        C.write_tsv(df, os.path.join(out, f"{name}.tsv"))
    return tabs, out


# ---------------------------------------------------------------------------
# Local draft from committed PUBLIC aggregates
# ---------------------------------------------------------------------------
def _num(x):
    return pd.to_numeric(x, errors="coerce")


def local_tables(report_root=REPORT_ROOT, panel_d_genes=None, top_n=3):
    """Builds the draft tables from committed PUBLIC files. Returns dict of DataFrames."""
    d51 = os.path.join(report_root, "51_kir_structure_qc")
    d43 = os.path.join(report_root, "43_kir_full_cohort")
    content = pd.read_csv(os.path.join(d51, "gene_content_haplotypes.tsv"), sep="\t", dtype=str)
    cb = pd.read_csv(os.path.join(d51, "cB_share_by_completeness.tsv"), sep="\t", dtype=str)
    allele = pd.read_csv(os.path.join(d43, "kir_allele_freq.tsv"), sep="\t", dtype=str)
    gba = pd.read_csv(os.path.join(d43, "kir_gene_by_ancestry.tsv"), sep="\t", dtype=str)

    comp = cb[cb["scope"] == "complete_haps"].set_index("ancestry")
    n_comp = {a: int(comp.loc[a, "n_haps"]) for a in ANC if a in comp.index and
              str(comp.loc[a, "n_haps"]).isdigit()}
    named = content[content["content_string"].notna() & (content["content_string"] != "")].copy()
    strings = named["content_string"].tolist()
    pct = {a: _num(named[f"{a}_pct"]).to_numpy() for a in ANC}
    shown_mass = {a: float(np.nansum(pct[a])) for a in ANC}
    tabs = {}

    # a. gene frequency on complete haplotypes, from content strings: lower bound = sum of the
    # shown % of strings carrying the gene; upper = lower + unexplained % (blank cells + pooled
    # "other" strings); estimate = lower / shown mass (assumes unexplained mass has the same gene
    # composition as the shown strings).
    rows = []
    for gene in KIR_GENES:
        gi = KIR_GENES.index(gene)
        carry = np.array([s[gi] == "1" for s in strings])
        for a in ANC:
            lo = float(np.nansum(pct[a][carry]))
            unexpl = max(0.0, 100.0 - shown_mass[a])
            est = 100.0 * lo / shown_mass[a] if shown_mass[a] > 0 else np.nan
            if gene in ("KIR3DL3", "KIR3DL2"):  # completeness rule => 100% by definition
                lo, est, unexpl = 100.0, 100.0, 0.0
            n = n_comp.get(a, 0)
            implied = est / 100.0 * n
            masked = n < SUPPRESS or (0 < round(implied) < SUPPRESS) or (0 < round(n - implied) < SUPPRESS)
            ref = gba[(gba["gene"] == gene) & (gba["ancestry"] == a)]
            rows.append({"gene": gene, "ancestry": a, "n_complete_haps": C.mask_n(n),
                         "hap_freq_pct": "" if masked else f"{est:.1f}",
                         "hap_freq_lower_bound": "" if masked else f"{lo:.1f}",
                         "hap_freq_upper_bound": "" if masked else f"{min(100.0, lo + unexpl):.1f}",
                         "all_haps_presence_pct_43": ref["presence_pct"].iloc[0] if len(ref) else "",
                         "basis": "draft_from_51_content_strings_ge20"})
    tabs["gene_freq_by_ancestry"] = pd.DataFrame(rows)

    # b. cA/cB on complete haplotypes, exact from 51
    rows = []
    for a in ANC:
        if a not in comp.index:
            continue
        r = comp.loc[a]
        rows.append({"ancestry": a, "n_complete_haps": r["n_haps"], "cB_pct": r["cB_pct"],
                     "cB_lo": r["lo"], "cB_hi": r["hi"],
                     "cB_pct_all_haps_51": cb[(cb["scope"] == "all_haps") &
                                              (cb["ancestry"] == a)]["cB_pct"].iloc[0],
                     "basis": "exact_from_51_cB_share_by_completeness"})
    tabs["content_by_ancestry"] = pd.DataFrame(rows)

    # c. genes per haplotype (renormalised over shown strings)
    ngenes = np.array([s.count("1") for s in strings])
    rows = []
    for a in ANC:
        n = n_comp.get(a, 0)
        sm = shown_mass[a]
        for g in range(int(ngenes.min()), int(ngenes.max()) + 1):
            sel = pct[a][ngenes == g]
            p = float(np.nansum(sel)) / sm * 100 if sm else np.nan
            implied = p / 100 * n
            # a 0 built only from blank (1-19, masked) cells is not a true zero -> blank it
            masked = (0 < round(implied) < SUPPRESS) or (p == 0 and np.isnan(sel).any())
            rows.append({"ancestry": a, "n_genes": g, "n_complete_haps": C.mask_n(n),
                         "pct": "" if masked else f"{p:.1f}",
                         "basis": "draft_from_51_content_strings_ge20_renormalised"})
        w = np.nan_to_num(pct[a])
        rows.append({"ancestry": a, "n_genes": "mean", "n_complete_haps": C.mask_n(n),
                     "pct": f"{(ngenes * w).sum() / w.sum():.2f}",
                     "basis": "draft_from_51_content_strings_ge20_renormalised"})
        rows.append({"ancestry": a, "n_genes": "unexplained_pct", "n_complete_haps": C.mask_n(n),
                     "pct": f"{100 - sm:.1f}", "basis": "100 - sum of shown content-string %"})
    tabs["genes_per_hap_by_ancestry"] = pd.DataFrame(rows)

    # d. pooled top alleles, protein field summed over named full-resolution rows (known calls,
    # all unrelated haplotypes; 43 does not publish per-ancestry allele frequencies)
    genes = panel_d_genes or DEFAULT_PANEL_D_GENES
    rows = []
    for gene in genes:
        g = allele[(allele["gene"] == gene) & ~allele["allele"].str.startswith("other")]
        if g.empty:
            continue
        total = int(g["n_known_calls_total_gene"].iloc[0])
        cnt = Counter()
        for lab, n in zip(g["allele"], g["n_haplotypes"]):
            if str(n).isdigit():
                cnt[protein_field(lab)] += int(n)
        cnt["other"] += total - sum(cnt.values())  # 43's pooled "<20 each" alleles
        rows += allele_rows(gene, "pooled", cnt, "protein", top_n,
                            basis="draft_from_43_kir_allele_freq_pooled_known_all_haps")
    tabs["top_alleles_by_ancestry"] = pd.DataFrame(rows)

    rows = [{"metric": "complete_haplotypes_unrelated_with_ancestry",
             "value": sum(n_comp.values())}]
    rows += [{"metric": f"complete_haplotypes_{a}", "value": C.mask_n(n_comp.get(a, 0))} for a in ANC]
    rows += [{"metric": f"content_string_mass_shown_pct_{a}", "value": f"{shown_mass[a]:.1f}"}
             for a in ANC]
    rows += [{"metric": "named_content_strings_ge20", "value": len(strings)},
             {"metric": "status", "value": "DRAFT from committed PUBLIC aggregates; VM run pending"}]
    tabs["run_summary"] = pd.DataFrame(rows)
    return tabs


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def load_tables(tables_dir):
    out = {}
    for name in ("gene_freq_by_ancestry", "content_by_ancestry", "genes_per_hap_by_ancestry",
                 "top_alleles_by_ancestry"):
        out[name] = pd.read_csv(os.path.join(tables_dir, f"{name}.tsv"), sep="\t", dtype=str)
    return out


def render_figure(tabs, out_stem, draft, resolution="protein", panel_d_genes=None, top_n=3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import _viz_common as vc

    cmap = LinearSegmentedColormap.from_list("seq", ["#F4F6FA", vc.ACCENT_COLOR])
    genes_d = panel_d_genes or DEFAULT_PANEL_D_GENES
    with vc.nature_style():
        fig = plt.figure(figsize=(vc.mm(183), vc.mm(150)))
        gs = fig.add_gridspec(2, 2, height_ratios=[2.3, 1], width_ratios=[1, 1.15],
                              hspace=0.42, wspace=0.55, left=0.09, right=0.97, top=0.92,
                              bottom=0.08)
        title = "KIR gene content and alleles by ancestry, complete haplotypes"
        if draft:
            title = ("DRAFT (from committed aggregates; VM run pending) -- " + title)
        fig.suptitle(title, fontsize=7, x=0.09, ha="left")

        # a. gene frequency heatmap
        ax = fig.add_subplot(gs[0, 0])
        gf = tabs["gene_freq_by_ancestry"]
        M = np.full((len(KIR_GENES), len(ANC)), np.nan)
        for i, g in enumerate(KIR_GENES):
            for j, a in enumerate(ANC):
                v = gf[(gf["gene"] == g) & (gf["ancestry"] == a)]["hap_freq_pct"]
                if len(v):
                    M[i, j] = _num(v.iloc[0])
        ax.imshow(np.nan_to_num(M, nan=0), cmap=cmap, vmin=0, vmax=100, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    vc.hatch_suppressed(ax, j - 0.5, i - 0.5, 1, 1)
                else:
                    ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=4.5,
                            color="white" if M[i, j] > 60 else "#222222")
        ax.set_xticks(range(len(ANC)))
        ax.set_xticklabels(ANC)
        ax.set_yticks(range(len(KIR_GENES)))
        ax.set_yticklabels([g.replace("KIR", "") for g in KIR_GENES])
        for t, g in zip(ax.get_yticklabels(), KIR_GENES):
            if g in FRAMEWORK:
                t.set_fontweight("bold")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("% of complete haplotypes carrying the gene"
                     + (" (estimated)" if draft else ""), loc="left")
        vc.panel_letter(ax, "a", dx=-0.2, dy=1.02)

        # b. top alleles
        ax = fig.add_subplot(gs[0, 1])
        ta = tabs["top_alleles_by_ancestry"]
        ta = ta[ta["resolution"] == resolution] if "resolution" in ta.columns else ta
        cols = ["pooled"] if draft else ANC
        rows = []
        for g in genes_d:
            sub = ta[(ta["gene"] == g)]
            ref = sub[sub["ancestry"] == ("pooled" if draft else "ALL")]
            ref = ref[~ref["allele"].str.startswith("other") & (ref["allele"] != NOVEL_PROTEIN_LABEL)]
            for lab in ref["allele"].head(top_n):
                rows.append((g, lab))
        D = np.full((len(rows), len(cols)), np.nan)
        for i, (g, lab) in enumerate(rows):
            for j, a in enumerate(cols):
                v = ta[(ta["gene"] == g) & (ta["ancestry"] == a) & (ta["allele"] == lab)]["freq_pct"]
                D[i, j] = _num(v.iloc[0]) if len(v) else np.nan
        if len(rows):
            ax.imshow(np.nan_to_num(D, nan=0), cmap=cmap, vmin=0, vmax=max(60, np.nanmax(D) if
                      np.isfinite(D).any() else 60), aspect="auto")
            for i in range(D.shape[0]):
                for j in range(D.shape[1]):
                    if np.isnan(D[i, j]):
                        vc.hatch_suppressed(ax, j - 0.5, i - 0.5, 1, 1)
                    else:
                        ax.text(j, i, f"{D[i, j]:.0f}", ha="center", va="center", fontsize=4.5,
                                color="white" if D[i, j] > 36 else "#222222")
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels([lab.replace("KIR", "") for _, lab in rows])
            for k in range(1, len(rows)):
                if rows[k][0] != rows[k - 1][0]:
                    ax.axhline(k - 0.5, color="white", lw=1.2)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(["pooled, all haps" if c == "pooled" else c for c in cols])
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        lab_res = "protein field" if resolution == "protein" else "full name, known calls"
        ax.set_title(f"Top alleles, % of the gene's calls ({lab_res})"
                     + ("\nknown calls, all haplotypes, pooled ancestries" if draft else ""),
                     loc="left")
        vc.panel_letter(ax, "b", dx=-0.32, dy=1.02)

        # c. cB share
        ax = fig.add_subplot(gs[1, 0])
        cbt = tabs["content_by_ancestry"].set_index("ancestry")
        for j, a in enumerate(ANC):
            if a not in cbt.index or not str(cbt.loc[a, "cB_pct"]).strip() or \
                    pd.isna(_num(cbt.loc[a, "cB_pct"])):
                vc.hatch_suppressed(ax, j - 0.35, 0, 0.7, 100)
                continue
            p, lo, hi = (float(_num(cbt.loc[a, k])) for k in ("cB_pct", "cB_lo", "cB_hi"))
            ax.bar(j, p, width=0.7, color=vc.ANCESTRY_COLORS[a])
            ax.errorbar(j, p, yerr=[[p - lo], [hi - p]], color="#333333", lw=0.6, capsize=1.5)
            ax.text(j, hi + 2, f"{p:.0f}", ha="center", va="bottom", fontsize=5)
        ax.set_xticks(range(len(ANC)))
        ax.set_xticklabels(ANC)
        ax.set_ylim(0, 80)
        ax.set_ylabel("% cB haplotypes")
        ax.set_title("cB share (rest cA), 95% CI", loc="left")
        vc.panel_letter(ax, "c", dx=-0.2, dy=1.02)

        # d. genes per haplotype
        ax = fig.add_subplot(gs[1, 1])
        gp = tabs["genes_per_hap_by_ancestry"]
        gp = gp[gp["n_genes"].astype(str).str.isdigit()]
        for a in ANC:
            s = gp[gp["ancestry"] == a]
            x = s["n_genes"].astype(int).to_numpy()
            y = _num(s["pct"]).to_numpy()
            ax.plot(x, y, marker="o", ms=2, lw=0.8, color=vc.ANCESTRY_COLORS[a], label=a)
        ax.set_xlabel("distinct KIR genes per haplotype (of 17)")
        ax.set_ylabel("% of complete haplotypes")
        ax.set_title("Genes per haplotype" + (" (shown strings, renormalised)" if draft else ""),
                     loc="left")
        ax.legend(ncol=3, loc="upper right", **vc.LEGEND_KW)
        vc.panel_letter(ax, "d", dx=-0.2, dy=1.02)

        fig.text(0.09, 0.012, "Hatched cell / gap in line: masked (n<20) or not available. Bold genes: framework. "
                 "Complete haplotype = all KIR calls on one contig with 3DL3 and 3DL2 present.",
                 fontsize=5, ha="left")
        return vc.save_fig(fig, out_stem)


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["vm", "local", "figure"], required=True)
    ap.add_argument("--kir-outroot", default=C.DEFAULT_KIR_OUTROOT)
    ap.add_argument("--cohort-membership", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--exclude-people", default=None,
                    help="59's INTERNAL exclude_people.tsv (~/s05/internal/59/exclude_people.tsv); "
                         "people with exclude=True are dropped. Also accepts a person_id TSV or a "
                         "plain id list. Optional.")
    ap.add_argument("--all-people", dest="unrelated_only", action="store_false",
                    help="keep relatives too (default: unrelated only)")
    ap.add_argument("--report-root", default=REPORT_ROOT, help="local mode: committed reports root")
    ap.add_argument("--tables-dir", default=None, help="figure mode: where the vm tables are")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--top-n", type=int, default=5, help="alleles named per gene x ancestry (tables)")
    ap.add_argument("--panel-d-top", type=int, default=3)
    ap.add_argument("--panel-d-genes", default=",".join(DEFAULT_PANEL_D_GENES))
    ap.add_argument("--panel-d-resolution", choices=["protein", "full_known"], default="protein")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)
    genes_d = [g if g.startswith("KIR") else "KIR" + g for g in args.panel_d_genes.split(",") if g]

    if args.mode == "vm":
        args.out_dir = args.out_dir or DEFAULT_VM_OUT
        tabs, out = run_vm(args)
        if not args.no_figure:
            try:
                render_figure(tabs, os.path.join(out, "fig_kir_by_ancestry_panel"), draft=False,
                              resolution=args.panel_d_resolution, panel_d_genes=genes_d,
                              top_n=args.panel_d_top)
            except ImportError as e:
                C.log(f"[61] figure skipped ({e}); run --mode figure on the pulled tables")
    elif args.mode == "local":
        out = args.out_dir or DEFAULT_LOCAL_OUT
        tabs = local_tables(args.report_root, panel_d_genes=genes_d, top_n=args.panel_d_top)
        for name, df in tabs.items():
            C.write_tsv(df, os.path.join(out, f"draft_{name}.tsv"))
        if not args.no_figure:
            render_figure(tabs, os.path.join(out, "fig_kir_by_ancestry_panel_DRAFT"), draft=True,
                          resolution="protein", panel_d_genes=genes_d, top_n=args.panel_d_top)
    else:
        if not args.tables_dir:
            ap.error("--mode figure needs --tables-dir")
        tabs = load_tables(args.tables_dir)
        out = args.out_dir or args.tables_dir
        render_figure(tabs, os.path.join(out, "fig_kir_by_ancestry_panel"), draft=False,
                      resolution=args.panel_d_resolution, panel_d_genes=genes_d,
                      top_n=args.panel_d_top)
    C.log("[61] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
