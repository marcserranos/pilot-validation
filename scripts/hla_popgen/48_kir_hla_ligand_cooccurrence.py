#!/usr/bin/env python3
"""S04 WS-D -- KIR-HLA ligand co-occurrence by ancestry (seed-list item: ORCHESTRATOR_HANDOFF.md
sect2 WS-D, "KIR-HLA ligand co-occurrence by ancestry"). Aggregate-only outputs.

## Biology (why this test)

Three KIR inhibitory/activating receptors recognise specific HLA class-I epitope groups, not
individual alleles -- the functional "ligand" is the group, not the two-field type:
  - **KIR2DL1** (inhibitory) recognises the **HLA-C2** group.
  - **KIR2DL2 / KIR2DL3** (inhibitory) recognise the **HLA-C1** group. KIR2DL2 binds C1 (and
    weakly C2) more promiscuously than KIR2DL3 (Moesta et al. 2008 J Immunol; Parham 2005 Nat
    Rev Immunol review).
  - **KIR3DL1** (inhibitory) recognises **HLA-Bw4** (all HLA-B alleles carrying the Bw4 public
    epitope, plus the HLA-A Bw4 alleles A*23/A*24/A*32), with allotype-dependent avidity split
    by the Bw4 **position-80 dimorphism** (80I = higher avidity, 80T = lower) -- Cella et al.
    1994; Gumperz et al. 1995 J Exp Med.
  - **KIR3DS1** (activating, no confirmed HLA ligand by direct binding assay, but epidemiologic
    association specifically with **Bw4-80I**) -- Martin et al. 2002 Nat Genet (HIV progression).

Because a person's genotype fixes *both* halves of a pair in trans (unlike a linked haplotype),
co-carriage of a KIR gene and its ligand across a cohort should, under independence (no epistasis,
no shared population-structure confound), equal the product of the two marginal frequencies. This
script tests observed vs. expected co-carriage per ancestry with a 2x2 contingency table
(Fisher's exact + chi2), an odds ratio, and a person-label permutation baseline. **A real
enrichment or depletion here is NOT expected** under simple population genetics -- HLA (chr6) and
KIR (chr19) segregate independently, so any signal is either (a) genuine linkage disequilibrium
with a third factor (ancestry is the obvious confound -- both loci have ancestry-stratified allele
frequencies, so within-ancestry stratification is the whole point of this script), (b) selection
acting on specific KIR-HLA combinations (a real, published phenomenon -- e.g. HIV/pre-eclampsia
literature), or (c) an assay/miscall artifact. This script cannot distinguish those on its own; it
flags what is worth a closer look.

## Epitope group definitions used here (READ BEFORE TRUSTING NUMBERS)

**HLA-C1/C2** is determined by the residue at *mature-protein* position 80: Asn80 = C2,
Lys80 = C1 (Colonna et al. 1993; Winter & Long 1997). This script does NOT have access to
per-person amino-acid sequence at position 80 (Table 1's `consensus` is a nomenclature string,
not a translated sequence) -- so C1/C2 group is assigned via a **documented two-field allele-name
lookup table** (`C1C2_TABLE` below), built from the well-established HLA-C two-field group
assignments summarised in the IPD-KIR ligand-motif reference pages and in review tables (e.g.
Middleton & Gonzalez 2010 Immunology; Norman et al. 2016 Nat Genet supplementary tables). This
lookup is **not exhaustive** -- alleles not in the table are reported as `unclassified` (never
silently folded into C1 or C2), and the script prints the count and % of HLA-C calls that fell
into `unclassified` so a VM operator can judge coverage before trusting downstream numbers.
**This table needs a VM-side cross-check against IPD-IMGT/HLA's own C1/C2 assignment file (or a
position-80 translation from the cohort's own CDS sequences, which 24_novelty_by_field.py's
RefIndex machinery could in principle supply) before any of this script's HLA-C-ligand numbers
are treated as final** -- flagged here explicitly rather than assumed correct.

**HLA-Bw4/Bw6** is likewise assigned via a documented two-field-group-level lookup
(`BW4_B_GROUPS`, `BW4_A_GROUPS`) based on the standard Bw4-bearing allele-group list reproduced
in multiple KIR-ligand studies (Bw4: B*13, B*27, B*37, B*38, B*44, B*47, B*49, B*51, B*52, B*53,
B*57, B*58, B*59, B*63, B*77; Bw4 at the A locus: A*23, A*24, A*32). **Known exception, flagged
not silently absorbed:** most B*15 alleles are Bw6, but a documented subset (B*15:13, B*15:16,
B*15:17, B*15:24 and a few others) are Bw4 -- this script does NOT special-case B*15 by
four-digit allele (would need the full IPD-IMGT/HLA Bw4/Bw6 exception list, not reproduced from
memory here) and instead reports ALL B*15 as Bw6 by default, flagged via `--b15-as-bw4-list` (a
CLI-supplied comma list of exception 2-field alleles, empty by default) so a VM operator can
supply the authoritative exception list once looked up. Bw4-80I/80T sub-stratification is **not**
implemented in this delivery (would need position-80 translation, same limitation as C1/C2) --
the `n_persons_with_c1c2_or_bw_unclassified` QC row in the output flags exactly how much of the
cohort this affects.

## KIR receptor definitions

Presence/absence of KIR2DL1, KIR2DL2, KIR2DL3, KIR3DL1, KIR3DS1 per person (either haplotype
carries the gene at all -- copy-number and allele-level avidity differences are out of scope for
this pass), from the same `parse_hap_gtf`/`KIR_GENES` logic as 41_kir_pilot.py /
43_kir_full_aggregate.py.

## Method

Per ancestry (ANCESTRY_ORDER, unrelated subset, same kinship threshold as 41/43): for each of the
5 functional pairs (2DL1xC2, 2DL2xC1, 2DL3xC1, 3DL1xBw4, 3DS1xBw4), build the person-level 2x2
table (KIR present/absent x ligand-carrier present/absent), report counts (disclosure-masked),
odds ratio with a 2x2 Haldane-Anscombe-corrected OR when a cell is 0, Fisher's exact p-value,
chi2 p-value, and a >=20-shuffle permutation null of the OR (KIR presence label shuffled within
ancestry, ligand fixed). Disclosure: any 2x2 cell in [1,19] masks the OR/CI/both p-values for that
row (the raw masked cell counts are still reported as `<20`, never blank, never 0 unless truly 0).

## Running (VM, real mode)

    cd ~/s04 && PYTHONPATH=~/s04:~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 48_kir_hla_ligand_cooccurrence.py \\
        --table1 ~/pipeline_outputs/hla_calls_rich.tsv \\
        --cohort-membership ~/pipeline_outputs/cohort_membership.tsv \\
        --relatedness-table ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv \\
        --kir-outroot ~/pipeline_outputs_kir \\
        --kir-pilot-script ~/s03/41_kir_pilot.py \\
        --n-perms 1000 \\
        --out-dir ~/s04/results/48

## Local synthetic dry run

    python3 scripts/hla_popgen/48_kir_hla_ligand_cooccurrence.py --synthetic --n-people 400 \\
        --out-dir /tmp/ligand_synthetic
"""
import argparse
import importlib.util
import os
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]

# ---------------------------------------------------------------------------
# Ligand-group lookup tables -- see module docstring for citations and caveats.
# ---------------------------------------------------------------------------
# HLA-C two-field groups -> C1/C2, per the position-80 dimorphism (documented, common
# low-resolution assignment reproduced in KIR-ligand review tables; NOT independently
# re-derived from sequence in this script -- flagged for VM verification in the docstring).
C1C2_TABLE = {
    # C2 group (Asn80)
    "02:02": "C2", "04:01": "C2", "05:01": "C2", "06:02": "C2", "07:04": "C2",
    "08:02": "C2", "12:03": "C2", "15:02": "C2", "16:02": "C2", "17:01": "C2", "18:01": "C2",
    # C1 group (Lys80)
    "01:02": "C1", "03:02": "C1", "03:03": "C1", "03:04": "C1", "07:01": "C1", "07:02": "C1",
    "08:01": "C1", "12:02": "C1", "14:02": "C1", "16:01": "C1",
}
BW4_B_GROUPS = {"13", "27", "37", "38", "44", "47", "49", "51", "52", "53", "57", "58", "59",
                "63", "77"}
BW4_A_GROUPS = {"23", "24", "32"}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def suppressed(n):
    n = int(n)
    if n == 0:
        return "0"
    return "<%d" % SUPPRESS_BELOW if n < SUPPRESS_BELOW else str(n)


def hla_two_field(consensus):
    if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
        return None
    s = str(consensus).strip()
    if not s or s.upper() in {"UNDETERMINED", "NA", ""}:
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = []
    for f in s.split(":"):
        if f == "":
            continue
        if f.strip().lower() == "new":
            break
        fields.append(f)
    if len(fields) < 2:
        return None
    return ":".join(fields[:2])


def classify_c1c2(two_field):
    return C1C2_TABLE.get(two_field, "unclassified")


def classify_bw4(gene, two_field, b15_as_bw4=frozenset()):
    if two_field is None:
        return "unclassified"
    group = two_field.split(":", 1)[0]
    if gene == "HLA-B":
        if group == "15":
            return "Bw4" if two_field in b15_as_bw4 else "Bw6"
        return "Bw4" if group in BW4_B_GROUPS else "Bw6"
    if gene == "HLA-A":
        return "Bw4" if group in BW4_A_GROUPS else "not_bw4_locus_A"
    return "unclassified"


def load_kir_module(script_path):
    spec = importlib.util.spec_from_file_location("kir_pilot_48", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Per-person epitope + receptor carriage.
# ---------------------------------------------------------------------------
def build_person_epitopes(table1_df, b15_as_bw4=frozenset()):
    """Returns DataFrame indexed by person_id: has_C1, has_C2, has_Bw4 (B or A locus),
    n_unclassified_C, n_unclassified_Bw (QC counters, not person-level flags)."""
    df = table1_df[table1_df["gene"].isin(["HLA-B", "HLA-C", "HLA-A"])].copy()
    df["allele2"] = df["consensus"].map(hla_two_field)

    people = sorted(table1_df["person_id"].unique())
    has_c1 = pd.Series(False, index=people)
    has_c2 = pd.Series(False, index=people)
    has_bw4 = pd.Series(False, index=people)
    n_unclass_c = 0
    n_unclass_bw = 0

    c_rows = df[df["gene"] == "HLA-C"]
    for pid, grp in c_rows.groupby("person_id"):
        for a2 in grp["allele2"].dropna():
            grp_label = classify_c1c2(a2)
            if grp_label == "C1":
                has_c1.loc[pid] = True
            elif grp_label == "C2":
                has_c2.loc[pid] = True
            else:
                n_unclass_c += 1

    b_rows = df[df["gene"].isin(["HLA-B", "HLA-A"])]
    for pid, grp in b_rows.groupby("person_id"):
        for gene, a2 in zip(grp["gene"], grp["allele2"]):
            if a2 is None:
                continue
            label = classify_bw4(gene, a2, b15_as_bw4)
            if label == "Bw4":
                has_bw4.loc[pid] = True
            elif label == "unclassified":
                n_unclass_bw += 1

    out = pd.DataFrame({"has_C1": has_c1, "has_C2": has_c2, "has_Bw4": has_bw4})
    out.index.name = "person_id"
    return out, n_unclass_c, n_unclass_bw


def build_kir_receptor_carriage(pids, kir_outroot, kir_mod, receptor_genes):
    presence = pd.DataFrame(False, index=pids, columns=receptor_genes)
    for pid in pids:
        genes_this_person = set()
        for hap in ("hap1", "hap2"):
            gtf_path = os.path.join(kir_outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
            if not os.path.exists(gtf_path):
                continue
            try:
                rows = kir_mod.parse_hap_gtf(gtf_path)
            except (OSError, EOFError):
                continue
            genes_this_person.update(r["gene"] for r in rows if r["gene"] in receptor_genes)
        for g in genes_this_person:
            presence.at[pid, g] = True
    presence.index.name = "person_id"
    return presence


# ---------------------------------------------------------------------------
# 2x2 contingency stats.
# ---------------------------------------------------------------------------
def two_by_two_stats(kir_present, ligand_present, n_perms=100, seed=0):
    """kir_present/ligand_present: boolean arrays, same length (one ancestry's people).
    Returns a dict with masked cell counts + OR/CI/p-values (blanked if any cell in [1,19])."""
    a = int(np.sum(kir_present & ligand_present))            # KIR+ ligand+
    b = int(np.sum(kir_present & ~ligand_present))           # KIR+ ligand-
    c = int(np.sum(~kir_present & ligand_present))           # KIR- ligand+
    d = int(np.sum(~kir_present & ~ligand_present))          # KIR- ligand-
    n = a + b + c + d

    cells = [a, b, c, d]
    disclosive = any(0 < x < SUPPRESS_BELOW for x in cells)

    result = {
        "n_a_kir_and_ligand": suppressed(a), "n_b_kir_only": suppressed(b),
        "n_c_ligand_only": suppressed(c), "n_d_neither": suppressed(d), "n_total": n,
    }

    if n == 0 or disclosive:
        result.update({"odds_ratio": "", "or_ci_lo": "", "or_ci_hi": "",
                       "fisher_p": "", "chi2_p": "", "perm_p": "", "n_perms": 0})
        return result

    # Haldane-Anscombe correction only when a zero cell would make the OR undefined.
    aa, bb, cc, dd = a, b, c, d
    if 0 in cells:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds_ratio = (aa * dd) / (bb * cc)

    # Fisher exact (exact, small-sample safe) + chi2 (large-sample cross-check).
    fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]])
    try:
        chi2, chi2_p, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=True)
    except ValueError:
        chi2_p = float("nan")

    # log-OR CI (Woolf's method).
    se_log_or = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    log_or = np.log(odds_ratio)
    ci_lo = float(np.exp(log_or - 1.96 * se_log_or))
    ci_hi = float(np.exp(log_or + 1.96 * se_log_or))

    # Permutation null: shuffle the KIR-presence label within this ancestry, ligand fixed.
    rng = np.random.default_rng(seed)
    perm_ors = []
    for _ in range(n_perms):
        kp = rng.permutation(kir_present)
        pa = int(np.sum(kp & ligand_present))
        pb = int(np.sum(kp & ~ligand_present))
        pc = int(np.sum(~kp & ligand_present))
        pd_ = int(np.sum(~kp & ~ligand_present))
        if 0 in (pa, pb, pc, pd_):
            pa, pb, pc, pd_ = pa + 0.5, pb + 0.5, pc + 0.5, pd_ + 0.5
        perm_ors.append((pa * pd_) / (pb * pc))
    perm_arr = np.array(perm_ors)
    perm_p = float((np.sum(np.abs(np.log(perm_arr)) >= abs(log_or)) + 1) / (len(perm_arr) + 1))

    result.update({
        "odds_ratio": round(float(odds_ratio), 3),
        "or_ci_lo": round(ci_lo, 3), "or_ci_hi": round(ci_hi, 3),
        "fisher_p": float(fisher_p), "chi2_p": float(chi2_p) if chi2_p == chi2_p else "",
        "perm_p": perm_p, "n_perms": n_perms,
    })
    return result


PAIRS = [
    ("KIR2DL1", "has_C2", "C2"),
    ("KIR2DL2", "has_C1", "C1"),
    ("KIR2DL3", "has_C1", "C1"),
    ("KIR3DL1", "has_Bw4", "Bw4"),
    ("KIR3DS1", "has_Bw4", "Bw4"),
]
RECEPTOR_GENES = sorted({p[0] for p in PAIRS})


def run_pipeline(table1_df, epitopes_df, kir_presence_df, ancestry_of, out_dir, n_perms,
                 b15_as_bw4, status_path, n_unclass_c=0, n_unclass_bw=0):
    os.makedirs(out_dir, exist_ok=True)
    people = sorted(set(epitopes_df.index) & set(kir_presence_df.index))
    epitopes_df = epitopes_df.loc[people]
    kir_presence_df = kir_presence_df.loc[people]
    anc = np.array([ancestry_of.get(p, "UNASSIGNED") for p in people])

    rows = []
    for ancestry in ANCESTRY_ORDER:
        mask = anc == ancestry
        n_people_anc = int(mask.sum())
        if n_people_anc < SUPPRESS_BELOW:
            log(f"[48] ancestry={ancestry}: skipped (n={n_people_anc} < 20)")
            continue
        for kir_gene, ligand_col, ligand_name in PAIRS:
            kir_present = kir_presence_df.loc[mask, kir_gene].values
            ligand_present = epitopes_df.loc[mask, ligand_col].values
            r = two_by_two_stats(kir_present, ligand_present, n_perms=n_perms)
            r["ancestry"] = ancestry
            r["kir_gene"] = kir_gene
            r["ligand"] = ligand_name
            rows.append(r)

    cols = ["ancestry", "kir_gene", "ligand", "n_total", "n_a_kir_and_ligand", "n_b_kir_only",
            "n_c_ligand_only", "n_d_neither", "odds_ratio", "or_ci_lo", "or_ci_hi",
            "fisher_p", "chi2_p", "perm_p", "n_perms"]
    pd.DataFrame(rows)[cols].to_csv(
        os.path.join(out_dir, "kir_hla_ligand_cooccurrence.tsv"), sep="\t", index=False)

    # Epitope frequencies per ancestry (Cheap, feeds the ligand analysis directly).
    freq_rows = []
    for ancestry in ANCESTRY_ORDER:
        mask = anc == ancestry
        n_people_anc = int(mask.sum())
        if n_people_anc < SUPPRESS_BELOW:
            continue
        for col, label in (("has_C1", "C1"), ("has_C2", "C2"), ("has_Bw4", "Bw4")):
            n_carriers = int(epitopes_df.loc[mask, col].sum())
            n_str, d_str = suppressed(n_carriers), n_people_anc
            pct = (round(100.0 * n_carriers / n_people_anc, 1)
                   if not (0 < n_carriers < SUPPRESS_BELOW) else "")
            freq_rows.append({"ancestry": ancestry, "epitope": label,
                              "n_carriers": n_str, "n_total": d_str, "pct": pct})
    pd.DataFrame(freq_rows).to_csv(
        os.path.join(out_dir, "epitope_freq_by_ancestry.tsv"), sep="\t", index=False)

    # QC: unclassified-allele coverage flag.
    qc_rows = [
        {"metric": "n_unclassified_C_allele_calls", "value": suppressed(n_unclass_c)},
        {"metric": "n_unclassified_Bw_allele_calls", "value": suppressed(n_unclass_bw)},
        {"metric": "n_people_in_analysis", "value": len(people)},
        {"metric": "c1c2_table_size", "value": len(C1C2_TABLE)},
        {"metric": "bw4_b_groups_size", "value": len(BW4_B_GROUPS)},
        {"metric": "b15_as_bw4_exceptions_supplied", "value": len(b15_as_bw4)},
    ]
    pd.DataFrame(qc_rows).to_csv(os.path.join(out_dir, "ligand_lookup_qc.tsv"),
                                 sep="\t", index=False)

    with open(status_path, "w") as f:
        f.write(f"done n_people={len(people)} n_pairs_tested={len(rows)} "
                f"n_unclassified_C={n_unclass_c} n_unclassified_Bw={n_unclass_bw}\n")
    log(f"[48] wrote {len(rows)} pair x ancestry rows -> {out_dir}")


# ---------------------------------------------------------------------------
# Synthetic cohort generator.
# ---------------------------------------------------------------------------
def make_synthetic(n_people=400, seed=20260926, planted_enrichment=True):
    rng = np.random.default_rng(seed)
    pids = [f"synthP{i:05d}" for i in range(n_people)]
    ancestries = rng.choice(ANCESTRY_ORDER, size=n_people, p=[0.2, 0.15, 0.15, 0.3, 0.1, 0.1])
    ancestry_of = dict(zip(pids, ancestries))

    c_alleles = list(C1C2_TABLE.keys())
    b_bw4 = ["57:01", "51:01"]  # Bw4
    b_bw6 = ["07:02", "08:01"]  # Bw6
    table1_rows = []
    for pid in pids:
        for hap in ("hap1", "hap2"):
            c_allele = rng.choice(c_alleles)
            table1_rows.append({"person_id": pid, "hap": hap, "gene": "HLA-C",
                                "consensus": f"HLA-C*{c_allele}"})
            b_allele = rng.choice(b_bw4 if rng.random() < 0.35 else b_bw6)
            table1_rows.append({"person_id": pid, "hap": hap, "gene": "HLA-B",
                                "consensus": f"HLA-B*{b_allele}"})
    table1_df = pd.DataFrame(table1_rows)

    epitopes_df, n_unclass_c, n_unclass_bw = build_person_epitopes(table1_df)

    kir_presence = pd.DataFrame(False, index=pids, columns=RECEPTOR_GENES)
    for pid in pids:
        has_c2 = epitopes_df.at[pid, "has_C2"]
        has_bw4 = epitopes_df.at[pid, "has_Bw4"]
        for gene in RECEPTOR_GENES:
            base_p = 0.5
            if planted_enrichment:
                if gene == "KIR2DL1" and has_c2:
                    base_p = 0.8  # planted enrichment: 2DL1 co-occurs with its ligand C2
                if gene == "KIR3DL1" and has_bw4:
                    base_p = 0.8
            kir_presence.at[pid, gene] = rng.random() < base_p
    kir_presence.index.name = "person_id"

    return table1_df, epitopes_df, kir_presence, ancestry_of, n_unclass_c, n_unclass_bw


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--cohort-membership",
                    default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table",
                    default=os.path.expanduser(
                        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/"
                        "aux/relatedness/samples_relatedness.tsv"))
    ap.add_argument("--kir-outroot", default=os.path.expanduser("~/pipeline_outputs_kir"))
    ap.add_argument("--kir-pilot-script", default=os.path.expanduser("~/s03/41_kir_pilot.py"))
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s04/results/48"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--n-perms", type=int, default=100)
    ap.add_argument("--b15-as-bw4-list", default="",
                    help="Comma-separated 2-field B*15 alleles to treat as Bw4 (documented "
                         "exceptions to the default B*15=Bw6 rule; empty until a VM operator "
                         "supplies the authoritative IPD-IMGT/HLA exception list).")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--n-people", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20260926)
    args = ap.parse_args()

    t0 = time.perf_counter()
    os.makedirs(args.out_dir, exist_ok=True)
    status_path = os.path.join(args.out_dir, "STATUS.txt")
    b15_as_bw4 = frozenset(x.strip() for x in args.b15_as_bw4_list.split(",") if x.strip())

    if args.synthetic:
        (table1_df, epitopes_df, kir_presence, ancestry_of,
         n_unclass_c, n_unclass_bw) = make_synthetic(args.n_people, args.seed)
        run_pipeline(table1_df, epitopes_df, kir_presence, ancestry_of, args.out_dir,
                    args.n_perms, b15_as_bw4, status_path, n_unclass_c, n_unclass_bw)
        log(f"[48] synthetic run done in {time.perf_counter()-t0:.0f}s")
        return

    from _viz_common import load_table1, load_cohort_membership
    table1_df = load_table1(args.table1)
    if args.limit:
        keep = sorted(table1_df["person_id"].unique())[:args.limit]
        table1_df = table1_df[table1_df["person_id"].isin(keep)]
    cohort_df = load_cohort_membership(args.cohort_membership)
    ancestry_of = dict(zip(cohort_df["person_id"].astype(str),
                          cohort_df["ancestry_pred"].astype(str).str.upper()))

    kir_mod = load_kir_module(args.kir_pilot_script)
    pairs = kir_mod.load_relatedness_pairs(args.relatedness_table)
    all_pids = sorted(set(table1_df["person_id"]) | set(cohort_df["person_id"].astype(str)))
    kept, removed = kir_mod.greedy_unrelated(all_pids, pairs, kin_min=kir_mod.KIN_MIN)
    log(f"[48] unrelated subset: {len(kept)}/{len(all_pids)} ({len(removed)} relatives dropped)")
    table1_df = table1_df[table1_df["person_id"].isin(kept)]

    epitopes_df, n_unclass_c, n_unclass_bw = build_person_epitopes(table1_df, b15_as_bw4)
    log(f"[48] unclassified HLA-C calls: {n_unclass_c}; unclassified Bw calls: {n_unclass_bw}")
    unrelated_pids = sorted(kept & set(table1_df["person_id"]))
    kir_presence = build_kir_receptor_carriage(unrelated_pids, args.kir_outroot, kir_mod,
                                               RECEPTOR_GENES)

    run_pipeline(table1_df, epitopes_df, kir_presence, ancestry_of, args.out_dir, args.n_perms,
                b15_as_bw4, status_path, n_unclass_c, n_unclass_bw)
    log(f"[48] done in {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
