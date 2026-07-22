#!/usr/bin/env python3
"""Diagnose + cross-method-compare the 60-person Immuannot pilot (2026-07-22, Marc's questions
after the run: how do the failures distribute -- by person? by gene? by stage? -- how does the
time distribute, and how does Immuannot compare to SpecImmune-LR per gene and per ancestry).

Reads only files already on the VM (no re-run needed for any of this):
  ~/pipeline_outputs/experiment_d/cohort.tsv    -- person_id, ancestry (the attempted 60)
  ~/pipeline_outputs/immuannot_calls.tsv        -- person_id, gene(HLA-*), immuannot_1/2  (calls)
  ~/pipeline_outputs/immuannot_timing.tsv       -- per-hap stage timings; NULLs mark where a hap
                                                    stopped, so failures are locatable by stage
  ~/pipeline_outputs/<pid>/comparison_log.csv   -- Experiment D's SpecImmune-LR calls (+ AoU/SpecHLA)

METHODOLOGY -- NO GROUND TRUTH (Marc, 2026-07-22): this is an evaluation of two callers against
each other, not against truth. Immuannot (assembly-based) and SpecImmune-LR (read-based) both
derive from the SAME long-read sequencing, just different data products (a phased assembly vs the
raw reads) -- so they are MORE independent than AoU-native+SpecHLA were (which shared short-read
data, DECISIONS.md) but NOT fully independent. Therefore where they disagree we CANNOT say which is
right from concordance alone -- exactly the same adjudication limit DECISIONS.md already logs for
the 3-way bake-off. What concordance *can* do: show WHERE the two methods diverge (which genes,
which ancestries), which is diagnostic of where each struggles and cross-references what we already
know (SpecImmune's DQA1/DRB1 behavior; AoU's DQA1 outlier finding; the field-cascade results). The
ground-truth-anchored answer is Aleix's separate HLA-Resolve/HPRC/Lai arm, not this.

DB-VERSION CONFOUND (carried from EXPERIMENTS.md + aleix/README.md): Immuannot's IMGT reference
(Zenodo "Data-2024Feb02") and SpecImmune's (3.64.0 in Experiment D) are different releases, so some
apparent field-3/4 discordances are allele *renaming*, not real calling differences. This is why
the headline concordance metric here is FIELD 2 (protein-level, non-synonymous -- Marc's standing
priority), which is far more robust to renaming than fields 3/4.

Aggregate-only output (counts/rates, no genotypes) -- keep raw inputs on the VM (egress caveat,
DECISIONS.md); paste the markdown / bring it into reports/immuannot_pilot/.

Usage (any pixi env with pandas -- no matplotlib needed):
  python3 scripts/diagnose_immuannot_pilot.py [--outroot ~/pipeline_outputs]
"""
import argparse
import gzip
import os
import re
import sys

import pandas as pd

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
NULL = {"", "NA", "-", "nan", "None", ".", "na"}

# --- Field-cascade machinery, kept in sync with analyze_experiment_d_field_cascade.py.
# Copied (not imported) so this runs in the `specimmune` env without pulling in matplotlib. ---


def parse_fields(raw):
    if raw is None:
        return None
    s = str(raw).strip()
    if s in NULL or (hasattr(pd, "isna") and pd.isna(raw)):
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    return fields or None


def compare_allele(sr_fields, truth_fields):
    result, mismatched, unassessable = [], False, False
    for lvl in range(1, 5):
        if mismatched:
            result.append(False); continue
        if unassessable:
            result.append(None); continue
        if lvl > len(sr_fields) or lvl > len(truth_fields):
            result.append(None); unassessable = True; continue
        if sr_fields[lvl - 1] == truth_fields[lvl - 1]:
            result.append(True)
        else:
            result.append(False); mismatched = True
    return result


def pairing_score(res_a, res_b):
    return tuple((1 if res_a[i] is True else 0) + (1 if res_b[i] is True else 0) for i in range(4))


def compare_genotype(a_pair, b_pair):
    a1, a2 = parse_fields(a_pair[0]), parse_fields(a_pair[1])
    b1, b2 = parse_fields(b_pair[0]), parse_fields(b_pair[1])
    if a1 is None or a2 is None or b1 is None or b2 is None:
        return None
    pa = (compare_allele(a1, b1), compare_allele(a2, b2))
    pb = (compare_allele(a1, b2), compare_allele(a2, b1))
    return list(max([pa, pb], key=lambda p: pairing_score(p[0], p[1])))


def geno_status_at(cascade, lvl):
    """Genotype-level status at field lvl, from the two per-allele cascades (best pairing):
      'discordant'  -- either allele is a real mismatch at this field (a False)
      'na'          -- no mismatch, but either allele doesn't resolve this deep (a None) -- one
                       caller simply reports fewer fields, NOT a disagreement (critical for
                       fields 3-4, where Immuannot's full 4-field resolution meets SpecImmune's
                       often-shallower calls; counting these as discordant would understate deep
                       concordance -- same True/False/None discipline as the field-cascade script)
      'concordant'  -- both alleles match at this field (both True)"""
    a, b = cascade[0][lvl - 1], cascade[1][lvl - 1]
    if a is False or b is False:
        return "discordant"
    if a is None or b is None:
        return "na"
    return "concordant"


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr); sys.exit(1)


def stage_reached(row):
    """Furthest pipeline stage a hap attempt reached, from which timing columns are non-null.
    (immuannot_seconds set even if immuannot produced no gtf -- 'reached immuannot' != 'succeeded';
    the calls file tells success. A suspiciously short immuannot run is flagged separately.)"""
    if pd.notna(row.get("immuannot_seconds")):
        return "reached_immuannot"
    if pd.notna(row.get("trim_seconds")):
        return "failed_at_trim"
    if pd.notna(row.get("contig_lookup_seconds")):
        return "failed_after_lookup"
    return "failed_at_target_finding"  # 0 contigs overlap / samtools view failed


def load_specimmune(cohort, outroot):
    """{(person_id, gene_bare): full comparison_log.csv row (dict)} -- kept as the full row, not
    just the (si_1, si_2) tuple, so callers can also reach the confidence columns confirmed real
    in compare_hla_results.py's own CSV_FIELDS: si_identity_1/2, si_reads_1/2, si_ambig_n_1/2,
    si_tied_1/2, si_step1_1/2. Not guessed -- grepped from the actual writer."""
    out = {}
    for pid in cohort["person_id"].astype(str):
        log_path = os.path.join(outroot, pid, "comparison_log.csv")
        if not os.path.exists(log_path):
            continue
        df = pd.read_csv(log_path, dtype=str)
        df = df[df["run_label"] == "experiment_d"]
        if df.empty:
            continue
        df = df.sort_values("timestamp").drop_duplicates(["person_id", "gene"], keep="last")
        for _, r in df.iterrows():
            out[(str(r["person_id"]), str(r["gene"]))] = r.to_dict()
    return out


def si_alleles(row):
    return (row.get("specimmune_1"), row.get("specimmune_2")) if row else None


def si_confidence_label(row):
    """'confident' / 'ambiguous(tied)' / 'lower_identity(x.xxx)' / 'unknown' -- reuses the exact
    'confident' definition analyze_experiment_d.py already established (not tied, identity>=0.99
    on both haplotypes), so this doesn't invent a new confidence bar."""
    if not row:
        return "unknown"

    def truthy(x):
        return str(x).strip().lower() in {"true", "1", "yes"}
    try:
        tied = truthy(row.get("si_tied_1")) or truthy(row.get("si_tied_2"))
        id1, id2 = float(row.get("si_identity_1")), float(row.get("si_identity_2"))
    except (TypeError, ValueError):
        return "unknown"
    if tied:
        return "ambiguous(tied)"
    if id1 >= 0.99 and id2 >= 0.99:
        return "confident"
    return f"lower_identity({min(id1, id2):.3f})"


def is_na(v):
    return v is None or str(v).strip() in NULL or (hasattr(pd, "isna") and pd.isna(v))


def classify_pair(imm_pair, si_pair):
    """Reason taxonomy for why a (person, gene) comparison is or isn't assessable -- Marc,
    2026-07-22: 'why are unresolved cases unresolved, and does SpecImmune ever fail to resolve
    too, or does it just claim something with less certainty?' This answers the first half
    (which caller, if either, produced no call at all); si_confidence_label answers the second
    half (when SpecImmune DOES produce a call, how sure is it)."""
    imm_na1 = imm_pair is None or is_na(imm_pair[0])
    imm_na2 = imm_pair is None or is_na(imm_pair[1])
    si_na1 = si_pair is None or is_na(si_pair[0])
    si_na2 = si_pair is None or is_na(si_pair[1])
    imm_full_na, si_full_na = imm_na1 and imm_na2, si_na1 and si_na2
    imm_partial, si_partial = imm_na1 != imm_na2, si_na1 != si_na2

    if imm_full_na and si_full_na:
        return "both_na"
    if imm_full_na and not si_full_na:
        return "immuannot_na_specimmune_resolved"
    if si_full_na and not imm_full_na:
        return "specimmune_na_immuannot_resolved"
    if imm_partial or si_partial:
        return "one_haplotype_na"
    return "both_resolved"


# GTF attributes mix quoted strings ('template_allele "HLA-A*01:01"') AND unquoted bare values
# ('template_distance 0') -- confirmed in the real sanity-check output (parse_gtf's first-line
# print showed exactly 'template_distance 0', no quotes). A quoted-only regex silently drops
# every unquoted attribute -- caught by a synthetic test before this ever touched real data.
ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^";]+))')


def parse_gtf_rich(gtf_gz_path):
    """Merges EVERY attribute seen across every GTF row (gene + transcript feature rows alike)
    for the same gene, keyed by gene_name. Deliberately generic (doesn't assume template_distance
    lives on the 'gene' row and consensus on the 'transcript' row, though the docs say so) --
    UNVERIFIED against a real transcript-row attribute string as of 2026-07-22 (only gene-row
    samples have actually been seen printed from a real run); confirm with a real file before
    trusting template_distance/'new'-tag values pulled from here. Returns {gene_name: {attr: val}}."""
    per_gene, gid_to_name = {}, {}
    with gzip.open(gtf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = {k: (q if q else b.strip()) for k, q, b in ATTR_RE.findall(fields[8])}
            gid, gname = attrs.get("gene_id"), attrs.get("gene_name")
            if gid and gname:
                gid_to_name[gid] = gname
            key = gname or (gid_to_name.get(gid) if gid else None)
            if not key:
                continue
            slot = per_gene.setdefault(key, {})
            for k, v in attrs.items():
                slot.setdefault(k, v)
    return per_gene


def load_immuannot_confidence(cohort, outroot):
    """{(person_id, hap, gene_bare): rich-attr dict} by reading each person's own hap{1,2}.gtf.gz
    directly (NOT immuannot_calls.tsv, which only carries the final allele string) -- these files
    are still on disk (never pruned) so no re-run is needed. Missing files are skipped, not fatal
    -- a person whose gtf.gz got cleaned up just won't have a confidence row, same discipline as
    everywhere else in this project (report what's missing, don't guess it)."""
    out = {}
    for pid in cohort["person_id"].astype(str):
        for hap in ("hap1", "hap2"):
            gz = os.path.join(outroot, pid, "immuannot_output", f"{hap}.gtf.gz")
            if not os.path.exists(gz):
                continue
            try:
                per_gene = parse_gtf_rich(gz)
            except OSError:
                continue
            for gname, attrs in per_gene.items():
                bare = gname.replace("HLA-", "")
                out[(pid, hap, bare)] = attrs
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--cohort", default=None)
    ap.add_argument("--analysis-dir", default=None)
    args = ap.parse_args()

    cohort_path = args.cohort or os.path.join(args.outroot, "experiment_d", "cohort.tsv")
    calls_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    timing_path = os.path.join(args.outroot, "immuannot_timing.tsv")
    for p in (cohort_path, calls_path, timing_path):
        if not os.path.exists(p):
            die(f"not found: {p}")

    cohort = pd.read_csv(cohort_path, sep="\t", dtype=str)
    anc = dict(zip(cohort["person_id"].astype(str), cohort["ancestry"].astype(str)))
    attempted = set(cohort["person_id"].astype(str))
    # keep_default_na=False: run_immuannot writes literal "NA" for a missing allele; pandas would
    # otherwise turn that back into NaN and break the NA counting (bug caught in summarize script).
    calls = pd.read_csv(calls_path, sep="\t", dtype=str, keep_default_na=False)
    timing = pd.read_csv(timing_path, sep="\t", dtype=str)
    for c in ("contig_lookup_seconds", "trim_seconds", "immuannot_seconds", "hap_total_seconds"):
        if c in timing.columns:
            timing[c] = pd.to_numeric(timing[c], errors="coerce")

    calls["gene_bare"] = calls["gene"].str.replace("^HLA-", "", regex=True)
    completed = set(calls["person_id"].astype(str))

    parts = ["# Immuannot 60-person pilot -- diagnosis + cross-method comparison", "",
             "Aggregate-only. No ground truth: Immuannot vs SpecImmune is cross-method AGREEMENT, "
             "not accuracy -- see the script header for the adjudication limit and the DB-version "
             "confound (why FIELD 2 is the headline metric).", ""]

    # ============ SECTION 1: failure distribution ============
    parts += ["## 1. Failure distribution", ""]
    people_in_timing = set(timing["person_id"].astype(str))
    pre_loop_skipped = sorted(attempted - people_in_timing)
    parts += [
        f"- Attempted: {len(attempted)}. Produced >=1 call: {len(completed)}. "
        f"Absent from timing entirely (skipped before the hap loop -- almost certainly no assembly "
        f"data / not revio, per reports/lr_data_census/README.md): {len(pre_loop_skipped)}.",
        "",
        "### 1a. Where each haplotype attempt stopped (stage)", "",
        "| Stage reached | Count (hap attempts) |", "|---|---|",
    ]
    timing["stage"] = timing.apply(stage_reached, axis=1)
    for stage, n in timing["stage"].value_counts().items():
        parts.append(f"| {stage} | {n} |")

    # per-person: both haps fail vs one hap vs both ok (person-concentration)
    parts += ["", "### 1b. Person-concentration (are failures clustered on whole people?)", ""]
    reached = timing[timing["stage"] == "reached_immuannot"]
    haps_ok = reached.groupby("person_id").size()
    both_ok = (haps_ok == 2).sum()
    one_ok = (haps_ok == 1).sum()
    parts += [
        f"- People with BOTH haps reaching Immuannot: {both_ok}",
        f"- People with exactly ONE hap reaching Immuannot: {one_ok}",
        f"- People skipped pre-loop (no hap attempted): {len(pre_loop_skipped)}",
        "",
        "### 1c. Gene-concentration (are specific genes systematically missing?)", "",
        "Among people who produced calls, how often each classical gene comes back NA "
        "(per person x haplotype = 2 slots per completed person). A high, uniform-across-people "
        "rate at one gene = systematic; scattered = per-person.", "",
        "| Gene | NA slots | NA rate |", "|---|---|---|",
    ]
    n_completed = len(completed)
    for g in GENES:
        sub = calls[calls["gene_bare"] == g]
        present_people = set(sub["person_id"].astype(str))
        na_slots = 0
        for _, r in sub.iterrows():
            na_slots += (parse_fields(r["immuannot_1"]) is None) + (parse_fields(r["immuannot_2"]) is None)
        # people with NO row for this gene at all also count as 2 NA slots each
        na_slots += 2 * len(completed - present_people)
        denom = 2 * n_completed
        parts.append(f"| {g} | {na_slots}/{denom} | {100*na_slots/denom:.0f}% |" if denom else f"| {g} | - | - |")

    # ============ SECTION 2: timing distribution ============
    parts += ["", "## 2. Time distribution", ""]
    reached_t = timing[timing["stage"] == "reached_immuannot"]

    def stat(s, unit):
        s = s.dropna()
        if s.empty:
            return "n=0"
        return (f"n={len(s)}: median {s.median():.1f}{unit}, mean {s.mean():.1f}{unit}, "
                f"range {s.min():.1f}-{s.max():.1f}{unit}")

    per_person_total = timing.groupby("person_id")["hap_total_seconds"].sum() / 60
    parts += [
        f"- immuannot.sh per hap: {stat(reached_t['immuannot_seconds'] / 60, ' min')}",
        f"- trim (samtools faidx) per hap: {stat(reached_t['trim_seconds'], ' s')}",
        f"- per-person total (both haps): {stat(per_person_total, ' min')}",
        "",
        "Suspiciously fast immuannot runs (< 0.5 min) -- 'reached immuannot' but likely produced "
        f"little (the gene_infra positive-control ran in ~0.3 min): "
        f"{int((reached_t['immuannot_seconds'] < 30).sum())} hap attempts.",
    ]

    # ============ SECTION 3: Immuannot vs SpecImmune, per gene x ancestry ============
    parts += ["", "## 3. Immuannot vs SpecImmune-LR -- concordance per gene x ancestry, all 4 fields",
              "",
              "Cross-method agreement (not accuracy). Cell = concordant% (assessable N). Assessable "
              "at a field level = both methods actually resolve that deep AND agree on all shallower "
              "fields -- pairs where one caller simply reports fewer fields are EXCLUDED (na), not "
              "counted as disagreement (True/False/None discipline, same as "
              "analyze_experiment_d_field_cascade.py). This is why N shrinks at fields 3-4: fewer "
              "people have both callers resolving that deep, not more failures.",
              "",
              "**Read fields 3-4 with the DB-version confound front of mind:** Immuannot (IMGT "
              "Data-2024Feb02) and SpecImmune (3.64.0) ship different IMGT releases, and synonymous "
              "(field 3) / non-coding (field 4) differences are exactly where allele RENAMING "
              "between releases shows up -- so a chunk of field-3/4 discordance is naming, not "
              "biology. Field 2 (protein, non-synonymous) is the headline for this reason.", ""]
    si = load_specimmune(cohort, args.outroot)
    imm = {(str(r["person_id"]), r["gene_bare"]): (r["immuannot_1"], r["immuannot_2"])
           for _, r in calls.iterrows()}

    def cell(gene, group, lvl):
        """Returns (concordant, discordant). 'na' pairs (one caller shallower) excluded from both,
        so the rate below is concordant/(concordant+discordant) -- deep-resolution gaps don't count
        against agreement."""
        conc = disc = 0
        for pid in attempted:
            if anc.get(pid) != group:
                continue
            a = imm.get((pid, gene))
            b = si_alleles(si.get((pid, gene)))
            if a is None or b is None:
                continue
            casc = compare_genotype(a, b)
            if casc is None:
                continue
            st = geno_status_at(casc, lvl)
            if st == "concordant":
                conc += 1
            elif st == "discordant":
                disc += 1
        return conc, disc

    for lvl, label in [(1, "Field 1 (allele group)"),
                       (2, "Field 2 (protein, non-synonymous -- HEADLINE)"),
                       (3, "Field 3 (synonymous, coding -- DB-version-confounded)"),
                       (4, "Field 4 (non-coding -- DB-version-confounded)")]:
        parts += [f"### {label}", "", "| Gene | " + " | ".join(GROUPS) + " | overall |",
                  "|---|" + "---|" * (len(GROUPS) + 1)]
        for g in GENES:
            cells, tot_c, tot_d = [], 0, 0
            for grp in GROUPS:
                c, d = cell(g, grp, lvl)
                tot_c += c; tot_d += d
                n = c + d
                cells.append(f"{100*c/n:.0f}% ({n})" if n else "-")
            tot_n = tot_c + tot_d
            overall = f"{100*tot_c/tot_n:.0f}% ({tot_n})" if tot_n else "-"
            parts.append(f"| {g} | " + " | ".join(cells) + f" | {overall} |")
        parts.append("")

    parts += [
        "N in () = assessable at that field (concordant + discordant; na excluded). N shrinks with "
        "depth because fewer pairs both resolve that deep. Read low-N cells (MID/SAS are thin) as "
        "directional. A gene x ancestry cell that is low on "
        "Field 2 specifically -- not just Field 3/4 -- is a real protein-level divergence worth "
        "cross-referencing against the known SpecImmune DQA1/DRB1 and AoU DQA1 findings.",
    ]

    # ============ SECTION 4: why unresolved cases are unresolved ============
    # Marc, 2026-07-22: "study why the cases that are unresolved are unresolved. Has SpecImmune
    # also unresolved any case? or do they just claim something with less uncertainty?"
    parts += ["", "## 4. Why comparisons are unresolved -- which caller, and how confident is the "
              "one that DOES resolve", "",
              "Reason taxonomy per (person, gene): `both_na` (neither caller called it), "
              "`immuannot_na_specimmune_resolved` / `specimmune_na_immuannot_resolved` (one caller "
              "produced nothing while the other did), `one_haplotype_na` (one allele missing on "
              "one or both sides), `both_resolved` (proceeds to the field tables above).", "",
              "| Gene | both_na | immuannot_na, SI resolved | SI_na, immuannot resolved | "
              "one_hap_na | both_resolved |",
              "|---|---|---|---|---|---|"]
    reason_counts_overall = {}
    for g in GENES:
        counts = {"both_na": 0, "immuannot_na_specimmune_resolved": 0,
                  "specimmune_na_immuannot_resolved": 0, "one_haplotype_na": 0, "both_resolved": 0}
        for pid in attempted:
            a = imm.get((pid, g))
            b = si_alleles(si.get((pid, g)))
            r = classify_pair(a, b)
            counts[r] += 1
            reason_counts_overall[r] = reason_counts_overall.get(r, 0) + 1
        parts.append(f"| {g} | {counts['both_na']} | {counts['immuannot_na_specimmune_resolved']} "
                      f"| {counts['specimmune_na_immuannot_resolved']} | "
                      f"{counts['one_haplotype_na']} | {counts['both_resolved']} |")
    parts += ["", f"Overall: {reason_counts_overall}", "",
              "### When SpecImmune DOES resolve a gene, how confident is it?", "",
              "Answers 'do they just claim something with less certainty' directly: `confident` = "
              "not tied, both haplotype identities >= 0.99 (same bar analyze_experiment_d.py "
              "already used); `ambiguous(tied)` = SpecImmune itself flagged a tie among candidates; "
              "`lower_identity(x)` = resolved but below the identity bar.", "",
              "| Gene | confident | ambiguous(tied) | lower_identity | unknown |",
              "|---|---|---|---|---|"]
    for g in GENES:
        labels = {"confident": 0, "ambiguous(tied)": 0, "unknown": 0}
        lower_n = 0
        for pid in attempted:
            row = si.get((pid, g))
            if row is None:
                continue
            lbl = si_confidence_label(row)
            if lbl.startswith("lower_identity"):
                lower_n += 1
            else:
                labels[lbl] = labels.get(lbl, 0) + 1
        parts.append(f"| {g} | {labels['confident']} | {labels['ambiguous(tied)']} | {lower_n} | "
                      f"{labels['unknown']} |")

    parts += ["", "### Immuannot's own confidence signal when it DOES resolve a gene", "",
              "Verified 2026-07-22 against a real hap1.gtf.gz (person 1001871): the docs' 'new' "
              "tag does not actually appear in this output -- every transcript row instead carries "
              "`template_warning` ('NA' when clean, or a specific flag like 'partial_CDS' when "
              "there's a real structural issue -- e.g. HLA-S showed this in the real file). Using "
              "template_warning here instead, since it's the signal that's actually present. "
              "template_distance (edit distance to the nearest IMGT allele, from the gene row) was "
              "separately confirmed present in the very first real run.", "",
              "| Gene | n with template_distance | mean template_distance | n with a real warning |",
              "|---|---|---|---|"]
    imm_conf = load_immuannot_confidence(cohort, args.outroot)
    for g in GENES:
        dists, n_warned = [], 0
        for (pid, hap, gname), attrs in imm_conf.items():
            if gname != g:
                continue
            td = attrs.get("template_distance")
            if td is not None:
                try:
                    dists.append(float(td))
                except ValueError:
                    pass
            warn = str(attrs.get("template_warning", "NA")).strip()
            if warn and warn != "NA":
                n_warned += 1
        mean_d = f"{sum(dists)/len(dists):.2f}" if dists else "-"
        parts.append(f"| {g} | {len(dists)} | {mean_d} | {n_warned} |")

    # ============ SECTION 5: per-individual detail (VM-local CSV, NOT pasted -- bare ids) ============
    detail_rows = []
    for pid in sorted(attempted, key=lambda p: (anc.get(p, "?"), p)):
        for g in GENES:
            a = imm.get((pid, g))
            row = si.get((pid, g))
            b = si_alleles(row)
            reason = classify_pair(a, b)
            f2 = "-"
            if reason == "both_resolved":
                casc = compare_genotype(a, b)
                f2 = geno_status_at(casc, 2) if casc else "unresolved(allele-level)"
            detail_rows.append({
                "person_id": pid, "ancestry": anc.get(pid, "?"), "gene": g,
                "immuannot_1": a[0] if a else "NA", "immuannot_2": a[1] if a else "NA",
                "specimmune_1": b[0] if b else "NA", "specimmune_2": b[1] if b else "NA",
                "si_confidence": si_confidence_label(row) if row else "no_si_row",
                "reason": reason, "field2_status": f2,
            })
    adir = args.analysis_dir or os.path.join(args.outroot, "immuannot_pilot", "analysis")
    os.makedirs(adir, exist_ok=True)
    detail_path = os.path.join(adir, "immuannot_vs_specimmune_detail.csv")
    pd.DataFrame(detail_rows).to_csv(detail_path, index=False)

    parts += ["", "## 5. Per-individual detail (grouped by ancestry)", "",
              f"Written to `{detail_path}` -- **bare person_ids, VM-local only, do not paste in "
              f"full or commit** (same rule as cohort.tsv / immuannot_timing.tsv). One row per "
              f"(person, gene): both callers' alleles, SpecImmune's confidence label, the "
              f"unresolved-reason, and the Field 2 status. Open directly on the VM to see exactly "
              f"which individuals/genes drive each ancestry's numbers above -- "
              f"e.g. `column -s, -t {detail_path} | less -S`, or filter to just the discordant "
              f"rows: `awk -F, '$10==\"discordant\"' {detail_path}`."]

    md = "\n".join(parts) + "\n"
    out_md = os.path.join(adir, "immuannot_pilot_diagnosis.md")
    with open(out_md, "w") as f:
        f.write(md)
    print(md)
    print(f"\n(written to {out_md} -- aggregate-only, safe to paste back / bring into the repo)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
