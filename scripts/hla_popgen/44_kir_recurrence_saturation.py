#!/usr/bin/env python3
"""WS-A (Sprint S04) -- KIR recurrence + saturation, on the same footing as HLA (WS-B setup).

Repeats for KIR the recurrence-class and discovery/saturation analyses already done for HLA
(`03_novel_alleles.py`, `34_novel_protein_recurrence.py`, `24_novelty_by_field.py`,
`04_allele_saturation.py`, `18_discovery_rate.py`, and the corrected equal-N per-ancestry
methodology in `39_saturation_by_ancestry.py` -- its README documents the S03 fixes and is
followed here, not re-derived), and runs the IDENTICAL analysis on HLA at the same time, on the
same unrelated people (KIR ∩ HLA, then one shared `greedy_unrelated()` pass) and the same
ancestry scheme and subsampling seeds, so WS-B (script 46, KIR-vs-HLA catalogue coverage) can
compare the two regions on equal footing without re-deriving anything.

REUSE (imported via importlib, not re-derived):
  - `41_kir_pilot.py`: classify_novelty_tier, greedy_unrelated, load_relatedness_pairs, KIN_MIN,
    KIR_GENES, ANCESTRY_ORDER. (parse_hap_gtf itself is NOT reused for the identity-matching pass
    below -- see "KIR IDENTITY EXTRACTION" -- because it does not expose the contig/copy-index
    fields a `cds.fa.gz` join needs; `parse_hap_gtf_full()` in this file adapts its regex
    extraction, extended with those two fields, and explicitly reuses
    `classify_novelty_tier()` rather than re-deriving novelty classification.)
  - `43_kir_full_aggregate.py`: discover_people, suppressed(), rate_row(), wilson_ci(),
    load_ancestry_map().
  - `39_saturation_by_ancestry.py`: build_labeled_calls (the full HLA call-labeling pipeline --
    field_class/seq_class/prot_id/cds_id/keep_clean, ancestry + unrelated flags), its own imports
    of `24_novelty_by_field.py` (m24) and `04_allele_saturation.py` (m04: incidence_freqs, chao2,
    chao2_extrapolate, clench_curve/fit_clench_asymptote).
  - `24_novelty_by_field.py`: normalize_allele_name (genomic-level HLA identity), sha8,
    protein_info (translate + strip terminal stop) -- reused for BOTH species' CDS/protein hashing
    so the two species get identically-computed hashes over identically-extracted sequences.
  - `03_novel_alleles.py`: parse_cds_fasta -- the `<hap>/cds.fa.gz` reader HLA's own novel-allele
    pipeline already uses (reference/IMMUANNOT_GTF_SPEC.md part D/E), reused verbatim for KIR's
    `cds.fa.gz` too (same Immuannot codepath writes it for both HLA and KIR runs -- see "KIR
    IDENTITY EXTRACTION" below for why this file should exist per haplotype for KIR as well).

ALLELE IDENTITY -- THREE MATCHED GRANULARITIES, SAME DEFINITION FOR BOTH SPECIES
(this section replaces an earlier, mismatched design flagged by the S04 coordinator as an open
item before the VM run: KIR previously used a genomic-level identity for its "all"/"any_novel"
tracks while HLA used a CDS-level identity for the same tracks. Fixed here.):

  - `genomic` -- the full called allele string, including non-coding differences. HLA: the
    UNTRUNCATED normalized consensus name (`normalize_allele_name(consensus)[1]`, all colon
    fields kept, e.g. `"HLA-A*01:01:01:new"`) -- this is the SAME string
    `24_novelty_by_field.field_class_of()`/`add_call_labels()` already derives `field_class` from
    (the "new" suffix's field DEPTH is what makes a call `known`/`f4_noncoding`/`f3_synonymous`/
    `f2_protein`), so using the full name as the genomic-level identity is not a new definition,
    just using the existing name at its full, untruncated resolution instead of the 2-/3-field
    truncations (`prot_id`/`cds_id`) 24 also derives from it. KIR: the raw Immuannot `consensus`
    string for a call (matches `43_kir_full_aggregate.py`'s own convention).
  - `cds` -- the coding-sequence identity (a hash of the observed CDS nucleotide sequence). HLA:
    `cds_id` (`24_novelty_by_field.allele_ids()`, already computed -- CDS-level 3-field name for
    known/synonymous-tier calls, a `sha8` hash of the reconstructed CDS for the rest). KIR: a NEW
    hash (`<gene>_cds_<sha8>`) over the observed CDS sequence extracted from `<hap>/cds.fa.gz` and
    joined to the GTF call by `(contig, gene)` -- see "KIR IDENTITY EXTRACTION" below.
  - `protein` -- the translated-protein identity. HLA: `prot_id` (already computed -- 2-field
    IPD-IMGT name for known alleles, a `sha8` hash of the translated protein for novel-protein
    calls). KIR: a NEW hash (`<gene>_prot_<sha8>`) over `24_novelty_by_field.protein_info(cds_seq)
    ["protein"]` (translate + strip one terminal stop codon -- reused verbatim, not re-derived),
    from the SAME `cds.fa.gz`-extracted sequence as the `cds` level.

  "any_novel" (any-level novel) = novel at the GENOMIC level: KIR `novelty_tier != "known"`
  (excludes `undetermined`); HLA `field_class != "known"` (excludes `uncalled`). Identity used for
  this track = `genomic`.
  "protein_novel" (protein-level novel) = novel at the PROTEIN level: KIR `novelty_tier ==
  "novel_protein"`; HLA `field_class == "f2_protein" AND seq_class == "novel_protein"`. Identity
  used for this track = `protein`.

  Five tracks are exported per gene x ancestry x species: the three baseline granularities
  (`genomic`, `cds`, `protein` -- every clean call, no novelty filter -- richness/coverage at each
  resolution, directly answering WS-B's "which region is better represented in its reference, at
  which resolution") plus the two novelty-filtered tracks (`any_novel`, `protein_novel`).

KIR IDENTITY EXTRACTION -- WHAT FILES EXIST, WHAT WAS CHECKED, WHAT COULD NOT BE VERIFIED
  `41_kir_pilot.py` and `43_kir_full_aggregate.py` have only ever read
  `<pid>/immuannot_output/hap{1,2}.gtf.gz` (confirmed by reading both files in full -- neither
  references `cds.fa.gz`, `.trimmed.fa`, or any other Immuannot intermediate). Per
  `reference/IMMUANNOT_GTF_SPEC.md` parts D/E (a static read of Immuannot's upstream source code,
  NOT independently verified by running the pipeline or inspecting the real VM filesystem -- no
  VM access from this local-only task), Immuannot's own `searchTemplate.py` writes a per-haplotype
  `<pid>/immuannot_output/hap{1,2}/cds.fa.gz` (a DIRECTORY named `hap1`/`hap2`, sibling to the
  flat `hap1.gtf.gz`/`hap2.gtf.gz` files) containing the actual observed, contig-derived CDS
  nucleotide sequence for every gene copy detected -- in coding (5'->3') orientation, header
  `>{contig}_{gene}_{i}` -- and `annot.combine.sh`'s cleanup step does NOT delete this file. This
  project's `run_immuannot_person.py` runs the unmodified upstream `immuannot.sh` and never
  deletes the `<outpref>` folder itself, for EITHER the HLA run (chr6 default region) or the KIR
  run (chr19 region, `41_kir_pilot.py`'s `--region`/`--pad`) -- the same script, same cleanup
  behavior, differing only in `--region`. `03_novel_alleles.py` already depends on this exact file
  existing for HLA's own novel-allele sequence matching. This document's own caveat: "this specific
  claim... is inferred from source logic, not empirically verified (no minimap2 available in this
  environment to run the pipeline end-to-end)" -- i.e. the file's existence for the real
  `~/pipeline_outputs_kir` tree has NOT been confirmed by this task, only argued from the shared
  code path. `44_kir_recurrence_saturation.py`'s own code (`build_kir_identity_sets()`) checks
  for `cds.fa.gz` at runtime and does not assume it: if it is not found for ANY haplotype across
  the whole run, `cds`/`protein`/`protein_novel` are exported as the literal string `"NA"` for
  every KIR row (never silently zero, never silently compared against HLA's real numbers) and a
  loud warning is logged; the `genomic`/`any_novel` tracks (which only need `hap{1,2}.gtf.gz`,
  confirmed to exist by 41/43's own successful full-cohort run) are computed normally regardless.

  Ambiguity: for a gene detected in >1 copy on the same contig (a real, if rare, structural
  duplication case), the `(contig, gene)` join against `cds.fa.gz` is inherently ambiguous (the
  detection-order index `i` in `cds.fa.gz`'s header does not reliably map onto the final,
  coordinate-sorted GTF row order for >1 copy -- reference/IMMUANNOT_GTF_SPEC.md part D's own
  documented caveat). Such calls are EXCLUDED from the `cds`/`protein`/`protein_novel` tracks
  (counted in the exported QC, never guessed), following `03_novel_alleles.match_novel_rows()`'s
  own precedent exactly.

DISCLOSURE (hard rule, same as 39/43 -- see AGENT_PREAMBLE.md / ORCHESTRATOR_HANDOFF.md sect 3):
  - Every raw count 1-19 is written as the literal string "<20" (never 0, never omitted). A true
    0 stays "0". Applies to recurrence-class ALLELE counts too, even though those are allele
    counts rather than participant counts -- per this sprint's own explicit instruction ("still
    write 1-19 as '<20' to be safe"): an allele carried by exactly 1-19 people is disclosive about
    those specific people's rare genotype regardless of whether the number being reported is
    labeled "people" or "alleles".
  - No allele name/string (KIR consensus/hash or HLA prot_id/cds_id/genomic name) is ever written
    next to a carrier count < 20.
  - Any rate whose numerator or denominator is 1-19 is blanked (point estimate AND CI).
  - Good-Turing coverage / Chao2 are computed on-VM from f1/f2 (Q1/Q2), which can themselves be
    small; only the resulting coverage/Chao2/SE ESTIMATES are exported, never the raw f1/f2 pair,
    per this sprint's instruction ("compute Chao on VM, export only the estimate if f1/f2 are
    small"). f1/f2 are written to the export TSV only when both are >= SUPPRESS_BELOW.
  - The literal string "NA" (a level that could not be matched for a species -- see above) is
    distinct from both "0" (a true zero) and "<20" (a masked small count) and must never be
    conflated with either downstream (45's figure script must not plot an "NA" row as a zero).

RECURRENCE CLASSES (allele counts, not mutually exclusive by design -- documented, not a bug):
  eq1  -- alleles carried by EXACTLY 1 unrelated person.
  eq2  -- alleles carried by EXACTLY 2 unrelated people.
  gt2  -- alleles carried by >= 3 unrelated people (a superset that INCLUDES ge20).
  ge20 -- alleles carried by >= 20 unrelated people (a subset of gt2, called out separately
          because 20 is this project's own disclosure floor -- ORCHESTRATOR_HANDOFF.md WS-A).

STATISTICS:
  - Good-Turing sample coverage, incidence-data (Chao & Jost 2012, Ecol. 93:2533-2547) bias
    correction: C_hat = 1 - (Q1/n) * [(m-1)Q1 / ((m-1)Q1 + 2Q2)] (Q2>0 case), n = total number of
    incidences (sum of per-allele occurrence counts across sampling units), m = number of units.
  - Chao2 richness + its classic incidence-based SE (Chao 1987, Biometrics 43:783-791; formulas as
    reproduced in Colwell/EstimateS's user guide eq. 21b/23b) -- `chao2()` itself imported
    verbatim from 04_allele_saturation.py; the SE is new (04/39 never computed/exported it) and is
    unit-tested here for non-negativity and the Q2==0 bias-corrected branch.
  - Equal-N discovery slope: linear fit of the mean rarefaction curve over the window
    [0.9*N*, N*], exactly 39_saturation_by_ancestry's own "primary" cross-ancestry-comparable
    number (see that script's README) -- generalized here to ALSO run per recurrence class, not
    just the >=1-carrier curve.
  - Rarefaction/discovery curves: adapted from 39_saturation_by_ancestry.rarefaction_threshold_curves
    (itself already generic over an arbitrary `thresholds` list) to accept a PRECOMPUTED list of
    permutation orders instead of generating its own -- so KIR and HLA curves, within one ancestry,
    walk the EXACT SAME people in the EXACT SAME order across all `n_permutations` replicates (the
    brief's "same subsampling seeds" requirement for WS-B). thresholds=[1,2,3,20]; exact-class
    curves are derived by differencing (eq1=thr[1]-thr[2], eq2=thr[2]-thr[3], gt2=thr[3],
    ge20=thr[20]).
  - Optional iNEXT-style 2N extrapolation via 04's `chao2_extrapolate` (cheap -- same Q1/f0_hat
    already computed for Chao2. On by default (`--no-extrapolate-2n` to skip).

OUTPUTS (--out-dir, default ~/s04/results/44/; only aggregates, safe to pull off the VM):
  recurrence_classes.tsv   gene, ancestry, level, species, eq1, eq2, gt2, ge20  (masked <20, or
                           "NA" -- see ALLELE IDENTITY)
  saturation_curves.tsv    gene, ancestry, level, species, n, mean_distinct, lo2_5, hi97_5,
                           mean_eq1, mean_eq2, mean_gt2, mean_ge20 (rarefaction curve points;
                           n=1..N per (gene,ancestry,level,species), subsampled to
                           --curve-stride points to keep the file small)
  equal_n_slope.tsv        gene, ancestry, level, species, n_star, slope_per_1000,
                           mean_distinct_at_n_star
  coverage_chao2.tsv       gene, ancestry, level, species, n_people, s_obs, good_turing_coverage,
                           chao2, chao2_se, chao2_undetected_f0hat, chao_new_by_2n (if requested),
                           q1, q2 (blank unless both >= SUPPRESS_BELOW)
  kir_cds_match_qc.tsv     aggregate-only cds.fa.gz join QC: n_join_rows, n_matched,
                           n_ambiguous_copy, n_missing_cds_record, n_cds_fasta_missing,
                           n_hap_seen, n_hap_cds_fasta_present, cds_available (bool)
  kir_protein_catalogue_qc.tsv  (2026-09-27) one row per KIR gene: protein_catalogue_covered (bool),
                           reason in {ok, no_refdata, gene_not_in_catalogue,
                           no_catalogued_protein_entries} -- see kir_protein_catalogue_status()/
                           kir_gene_protein_covered(). A gene with covered=False has its
                           protein/protein_novel rows exported as "NA" everywhere in
                           recurrence_classes.tsv/coverage_chao2.tsv/etc, never a hash-based guess.
  STATUS.txt               aggregate-progress-only status file, rewritten as the run proceeds.

USAGE (VM; see the module's own --help for every flag):
  cd ~/s04 && PYTHONPATH=~/s04:~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 -u 44_kir_recurrence_saturation.py \\
      --kir-outroot ~/pipeline_outputs_kir \\
      --hla-table ~/pipeline_outputs/hla_calls_rich.tsv \\
      --hla-people-outroot ~/pipeline_outputs/people \\
      --cohort-membership ~/pipeline_outputs/cohort_membership.tsv \\
      --relatedness-table ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv \\
      --refdata ~/tools/Immuannot_refdata \\
      --out-dir ~/s04/results/44 --workers 4 2>&1 | tee -a ~/s04/results/44/run.log

Never writes anywhere under ~/pipeline_outputs* (read-only inputs; verify with a pre/post md5
check of the 14 HLA production tables per this project's standing rule, same as 43's run did).
"""
import argparse
import gzip
import importlib.util
import math
import multiprocessing as mp
import os
import re
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
LEVELS = ["genomic", "cds", "protein", "any_novel", "protein_novel"]
KIR_LEVELS_NEEDING_CDS_FASTA = {"cds", "protein", "protein_novel"}
SPECIES = ["kir", "hla"]
KIR_THRESHOLDS = [1, 2, 3, 20]
# 39_saturation_by_ancestry.py's own DEFAULT_PEOPLE_OUTROOT convention: cds.fa.gz lives under
# <outroot>/people/<pid>/immuannot_output/hap{1,2}/cds.fa.gz, one level below --hla-table's own
# default parent dir. BUG FIX 2026-09-26: this module previously hardcoded plain
# "~/pipeline_outputs" (no "/people") as build_labeled_calls()'s `outroot` arg, so
# match_sequences() found zero cds.fa.gz files for every depth-2/3 HLA call and every HLA
# protein_novel s_obs came out 0 -- see sprints/S04_kir_recurrence_style_share/LOG.md.
DEFAULT_HLA_PEOPLE_OUTROOT = os.path.expanduser("~/pipeline_outputs/people")
N_PERMUTATIONS_DEFAULT = 25       # Pakistan Fig 3e convention, matches 39's default
MIN_PEOPLE_PER_ANCESTRY = 100      # "well-powered" floor for N*, matches 39
CURVE_STRIDE_DEFAULT = 25          # export every Nth curve point (keeps saturation_curves.tsv small)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_m41 = _m43 = _m39 = _m04 = _m24 = _m03 = None


def m41():
    global _m41
    if _m41 is None:
        _m41 = _load_module("41_kir_pilot.py", "kir_pilot_44")
    return _m41


def m43():
    global _m43
    if _m43 is None:
        _m43 = _load_module("43_kir_full_aggregate.py", "kir_full_aggregate_44")
        _m43._kir = m41()
    return _m43


def m39():
    global _m39
    if _m39 is None:
        _m39 = _load_module("39_saturation_by_ancestry.py", "saturation_by_ancestry_44")
    return _m39


def m04():
    global _m04
    if _m04 is None:
        _m04 = _load_module("04_allele_saturation.py", "allele_saturation_44")
    return _m04


def m24():
    global _m24
    if _m24 is None:
        _m24 = _load_module("24_novelty_by_field.py", "novelty_by_field_44")
    return _m24


def m03():
    global _m03
    if _m03 is None:
        _m03 = _load_module("03_novel_alleles.py", "novel_alleles_44")
    return _m03


# ---------------------------------------------------------------------------
# Disclosure helpers (reused verbatim from 43, not re-derived).
# ---------------------------------------------------------------------------
def suppressed(n):
    return m43().suppressed(n)


def rate_row(n, d, pct_decimals=1):
    return m43().rate_row(n, d, pct_decimals)


def _cov_s_obs_numeric(row):
    """s_obs is exported as int for real rows or the string 'NA' for an unmatched level (see
    build_na_rows) -- coerce to float('nan') for 'NA' so callers can uniformly use pd.isna/math."""
    v = row.get("s_obs")
    if v == "NA" or v is None:
        return float("nan")
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def sanity_check_coverage(cov_rows):
    """Runtime sanity gate on `coverage_chao2.tsv` rows (list of dicts, as built by
    build_coverage_chao2_row/build_na_rows), run once per full 44 invocation right before export.
    Catches, LOUDLY (never silently), the two failure modes this module has actually shipped with:
      (1) protein_novel s_obs == 0 for every gene of a species that otherwise has real (>0)
          any_novel novelty -- the HLA `--hla-people-outroot` path bug (match_sequences() silently
          matching zero calls, so seq_class/prot_id never populate for any depth-2/3 call).
      (2) any gene x ancestry's pct_novel_protein > 100% -- the KIR protein-hash identity-
          granularity mismatch (see build_person_kir_identity's docstring). 2026-09-27: the
          denominator here is `protein` s_obs (the SAME level's baseline richness), not `genomic` --
          matching 46_kir_vs_hla_coverage.py's own fixed `pct_novel_protein` definition
          (protein_novel/protein, not protein_novel/genomic -- comparing a novelty COUNT against a
          richness figure from a DIFFERENT identity granularity was itself part of the original bug;
          `protein_novel` is by construction a subset of `protein` at the same granularity, so this
          ratio can only exceed 100% from a real remaining bug, never a granularity mismatch).
    Does not raise (a warning, not a hard failure -- a real biological zero or a small-N NA/blank
    row from disclosure masking are both possible and not necessarily bugs) but always prints to
    stderr so it can't be missed in run.log. Returns the list of warning strings emitted (empty if
    clean), so tests can assert on content instead of stderr capture."""
    warnings = []
    by_key = {}   # (species, gene, ancestry) -> {level: s_obs}
    for r in cov_rows:
        key = (r.get("species"), r.get("gene"), r.get("ancestry"))
        by_key.setdefault(key, {})[r.get("level")] = _cov_s_obs_numeric(r)

    any_novel_total = defaultdict(float)
    protein_novel_total = defaultdict(float)
    protein_novel_real_seen = defaultdict(bool)  # True once a species has >=1 non-"NA" row
    for (species, gene, ancestry), levels in by_key.items():
        any_novel = levels.get("any_novel")
        protein_novel = levels.get("protein_novel")
        protein = levels.get("protein")
        if any_novel is not None and not math.isnan(any_novel):
            any_novel_total[species] += any_novel
        if protein_novel is not None and not math.isnan(protein_novel):
            protein_novel_total[species] += protein_novel
            protein_novel_real_seen[species] = True
        if (protein is not None and protein_novel is not None
                and not math.isnan(protein) and not math.isnan(protein_novel) and protein > 0):
            pct = 100.0 * protein_novel / protein
            if pct > 100.0:
                msg = (f"[44] SANITY WARNING: pct_novel_protein > 100% for species={species} "
                       f"gene={gene} ancestry={ancestry} (protein_novel s_obs={protein_novel:g} > "
                       f"protein s_obs={protein:g}, {pct:.1f}%) -- protein_novel should be a subset "
                       f"of protein at the SAME granularity; see build_person_kir_identity's "
                       f"docstring / kir_gene_protein_covered() for the known KIR failure modes.")
                warnings.append(msg)
                log(msg)
    for species in any_novel_total:
        if (any_novel_total[species] > 0 and protein_novel_real_seen.get(species, False)
                and protein_novel_total.get(species, 0.0) == 0.0):
            msg = (f"[44] SANITY WARNING: protein_novel s_obs == 0 for EVERY gene/ancestry of "
                   f"species={species}, despite any_novel s_obs summing to "
                   f"{any_novel_total[species]:g} > 0 -- protein-level novelty pipeline likely "
                   f"broken for this species (e.g. the HLA --hla-people-outroot cds.fa.gz path "
                   f"bug), not a real biological zero.")
            warnings.append(msg)
            log(msg)
    return warnings


def write_status(out_dir, msg):
    try:
        with open(os.path.join(out_dir, "STATUS.txt"), "a") as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
    except OSError as e:
        log(f"WARNING: could not write STATUS.txt: {e}")


# ---------------------------------------------------------------------------
# Core, pure, unit-tested statistics -- shared by both species.
# ---------------------------------------------------------------------------
def carrier_counts(unit_sets):
    """unit_sets: dict person_id -> set(allele_id). Returns Counter(allele_id -> n_carriers)."""
    c = Counter()
    for s in unit_sets.values():
        for a in s:
            c[a] += 1
    return c


def classify_recurrence(counts):
    """counts: Counter/dict allele_id -> n_carriers (unrelated people). Returns raw (unmasked)
    dict eq1/eq2/gt2/ge20 -> n_alleles. gt2 and ge20 overlap by design -- see module docstring."""
    eq1 = eq2 = gt2 = ge20 = 0
    for n in counts.values():
        if n == 1:
            eq1 += 1
        elif n == 2:
            eq2 += 1
        else:
            gt2 += 1  # n >= 3, includes n >= 20
        if n >= 20:
            ge20 += 1
    return {"eq1": eq1, "eq2": eq2, "gt2": gt2, "ge20": ge20}


def incidence_freqs(unit_sets):
    """dict person_id -> set(allele_id) form of 04_allele_saturation.incidence_freqs (which takes
    a list of sets) -- thin adapter, not a re-derivation."""
    return m04().incidence_freqs(list(unit_sets.values()))


def good_turing_coverage(Q, m, n_incidences):
    """Chao & Jost (2012), Ecology 93:2533-2547 -- incidence-data sample coverage estimator.
    Q: Counter mapping occurrence-count -> number of alleles with that occurrence count (Q[1]=f1
    uniques, Q[2]=f2 duplicates). m: number of sampling units. n_incidences: total number of
    (unit, allele) incidences = sum(occ.values()) in 04's incidence_freqs() terms.
    C_hat = 1 - (Q1/n) * [(m-1)Q1 / ((m-1)Q1 + 2Q2)]           if Q2 > 0
          = 1 - (Q1/n) * [(m-1)(Q1-1) / ((m-1)(Q1-1) + 2)]     if Q2 == 0, Q1 > 1
          = 1                                                   if Q2 == 0, Q1 <= 1 (no unseen mass)
    """
    Q1, Q2 = Q.get(1, 0), Q.get(2, 0)
    if n_incidences <= 0:
        return float("nan")
    if Q2 > 0:
        corr = ((m - 1) * Q1) / ((m - 1) * Q1 + 2 * Q2)
    elif Q1 > 1:
        corr = ((m - 1) * (Q1 - 1)) / ((m - 1) * (Q1 - 1) + 2)
    else:
        corr = 0.0
    return 1.0 - (Q1 / n_incidences) * corr


def chao2_se(S_obs, Q, m):
    """Incidence-based Chao2 standard error (Chao 1987, Biometrics 43:783-791; formulas as
    reproduced in Colwell's EstimateS User's Guide sect. "Variance of Chao2"). Two branches,
    matching 04_allele_saturation.chao2()'s own Q2>0 / Q2==0 split so the SE is always consistent
    with whichever point-estimate branch chao2() used."""
    Q1, Q2 = Q.get(1, 0), Q.get(2, 0)
    if m <= 1:
        return 0.0
    A = (m - 1) / m
    if Q2 > 0:
        r = Q1 / Q2
        var = Q2 * (0.5 * A * r ** 2 + A ** 2 * r ** 3 + 0.25 * A ** 2 * r ** 4)
    else:
        s_c = S_obs + A * Q1 * (Q1 - 1) / 2.0
        var = (A * Q1 * (Q1 - 1) / 2.0
               + A ** 2 * Q1 * (2 * Q1 - 1) ** 2 / 4.0
               - (A ** 2 * Q1 ** 4 / (4.0 * s_c) if s_c > 0 else 0.0))
    return math.sqrt(max(var, 0.0))


def make_orders(n_units, n_permutations, seed):
    """Precomputed list of permutation index-arrays over range(n_units), one draw per replicate.
    Sharing ONE such list between the KIR and HLA extraction for a given ancestry (same n_units,
    same pids order feeding unit_sets_list) is what guarantees identical subsampling paths across
    species -- the brief's WS-B requirement."""
    rng = np.random.default_rng(seed)
    return [rng.permutation(n_units) for _ in range(n_permutations)]


def rarefaction_curves_fixed_orders(unit_sets_list, orders, thresholds=KIR_THRESHOLDS):
    """Adapted from 39_saturation_by_ancestry.rarefaction_threshold_curves: identical inner loop,
    generalized to take PRECOMPUTED permutation orders (see make_orders()) instead of drawing its
    own, so two calls with the same `orders` (one per species) walk exactly the same people in
    exactly the same sequence. unit_sets_list: list of sets, index-aligned with `orders`' values
    (i.e. unit_sets_list[i] is person i's allele set for whichever gene/level/species is being
    curved). Returns {'distinct': (n_perm, m) array, thresholds[k]: (n_perm, m) array}."""
    m = len(unit_sets_list)
    n_perm = len(orders)
    distinct = np.zeros((n_perm, m))
    thr = {k: np.zeros((n_perm, m)) for k in thresholds}
    for p, order in enumerate(orders):
        carrier_count = {}
        n_distinct = 0
        n_reached = {k: 0 for k in thresholds}
        for step, i in enumerate(order):
            for allele in unit_sets_list[i]:
                c = carrier_count.get(allele, 0)
                if c == 0:
                    n_distinct += 1
                c += 1
                carrier_count[allele] = c
                if c in n_reached:
                    n_reached[c] += 1
            distinct[p, step] = n_distinct
            for k in thresholds:
                thr[k][p, step] = n_reached[k]
    out = {"distinct": distinct}
    out.update(thr)
    return out


def summarize_curve(arr):
    """arr: (n_permutations, m). Returns steps(1..m), mean, lo2.5, hi97.5 -- identical convention
    to 39_saturation_by_ancestry.summarize_curve."""
    m = arr.shape[1]
    steps = np.arange(1, m + 1)
    return steps, arr.mean(axis=0), np.percentile(arr, 2.5, axis=0), np.percentile(arr, 97.5, axis=0)


def class_curves_from_thresholds(curves):
    """curves: output of rarefaction_curves_fixed_orders (KIR_THRESHOLDS=[1,2,3,20]). Returns
    {'eq1':arr,'eq2':arr,'gt2':arr,'ge20':arr} via differencing of the >=k "reached" curves --
    e.g. eq1(N) = (#alleles with >=1 carrier at N) - (#alleles with >=2 carriers at N)."""
    return {
        "eq1": curves[1] - curves[2],
        "eq2": curves[2] - curves[3],
        "gt2": curves[3],          # >=3 carriers, includes >=20
        "ge20": curves[20],
    }


def equal_n_slope(steps, mean, n_star, window_frac=0.9):
    """Linear fit of the mean rarefaction curve over [window_frac*n_star, n_star], reported as new
    alleles per next 1,000 people -- 39_saturation_by_ancestry's own "primary, cross-ancestry-
    comparable" equal-N slope, generalized here to any curve (overall or one recurrence class)."""
    if n_star is None or n_star < 2 or n_star > steps[-1]:
        return float("nan")
    lo = max(1, int(round(window_frac * n_star)))
    mask = (steps >= lo) & (steps <= n_star)
    if mask.sum() < 2:
        return float("nan")
    slope = np.polyfit(steps[mask], mean[mask], 1)[0]
    return float(slope * 1000)


def pick_n_star(n_people_by_ancestry, min_people=MIN_PEOPLE_PER_ANCESTRY, exclude=("MID",)):
    """39's own convention: the largest N every well-powered (>=min_people), non-excluded
    ancestry reaches -- i.e. min(count) over those ancestries. Returns None if none qualify."""
    powered = {a: n for a, n in n_people_by_ancestry.items()
               if a not in exclude and n >= min_people}
    return min(powered.values()) if powered else None


# ---------------------------------------------------------------------------
# KIR identity extraction: genomic (raw consensus) needs only hap{1,2}.gtf.gz. cds/protein need
# a (contig, gene) join against hap{1,2}/cds.fa.gz -- see module docstring "KIR IDENTITY
# EXTRACTION" for exactly what was checked and what could not be verified without VM access.
# ---------------------------------------------------------------------------
_KIR_GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
_KIR_CONSENSUS_RE = re.compile(r'consensus "([^"]+)"')
_KIR_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
_KIR_CDS_DIST_RE = re.compile(r'cds_distance (\d+)')
_KIR_CDS_MUT_RE = re.compile(r'cds_mut "([^"]*)"')
_KIR_COPY_SUFFIX_RE = re.compile(r'\.(\d+)$')


def parse_hap_gtf_full(gtf_path, kir):
    """One-pass parse of a KIR hap GTF, combining 41_kir_pilot.parse_hap_gtf()'s own extraction
    (gene, consensus, cds_distance, cds_mut, is_novel, is_undetermined) -- reusing its
    `classify_novelty_tier()` verbatim for the novelty_tier field, not re-deriving that logic --
    with two EXTRA fields (`contig` = GTF column 1, `copy_index` = the gene_id's ".N" suffix per
    reference/IMMUANNOT_GTF_SPEC.md part A) that 41's own parser has no reason to expose but a
    cds.fa.gz join (03_novel_alleles.py's own convention) needs. A second regex pass over the same
    9-column transcript rows, not a re-derivation of 41's novelty CLASSIFICATION (that part is a
    direct call to `kir.classify_novelty_tier`)."""
    rows = []
    if not gtf_path or not os.path.exists(gtf_path):
        return rows
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            gene_m = _KIR_GENE_NAME_RE.search(attrs)
            cons_m = _KIR_CONSENSUS_RE.search(attrs)
            if not gene_m or not cons_m:
                continue
            consensus = cons_m.group(1)
            geneid_m = _KIR_GENE_ID_RE.search(attrs)
            copy_m = _KIR_COPY_SUFFIX_RE.search(geneid_m.group(1)) if geneid_m else None
            cds_dist_m = _KIR_CDS_DIST_RE.search(attrs)
            cds_mut_m = _KIR_CDS_MUT_RE.search(attrs)
            row = {
                "gene": gene_m.group(1), "consensus": consensus, "contig": fields[0],
                "copy_index": int(copy_m.group(1)) if copy_m else 1,
                "cds_distance": int(cds_dist_m.group(1)) if cds_dist_m else None,
                "cds_mut": cds_mut_m.group(1) if cds_mut_m else None,
                "is_novel": consensus.rstrip('"').endswith("new") and consensus != "undetermined",
                "is_undetermined": consensus == "undetermined",
            }
            row["novelty_tier"] = kir.classify_novelty_tier(row)
            rows.append(row)
    return rows


def resolve_kir_cds_sequences(rows, cds_fa_path, m03mod):
    """Joins already-parsed KIR transcript rows (from parse_hap_gtf_full, carrying contig/gene/
    consensus/copy_index) to one hap's `cds.fa.gz` (m03mod.parse_cds_fasta, reused verbatim),
    following 03_novel_alleles.match_novel_rows()'s own unambiguous-join contract: a gene with
    exactly 1 row for a given (contig, gene) on this hap matches unambiguously against
    cds.fa.gz's exactly-1-candidate list for that key; anything else (copy_index>1 on this contig,
    or 0/>1 candidates in cds.fa.gz) is left UNRESOLVED and counted, never guessed.
    Returns (seq_by_key: {(gene, consensus): seq}, stats: Counter -- aggregate-only:
    n_join_rows, n_matched, n_ambiguous_copy, n_missing_cds_record, n_cds_fasta_missing)."""
    stats = Counter()
    seq_by_key = {}
    stats["n_join_rows"] = len(rows)
    if not rows:
        return seq_by_key, stats
    if not cds_fa_path or not os.path.exists(cds_fa_path):
        stats["n_cds_fasta_missing"] = len(rows)
        return seq_by_key, stats
    group_sizes = Counter((r["contig"], r["gene"]) for r in rows)
    fasta_index = m03mod.parse_cds_fasta(cds_fa_path)
    for r in rows:
        gkey = (r["contig"], r["gene"])
        if group_sizes[gkey] > 1:
            stats["n_ambiguous_copy"] += 1
            continue
        rest = f"{r['contig']}_{r['gene']}"
        candidates = fasta_index.get(rest, [])
        if len(candidates) == 0:
            stats["n_missing_cds_record"] += 1
            continue
        if len(candidates) > 1:
            stats["n_ambiguous_copy"] += 1
            continue
        _i, seq = candidates[0]
        if not seq:
            stats["n_missing_cds_record"] += 1
            continue
        seq_by_key[(r["gene"], r["consensus"])] = seq
        stats["n_matched"] += 1
    return seq_by_key, stats


def kir_gene_protein_covered(kir_ref, gene):
    """True iff `kir_ref` (an `m24().RefIndex` built from IPD-KIR's CDSseq/*.fa.gz) has at least
    one usable, in-catalogue protein record for `gene`. Two distinct ways a gene can fail this that
    a bare `gene in kir_ref.genes()` check would miss (2026-09-27 follow-up fix, S04 coordinator
    item -- VM smoke test flagged KIR2DL2/KIR2DL5B/KIR2DP1 still exceeding 100%):
      - `gene in kir_ref.genes()` only means the catalogue has >=1 CDS record for it (RefIndex.cds),
        NOT a usable PROTEIN record -- a pseudogene like KIR2DP1, whose reference alleles are all
        frameshifted/premature-stop, adds every record to `.cds` but NEVER to `.prot`
        (`RefIndex.add()`'s own gate), so `.genes()` says "covered" while `.prot[gene]` stays empty
        forever -- every observed translation would then look "not catalogued", including
        genuinely-known ones, which is exactly the kind of runaway protein_novel count this
        function exists to prevent.
      - gene attribution inside `kir_ref` always comes from each FASTA record's own header
        (`RefIndex.add()`'s `g` argument, populated by `iter_fasta()`/`parse_ref_header()` from the
        `*`-bearing token in the header line), never from the CDSseq file's own name -- so this
        function is filename-agnostic by construction and does not care whether IPD-KIR happens to
        bundle e.g. KIR2DL5A/KIR2DL5B or KIR2DL2/KIR2DL3 alleles into one file or several. If a
        gene still comes back uncovered, the header-level allele nomenclature genuinely doesn't
        carry that exact KIR_GENES spelling (e.g. undifferentiated "KIR2DL5" headers instead of
        "KIR2DL5A"/"KIR2DL5B") -- a real catalogue-naming gap, not a file-skip artifact (see
        `load_refdata(..., genes_needed=None)` in `_identity_worker_init`/`main()`, which no longer
        pre-filters CDSseq files by filename stem for KIR, precisely so a bundled file can't hide a
        gene's own headers from this check)."""
    return kir_ref is not None and gene in kir_ref.genes() and bool(kir_ref.prot.get(gene))


def kir_protein_catalogue_status(kir_ref, kir_genes):
    """Per-gene {gene: (covered: bool, reason: str)} snapshot of whether `kir_ref` gives this KIR
    gene a usable protein catalogue -- computed ONCE in the main process (not per-worker; this is
    bookkeeping for NA-row/QC decisions in main(), not per-call identity assignment, so it doesn't
    need to survive a multiprocessing pickle). reason in {"ok", "no_refdata",
    "gene_not_in_catalogue", "no_catalogued_protein_entries"} -- see kir_gene_protein_covered's
    docstring for what each uncovered reason means."""
    status = {}
    for gene in kir_genes:
        if kir_ref is None:
            status[gene] = (False, "no_refdata")
        elif gene not in kir_ref.genes():
            status[gene] = (False, "gene_not_in_catalogue")
        elif not kir_ref.prot.get(gene):
            status[gene] = (False, "no_catalogued_protein_entries")
        else:
            status[gene] = (True, "ok")
    return status


def build_person_kir_identity(pid, kir_outroot, kir, m03mod, m24mod, kir_ref=None):
    """One person's contribution across all 5 KIR identity levels (see module docstring).
    Returns (pid, per_level: {level: {gene: set(id)}}, qc: Counter, n_hap_seen: int,
    n_hap_cds_present: int) -- picklable (plain dicts/sets/Counter), safe for multiprocessing.

    `kir_ref`: optional `m24().RefIndex` built from IPD-KIR's own CDSseq/*.fa.gz snapshot (same
    directory 24_novelty_by_field.py's HLA RefIndex globs -- confirmed present for all 17 KIR genes,
    see 41_kir_pilot.classify_novelty_tier's own docstring). BUG FIX 2026-09-26/27
    (identity-granularity mismatch, S04 coordinator item): the `protein`/`protein_novel` tracks
    previously hashed the raw observed CDS translation for EVERY call, known or not, so two people
    carrying the textbook-known same allele could land on different hashes from ordinary per-sample
    sequencing noise -- richer (finer-grained) than the curated `genomic` identity, which is exactly
    backwards for a protein-vs-genomic comparison and is why `pct_novel_protein` could exceed 100%.
    Mirroring HLA's own `prot_id` scheme (a curated name for a known/catalogued protein, a hash only
    for a genuinely uncatalogued one) fixes this for genes the catalogue actually covers
    (`kir_gene_protein_covered()` -- checked per-call here, not just per-gene, since it's cheap and
    keeps this function self-contained): a translated protein already in `kir_ref.prot[gene]`
    collapses to ONE catalogue-name-derived id across everyone; a genuinely uncatalogued one gets a
    hash id AND counts as `protein_novel`. 2026-09-27 tightening: for a gene the catalogue does NOT
    cover at all (`kir_gene_protein_covered()` False -- no refdata loaded, the gene has no catalogue
    records under ANY header spelling, or every catalogue record for it is a nonfunctional
    pseudogene entry never added to `.prot`), this function now adds NOTHING to `protein`/
    `protein_novel` for that gene -- no hash fallback, ever, for an uncovered gene (the previous
    41-novelty_tier-heuristic-plus-hash fallback was itself the source of the >100% bug for
    KIR2DL2/KIR2DL5B/KIR2DP1 on the 2026-09-26 VM smoke test: it silently produced hash ids finer
    than `genomic` for exactly the genes whose catalogue coverage was incomplete). Callers (`main()`)
    are responsible for exporting an explicit `"NA"` row (with a reason, from
    `kir_protein_catalogue_status()`) for any (gene, protein/protein_novel) pair left unpopulated
    this way, never a silent zero."""
    per_level = {lvl: defaultdict(set) for lvl in LEVELS}
    qc = Counter()
    n_hap_seen = n_hap_cds_present = 0
    for hap in ("hap1", "hap2"):
        gtf_path = os.path.join(kir_outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
        rows = parse_hap_gtf_full(gtf_path, kir)
        kir_rows = [r for r in rows if r["gene"] in kir.KIR_GENES]
        if not kir_rows:
            continue
        n_hap_seen += 1
        cds_fa_path = os.path.join(kir_outroot, str(pid), "immuannot_output", hap, "cds.fa.gz")
        if os.path.exists(cds_fa_path):
            n_hap_cds_present += 1
        seq_by_key, stats = resolve_kir_cds_sequences(kir_rows, cds_fa_path, m03mod)
        qc.update(stats)
        for r in kir_rows:
            gene, tier, consensus = r["gene"], r["novelty_tier"], r["consensus"]
            if tier != "undetermined":
                per_level["genomic"][gene].add(consensus)
            if tier not in ("known", "undetermined"):
                per_level["any_novel"][gene].add(consensus)
            seq = seq_by_key.get((gene, consensus))
            if seq and tier != "undetermined":
                cds_id = f"{gene}_cds_{m24mod.sha8(seq)}"
                per_level["cds"][gene].add(cds_id)
                if not kir_gene_protein_covered(kir_ref, gene):
                    qc["n_protein_gene_uncovered_calls"] += 1
                    continue  # never guess a finer-than-genomic id for an uncovered gene
                protein = m24mod.protein_info(seq)["protein"]
                catalog_names = kir_ref.prot[gene].get(protein)
                if catalog_names:
                    prot_id = f"{gene}_{sorted(catalog_names)[0]}"
                else:
                    prot_id = f"{gene}_prot_{m24mod.sha8(protein)}"
                per_level["protein"][gene].add(prot_id)
                if not catalog_names:
                    per_level["protein_novel"][gene].add(prot_id)
    per_level = {lvl: dict(d) for lvl, d in per_level.items()}
    return pid, per_level, qc, n_hap_seen, n_hap_cds_present


_id_kir = _id_m03 = _id_m24 = _id_kir_ref = None


def _identity_worker_init(refdata=None):
    global _id_kir, _id_m03, _id_m24, _id_kir_ref
    _id_kir, _id_m03, _id_m24 = m41(), m03(), m24()
    _id_kir_ref = None
    if refdata is not None:
        try:
            # genes_needed=None (2026-09-27 fix): load_refdata's file-level pre-filter matches a
            # CDSseq/*.fa.gz file's own FILENAME stem against genes_needed and skips the whole file
            # if it doesn't match -- but gene attribution for every record inside a file comes from
            # that record's HEADER (RefIndex.add() via iter_fasta/parse_ref_header), not the
            # filename. If IPD-KIR bundles e.g. KIR2DL2+KIR2DL3, or KIR2DL5A+KIR2DL5B, alleles into
            # one file whose stem matches only ONE (or neither) of those exact KIR_GENES spellings,
            # a genes_needed={KIR_GENES} filter would skip that file entirely and silently starve
            # BOTH genes of catalogue coverage -- confirmed as one of the two mechanisms behind the
            # 2026-09-26 VM smoke test's KIR2DL2/KIR2DL5B/KIR2DP1 pct_novel_protein>100% warnings.
            # The 17 KIR genes' CDSseq files are small; loading all of them (instead of filtering by
            # filename) costs a negligible amount of I/O and removes this failure mode entirely.
            _id_kir_ref = _id_m24.load_refdata(refdata, genes_needed=None)
        except SystemExit as e:
            log(f"WARNING: could not load KIR protein catalogue from {refdata!r} ({e}); "
                f"'protein'/'protein_novel' will be exported as 'NA' (reason=no_refdata) for "
                f"every KIR gene -- see kir_protein_catalogue_status()/kir_gene_protein_covered().")
            _id_kir_ref = None


def _identity_worker_task(args):
    pid, kir_outroot = args
    return build_person_kir_identity(pid, kir_outroot, _id_kir, _id_m03, _id_m24, _id_kir_ref)


def build_kir_identity_sets(pids, kir_outroot, workers, refdata=None):
    """Orchestrates build_person_kir_identity over all `pids` (optionally multiprocessed, mirroring
    43_kir_full_aggregate.parse_all's Pool(initializer=...) pattern). Returns
    (sets_by_level: {level: {gene: {pid: set(id)}}}, qc: Counter, cds_available: bool,
    n_hap_seen: int, n_hap_cds_present: int). `cds_available` is False iff cds.fa.gz was found for
    ZERO haplotypes across the entire cohort (a hard failure mode -- e.g. the file was cleaned up
    or never written for the KIR run -- not ordinary per-person sparsity).

    `refdata`: passed through to each worker's own `_identity_worker_init` (not shared as a live
    object across process boundaries -- each worker loads its own KIR RefIndex, same pattern this
    file already uses for m41/m03/m24). None (the default) reproduces the pre-fix hash-only /
    novelty_tier-heuristic behavior -- see build_person_kir_identity's docstring."""
    sets_by_level = {lvl: defaultdict(dict) for lvl in LEVELS}
    qc_total = Counter()
    n_hap_seen_total = n_hap_cds_present_total = 0
    if workers <= 1:
        _identity_worker_init(refdata)
        results = (_identity_worker_task((pid, kir_outroot)) for pid in pids)
    else:
        pool = mp.Pool(processes=workers, initializer=_identity_worker_init, initargs=(refdata,))
        results = pool.imap_unordered(_identity_worker_task,
                                       [(pid, kir_outroot) for pid in pids], chunksize=32)
    for pid, per_level, qc, n_seen, n_cds in results:
        for lvl, gene_sets in per_level.items():
            for gene, ids in gene_sets.items():
                sets_by_level[lvl][gene][pid] = ids
        qc_total.update(qc)
        n_hap_seen_total += n_seen
        n_hap_cds_present_total += n_cds
    if workers > 1:
        pool.close()
        pool.join()
    cds_available = n_hap_cds_present_total > 0
    return sets_by_level, qc_total, cds_available, n_hap_seen_total, n_hap_cds_present_total


# ---------------------------------------------------------------------------
# HLA-side extraction: reuses 39_saturation_by_ancestry.build_labeled_calls() (the full
# field_class/seq_class/prot_id/cds_id pipeline) -- no re-derivation of sequence matching. The
# genomic-level identity (full, untruncated normalized consensus name) is the one new column added
# here, via 24_novelty_by_field.normalize_allele_name (reused, not re-derived).
# ---------------------------------------------------------------------------
def add_hla_genomic_id(calls, m24mod):
    """Adds a 'genomic_id' column: the full normalized consensus name (all colon fields kept),
    cached per distinct consensus string (there are far fewer distinct strings than rows)."""
    cache = {}

    def f(cons):
        if not isinstance(cons, str):
            return None
        if cons not in cache:
            _, norm = m24mod.normalize_allele_name(cons)
            cache[cons] = norm
        return cache[cons]

    calls = calls.copy()
    calls["genomic_id"] = calls["consensus"].map(f) if "consensus" in calls.columns else None
    return calls


def hla_level_mask(calls, level):
    if level in ("genomic", "cds", "protein"):
        return calls["keep_clean"]
    if level == "any_novel":
        return calls["keep_clean"] & (calls["field_class"] != "known")
    if level == "protein_novel":
        return (calls["keep_clean"] & (calls["field_class"] == "f2_protein")
                & (calls["seq_class"] == "novel_protein"))
    raise ValueError(level)


def hla_id_col(level):
    return {"genomic": "genomic_id", "cds": "cds_id", "protein": "prot_id",
            "any_novel": "genomic_id", "protein_novel": "prot_id"}[level]


def hla_person_gene_sets(calls, level, genes_bare):
    """Returns {gene: {person_id: set(allele_id)}}, restricted to `genes_bare` and `level`."""
    id_col = hla_id_col(level)
    mask = hla_level_mask(calls, level) & calls["gene_b"].isin(genes_bare) & calls[id_col].notna()
    sub = calls.loc[mask, ["person_id", "gene_b", id_col]]
    out = defaultdict(lambda: defaultdict(set))
    for pid, gene, aid in zip(sub["person_id"], sub["gene_b"], sub[id_col]):
        out[gene][pid].add(aid)
    return out


# ---------------------------------------------------------------------------
# Joint unrelated set (WS-B "same unrelated people who have both").
# ---------------------------------------------------------------------------
def joint_unrelated_set(kir_pids, hla_pids, relatedness_path, kin_min):
    """Intersect the two species' available-person pools FIRST, then run ONE greedy_unrelated()
    pass over that intersection -- so both species get the exact same unrelated set (running
    greedy_unrelated separately on each species' own full pool would not guarantee this, since the
    kinship graph is restricted to whichever id set is passed in)."""
    kir = m41()
    common = sorted(set(map(str, kir_pids)) & set(map(str, hla_pids)))
    pairs = kir.load_relatedness_pairs(relatedness_path)
    kept, removed = kir.greedy_unrelated(common, pairs, kin_min=kin_min)
    return sorted(kept), len(removed), len(common)


# ---------------------------------------------------------------------------
# Table builders (apply disclosure masking here, at the export boundary).
# ---------------------------------------------------------------------------
def build_recurrence_row(gene, ancestry, level, species, unit_sets):
    counts = carrier_counts(unit_sets)
    cls = classify_recurrence(counts)
    return {"gene": gene, "ancestry": ancestry, "level": level, "species": species,
            "n_people": len(unit_sets), "n_distinct_alleles": suppressed(len(counts)),
            "eq1": suppressed(cls["eq1"]), "eq2": suppressed(cls["eq2"]),
            "gt2": suppressed(cls["gt2"]), "ge20": suppressed(cls["ge20"])}


def build_coverage_chao2_row(gene, ancestry, level, species, unit_sets, extrapolate_2n=True):
    S_obs, m, Q = incidence_freqs(unit_sets)
    n_incidences = sum(k * v for k, v in Q.items())
    coverage = good_turing_coverage(Q, m, n_incidences)
    chao = m04().chao2(S_obs, Q, m)
    se = chao2_se(S_obs, Q, m)
    row = {"gene": gene, "ancestry": ancestry, "level": level, "species": species,
           "n_people": m, "s_obs": S_obs,
           "good_turing_coverage": round(coverage, 4) if coverage == coverage else "",
           "chao2": round(chao, 1), "chao2_se": round(se, 1),
           "chao2_undetected_f0hat": round(chao - S_obs, 1)}
    q1, q2 = Q.get(1, 0), Q.get(2, 0)
    row["q1"] = suppressed(q1) if q1 < SUPPRESS_BELOW else str(q1)
    row["q2"] = suppressed(q2) if q2 < SUPPRESS_BELOW else str(q2)
    if q1 < SUPPRESS_BELOW or q2 < SUPPRESS_BELOW:
        row["q1"] = row["q2"] = ""  # blank both together, never one alone (back-reveal risk)
    if extrapolate_2n:
        row["chao_new_by_2n"] = round(
            m04().chao2_extrapolate(S_obs, Q, m, 2 * m) - S_obs, 1) if m > 0 else ""
    return row


def build_curve_rows(gene, ancestry, level, species, unit_sets_list, orders, n_star, stride):
    curves = rarefaction_curves_fixed_orders(unit_sets_list, orders, KIR_THRESHOLDS)
    class_curves = class_curves_from_thresholds(curves)
    steps, mean_d, lo_d, hi_d = summarize_curve(curves["distinct"])
    means = {k: summarize_curve(v)[1] for k, v in class_curves.items()}
    rows = []
    idxs = list(range(0, len(steps), max(1, stride)))
    if idxs[-1] != len(steps) - 1:
        idxs.append(len(steps) - 1)
    for i in idxs:
        rows.append({
            "gene": gene, "ancestry": ancestry, "level": level, "species": species,
            "n": int(steps[i]), "mean_distinct": round(float(mean_d[i]), 2),
            "lo2_5": round(float(lo_d[i]), 2), "hi97_5": round(float(hi_d[i]), 2),
            "mean_eq1": round(float(means["eq1"][i]), 2),
            "mean_eq2": round(float(means["eq2"][i]), 2),
            "mean_gt2": round(float(means["gt2"][i]), 2),
            "mean_ge20": round(float(means["ge20"][i]), 2),
        })
    slope_rows = []
    if n_star is not None:
        slope_rows.append({
            "gene": gene, "ancestry": ancestry, "level": level, "species": species,
            "n_star": n_star, "curve": "distinct",
            "slope_per_1000": round(equal_n_slope(steps, mean_d, n_star), 2),
            "mean_at_n_star": round(float(np.interp(n_star, steps, mean_d)), 2),
        })
        for cname, arr in means.items():
            slope_rows.append({
                "gene": gene, "ancestry": ancestry, "level": level, "species": species,
                "n_star": n_star, "curve": cname,
                "slope_per_1000": round(equal_n_slope(steps, arr, n_star), 2),
                "mean_at_n_star": round(float(np.interp(n_star, steps, arr)), 2),
            })
    return rows, slope_rows


def build_na_rows(gene, level, species, pids_by_ancestry, n_star_by_ancestry):
    """A level that could not be matched for a species (see module docstring "KIR IDENTITY
    EXTRACTION") is exported as the literal string "NA" everywhere a real value would go -- never
    blank, never 0, never silently omitted, and never a censored "<20" (a different, disjoint
    meaning). `n_people` is still reported (a real, non-disclosive ancestry-group size, not
    derived from the unavailable sequence data)."""
    rec_rows, curve_rows, slope_rows, cov_rows = [], [], [], []
    for anc, pids in pids_by_ancestry.items():
        rec_rows.append({"gene": gene, "ancestry": anc, "level": level, "species": species,
                          "n_people": len(pids), "n_distinct_alleles": "NA",
                          "eq1": "NA", "eq2": "NA", "gt2": "NA", "ge20": "NA"})
        cov_rows.append({"gene": gene, "ancestry": anc, "level": level, "species": species,
                          "n_people": len(pids), "s_obs": "NA", "good_turing_coverage": "NA",
                          "chao2": "NA", "chao2_se": "NA", "chao2_undetected_f0hat": "NA",
                          "q1": "NA", "q2": "NA", "chao_new_by_2n": "NA"})
        curve_rows.append({"gene": gene, "ancestry": anc, "level": level, "species": species,
                            "n": len(pids), "mean_distinct": "NA", "lo2_5": "NA", "hi97_5": "NA",
                            "mean_eq1": "NA", "mean_eq2": "NA", "mean_gt2": "NA",
                            "mean_ge20": "NA"})
        n_star = n_star_by_ancestry.get(anc)
        for curve_name in ("distinct", "eq1", "eq2", "gt2", "ge20"):
            slope_rows.append({"gene": gene, "ancestry": anc, "level": level, "species": species,
                                "n_star": n_star if n_star is not None else "NA",
                                "curve": curve_name, "slope_per_1000": "NA",
                                "mean_at_n_star": "NA"})
    return rec_rows, curve_rows, slope_rows, cov_rows


# ---------------------------------------------------------------------------
# Orchestration: for one gene x level, build recurrence/curve/slope/coverage rows for ALL + each
# ancestry, using a SHARED pid ordering/permutation-order set per ancestry (species-agnostic).
# ---------------------------------------------------------------------------
def run_gene_level(gene, level, species, full_gene_sets, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, stride, extrapolate_2n):
    """full_gene_sets: {person_id: set(allele_id)} already restricted to `gene`/`level` (ALL
    people, not yet sliced by ancestry). Returns the four output-row lists for this
    (gene, level, species) across {ALL + each ancestry}."""
    rec_rows, curve_rows, slope_rows, cov_rows = [], [], [], []
    for anc, pids in pids_by_ancestry.items():   # "ALL" is included as a pseudo-ancestry
        orders = orders_by_ancestry[anc]
        unit_sets = {p: full_gene_sets.get(p, set()) for p in pids}
        rec_rows.append(build_recurrence_row(gene, anc, level, species, unit_sets))
        cov_rows.append(build_coverage_chao2_row(gene, anc, level, species, unit_sets,
                                                  extrapolate_2n))
        unit_sets_list = [unit_sets[p] for p in pids]
        n_star = n_star_by_ancestry.get(anc)
        c_rows, s_rows = build_curve_rows(gene, anc, level, species, unit_sets_list, orders,
                                           n_star, stride)
        curve_rows.extend(c_rows)
        slope_rows.extend(s_rows)
    return rec_rows, curve_rows, slope_rows, cov_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kir-outroot", default=os.path.expanduser("~/pipeline_outputs_kir"))
    ap.add_argument("--hla-table", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--hla-people-outroot", default=DEFAULT_HLA_PEOPLE_OUTROOT,
                     help="per-person dir containing <pid>/immuannot_output/hap{1,2}/cds.fa.gz for "
                          "HLA sequence matching (39_saturation_by_ancestry.build_labeled_calls's "
                          "own --outroot / DEFAULT_PEOPLE_OUTROOT convention -- NOT the same as "
                          "--hla-table's parent dir; BUG FIX 2026-09-26: this file previously "
                          "hardcoded '~/pipeline_outputs' here, one directory short of where "
                          "cds.fa.gz actually lives, so match_sequences() silently matched zero "
                          "depth-2/3 calls and every HLA protein_novel S_obs came out 0 -- see "
                          "sprints/S04_kir_recurrence_style_share/LOG.md).")
    ap.add_argument("--cohort-membership", default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table", default=os.path.expanduser(
        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
        "samples_relatedness.tsv"))
    ap.add_argument("--refdata", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s04/results/44"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=4, help="threads for HLA sequence matching (39's match_sequences)")
    ap.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS_DEFAULT)
    ap.add_argument("--curve-stride", type=int, default=CURVE_STRIDE_DEFAULT)
    ap.add_argument("--seed", type=int, default=44)
    ap.add_argument("--kin-min", type=float, default=0.0442)
    ap.add_argument("--strict-threshold", type=float, default=0.9)
    ap.add_argument("--limit", type=int, default=None, help="smoke test: cap discovered KIR people")
    ap.add_argument("--no-extrapolate-2n", dest="extrapolate_2n", action="store_false")
    ap.add_argument("--skip-relatedness", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    write_status(args.out_dir, f"START outroot={args.kir_outroot} out_dir={args.out_dir}")
    t0 = time.perf_counter()

    kir = m41()
    kir_agg = m43()
    m39mod = m39()
    m24mod = m24()

    log("[44] discovering KIR people ...")
    kir_pids = kir_agg.discover_people(args.kir_outroot, args.limit)
    write_status(args.out_dir, f"discovered {len(kir_pids)} KIR people")

    log("[44] loading HLA labeled calls (39's pipeline) ...")
    hla_args = argparse.Namespace(
        table1=args.hla_table, limit=args.limit, cohort_membership=args.cohort_membership,
        relatedness_table=args.relatedness_table, kin_min=args.kin_min,
        strict_threshold=args.strict_threshold, skip_relatedness=args.skip_relatedness,
        refdata=args.refdata, outroot=args.hla_people_outroot, threads=args.threads)
    calls, ref_catalogue_size, mstats, n_removed = m39mod.build_labeled_calls(hla_args)
    calls = add_hla_genomic_id(calls, m24mod)
    hla_pids = sorted(calls["person_id"].unique())
    write_status(args.out_dir, f"HLA calls loaded: {len(hla_pids)} people")

    log("[44] computing joint unrelated set (KIR ∩ HLA, one greedy_unrelated pass) ...")
    unrelated_pids, n_removed_joint, n_common = joint_unrelated_set(
        kir_pids, hla_pids, args.relatedness_table, args.kin_min)
    log(f"[44] joint unrelated set: {len(unrelated_pids)}/{n_common} common people "
        f"({n_removed_joint} relatives dropped)")
    write_status(args.out_dir, f"joint unrelated set: {len(unrelated_pids)}")

    ancestry_of = kir_agg.load_ancestry_map(args.cohort_membership)
    pids_by_ancestry_base = defaultdict(list)
    for p in unrelated_pids:
        anc = ancestry_of.get(str(p), "UNASSIGNED")
        if anc in ANCESTRY_ORDER:
            pids_by_ancestry_base[anc].append(p)
    pids_by_ancestry = {"ALL": unrelated_pids, **pids_by_ancestry_base}
    n_people_by_ancestry = {a: len(ps) for a, ps in pids_by_ancestry_base.items()}
    n_star = pick_n_star(n_people_by_ancestry)
    n_star_by_ancestry = {a: (n_star if a != "ALL" else len(unrelated_pids)) for a in pids_by_ancestry}
    log(f"[44] N* (equal-N reference, excl. MID) = {n_star}")

    orders_by_ancestry = {a: make_orders(len(ps), args.n_permutations, args.seed + hash(a) % 10_000)
                          for a, ps in pids_by_ancestry.items()}

    kir_genes = kir.KIR_GENES
    hla_genes = m39mod.CLASSICAL_GENES_BARE

    log("[44] building KIR identity sets (all 5 levels, one pass per person) ...")
    (kir_sets_by_level, kir_qc, kir_cds_available,
     n_hap_seen, n_hap_cds_present) = build_kir_identity_sets(
        unrelated_pids, args.kir_outroot, args.workers, refdata=args.refdata)
    if not kir_cds_available:
        log("[44] WARNING: cds.fa.gz was not found for ANY KIR haplotype -- 'cds'/'protein'/"
            "'protein_novel' will be exported as 'NA' for every KIR row. See module docstring "
            "'KIR IDENTITY EXTRACTION' for what this project's other scripts (41/43) actually "
            "read (only hap{1,2}.gtf.gz) and what could not be verified without VM access.")
    write_status(args.out_dir, f"KIR identity sets built; cds_available={kir_cds_available} "
                                f"(n_hap_seen={n_hap_seen}, n_hap_cds_present={n_hap_cds_present})")

    # 2026-09-27 fix: per-gene KIR protein-catalogue coverage (see kir_gene_protein_covered's/
    # kir_protein_catalogue_status's docstrings) -- computed once in the main process, purely for
    # NA-row/QC bookkeeping (the actual per-call identity decision already happened inside each
    # worker via its OWN kir_ref, built the same way). Loaded a second time here rather than passed
    # from a worker to avoid pickling a RefIndex (its nested defaultdict(lambda: ...) isn't
    # picklable) across a process boundary for what is a cheap, one-time reference-data parse.
    log("[44] loading KIR protein catalogue for per-gene coverage status ...")
    try:
        kir_ref_for_status = m24mod.load_refdata(args.refdata, genes_needed=None)
    except SystemExit as e:
        log(f"[44] WARNING: could not load KIR protein catalogue from {args.refdata!r} ({e}); "
            f"every KIR gene's protein/protein_novel level will be NA (reason=no_refdata).")
        kir_ref_for_status = None
    kir_protein_status = kir_protein_catalogue_status(kir_ref_for_status, kir.KIR_GENES)
    uncovered = {g: reason for g, (ok, reason) in kir_protein_status.items() if not ok}
    if uncovered:
        log(f"[44] WARNING: {len(uncovered)}/{len(kir.KIR_GENES)} KIR genes have no usable protein "
            f"catalogue entry -- their 'protein'/'protein_novel' rows will be exported as 'NA', "
            f"never a hash-based guess: {uncovered}")
    pd.DataFrame([{"gene": g, "protein_catalogue_covered": ok, "reason": reason}
                  for g, (ok, reason) in sorted(kir_protein_status.items())]).to_csv(
        os.path.join(args.out_dir, "kir_protein_catalogue_qc.tsv"), sep="\t", index=False)

    log("[44] building HLA identity sets (all 5 levels) ...")
    hla_sets_by_level = {lvl: hla_person_gene_sets(calls, lvl, hla_genes) for lvl in LEVELS}

    all_rec, all_curve, all_slope, all_cov = [], [], [], []
    for level in LEVELS:
        for gene in kir_genes:
            if level in KIR_LEVELS_NEEDING_CDS_FASTA and not kir_cds_available:
                rec, curve, slope, cov = build_na_rows(gene, level, "kir", pids_by_ancestry,
                                                        n_star_by_ancestry)
            elif level in ("protein", "protein_novel") and not kir_protein_status[gene][0]:
                # Per-gene catalogue-coverage NA (2026-09-27 fix) -- distinct from the cds_available
                # gate above: cds.fa.gz can be present and joined fine (so 'cds' is real) while the
                # IPD-KIR protein catalogue still doesn't cover this specific gene (pseudogene with
                # no functional reference protein, or a naming gap) -- see
                # kir_protein_catalogue_status()/kir_gene_protein_covered().
                rec, curve, slope, cov = build_na_rows(gene, level, "kir", pids_by_ancestry,
                                                        n_star_by_ancestry)
            else:
                full_gene_sets = kir_sets_by_level[level].get(gene, {})
                rec, curve, slope, cov = run_gene_level(
                    gene, level, "kir", full_gene_sets, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, args.curve_stride, args.extrapolate_2n)
            all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
        for gene in hla_genes:
            full_gene_sets = hla_sets_by_level[level].get(gene, {})
            rec, curve, slope, cov = run_gene_level(
                gene, level, "hla", full_gene_sets, pids_by_ancestry, orders_by_ancestry,
                n_star_by_ancestry, args.curve_stride, args.extrapolate_2n)
            all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
        write_status(args.out_dir, f"level={level} done ({len(all_rec)} recurrence rows so far)")

    sanity_check_coverage(all_cov)

    pd.DataFrame(all_rec).to_csv(os.path.join(args.out_dir, "recurrence_classes.tsv"),
                                  sep="\t", index=False)
    pd.DataFrame(all_curve).to_csv(os.path.join(args.out_dir, "saturation_curves.tsv"),
                                    sep="\t", index=False)
    pd.DataFrame(all_slope).to_csv(os.path.join(args.out_dir, "equal_n_slope.tsv"),
                                    sep="\t", index=False)
    pd.DataFrame(all_cov).to_csv(os.path.join(args.out_dir, "coverage_chao2.tsv"),
                                  sep="\t", index=False)
    pd.DataFrame([{**dict(kir_qc), "n_hap_seen": n_hap_seen,
                   "n_hap_cds_fasta_present": n_hap_cds_present,
                   "cds_available": kir_cds_available}]).to_csv(
        os.path.join(args.out_dir, "kir_cds_match_qc.tsv"), sep="\t", index=False)
    write_status(args.out_dir, f"DONE in {time.perf_counter()-t0:.0f}s -- 5 tables written")
    log(f"[44] done in {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
