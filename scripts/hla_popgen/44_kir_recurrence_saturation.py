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

ALLELE IDENTITY -- THREE MATCHED GRANULARITIES, SAME DEFINITION FOR BOTH SPECIES, ALL
SEQUENCE-BASED (2026-09-27 identity-impossibility fix; see build_person_kir_identity's docstring
for the full root-cause writeup and hla_id_col's for the mirrored HLA-side fix). This is the
SECOND identity fix in this file's history: the first (2026-09-26) matched KIR's genomic-level
identity to HLA's for the "all"/"any_novel" tracks but left both NAME-based; this one makes ALL
THREE granularities sequence-based for both species, closing the impossibility that let v3 report
more distinct KIR proteins (6,300) than distinct KIR genomic alleles (1,460) -- protein identity is
a coarsening of full-sequence identity and can never be richer than it.

  - `genomic` -- the best-available sequence-based identity for the FULL called allele. Neither
    species has a verified, already-used extraction path for a haplotype's full genomic sequence
    (introns/UTRs) -- only `cds.fa.gz` is confirmed extracted/joined per haplotype
    (reference/IMMUANNOT_GTF_SPEC.md; `hap{N}.trimmed.fa` is kept on disk but a gene-span
    coordinate extraction from it is new, unverified logic, out of scope here) -- so `genomic` is
    DEFINED AS THE SAME IDENTITY AS `cds` for BOTH species: a hash of the observed CDS nucleotide
    sequence, `<gene>_cds_<sha8>` for KIR, `cds_id` (24_novelty_by_field.allele_ids(), already
    hash-based for novel calls) for HLA. This is an intentional, documented equivalence (CDS as the
    best common denominator for "genomic," not a redundant duplicate column by oversight): it makes
    `S_obs(genomic) == S_obs(cds)` hold by construction, which is what
    `check_identity_invariants()` enforces at export time. Previously (both the original design and
    the 2026-09-26 partial fix) `genomic` was NAME-based -- HLA's untruncated normalized consensus
    name, KIR's raw Immuannot `consensus` string -- which is the root cause of the impossibility
    (hypothesis (a)): a NOVEL call's consensus name is just "<nearest-known-template>...new", not
    the actual novel sequence, so two genuinely different novel sequences sharing the same nearest
    template collapsed onto one "genomic allele" while the hash-based `cds`/`protein` levels
    correctly kept them apart.
  - `cds` -- the coding-sequence identity (a hash of the observed CDS nucleotide sequence). HLA:
    `cds_id` (`24_novelty_by_field.allele_ids()`, already computed -- CDS-level 3-field name for
    known/synonymous-tier calls, a `sha8` hash of the reconstructed CDS for the rest -- the 3-field
    IPD-IMGT nomenclature is bijective with CDS sequence by definition, so the name IS a legitimate
    sequence identity for known alleles, not a name-based shortcut). KIR: a hash
    (`<gene>_cds_<sha8>`) over the observed CDS sequence extracted from `<hap>/cds.fa.gz` and
    joined to the GTF call by `(contig, gene)` -- see "KIR IDENTITY EXTRACTION" below. BOTH species
    now additionally require the call to pass the SAME artifact filter before being hashed
    (hypothesis (b) -- see build_person_kir_identity's docstring fix (b)): HLA already excluded
    `partial_cds`/`inframe_stop`/`homopolymer_indel`-flagged calls and CDS sequences that don't
    translate cleanly (`frameshift`/`premature_stop`) via `keep_clean`; KIR's `cds.fa.gz`
    extraction never checked any of this before this fix, so an assembly/annotation artifact could
    become its own spurious "novel" identity. `artifact_qc.tsv` reports the per-species,
    per-gene artifact counts (task item 2).
  - `protein` -- the translated-protein identity. HLA: `prot_id` (already computed -- 2-field
    IPD-IMGT name for known alleles, a `sha8` hash of the translated protein for novel-protein
    calls). KIR: a hash (`<gene>_prot_<sha8>`) over `24_novelty_by_field.protein_info(cds_seq)
    ["protein"]` (translate + strip one terminal stop codon -- reused verbatim, not re-derived),
    from the SAME `cds.fa.gz`-extracted, artifact-filtered sequence as the `cds` level.

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

TRUE GENOMIC IDENTITY (2026-09-27b, coordinator follow-up to the same-day genomic==cds alias):
Marc's explicit ask is "general" novelty (ANY sequence difference, including non-coding) shown
side by side with protein novelty -- the genomic==cds alias lost this, since CDS excludes introns/
UTR and non-coding differences are the bulk of KIR's novelty (S03). `genomic`/`any_novel` are now
built by `build_person_true_genomic_identity()` (species-agnostic): the GTF `gene` feature row's
own (contig, start, end, strand) is used to slice `hap{N}.trimmed.fa` (kept, not gzipped, per
reference/IMMUANNOT_GTF_SPEC.md part D/E and `00_recon_vm.py`'s HAP_ROOT_FILES), reverse-
complemented on the `-` strand, then hashed. Novelty is decided by exact-or-containment comparison
against `<refdata>/gen.fa.gz` (IPD's own genomic-allele FASTA -- see `genomic_span_novel()`'s
docstring for the exact rule and why exact-length equality is too strict). The SAME
`artifact_label_of()` gate `cds`/`protein` use applies here too, PLUS a genomic-specific check: a
span touching a trimmed contig's own edge (position 1 or the contig's last base) is excluded as a
likely truncation, never silently kept. Both `hap{N}.trimmed.fa` and `gen.fa.gz` availability are
checked explicitly per species/globally; either missing degrades `genomic`/`any_novel` to the
literal `"NA"` for that species (see LEVELS_NEEDING_TRIMMED_FASTA), never a guess or a silent
fallback to the CDS-based identity. Runtime: reading one extra plain-text FASTA (typically tens to
a few hundred KB, region-limited) and doing string-slice+hash work per haplotype, PLUS the one-time
`gen.fa.gz` catalogue load (comparable in size/cost to the existing CDSseq catalogue load) and a
per-distinct-(gene,sequence) containment scan (memoized via `novelty_cache`, so repeat known
alleles cost one scan, not one per carrier) -- expected to add a similar order of magnitude to the
existing `cds.fa.gz` pass (~10-20 min at 4 workers per the original identity-extraction estimate
below), i.e. roughly +10-25 min on top of the current total, not a new dominant cost; not yet
measured on real VM data.

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
  artifact_qc.tsv          (2026-09-27, task item 2) species, gene, artifact_label, n_calls --
                           aggregate-only counts of {clean, homopolymer_indel, partial_cds,
                           inframe_stop, frameshift_or_stop} calls, SAME artifact_label_of()
                           definition for both species (see build_artifact_qc_rows()).
  diagnostics_identity.tsv (2026-09-27, task item 4) species, gene, n_calls, n_distinct_names
                           (the OLD name-based identity), n_distinct_genomic/cds/protein (the NEW
                           S_obs, from coverage_chao2.tsv's own s_obs at ancestry=ALL -- never
                           masked), n_artifact_<label> per type, and
                           median_distinct_cds_hashes_per_name -- the direct evidence for how much
                           the old name-based identity collapsed distinct sequences (hypothesis
                           (a)), per gene, for both species (see build_diagnostics_identity_rows()).
  genomic_artifact_qc.tsv  (2026-09-27b) species, gene, artifact_label, n_calls -- genomic-SPAN-
                           level artifact counts (contig_edge_truncated, partial_cds, inframe_stop,
                           homopolymer_indel, unresolved, no_catalogue_for_gene, clean), separate
                           from artifact_qc.tsv's CDS-level counts.
  genomic_identity_qc.tsv (2026-09-27b) one row: gen_catalogue_status, kir_genomic_available,
                           hla_genomic_available (bools), n_hap_seen/n_hap_trimmed_present per
                           species -- read this before trusting any 'genomic'/'any_novel' row is
                           real rather than 'NA'.
  STATUS.txt               aggregate-progress-only status file, rewritten as the run proceeds.

CHECKPOINTING / RESUME (2026-09-27c, after a ~1h VM run died at ~33min and lost everything --
Marc's rule: save gradually): every expensive stage (KIR CDS/protein identity, KIR true-genomic
identity, HLA labeled calls, HLA true-genomic identity) and the per-gene-level output-row loop
writes its result ATOMICALLY (tmp file + os.replace) to `<out-dir>/_checkpoints/<stage>.pkl` as
soon as it finishes; the two per-person identity loops ALSO checkpoint every `--checkpoint-every`
(default 1000) people within the loop itself, so a crash loses at most one chunk, not the whole
stage. On start (`--resume`, on by default), each checkpoint is validated against a small header
{script md5, args hash} written alongside it -- a checkpoint from a different code version or a
different `--kir-outroot`/`--hla-table`/`--limit`/etc. is a MISMATCH, logged and ignored (never
silently reused), and that stage/chunk is recomputed from scratch. `--no-resume` ignores every
checkpoint unconditionally (fresh run).

`_checkpoints/` holds PER-PERSON IDENTITY SETS (the exact `~/pipeline_outputs*`-derived per-person
allele sets this whole module exists to keep off the VM's shared results tree) -- it MUST NOT be
pulled off the VM, ever, disclosure rules or not: it is intermediate working state, not an
aggregate. It lives under `~/s04/results/44/_checkpoints` (i.e. under `~/s04/results/`, never
under `~/pipeline_outputs*`), one single, easily-excluded directory name -- a pull step should use
`rsync --exclude='_checkpoints'` or `find <out-dir> -maxdepth 1 -name '_checkpoints' -prune -o
-print` rather than a blanket recursive copy of `--out-dir`. Safe to delete once the run's 9
aggregate TSVs are written and pulled; a fresh `--no-resume` run also ignores (does not delete) any
stale `_checkpoints/` left over from a previous, differently-configured run.

INVARIANT (task item 3, `check_identity_invariants()`, raises ValueError -- never just warns): for
every (species, gene, ancestry) with unmasked, non-"NA" values, S_obs(protein) <= S_obs(cds) ==
S_obs(genomic), and any_novel/protein_novel n_distinct_alleles never exceed their own baseline
track's. Run once, right before the TSVs are written, on the real numbers about to be exported --
this is the hard gate that would have caught the v3 impossibility (6,300 KIR proteins > 1,460 KIR
genomic alleles) before it ever reached a report.

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
import glob
import gzip
import hashlib
import importlib.util
import json
import math
import multiprocessing as mp
import os
import pickle
import re
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _disclosure as _disc

SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
# Set once in main() from --disclosure; also propagated to m43()'s own module-level flag so
# suppressed()/rate_row() below (which delegate to m43()) stay in sync -- see set_disclosure_mode().
_DISCLOSURE_MODE = _disc.PUBLIC


def set_disclosure_mode(mode):
    global _DISCLOSURE_MODE
    _DISCLOSURE_MODE = mode
    m43()._DISCLOSURE_MODE = mode
LEVELS = ["genomic", "cds", "protein", "any_novel", "protein_novel"]
# 2026-09-27b TRUE-GENOMIC fix (supersedes the same-day genomic==cds alias): "genomic"/"any_novel"
# are back to their OWN identity, now a hash of the ACTUAL genomic span (gene start-end, including
# introns/UTR) read from hap{N}.trimmed.fa by GTF gene-row coordinates -- see
# build_person_true_genomic_identity()'s docstring. This restores non-coding novelty (the bulk of
# KIR's novelty signal, per S03) which the genomic==cds alias had collapsed away. "cds"/"protein"/
# "protein_novel" are UNCHANGED (still need cds.fa.gz); "genomic"/"any_novel" now need
# hap{N}.trimmed.fa + the genomic reference catalogue (gen.fa.gz) instead.
KIR_LEVELS_NEEDING_CDS_FASTA = {"cds", "protein", "protein_novel"}
LEVELS_NEEDING_TRIMMED_FASTA = {"genomic", "any_novel"}
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
CHECKPOINT_EVERY_DEFAULT = 1000    # people per within-loop checkpoint chunk (2026-09-27c resume fix)
CHECKPOINT_DIRNAME = "_checkpoints"  # under --out-dir; per-person state, NEVER pulled off the VM


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


def _rec_int_or_none(v):
    """Parse a recurrence_classes.tsv-style cell (int, '<20', '0', or 'NA') to an int, or None if
    it cannot be compared exactly (masked '<20' or missing 'NA' -- both hide the true value)."""
    if v is None:
        return None
    if isinstance(v, (int, np.integer)):
        return int(v)
    s = str(v)
    if s in ("NA", "<20", ""):
        return None
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return None


def check_identity_invariants(cov_rows, rec_rows):
    """HARD invariant check on `coverage_chao2.tsv`/`recurrence_classes.tsv` rows -- raises
    ValueError (never just warns, unlike `sanity_check_coverage`) if the identity hierarchy this
    whole module depends on is violated for any (species, gene, ancestry):
      1. s_obs(protein) <= s_obs(cds) <= s_obs(genomic) -- protein identity is a coarsening of CDS
         identity (translation), and CDS/genomic are now the SAME identity by construction (see
         build_person_kir_identity's/hla_id_col's 2026-09-27 fix docstrings) -- so genomic and cds
         must be EQUAL and protein must never exceed either. This is the exact invariant whose
         violation (KIR: 6,300 distinct v3 "proteins" > 1,460 distinct "genomic" alleles) triggered
         this fix.
      2. any_novel n_distinct_alleles <= genomic n_distinct_alleles, and protein_novel
         n_distinct_alleles <= protein n_distinct_alleles -- a novelty-filtered track is always a
         subset of its own baseline track.
    Only rows where BOTH sides are exact, unmasked, non-"NA" values are compared -- a masked '<20'
    or missing 'NA' cell hides the true count, so it is skipped, never treated as a violation or a
    pass. Returns nothing on success; raises ValueError with every violation found (not just the
    first) on failure."""
    cov_by_key = {}  # (species, gene, ancestry, level) -> s_obs (int) or None if uncomparable
    for r in cov_rows:
        key = (r.get("species"), r.get("gene"), r.get("ancestry"), r.get("level"))
        cov_by_key[key] = _rec_int_or_none(r.get("s_obs"))
    rec_by_key = {}  # (species, gene, ancestry, level) -> n_distinct_alleles (int) or None
    for r in rec_rows:
        key = (r.get("species"), r.get("gene"), r.get("ancestry"), r.get("level"))
        rec_by_key[key] = _rec_int_or_none(r.get("n_distinct_alleles"))

    violations = []
    triples = ({(s, g, a) for (s, g, a, _lvl) in cov_by_key}
               | {(s, g, a) for (s, g, a, _lvl) in rec_by_key})
    for species, gene, ancestry in sorted(triples, key=lambda t: (t[0] or "", t[1] or "", t[2] or "")):
        genomic = cov_by_key.get((species, gene, ancestry, "genomic"))
        cds = cov_by_key.get((species, gene, ancestry, "cds"))
        protein = cov_by_key.get((species, gene, ancestry, "protein"))
        if cds is not None and genomic is not None and cds > genomic:
            violations.append(f"{species}/{gene}/{ancestry}: s_obs(cds)={cds} > "
                               f"s_obs(genomic)={genomic}")
        if protein is not None and cds is not None and protein > cds:
            violations.append(f"{species}/{gene}/{ancestry}: s_obs(protein)={protein} > "
                               f"s_obs(cds)={cds}")

        any_novel = rec_by_key.get((species, gene, ancestry, "any_novel"))
        genomic_n = rec_by_key.get((species, gene, ancestry, "genomic"))
        if any_novel is not None and genomic_n is not None and any_novel > genomic_n:
            violations.append(f"{species}/{gene}/{ancestry}: any_novel n_distinct_alleles="
                               f"{any_novel} > genomic n_distinct_alleles={genomic_n}")
        protein_novel = rec_by_key.get((species, gene, ancestry, "protein_novel"))
        protein_n = rec_by_key.get((species, gene, ancestry, "protein"))
        if protein_novel is not None and protein_n is not None and protein_novel > protein_n:
            violations.append(f"{species}/{gene}/{ancestry}: protein_novel n_distinct_alleles="
                               f"{protein_novel} > protein n_distinct_alleles={protein_n}")

    if violations:
        raise ValueError(
            f"[44] IDENTITY INVARIANT VIOLATED ({len(violations)} case(s)) -- the identity "
            f"hierarchy (protein <= cds == genomic; novel <= baseline) does not hold on real "
            f"exported numbers. This is the exact class of bug the 2026-09-27 identity-"
            f"impossibility fix was meant to eliminate; a violation here means it is NOT fully "
            f"fixed. First few: " + "; ".join(violations[:10]))


def write_status(out_dir, msg):
    try:
        with open(os.path.join(out_dir, "STATUS.txt"), "a") as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
    except OSError as e:
        log(f"WARNING: could not write STATUS.txt: {e}")


class Timer:
    """Context manager: logs + writes to STATUS.txt a stage's start/end/elapsed. Used for every
    checkpointed stage below so a killed/resumed run's STATUS.txt shows exactly which stage was
    running when it died and how long each completed stage actually took."""
    def __init__(self, out_dir, stage_name):
        self.out_dir, self.stage_name = out_dir, stage_name

    def __enter__(self):
        self.t0 = time.perf_counter()
        write_status(self.out_dir, f"STAGE START {self.stage_name}")
        return self

    def __exit__(self, exc_type, exc, tb):
        elapsed = time.perf_counter() - self.t0
        if exc_type is None:
            write_status(self.out_dir, f"STAGE END {self.stage_name} elapsed={elapsed:.1f}s")
        else:
            write_status(self.out_dir,
                         f"STAGE FAILED {self.stage_name} elapsed={elapsed:.1f}s ({exc_type.__name__}: {exc})")
        return False


# ---------------------------------------------------------------------------
# Checkpointing (2026-09-27c): a ~1h VM run died at ~33min and lost EVERYTHING -- Marc's rule is
# "save gradually". Every expensive stage writes its result atomically (tmp file + os.replace)
# under `<out-dir>/_checkpoints/`, validated on load by a small header (script md5 + args hash) so
# a checkpoint from a stale code version or a different input configuration is ignored (logged,
# never silently reused) rather than trusted. See module docstring "CHECKPOINTING / RESUME" for
# the full contract, including why `_checkpoints/` must never be pulled off the VM.
# ---------------------------------------------------------------------------
def script_md5():
    with open(os.path.abspath(__file__), "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


# Args that affect NOTHING about the computed result (only speed, I/O location, or the
# resume/checkpoint machinery itself) -- excluded from the fingerprint so e.g. changing --workers
# or --out-dir doesn't spuriously invalidate every checkpoint.
# "disclosure" is excluded too: the export/masking layer never changes the underlying computed
# per-person identity sets that _checkpoints/ hold, only how the final tables are rendered -- so a
# checkpoint from a --disclosure public run is exactly as valid for a subsequent --disclosure
# internal run (and vice versa), letting the VM re-export both versions without recomputation
# (sprints/S04_kir_recurrence_style_share/VM_OPERATOR_HANDOFF.md's "44 INTERNAL re-export" plan).
_ARGS_EXCLUDED_FROM_FINGERPRINT = {"workers", "out_dir", "resume", "checkpoint_every", "disclosure"}


def args_fingerprint(args):
    d = {k: v for k, v in sorted(vars(args).items()) if k not in _ARGS_EXCLUDED_FROM_FINGERPRINT}
    blob = json.dumps(d, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def checkpoint_header(args):
    """The small validity header stored alongside every checkpoint -- see module docstring. One
    header per run, shared by every stage/chunk checkpoint (simpler than a per-stage subset of
    args; any input change invalidates every checkpoint, which is the safe direction to err in)."""
    return {"script_md5": script_md5(), "args_hash": args_fingerprint(args)}


def checkpoint_dir(out_dir):
    d = os.path.join(out_dir, CHECKPOINT_DIRNAME)
    os.makedirs(d, exist_ok=True)
    return d


def save_checkpoint(out_dir, name, payload, header):
    """Atomic write: a tmp file (unique per-process, so concurrent runs in different --out-dirs
    never collide) written in full, then os.replace()'d over the final path -- a reader never sees
    a partially-written checkpoint, and a crash mid-write leaves the OLD checkpoint (or none)
    intact, never a corrupt one."""
    d = checkpoint_dir(out_dir)
    path = os.path.join(d, f"{name}.pkl")
    tmp = os.path.join(d, f".{name}.tmp.{os.getpid()}")
    with open(tmp, "wb") as f:
        pickle.dump({"header": header, "payload": payload}, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)
    log(f"[44][checkpoint] saved {name} -> {path}")


def load_checkpoint(out_dir, name, header):
    """Returns the checkpointed payload if a valid (header-matching) checkpoint exists, else None
    (with a log line explaining why -- missing, unreadable, or stale). Never raises: a corrupt or
    stale checkpoint is treated exactly like a missing one -- recompute, don't crash the run."""
    path = os.path.join(checkpoint_dir(out_dir), f"{name}.pkl")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as f:
            data = pickle.load(f)
    except (OSError, EOFError, pickle.UnpicklingError, AttributeError) as e:
        log(f"[44][checkpoint] {name}: unreadable checkpoint ({e}) -- ignoring, recomputing.")
        return None
    if data.get("header") != header:
        log(f"[44][checkpoint] {name}: STALE checkpoint (script/args changed) -- ignoring, "
            f"recomputing. old={data.get('header')} new={header}")
        return None
    log(f"[44][checkpoint] {name}: valid checkpoint found -- resuming from it.")
    return data["payload"]


def _process_people_with_checkpoint(name, pids, make_arg, worker_fn, init_fn, init_args,
                                     workers, merge_fn, empty_state, out_dir=None, header=None,
                                     resume=True, checkpoint_every=CHECKPOINT_EVERY_DEFAULT):
    """Generic chunked/checkpointed per-person processing driver, shared by
    build_kir_identity_sets and build_true_genomic_identity_sets (2026-09-27c). Runs `worker_fn`
    (optionally via a `workers`-way multiprocessing.Pool, `init_fn`/`init_args` as its
    initializer) over `make_arg(pid)` for every pid in `pids` not already recorded as done in a
    resumed checkpoint, merging each result into an accumulator via `merge_fn(state, result)`
    (`result[0]` must be the pid). Every `checkpoint_every` completed people, the accumulator +
    the set of done pids is saved atomically to `<out_dir>/_checkpoints/<name>_progress.pkl` (see
    save_checkpoint) -- so a crash loses at most one chunk's worth of people, not the whole stage.

    `out_dir`/`header` of None disables checkpointing entirely (plain in-memory processing, no
    disk I/O, no resume) -- used by unit tests and any caller that doesn't want VM-only disk
    behavior; `resume=False` still writes checkpoints as it goes (so a LATER run with `--resume`
    can pick them up) but ignores any existing one on entry."""
    checkpointing = out_dir is not None and header is not None
    state, done_pids = None, set()
    if checkpointing and resume:
        ckpt = load_checkpoint(out_dir, f"{name}_progress", header)
        if ckpt is not None:
            state, done_pids = ckpt["state"], set(ckpt["done_pids"])
            log(f"[44][checkpoint] {name}: resuming with {len(done_pids)}/{len(pids)} people "
                f"already done.")
    if state is None:
        state = empty_state()

    remaining = [p for p in pids if p not in done_pids]
    if not remaining:
        log(f"[44] {name}: all {len(pids)} people already done (from checkpoint) -- nothing to "
            f"compute.")
        return state

    if len(done_pids):
        log(f"[44] {name}: processing {len(remaining)}/{len(pids)} remaining people "
            f"({len(done_pids)} resumed from checkpoint).")
    else:
        log(f"[44] {name}: processing {len(remaining)} people.")

    pool = None
    if workers <= 1:
        init_fn(*init_args)
        results_iter = (worker_fn(make_arg(p)) for p in remaining)
    else:
        pool = mp.Pool(processes=workers, initializer=init_fn, initargs=init_args)
        results_iter = pool.imap_unordered(worker_fn, (make_arg(p) for p in remaining),
                                            chunksize=32)

    since_checkpoint = 0
    try:
        for result in results_iter:
            merge_fn(state, result)
            done_pids.add(result[0])
            since_checkpoint += 1
            if checkpointing and since_checkpoint >= checkpoint_every:
                save_checkpoint(out_dir, f"{name}_progress",
                                 {"state": state, "done_pids": done_pids}, header)
                write_status(out_dir,
                             f"{name}: checkpoint at {len(done_pids)}/{len(pids)} people")
                since_checkpoint = 0
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    if checkpointing:
        save_checkpoint(out_dir, f"{name}_progress", {"state": state, "done_pids": done_pids},
                         header)
    return state


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
_KIR_TEMPLATE_WARNING_RE = re.compile(r'template_warning "([^"]*)"')
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
            warn_m = _KIR_TEMPLATE_WARNING_RE.search(attrs)
            row = {
                "gene": gene_m.group(1), "consensus": consensus, "contig": fields[0],
                "gene_id": geneid_m.group(1) if geneid_m else None,
                "copy_index": int(copy_m.group(1)) if copy_m else 1,
                "cds_distance": int(cds_dist_m.group(1)) if cds_dist_m else None,
                "cds_mut": cds_mut_m.group(1) if cds_mut_m else None,
                # 2026-09-27 identity-impossibility fix: `template_warning` is an ORDINARY
                # Immuannot transcript attribute (reference/IMMUANNOT_GTF_SPEC.md part A) written
                # by the SAME searchTemplate.py codepath for every gene family it annotates, HLA or
                # KIR (41_kir_pilot.py's own docstring: KIR calling is "the only change needed... a
                # different --region value", same tool, same attribute schema) -- it was simply
                # never extracted here before this fix, so KIR calls with a truncated/broken CDS
                # reconstruction (partial_CDS, inframe_stop) were never excluded from cds/protein
                # hashing the way 24_novelty_by_field.artifact_label_of()/keep_clean already
                # exclude the analogous HLA calls. See build_person_kir_identity's docstring.
                "template_warning": warn_m.group(1) if warn_m else None,
                "is_novel": consensus.rstrip('"').endswith("new") and consensus != "undetermined",
                "is_undetermined": consensus == "undetermined",
            }
            row["novelty_tier"] = kir.classify_novelty_tier(row)
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# TRUE GENOMIC identity (2026-09-27b, coordinator follow-up): restores non-coding novelty, which
# the same-day genomic==cds alias collapsed away -- Marc's explicit ask is "general" novelty (ANY
# sequence difference, including non-coding) side by side with protein novelty, and KIR's novelty
# is mostly non-coding (S03). Species-agnostic (the SAME functions drive both HLA and KIR, unlike
# the KIR-only cds.fa.gz-join code above) -- both species' Immuannot output shares the exact same
# GTF schema (reference/IMMUANNOT_GTF_SPEC.md part A) and the same per-person tree shape
# (<outroot>/<pid>/immuannot_output/hap{N}.gtf.gz + hap{N}.trimmed.fa, confirmed by
# 00_recon_vm.py's own HAP_ROOT_FILES convention).
# ---------------------------------------------------------------------------
_COMPLEMENT_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def load_trimmed_fasta(path):
    """hap{N}.trimmed.fa is PLAIN TEXT (NOT gzipped), per reference/IMMUANNOT_GTF_SPEC.md part D/E
    and 00_recon_vm.py's own HAP_ROOT_FILES list (unlike cds.fa.gz, which is gzipped). Returns
    {contig_id: seq} -- an empty dict, never an exception, if the file is missing; callers gate on
    this explicitly (this task's "never a silent fallback" requirement)."""
    seqs = {}
    if not path or not os.path.exists(path):
        return seqs
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def parse_hap_gtf_transcript_minimal(gtf_path):
    """Species-agnostic parse of Immuannot 'transcript' feature rows, for TRUE-GENOMIC identity
    only: gene, gene_id, contig, consensus, cds_mut, template_warning. Deliberately does NOT call
    41_kir_pilot.classify_novelty_tier (KIR-naming-convention-specific, e.g. its no-colon 'new'
    suffix rule) -- true-genomic novelty here is decided by sequence containment against
    gen.fa.gz (genomic_span_novel below), not by the consensus 'new'-suffix heuristic, so no
    species-specific novelty classification is needed at this layer for either species."""
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
            geneid_m = _KIR_GENE_ID_RE.search(attrs)
            if not gene_m or not cons_m or not geneid_m:
                continue
            cds_mut_m = _KIR_CDS_MUT_RE.search(attrs)
            warn_m = _KIR_TEMPLATE_WARNING_RE.search(attrs)
            rows.append({
                "gene": gene_m.group(1), "gene_id": geneid_m.group(1), "contig": fields[0],
                "consensus": cons_m.group(1),
                "cds_mut": cds_mut_m.group(1) if cds_mut_m else None,
                "template_warning": warn_m.group(1) if warn_m else None,
            })
    return rows


def parse_hap_gtf_gene_rows(gtf_path):
    """Species-agnostic parse of Immuannot 'gene' feature rows (reference/IMMUANNOT_GTF_SPEC.md
    part A): gene_id -> {contig, start, end, strand}. Columns 4/5 are 1-based inclusive,
    contig-relative (part C) -- directly usable as a trimmed.fa slice, no off-by-one adjustment
    (part D's "indirect route", step 2)."""
    spans = {}
    if not gtf_path or not os.path.exists(gtf_path):
        return spans
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gid_m = _KIR_GENE_ID_RE.search(fields[8])
            if not gid_m:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            spans[gid_m.group(1)] = {"contig": fields[0], "start": start, "end": end,
                                     "strand": fields[6]}
    return spans


def extract_genomic_span(trimmed_seqs, contig, start, end, strand):
    """1-based inclusive GTF gene-row coords -> the observed genomic span sequence (introns/UTR
    included), reverse-complemented to coding-sense orientation on the '-' strand (reference/
    IMMUANNOT_GTF_SPEC.md part D, "indirect route" step 3). Returns (seq_or_None,
    touches_contig_edge: bool). `touches_contig_edge` is the coordinator's flagged artifact: a gene
    mapped right up to a trimmed contig's own boundary (position 1, or the contig's last base) is a
    truncation candidate (the assembly/trim may have cut off real flanking sequence), not a
    complete genomic call -- callers exclude these from genomic identity entirely, same treatment
    as any other artifact class."""
    seq_full = trimmed_seqs.get(contig)
    if not seq_full:
        return None, False
    n = len(seq_full)
    if start < 1 or end > n or start > end:
        return None, False
    touches_edge = (start == 1) or (end == n)
    span = seq_full[start - 1:end]
    if strand == "-":
        span = reverse_complement(span)
    return span, touches_edge


def load_genomic_catalogue(refdata_dir, m24mod):
    """Loads the IPD genomic-allele reference FASTA into {gene_bare: set(seq)} for the
    true-genomic novelty check. Path: `<refdata>/gen.fa.gz` (the coordinator's specified location);
    falls back to a recursive glob (`**/gen.fa.gz`) in case it's nested like `CDSseq/*.fa.gz` is.
    Reuses `24_novelty_by_field.iter_fasta`/`parse_ref_header` verbatim -- the SAME header-parsing
    convention already used for CDSseq, so this covers HLA and KIR from one file/one parse if IPD
    bundles both gene families together (as CDSseq does). GUARDED: returns (None, reason) -- never
    raises, never fabricates an empty-but-'ok' catalogue -- if the file can't be found or yields no
    records. Callers must export genomic-level NOVELTY as 'NA' with that reason when this is None
    (the genomic BASELINE hash doesn't need this catalogue and is unaffected)."""
    path = os.path.join(refdata_dir, "gen.fa.gz")
    if not os.path.exists(path):
        cands = sorted(glob.glob(os.path.join(refdata_dir, "**", "gen.fa.gz"), recursive=True))
        path = cands[0] if cands else None
    if not path:
        return None, "no_gen_fasta"
    catalogue = defaultdict(set)
    for g, name, seq, frame in m24mod.iter_fasta(path):
        if g and seq:
            catalogue[_bare_gene(g)].add(seq.upper())
    if not catalogue:
        return None, "gen_fasta_empty"
    return dict(catalogue), "ok"


def _bare_gene(g):
    """'HLA-A' -> 'A'; KIR names and bare names pass through unchanged."""
    return g[4:] if isinstance(g, str) and g.startswith("HLA-") else g


def genomic_span_novel(obs_seq, ref_seqs):
    """Novelty rule for TRUE-GENOMIC identity (defensible choice, documented per the coordinator's
    ask -- "choose a defensible rule, document it, and test it"): IPD's own genomic reference
    records do not necessarily span the identical coordinates our own gene-row extraction does (UTR
    extent conventions differ between an assembly-derived trim and IPD's own record boundary), so
    requiring an EXACT match would misclassify almost every textbook-known allele as "novel" purely
    over a UTR-length mismatch. The precise fix (trimming to the shared aligned extent via
    `mm2.ipd.gen.paf.gz`, when present) needs per-haplotype PAF re-parsing not implemented in this
    pass -- flagged as a follow-up, not silently approximated as something stronger than it is.
    RULE USED: NOT novel iff the observed span EQUALS a reference genomic sequence for this gene,
    OR is CONTAINED IN one (our trim is the tighter of the two), OR CONTAINS one (the reference
    record is the tighter one) -- i.e. "novel" means no reference genomic sequence for this gene
    shares a containment relationship with the observed span in EITHER direction. Returns None
    (never True/False) when `ref_seqs` is empty/unavailable for this gene -- novelty is
    UNDETERMINED then, not silently "not novel"."""
    if not ref_seqs:
        return None
    for ref in ref_seqs:
        if obs_seq == ref or obs_seq in ref or ref in obs_seq:
            return False
    return True


def build_person_true_genomic_identity(pid, hap_root, gene_names, gen_catalogue, m24mod,
                                        novelty_cache=None):
    """Species-agnostic (SAME function drives HLA and KIR -- 2026-09-27b true-genomic fix).
    `hap_root`: the <outroot>/<pid>/immuannot_output tree's PARENT (KIR: --kir-outroot; HLA:
    --hla-people-outroot, the SAME per-person convention 39/03 already use for cds.fa.gz).
    `gene_names`: the species' own bare gene set (KIR_GENES / CLASSICAL_GENES_BARE).
    `gen_catalogue`: {gene: set(seq)} from load_genomic_catalogue(), or None if unavailable.
    `novelty_cache`: optional shared dict {(gene, seq): bool_or_None} -- many calls across people
    share the SAME textbook-known genomic sequence, so caching genomic_span_novel()'s O(catalogue
    size) containment scan per DISTINCT (gene, seq) pair (not per call) is the difference between
    "cheap" and "re-scanning a multi-thousand-allele catalogue per haplotype."

    Returns (pid, per_gene: {gene: {'genomic': set(id), 'any_novel': set(id)}}, qc: Counter,
    n_hap_seen: int, n_hap_trimmed_present: int).

    GATE per call (task: "apply the same artifact gate to genomic spans that touch a contig end
    (truncated)" plus the pre-existing artifact_label_of() gate cds/protein already use):
      1. `consensus == "undetermined"` -> skipped (shared HLA/KIR semantics, callIPDallele.py:140 --
         Immuannot could not resolve this call at all).
      2. `artifact_label_of(template_warning, cds_mut) != "clean"` (partial_cds/inframe_stop/
         homopolymer_indel, from the SAME attributes cds/protein-level processing already reads) ->
         excluded, `qc[("artifact_genomic", gene, label)]`.
      3. the gene-row span can't be resolved (no matching gene_id, or hap{N}.trimmed.fa missing/
         doesn't contain the contig) -> excluded, `qc[("artifact_genomic", gene, "unresolved")]`.
      4. the extracted span TOUCHES A CONTIG END (`extract_genomic_span`'s `touches_contig_edge`) ->
         excluded as a truncation artifact, `qc[("artifact_genomic", gene, "contig_edge_truncated")]`.
    A call surviving all four gets `genomic_id = f"{gene}_gen_{sha8(span)}"`, added to 'genomic'.
    Novelty (added to 'any_novel') is `genomic_span_novel(span, gen_catalogue.get(gene))` -- if that
    is None (no catalogue entries for this gene), the call still counts toward 'genomic' baseline
    richness but NEVER toward 'any_novel' either way (novelty is genuinely undetermined, not
    guessed), and `qc[("artifact_genomic", gene, "no_catalogue_for_gene")]` records it."""
    per_gene = defaultdict(lambda: {"genomic": set(), "any_novel": set()})
    qc = Counter()
    n_hap_seen = n_hap_trimmed_present = 0
    if novelty_cache is None:
        novelty_cache = {}
    for hap in ("hap1", "hap2"):
        person_dir = os.path.join(hap_root, str(pid), "immuannot_output")
        gtf_path = os.path.join(person_dir, f"{hap}.gtf.gz")
        # HLA GTFs name genes "HLA-A" while CLASSICAL_GENES_BARE is "A": compare bare names
        # (bug found on the VM smoke test: every HLA row was dropped -> genomic level NA).
        wanted = {_bare_gene(g) for g in gene_names}
        rows = []
        for r in parse_hap_gtf_transcript_minimal(gtf_path):
            if _bare_gene(r["gene"]) in wanted:
                r = dict(r, gene=_bare_gene(r["gene"]))
                rows.append(r)
        if not rows:
            continue
        n_hap_seen += 1
        trimmed_seqs = load_trimmed_fasta(os.path.join(person_dir, f"{hap}.trimmed.fa"))
        if trimmed_seqs:
            n_hap_trimmed_present += 1
        gene_spans = parse_hap_gtf_gene_rows(gtf_path)
        for r in rows:
            gene, consensus = r["gene"], r["consensus"]
            if consensus == "undetermined":
                continue
            artifact_label = m24mod.artifact_label_of(r.get("template_warning"), r["cds_mut"])
            if artifact_label != "clean":
                qc[("artifact_genomic", gene, artifact_label)] += 1
                continue
            span_info = gene_spans.get(r["gene_id"]) if r["gene_id"] else None
            seq = touches_edge = None
            if span_info and trimmed_seqs:
                seq, touches_edge = extract_genomic_span(
                    trimmed_seqs, span_info["contig"], span_info["start"], span_info["end"],
                    span_info["strand"])
            if not seq:
                qc[("artifact_genomic", gene, "unresolved")] += 1
                continue
            if touches_edge:
                qc[("artifact_genomic", gene, "contig_edge_truncated")] += 1
                continue
            qc[("artifact_genomic", gene, "clean")] += 1
            genomic_id = f"{gene}_gen_{m24mod.sha8(seq)}"
            per_gene[gene]["genomic"].add(genomic_id)
            cache_key = (gene, seq)
            if cache_key in novelty_cache:
                novel = novelty_cache[cache_key]
            else:
                ref_seqs = gen_catalogue.get(gene) if gen_catalogue else None
                novel = genomic_span_novel(seq, ref_seqs)
                novelty_cache[cache_key] = novel
            if novel is None:
                qc[("artifact_genomic", gene, "no_catalogue_for_gene")] += 1
            elif novel:
                per_gene[gene]["any_novel"].add(genomic_id)
    per_gene = {g: {"genomic": d["genomic"], "any_novel": d["any_novel"]}
                for g, d in per_gene.items()}
    return pid, per_gene, qc, n_hap_seen, n_hap_trimmed_present


_id_gen_catalogue = None
_id_gen_cache = None


def _genomic_worker_init(gen_catalogue):
    global _id_gen_catalogue, _id_gen_cache
    _id_gen_catalogue = gen_catalogue
    _id_gen_cache = {}


def _genomic_worker_task(args):
    pid, hap_root, gene_names = args
    return build_person_true_genomic_identity(pid, hap_root, gene_names, _id_gen_catalogue, m24(),
                                              _id_gen_cache)


def _true_genomic_empty_state():
    return {"sets_by_level": {"genomic": defaultdict(dict), "any_novel": defaultdict(dict)},
            "qc": Counter(), "n_hap_seen": 0, "n_hap_trimmed_present": 0}


def _true_genomic_merge(state, result):
    pid, per_gene, qc, n_seen, n_trim = result
    for gene, levels in per_gene.items():
        for lvl in ("genomic", "any_novel"):
            state["sets_by_level"][lvl][gene][pid] = levels[lvl]
    state["qc"].update(qc)
    state["n_hap_seen"] += n_seen
    state["n_hap_trimmed_present"] += n_trim


def build_true_genomic_identity_sets(pids, hap_root, gene_names, gen_catalogue, workers,
                                      out_dir=None, header=None, resume=True,
                                      checkpoint_every=CHECKPOINT_EVERY_DEFAULT,
                                      checkpoint_name="true_genomic"):
    """Orchestrates build_person_true_genomic_identity over all `pids` -- species-agnostic; the
    SAME function drives KIR's and HLA's own per-person trees (see main()). Returns
    (sets_by_level: {'genomic': {gene: {pid: set}}, 'any_novel': {gene: {pid: set}}},
    qc: Counter, trimmed_available: bool, n_hap_seen: int, n_hap_trimmed_present: int).
    `trimmed_available` is False iff hap{N}.trimmed.fa was found for ZERO haplotypes across the
    whole cohort -- main() exports 'genomic'/'any_novel' as 'NA' (reason='no_trimmed_fasta') in
    that case for this species, never a silent fallback.

    2026-09-27c: this per-person loop is one of the checkpointed stages (module docstring
    "CHECKPOINTING / RESUME") -- `out_dir`/`header` given (main() passes both, with a species-
    specific `checkpoint_name` so KIR's and HLA's own true-genomic passes don't collide) turns on
    resumable, chunked (`checkpoint_every` people) checkpointing via
    `_process_people_with_checkpoint`; left as None (the default, used by every existing caller/
    test) reproduces the exact prior in-memory-only behavior."""
    # _genomic_worker_init/_genomic_worker_task (module-level globals) are used for BOTH the
    # workers<=1 and workers>1 paths so a single shared novelty_cache persists across the whole
    # cohort either way -- matching the pre-checkpointing behavior exactly (previously the
    # workers<=1 branch built one local `cache` dict and reused it across all `pids`).
    state = _process_people_with_checkpoint(
        checkpoint_name, pids,
        make_arg=lambda pid: (pid, hap_root, gene_names),
        worker_fn=_genomic_worker_task,
        init_fn=_genomic_worker_init, init_args=(gen_catalogue,),
        workers=workers, merge_fn=_true_genomic_merge, empty_state=_true_genomic_empty_state,
        out_dir=out_dir, header=header, resume=resume, checkpoint_every=checkpoint_every)
    trimmed_available = state["n_hap_trimmed_present"] > 0
    return (state["sets_by_level"], state["qc"], trimmed_available, state["n_hap_seen"],
            state["n_hap_trimmed_present"])


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

    2026-09-27 IDENTITY-IMPOSSIBILITY FIX (S04 coordinator item -- see LOG.md "Orchestrator
    sanity check REJECTS the KIR protein numbers"): the v3 pipeline could report MORE distinct
    KIR proteins (6,300) than distinct KIR "genomic" alleles (1,460), which is a logical
    impossibility -- protein identity is a coarsening of the full sequence, so it can never be
    finer. Root-caused to TWO independent bugs, both fixed here:

    (a) NAME-based, not sequence-based, "genomic"/"any_novel" identity. `genomic` used to be the
        raw Immuannot `consensus` string (and HLA's mirrored it with the untruncated normalized
        name -- see `hla_id_col` below). For a NOVEL call, Immuannot's consensus is just
        "<nearest-known-prefix>...new" -- it names WHICH known template the call is closest to
        and THAT a novel difference exists, not WHAT the actual novel sequence is. Two people
        with two genuinely DIFFERENT novel sequences that happen to diverge from the same known
        template at the same depth get the IDENTICAL consensus string -- collapsing distinct
        sequences onto one "genomic allele" while `cds`/`protein` (already hash-based for novel
        calls) correctly kept them apart. This is exactly backwards for a genomic-vs-protein
        richness comparison. FIX: `genomic` (and `any_novel`, its novelty-filtered counterpart)
        are now the SAME sequence-hash identity as `cds` (`<gene>_cds_<sha8>`, a hash of the
        observed CDS nucleotide sequence). This project has no code path that extracts a
        haplotype's FULL genomic sequence (introns/UTRs) for either species -- only `cds.fa.gz`
        is confirmed extracted/joined (reference/IMMUANNOT_GTF_SPEC.md; `hap{N}.trimmed.fa` is
        kept on disk but extracting a gene's genomic span from it needs new, VM-unverified
        coordinate logic) -- so CDS is used as the best common denominator for "genomic" on BOTH
        species (HLA's mirrored fix is in `hla_id_col`/`add_hla_genomic_id` below). This means
        `genomic` and `cds` are IDENTICAL by construction now (documented, not a redundancy bug):
        `S_obs(genomic) == S_obs(cds)` always holds, satisfying the invariant
        `S_obs(protein) <= S_obs(cds) <= S_obs(genomic)` (`check_identity_invariants()` below
        enforces this at export time, raising rather than warning).
    (b) NO artifact filtering for KIR before hashing. HLA's own pipeline
        (`24_novelty_by_field.artifact_label_of`/`keep_clean`) excludes `partial_cds`,
        `inframe_stop` (from Immuannot's `template_warning` attribute) and `homopolymer_indel`
        (from `cds_mut`), plus any CDS whose length isn't a multiple of 3 or that translates with
        an internal stop (`seq_class == "frameshift_or_stop"`), before it ever hashes a CDS/
        protein. KIR's `cds.fa.gz`/`cds_mut` extraction never read `template_warning` at all and
        never checked for frameshift/premature-stop after translating -- every assembly/annotation
        artifact (a truncated CDS reconstruction, a broken reading frame, an internal stop) became
        its own unique "novel protein" hash, inflating KIR's protein-level `S_obs`/novelty/
        singleton-share numbers relative to HLA's already-filtered ones. FIX: `parse_hap_gtf_full`
        now also extracts `template_warning`; this function computes the SAME
        `m24mod.artifact_label_of(template_warning, cds_mut)` HLA uses and excludes any non-clean
        call, PLUS excludes any resolved CDS whose `protein_info()` reports `frameshift` or
        `premature_stop` -- the identical two-part gate as HLA's `keep_clean`
        (`artifact_label == "clean" and seq_class != "frameshift_or_stop"`), applied symmetrically.
        Per-gene artifact counts are accumulated into `qc` (tuple keys `("artifact", gene, label)`)
        for the new `artifact_qc.tsv` aggregate (task item 2) and diagnostic tuple keys
        `("namehash", gene, consensus, cds_hash_or_"unresolved")` / `("n_calls", gene)` for the new
        `diagnostics_identity.tsv` (task item 4) -- see `main()`'s aggregation of these.

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
            qc[("n_calls", gene)] += 1
            if tier == "undetermined":
                continue
            # Same artifact gate HLA's keep_clean already applies (fix (b) above) -- computed
            # regardless of whether the cds.fa.gz join succeeds, since template_warning/cds_mut
            # come straight off the GTF row.
            artifact_label = m24mod.artifact_label_of(r.get("template_warning"), r["cds_mut"])
            qc[("artifact", gene, artifact_label)] += 1
            if artifact_label != "clean":
                continue
            seq = seq_by_key.get((gene, consensus))
            if not seq:
                qc[("namehash", gene, consensus, "unresolved")] += 1
                continue
            pinfo = m24mod.protein_info(seq)
            if pinfo["frameshift"] or pinfo["premature_stop"]:
                qc[("artifact", gene, "frameshift_or_stop")] += 1
                qc[("namehash", gene, consensus, "unresolved")] += 1
                continue
            # NOTE 2026-09-27b: this function no longer populates 'genomic'/'any_novel' -- those
            # are now built separately by build_person_true_genomic_identity() from
            # hap{N}.trimmed.fa gene spans (restores non-coding novelty; see LEVELS_NEEDING_
            # TRIMMED_FASTA / that function's docstring). 'cds'/'protein'/'protein_novel' below are
            # unchanged.
            cds_id = f"{gene}_cds_{m24mod.sha8(seq)}"
            per_level["cds"][gene].add(cds_id)
            qc[("namehash", gene, consensus, cds_id)] += 1
            if not kir_gene_protein_covered(kir_ref, gene):
                qc["n_protein_gene_uncovered_calls"] += 1
                continue  # never guess a finer-than-genomic id for an uncovered gene
            protein = pinfo["protein"]
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


def _kir_identity_empty_state():
    return {"sets_by_level": {lvl: defaultdict(dict) for lvl in LEVELS}, "qc": Counter(),
            "n_hap_seen": 0, "n_hap_cds_present": 0}


def _kir_identity_merge(state, result):
    pid, per_level, qc, n_seen, n_cds = result
    for lvl, gene_sets in per_level.items():
        for gene, ids in gene_sets.items():
            state["sets_by_level"][lvl][gene][pid] = ids
    state["qc"].update(qc)
    state["n_hap_seen"] += n_seen
    state["n_hap_cds_present"] += n_cds


def build_kir_identity_sets(pids, kir_outroot, workers, refdata=None, out_dir=None, header=None,
                             resume=True, checkpoint_every=CHECKPOINT_EVERY_DEFAULT,
                             checkpoint_name="kir_identity"):
    """Orchestrates build_person_kir_identity over all `pids` (optionally multiprocessed, mirroring
    43_kir_full_aggregate.parse_all's Pool(initializer=...) pattern). Returns
    (sets_by_level: {level: {gene: {pid: set(id)}}}, qc: Counter, cds_available: bool,
    n_hap_seen: int, n_hap_cds_present: int). `cds_available` is False iff cds.fa.gz was found for
    ZERO haplotypes across the entire cohort (a hard failure mode -- e.g. the file was cleaned up
    or never written for the KIR run -- not ordinary per-person sparsity).

    `refdata`: passed through to each worker's own `_identity_worker_init` (not shared as a live
    object across process boundaries -- each worker loads its own KIR RefIndex, same pattern this
    file already uses for m41/m03/m24). None (the default) reproduces the pre-fix hash-only /
    novelty_tier-heuristic behavior -- see build_person_kir_identity's docstring.

    2026-09-27c: checkpointed/resumable exactly like build_true_genomic_identity_sets -- see that
    function's docstring and the module docstring "CHECKPOINTING / RESUME". `out_dir`/`header` of
    None (the default, used by every existing caller/test) disables checkpointing entirely."""
    state = _process_people_with_checkpoint(
        checkpoint_name, pids,
        make_arg=lambda pid: (pid, kir_outroot),
        worker_fn=_identity_worker_task,
        init_fn=_identity_worker_init, init_args=(refdata,),
        workers=workers, merge_fn=_kir_identity_merge, empty_state=_kir_identity_empty_state,
        out_dir=out_dir, header=header, resume=resume, checkpoint_every=checkpoint_every)
    cds_available = state["n_hap_cds_present"] > 0
    return (state["sets_by_level"], state["qc"], cds_available, state["n_hap_seen"],
            state["n_hap_cds_present"])


# ---------------------------------------------------------------------------
# HLA-side extraction: reuses 39_saturation_by_ancestry.build_labeled_calls() (the full
# field_class/seq_class/prot_id/cds_id pipeline) -- no re-derivation of sequence matching.
#
# 2026-09-27 identity-impossibility fix (mirrors the KIR-side fix in build_person_kir_identity's
# docstring, fix (a)): "genomic_id" (the full, untruncated normalized consensus name) is KEPT here
# as a diagnostic-only column -- it is exactly the NAME-based identity that collapses distinct
# novel sequences sharing a known prefix+depth onto one string (a novel HLA consensus is also just
# "<known-prefix>:new", not the actual novel sequence) -- but it is NO LONGER used as the
# `genomic`/`any_novel` IDENTITY (see hla_id_col below, now aliased to `cds_id`). `diagnostics_
# identity.tsv`'s "median distinct CDS hashes per reference name" statistic groups by this column
# specifically to quantify how much collapse the OLD name-based identity caused, for both species.
# ---------------------------------------------------------------------------
def add_hla_genomic_id(calls, m24mod):
    """Adds a 'genomic_id' column: the full normalized consensus name (all colon fields kept),
    cached per distinct consensus string (there are far fewer distinct strings than rows).
    DIAGNOSTIC-ONLY as of the 2026-09-27 fix -- see module docstring above; not used for
    `genomic`/`any_novel` identity any more (hla_id_col aliases those to `cds_id`)."""
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
    """2026-09-27 fix: 'genomic'/'any_novel' now alias 'cds_id' (a real sequence-hash identity for
    novel calls, a nomenclature-bijective 3-field name for known/cds_known/f4_noncoding calls --
    see 24_novelty_by_field.allele_ids()) instead of the name-based 'genomic_id', which collapsed
    distinct novel sequences sharing a known prefix+depth (see build_person_kir_identity's
    docstring fix (a) for the full rationale, mirrored here). This makes genomic == cds for HLA by
    construction too, matching KIR's own fix and satisfying check_identity_invariants()."""
    return {"genomic": "cds_id", "cds": "cds_id", "protein": "prot_id",
            "any_novel": "cds_id", "protein_novel": "prot_id"}[level]


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
    """n_distinct_alleles: count_type=ALLELE_DISTINCT (richness) -- exact in both disclosure
    modes. eq1/eq2/gt2/ge20: count_type=RECURRENCE_CLASS -- exact in both modes too (the S04
    two-disclosure-versions rule loosens recurrence-class counts even under PUBLIC), but PUBLIC
    additionally carries a `review_flag` column (True when this gene/ancestry's total number of
    unrelated carriers is below _disclosure.REVIEW_TOTAL_CARRIERS_FLOOR, i.e. even a handful of
    eq1/eq2 alleles could plausibly account for most of a thin stratum's carriers -- see
    AOU_SMALL_CELL_POLICY.md answer (2)) for a supervisor to check before treating the cell as
    publishable as-is. INTERNAL mode never sets review_flag (not a PUBLIC concept) but does carry
    per-class lt20 flags instead."""
    counts = carrier_counts(unit_sets)
    cls = classify_recurrence(counts)
    total_carriers = sum(1 for s in unit_sets.values() if len(s) > 0)
    row = {"gene": gene, "ancestry": ancestry, "level": level, "species": species,
           "n_people": len(unit_sets),
           "n_distinct_alleles": _disc.mask(len(counts), _disc.ALLELE_DISTINCT,
                                             _DISCLOSURE_MODE).display}
    rec_results = {}
    for key in ("eq1", "eq2", "gt2", "ge20"):
        rec_results[key] = _disc.mask(cls[key], _disc.RECURRENCE_CLASS, _DISCLOSURE_MODE,
                                       total_carriers=total_carriers)
        row[key] = rec_results[key].display
    if _DISCLOSURE_MODE == _disc.PUBLIC:
        row["review_flag"] = any(r.review_flag for r in rec_results.values())
    else:
        for key, r in rec_results.items():
            row[f"{key}_lt20"] = r.lt20
    return row


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
    # s_obs above is richness (ALLELE_DISTINCT/RICHNESS) -- already exported exact unconditionally,
    # correct in both disclosure modes. q1/q2 (raw Good-Turing f1/f2) are a stricter, PROJECT-
    # SPECIFIC caveat on top of the general policy (module docstring "DISCLOSURE": "export only
    # the estimate if f1/f2 are small" -- these raw incidence counts are closer to recurrence
    # structure than a plain richness tally, so PUBLIC keeps the existing blank-both-together rule
    # regardless of _DISCLOSURE_MODE's general looseness elsewhere). INTERNAL mode shows them
    # exact (with lt20 flags), matching every other cell.
    q1, q2 = Q.get(1, 0), Q.get(2, 0)
    if _DISCLOSURE_MODE == _disc.INTERNAL:
        row["q1"], row["q2"] = str(q1), str(q2)
        row["q1_lt20"] = 0 < q1 < SUPPRESS_BELOW
        row["q2_lt20"] = 0 < q2 < SUPPRESS_BELOW
    else:
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


# ---------------------------------------------------------------------------
# Diagnostic aggregates (task item 4, 2026-09-27 identity-impossibility fix): artifact_qc.tsv and
# diagnostics_identity.tsv, built symmetrically for both species from data already collected above
# (kir_qc's tuple-keyed diagnostic entries; HLA's own `calls` frame, which already carries
# artifact_label/genomic_id/cds_id/prot_id per call).
# ---------------------------------------------------------------------------
def split_kir_qc(kir_qc):
    """Separates kir_qc's SCALAR (string-keyed, aggregate) entries -- the pre-existing
    kir_cds_match_qc.tsv fields (n_matched, n_ambiguous_copy, etc.) -- from its new tuple-keyed
    diagnostic entries added by build_person_kir_identity's 2026-09-27 fix:
      ("n_calls", gene) -> count of every KIR transcript row seen for that gene (pre-filter).
      ("artifact", gene, label) -> count of calls with that artifact_label (label in {"clean",
        "homopolymer_indel", "partial_cds", "inframe_stop", "frameshift_or_stop"}).
      ("namehash", gene, consensus_name, cds_hash_or_"unresolved") -> count of calls where this
        (pre-fix) consensus NAME mapped to this (post-fix) CDS HASH -- the direct evidence for
        name-collapse magnitude (task item 1/hypothesis (a)).
    Returns (scalar_qc: dict, n_calls_by_gene: dict, artifact_by_gene_label: dict,
    name_to_hashes_by_gene: {gene: {name: set(hash)}})."""
    scalar_qc = {}
    n_calls_by_gene = Counter()
    artifact_by_gene_label = Counter()
    name_to_hashes_by_gene = defaultdict(lambda: defaultdict(set))
    for k, v in kir_qc.items():
        if isinstance(k, tuple):
            if k[0] == "n_calls":
                n_calls_by_gene[k[1]] += v
            elif k[0] == "artifact":
                artifact_by_gene_label[(k[1], k[2])] += v
            elif k[0] == "namehash":
                _, gene, name, h = k
                if h != "unresolved":
                    name_to_hashes_by_gene[gene][name].add(h)
        else:
            scalar_qc[k] = v
    return scalar_qc, dict(n_calls_by_gene), dict(artifact_by_gene_label), dict(name_to_hashes_by_gene)


def build_artifact_qc_rows(kir_artifact_by_gene_label, hla_calls):
    """artifact_qc.tsv: gene, species, artifact_label, n_calls -- aggregate-only, symmetric across
    species (task item 2). HLA counts are read straight from `calls['artifact_label']`
    (24_novelty_by_field.add_call_labels, already computed for every HLA call by
    39.build_labeled_calls); KIR counts come from the SAME artifact_label_of() function, newly
    applied per build_person_kir_identity's 2026-09-27 fix."""
    rows = []
    for (gene, label), n in sorted(kir_artifact_by_gene_label.items()):
        rows.append({"species": "kir", "gene": gene, "artifact_label": label, "n_calls": int(n)})
    if hla_calls is not None and len(hla_calls):
        tab = hla_calls.groupby(["gene_b", "artifact_label"]).size()
        for (gene, label), n in tab.items():
            rows.append({"species": "hla", "gene": gene, "artifact_label": label, "n_calls": int(n)})
    return rows


def _median(values):
    if not values:
        return float("nan")
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return float(s[mid]) if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def build_diagnostics_identity_rows(all_cov, kir_n_calls_by_gene, kir_artifact_by_gene_label,
                                     kir_name_to_hashes, hla_calls):
    """diagnostics_identity.tsv, one row per (species, gene) -- pooled ALL ancestry (task item 4):
      n_calls, n_distinct_names (the OLD, pre-fix name-based identity), n_distinct_genomic/cds/
      protein (the NEW, post-fix sequence-hash S_obs, read from coverage_chao2.tsv's own `s_obs`
      at ancestry=ALL -- s_obs is never disclosure-masked, unlike n_distinct_alleles in
      recurrence_classes.tsv, so this is always a real number here), n_artifact_<label> per
      artifact type, and median_distinct_cds_hashes_per_name -- directly answers "does one NAME
      collapse multiple distinct sequences" (hypothesis (a)), for both species side by side."""
    s_obs_by = {}  # (species, gene, level) -> s_obs at ancestry=ALL
    for r in all_cov:
        if r.get("ancestry") != "ALL":
            continue
        s_obs_by[(r.get("species"), r.get("gene"), r.get("level"))] = _rec_int_or_none(r.get("s_obs"))

    artifact_labels = sorted({label for (_, label) in kir_artifact_by_gene_label})
    if hla_calls is not None and "artifact_label" in hla_calls.columns:
        artifact_labels = sorted(set(artifact_labels) | set(hla_calls["artifact_label"].dropna().unique()))

    rows = []

    def add_row(species, gene, n_calls, name_to_hashes, artifact_counts):
        row = {
            "species": species, "gene": gene, "n_calls": n_calls,
            "n_distinct_names": len(name_to_hashes),
            "n_distinct_genomic": s_obs_by.get((species, gene, "genomic"), ""),
            "n_distinct_cds": s_obs_by.get((species, gene, "cds"), ""),
            "n_distinct_protein": s_obs_by.get((species, gene, "protein"), ""),
            "median_distinct_cds_hashes_per_name": round(
                _median([len(hs) for hs in name_to_hashes.values()]), 2),
        }
        for label in artifact_labels:
            row[f"n_artifact_{label}"] = int(artifact_counts.get(label, 0))
        rows.append(row)

    for gene in sorted(kir_n_calls_by_gene):
        name_to_hashes = kir_name_to_hashes.get(gene, {})
        artifact_counts = {label: n for (g, label), n in kir_artifact_by_gene_label.items()
                            if g == gene}
        add_row("kir", gene, int(kir_n_calls_by_gene[gene]), name_to_hashes, artifact_counts)

    if hla_calls is not None and len(hla_calls):
        for gene, gdf in hla_calls.groupby("gene_b"):
            name_to_hashes = defaultdict(set)
            for name, cds_id in zip(gdf["genomic_id"], gdf["cds_id"]):
                if isinstance(name, str) and isinstance(cds_id, str):
                    name_to_hashes[name].add(cds_id)
            artifact_counts = gdf["artifact_label"].value_counts().to_dict()
            add_row("hla", gene, int(len(gdf)), dict(name_to_hashes), artifact_counts)

    return rows


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
    ap.add_argument("--resume", dest="resume", action="store_true", default=True,
                     help="(default) resume from any valid checkpoint under "
                          "<out-dir>/_checkpoints/ instead of recomputing a finished stage/chunk; "
                          "a checkpoint whose script-md5+args-hash header doesn't match this run "
                          "is ignored (logged), never silently reused.")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                     help="ignore any existing checkpoint and recompute everything from scratch "
                          "(still WRITES fresh checkpoints as it goes, unless --out-dir is shared "
                          "with a run you don't want to disturb).")
    ap.add_argument("--checkpoint-every", type=int, default=CHECKPOINT_EVERY_DEFAULT,
                     help="checkpoint the per-person identity loops every N people (default "
                          f"{CHECKPOINT_EVERY_DEFAULT}) -- a crash loses at most one chunk.")
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args()

    if args.disclosure == _disc.INTERNAL and args.out_dir == ap.get_default("out_dir"):
        args.out_dir = _disc.default_out_dir(_disc.INTERNAL, "44", args.out_dir)
    if args.disclosure == _disc.INTERNAL:
        _disc.assert_internal_path_allowed(args.out_dir)
    set_disclosure_mode(args.disclosure)

    os.makedirs(args.out_dir, exist_ok=True)
    header = checkpoint_header(args)
    write_status(args.out_dir, f"START outroot={args.kir_outroot} out_dir={args.out_dir} "
                                f"resume={args.resume} checkpoint_header={header}")
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
    # 2026-09-27c: whole-stage checkpoint (not chunked internally -- build_labeled_calls() is
    # 39_saturation_by_ancestry's own pipeline, reused verbatim per this module's own convention,
    # so it has no per-person loop exposed here to checkpoint WITHIN; see module docstring
    # "CHECKPOINTING / RESUME"). A crash during this stage still loses the whole stage's progress,
    # same as before this fix -- only a crash AFTER it completes is now protected.
    hla_ckpt = load_checkpoint(args.out_dir, "hla_labeled_calls", header) if args.resume else None
    if hla_ckpt is not None:
        calls, ref_catalogue_size, mstats, n_removed = hla_ckpt
    else:
        with Timer(args.out_dir, "hla_labeled_calls"):
            calls, ref_catalogue_size, mstats, n_removed = m39mod.build_labeled_calls(hla_args)
        save_checkpoint(args.out_dir, "hla_labeled_calls",
                         (calls, ref_catalogue_size, mstats, n_removed), header)
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
    with Timer(args.out_dir, "kir_identity"):
        (kir_sets_by_level, kir_qc, kir_cds_available,
         n_hap_seen, n_hap_cds_present) = build_kir_identity_sets(
            unrelated_pids, args.kir_outroot, args.workers, refdata=args.refdata,
            out_dir=args.out_dir, header=header, resume=args.resume,
            checkpoint_every=args.checkpoint_every)
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
    _kir_protein_qc_rows = [{"gene": g, "protein_catalogue_covered": ok, "reason": reason}
                             for g, (ok, reason) in sorted(kir_protein_status.items())]
    _disc.write_tsv(pd.DataFrame(_kir_protein_qc_rows),
                     os.path.join(args.out_dir, "kir_protein_catalogue_qc.tsv"), _DISCLOSURE_MODE)

    log("[44] building HLA identity sets (cds/protein/protein_novel levels) ...")
    hla_sets_by_level = {lvl: hla_person_gene_sets(calls, lvl, hla_genes)
                         for lvl in LEVELS if lvl not in LEVELS_NEEDING_TRIMMED_FASTA}

    # 2026-09-27b TRUE-GENOMIC fix (coordinator follow-up): 'genomic'/'any_novel' are no longer
    # aliased to 'cds' -- they are now built from hap{N}.trimmed.fa gene spans + gen.fa.gz, the SAME
    # species-agnostic code path for BOTH species (see build_person_true_genomic_identity's
    # docstring). Loaded ONCE and shared: IPD's gen.fa.gz is expected to carry both gene families
    # (mirrors how CDSseq/*.fa.gz is organized), so one load covers KIR and HLA.
    log("[44] loading genomic reference catalogue (gen.fa.gz) ...")
    gen_catalogue, gen_catalogue_reason = load_genomic_catalogue(args.refdata, m24mod)
    if gen_catalogue is None:
        log(f"[44] WARNING: genomic reference catalogue unavailable ({gen_catalogue_reason}) -- "
            f"'genomic'/'any_novel' baseline richness may still be real (from trimmed.fa alone) "
            f"but 'any_novel' (genomic-level NOVELTY) cannot be determined for any gene without a "
            f"catalogue to compare against; see genomic_span_novel()'s None return.")
    else:
        kir_genes_in_cat = sum(1 for g in kir_genes if g in gen_catalogue)
        log(f"[44] genomic catalogue loaded: {len(gen_catalogue)} genes total, "
            f"{kir_genes_in_cat}/{len(kir_genes)} KIR genes present.")

    log("[44] building KIR true-genomic identity sets (hap{N}.trimmed.fa + gen.fa.gz) ...")
    with Timer(args.out_dir, "kir_true_genomic"):
        (kir_gen_sets, kir_gen_qc, kir_trimmed_available,
         n_hap_seen_gk, n_hap_trim_gk) = build_true_genomic_identity_sets(
            unrelated_pids, args.kir_outroot, set(kir_genes), gen_catalogue, args.workers,
            out_dir=args.out_dir, header=header, resume=args.resume,
            checkpoint_every=args.checkpoint_every, checkpoint_name="kir_true_genomic")
    kir_genomic_available = kir_trimmed_available and gen_catalogue is not None
    if not kir_trimmed_available:
        log("[44] WARNING: hap{N}.trimmed.fa was not found for ANY KIR haplotype -- KIR "
            "'genomic'/'any_novel' will be exported as 'NA' (reason=no_trimmed_fasta).")
    elif gen_catalogue is None:
        log(f"[44] WARNING: no genomic catalogue ({gen_catalogue_reason}) -- KIR "
            f"'genomic'/'any_novel' will be exported as 'NA' (reason={gen_catalogue_reason}).")
    write_status(args.out_dir, f"KIR true-genomic sets built; available={kir_genomic_available} "
                                f"(n_hap_seen={n_hap_seen_gk}, n_hap_trimmed={n_hap_trim_gk})")
    if kir_genomic_available:
        kir_sets_by_level["genomic"] = kir_gen_sets["genomic"]
        kir_sets_by_level["any_novel"] = kir_gen_sets["any_novel"]

    log("[44] building HLA true-genomic identity sets (hap{N}.trimmed.fa + gen.fa.gz) ...")
    with Timer(args.out_dir, "hla_true_genomic"):
        (hla_gen_sets, hla_gen_qc, hla_trimmed_available,
         n_hap_seen_gh, n_hap_trim_gh) = build_true_genomic_identity_sets(
            unrelated_pids, args.hla_people_outroot, set(hla_genes), gen_catalogue, args.workers,
            out_dir=args.out_dir, header=header, resume=args.resume,
            checkpoint_every=args.checkpoint_every, checkpoint_name="hla_true_genomic")
    hla_genomic_available = hla_trimmed_available and gen_catalogue is not None
    if not hla_trimmed_available:
        log("[44] WARNING: hap{N}.trimmed.fa was not found for ANY HLA haplotype -- HLA "
            "'genomic'/'any_novel' will be exported as 'NA' (reason=no_trimmed_fasta).")
    write_status(args.out_dir, f"HLA true-genomic sets built; available={hla_genomic_available} "
                                f"(n_hap_seen={n_hap_seen_gh}, n_hap_trimmed={n_hap_trim_gh})")
    if hla_genomic_available:
        hla_sets_by_level["genomic"] = hla_gen_sets["genomic"]
        hla_sets_by_level["any_novel"] = hla_gen_sets["any_novel"]

    # 2026-09-27c: per-LEVEL checkpoint of the accumulated output rows (module docstring
    # "CHECKPOINTING / RESUME", stage "each finished output table") -- this loop's per-gene work is
    # in-memory-only (the expensive per-person disk work already happened, and was already
    # checkpointed, in the three build_*_identity_sets() stages above), but it can still run for a
    # while across 17+ KIR genes x ~7 HLA genes x 5 levels x up to 25 permutations each, so a level
    # boundary is a natural, cheap place to persist progress: a crash mid-loop resumes at the next
    # un-finished level instead of re-running every earlier one.
    level_ckpt = load_checkpoint(args.out_dir, "gene_level_rows", header) if args.resume else None
    if level_ckpt is not None:
        levels_done = set(level_ckpt["levels_done"])
        all_rec, all_curve, all_slope, all_cov = (
            level_ckpt["all_rec"], level_ckpt["all_curve"], level_ckpt["all_slope"],
            level_ckpt["all_cov"])
        log(f"[44][checkpoint] gene_level_rows: resuming with levels already done = "
            f"{sorted(levels_done)}")
    else:
        levels_done = set()
        all_rec, all_curve, all_slope, all_cov = [], [], [], []
    for level in LEVELS:
        if level in levels_done:
            log(f"[44] level={level} already done (from checkpoint) -- skipping.")
            continue
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
            elif level in LEVELS_NEEDING_TRIMMED_FASTA and not kir_genomic_available:
                rec, curve, slope, cov = build_na_rows(gene, level, "kir", pids_by_ancestry,
                                                        n_star_by_ancestry)
            else:
                full_gene_sets = kir_sets_by_level[level].get(gene, {})
                rec, curve, slope, cov = run_gene_level(
                    gene, level, "kir", full_gene_sets, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, args.curve_stride, args.extrapolate_2n)
            all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
        for gene in hla_genes:
            if level in LEVELS_NEEDING_TRIMMED_FASTA and not hla_genomic_available:
                rec, curve, slope, cov = build_na_rows(gene, level, "hla", pids_by_ancestry,
                                                        n_star_by_ancestry)
                all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
                continue
            full_gene_sets = hla_sets_by_level[level].get(gene, {})
            rec, curve, slope, cov = run_gene_level(
                gene, level, "hla", full_gene_sets, pids_by_ancestry, orders_by_ancestry,
                n_star_by_ancestry, args.curve_stride, args.extrapolate_2n)
            all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
        levels_done.add(level)
        save_checkpoint(args.out_dir, "gene_level_rows",
                         {"levels_done": sorted(levels_done), "all_rec": all_rec,
                          "all_curve": all_curve, "all_slope": all_slope, "all_cov": all_cov},
                         header)
        write_status(args.out_dir, f"level={level} done ({len(all_rec)} recurrence rows so far)")

    sanity_check_coverage(all_cov)
    # HARD gate (task item 3): unlike sanity_check_coverage (warns on known bug SIGNATURES), this
    # RAISES if the identity hierarchy (protein <= cds == genomic; novel <= baseline) is violated
    # on the actual exported numbers -- the exact class of impossibility that triggered the
    # 2026-09-27 identity fix (KIR: 6,300 "proteins" > 1,460 "genomic" alleles).
    check_identity_invariants(all_cov, all_rec)

    _disc.write_tsv(pd.DataFrame(all_rec), os.path.join(args.out_dir, "recurrence_classes.tsv"),
                     _DISCLOSURE_MODE)
    _disc.write_tsv(pd.DataFrame(all_curve), os.path.join(args.out_dir, "saturation_curves.tsv"),
                     _DISCLOSURE_MODE)
    _disc.write_tsv(pd.DataFrame(all_slope), os.path.join(args.out_dir, "equal_n_slope.tsv"),
                     _DISCLOSURE_MODE)
    _disc.write_tsv(pd.DataFrame(all_cov), os.path.join(args.out_dir, "coverage_chao2.tsv"),
                     _DISCLOSURE_MODE)

    (kir_qc_scalar, kir_n_calls_by_gene, kir_artifact_by_gene_label,
     kir_name_to_hashes) = split_kir_qc(kir_qc)
    _disc.write_tsv(
        pd.DataFrame([{**kir_qc_scalar, "n_hap_seen": n_hap_seen,
                       "n_hap_cds_fasta_present": n_hap_cds_present,
                       "cds_available": kir_cds_available}]),
        os.path.join(args.out_dir, "kir_cds_match_qc.tsv"), _DISCLOSURE_MODE)

    # New (2026-09-27 identity-impossibility fix, task items 2/4): artifact counts and name/hash
    # collapse diagnostics, symmetric across species.
    artifact_rows = build_artifact_qc_rows(kir_artifact_by_gene_label, calls)
    _disc.write_tsv(pd.DataFrame(artifact_rows), os.path.join(args.out_dir, "artifact_qc.tsv"),
                     _DISCLOSURE_MODE)
    diag_rows = build_diagnostics_identity_rows(
        all_cov, kir_n_calls_by_gene, kir_artifact_by_gene_label, kir_name_to_hashes, calls)
    _disc.write_tsv(pd.DataFrame(diag_rows), os.path.join(args.out_dir, "diagnostics_identity.tsv"),
                     _DISCLOSURE_MODE)

    # 2026-09-27b: genomic-span-level artifact counts (contig_edge_truncated, partial_cds,
    # inframe_stop, homopolymer_indel, unresolved, no_catalogue_for_gene, clean) -- separate from
    # artifact_qc.tsv's CDS-level counts since they're a different identity level's gate.
    genomic_artifact_rows = [
        {"species": species, "gene": gene, "artifact_label": label, "n_calls": int(n)}
        for qc, species in ((kir_gen_qc, "kir"), (hla_gen_qc, "hla"))
        for (kind, gene, label), n in qc.items() if kind == "artifact_genomic"
    ]
    _disc.write_tsv(pd.DataFrame(genomic_artifact_rows),
                     os.path.join(args.out_dir, "genomic_artifact_qc.tsv"), _DISCLOSURE_MODE)
    _disc.write_tsv(
        pd.DataFrame([{"gen_catalogue_status": gen_catalogue_reason,
                       "kir_genomic_available": kir_genomic_available,
                       "hla_genomic_available": hla_genomic_available,
                       "n_hap_seen_kir": n_hap_seen_gk, "n_hap_trimmed_present_kir": n_hap_trim_gk,
                       "n_hap_seen_hla": n_hap_seen_gh, "n_hap_trimmed_present_hla": n_hap_trim_gh}]),
        os.path.join(args.out_dir, "genomic_identity_qc.tsv"), _DISCLOSURE_MODE)

    write_status(args.out_dir, f"DONE in {time.perf_counter()-t0:.0f}s -- 9 tables written")
    log(f"[44] done in {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
