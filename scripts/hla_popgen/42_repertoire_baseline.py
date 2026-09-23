#!/usr/bin/env python3
"""WS6 — BenchRep-T-style repertoire-only disease-classification baseline (VJ-kmer).

Faithful replication of the VJ-kmer track of Im, Cohen-Lavi, Buendia, Kundaje, Boyd
("BenchRep-T", bioRxiv 10.64898/2026.06.09.727013v1) on AoU whole-blood RNA-seq TRUST4
repertoires. Full protocol extraction, with section refs: sprints/S03_call8_figures_kir_prediction/
WS6_benchrept_protocol.md -- read that first; this docstring only covers what's AoU-specific.

## Pipeline

1. Load per-sample TRUST4 `*_report.tsv` (columns: read_count, frequency, CDR3_dna,
   CDR3_amino_acids, V, D, J, C, consensus_id, consensus_id_complete_vdj), productive-only
   (drop `_`/`?` in CDR3_amino_acids), TRB and IGH kept separately and combined.
2. VJ-kmer features per BenchRep-T: V usage and J usage each normalized to a separate relative-
   frequency simplex; CDR3 (first/last residue trimmed) 4-mers + single-position-gapped 4-mers,
   normalized by total k-mer count per sample. Clone-level (unique CDR3+V+J+C row), not
   read-weighted -- report.tsv is already one row per clone with its own read_count/frequency,
   so "unique sequences" here means "one vote per report.tsv row," matching BenchRep-T's
   unique-sequence (not abundance-weighted) VJ-kmer convention. Because the exact BenchRep-T
   k-mer alphabet (4-mers + 4 gapped variants over ~20 amino acids) is combinatorially large
   (should be full: 20^4 * 5 ~= 800K columns), we cap the k-mer vocabulary to the top
   `--max-kmers` (default 4000) most-variable columns across the cohort being modeled --
   this is a pragmatic, documented deviation from the paper (which appears to use the full
   alphabet on much larger reference panels); V/J usage columns are never capped (small,
   fixed gene lists).
3. Phenotype labels via BigQuery over the AoU CDR (chronic/autoimmune per Cole's 2026-09-22
   guidance: HIV, B-cell disorders, RA/SLE/other rheum autoimmune, T1D/autoimmune thyroid, MS,
   IBD -- reusing 12_disease_phenotypes.py's ICD3+SNOMED-substring pattern and, where
   overlapping, Aleix's aleix/RNA-seq/scripts/query_overlap_phenotypes.py CDR-access pattern).
   Cases = >=2 condition_occurrence rows matching a phenotype's codes; controls = zero rows
   matching ANY of the immune-disease phenotypes below (a clean-control policy, stricter than
   "just not this one disease"). Only phenotypes with >=100 cases are modeled.
4. Four model variants per phenotype: covariate-only, repertoire-only, repertoire+covariate,
   and a depth-rarefied repertoire-only variant (each person's clones downsampled by read_count
   to a common target before recomputing features, nested-subset draws per BenchRep-T's
   Appendix B.4 fixed-seed nesting).
5. 3-fold stratified CV, SAME fold assignment shared across every model variant (BenchRep-T
   Section 2.2), repeated over `--n-repeats` (default 3) different fixed seeds, pooled
   out-of-fold AUROC/AUPRC per repeat with bootstrap CIs (1000 person-level resamples).
6. Outputs: aggregate metrics TSV, top-feature TSV (V/J genes and k-mers by |coefficient| /
   XGBoost importance), and a Nature-style figure (AUROC per phenotype x model with CIs,
   covariate-only baseline as a reference line). NOTHING per-person leaves the VM in real mode
   -- only these three aggregate files.

## Running on the VM (real mode)

    cd ~/s03 && PYTHONPATH=~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 42_repertoire_baseline.py \\
        --cohort-tsv ~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv \\
        --results-bucket gs://aleix-rnaseq-wb-cordial-leechee-9743/repertoire_results \\
        --cdr "$WORKSPACE_CDR" --project "$GOOGLE_PROJECT" \\
        --scratch ~/pipeline_outputs/rnaseq/ws6_scratch \\
        --out-dir ~/s03/results/42

`--cohort-tsv` needs a `research_id` column (and ideally `ancestry`, `sex_at_birth`,
`year_of_birth` if already joined -- otherwise this script queries `person` itself).
Downloads each person's `*_report.tsv` from `--results-bucket/<research_id>/` via
`gcloud storage cp` into `--scratch` (falls back to `gsutil cp`), deletes the local copy once
that person's features are extracted (keeps disk bounded the same way run_rnaseq_batch_local.sh
does), builds phenotype labels via BigQuery (needs `--cdr`, defaults to `$WORKSPACE_CDR`), and
writes ONLY the three aggregate files above under `--out-dir`.

## Local dry run (no VM, no data, no network)

    python3 scripts/hla_popgen/42_repertoire_baseline.py --synthetic --n-people 300 \\
        --out-dir /tmp/ws6_synthetic

Generates a synthetic cohort of TRUST4-shaped report.tsv rows in memory (fixed seed), synthetic
phenotype labels with a small planted repertoire signal, synthetic covariates, and runs the
identical pipeline end to end -- the only difference from real mode is the data source.
Requires scikit-learn (and, optionally, xgboost -- skipped gracefully with a warning if absent).
"""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

try:
    from _viz_common import (nature_style, mm, save_fig, NATURE_DOUBLE_COL_MM,
                              load_cohort_membership, ANCESTRY_ORDER)
    HAVE_VIZ_COMMON = True
except ImportError:
    HAVE_VIZ_COMMON = False
    ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


# Reuse 24's greedy_unrelated / relatedness loader -- lazy (24 pulls in matplotlib etc.), only
# touched in real-mode cohort building, same pattern as 29_hla_ld_by_ancestry.py's m24().
_m24 = None


def m24():
    global _m24
    if _m24 is None:
        _m24 = _load_module("24_novelty_by_field.py", "novelty_by_field_for42")
    return _m24


DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
DEFAULT_COHORT_MEMBERSHIP = os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv")
KIN_MIN = 0.0442  # third degree or closer, same threshold as 24/29
STRICT_MIN = 0.9

try:
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score, average_precision_score
    from sklearn.preprocessing import StandardScaler
    HAVE_SKLEARN = True
except ImportError:
    HAVE_SKLEARN = False

try:
    import xgboost as xgb
    HAVE_XGBOOST = True
except Exception:
    # ImportError if the package isn't installed; on some platforms (e.g. macOS without
    # libomp) the import raises xgboost's own XGBoostError instead -- either way, this
    # script must degrade gracefully rather than crash, since XGBoost is optional
    # (BenchRep-T's "VJ-kmer XG" track) and the L1-LR track ("VJ-kmer Reg") is the
    # faithfully-replicated baseline this script always runs.
    HAVE_XGBOOST = False


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Phenotype definitions -- Cole's 2026-09-22 guidance: chronic/autoimmune, avoid cancer.
# Same (label, icd3_codes, name_substrings) shape as 12_disease_phenotypes.py's HLA_LINKED,
# reused for consistency; a person matches a phenotype if EITHER an ICD-10 3-char code OR a
# SNOMED concept-name substring matches. "Controls" for every phenotype modeled = people who
# match ZERO of the phenotypes below (a clean immune-disease-free control pool), not just
# "doesn't have this one disease" -- avoids a control group secretly enriched for a sibling
# autoimmune condition.
# ---------------------------------------------------------------------------
IMMUNE_PHENOTYPES = [
    ("HIV", ["B20"], ["human immunodeficiency virus", "hiv disease"]),
    ("B-cell disorder", ["D89", "C90"], ["monoclonal gammopathy", "common variable immunodeficiency"]),
    ("Rheumatoid arthritis", ["M05", "M06"], ["rheumatoid arthritis"]),
    ("Systemic lupus erythematosus", ["M32"], ["systemic lupus", "lupus erythematosus"]),
    ("Other rheum. autoimmune", ["M35"], ["sjogren", "sicca syndrome", "systemic sclerosis", "scleroderma"]),
    ("Type 1 diabetes", ["E10"], ["type 1 diabetes"]),
    ("Autoimmune thyroid disease", ["E06"], ["hashimoto", "graves", "autoimmune thyroiditis"]),
    ("Multiple sclerosis", ["G35"], ["multiple sclerosis"]),
    ("Inflammatory bowel disease", ["K50", "K51"], ["crohn", "ulcerative colitis", "inflammatory bowel disease"]),
]
MIN_CASES = 100
MIN_CELL = 20  # never write a raw count below this; write "<20" instead


# ---------------------------------------------------------------------------
# TRUST4 report.tsv loading + VJ-kmer feature construction
# ---------------------------------------------------------------------------
AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")


def is_productive(cdr3_aa):
    if not isinstance(cdr3_aa, str) or len(cdr3_aa) < 2:
        return False
    return ("_" not in cdr3_aa) and ("?" not in cdr3_aa)


def chain_of(v_gene):
    """TRB / IGH / other, from the V gene name prefix (TRUST4's `V` column)."""
    if not isinstance(v_gene, str):
        return "other"
    if v_gene.startswith("TRB"):
        return "TRB"
    if v_gene.startswith("IGH"):
        return "IGH"
    return "other"


def trim_cdr3(seq):
    """Drop the flanking conserved Cys/Phe (first + last residue) -- BenchRep-T Appendix B.1,
    matching Mal-ID formatting. No-op (returns "") if too short to trim meaningfully."""
    if not isinstance(seq, str) or len(seq) < 3:
        return ""
    return seq[1:-1]


def gapped_4mers(trimmed):
    """4-mers plus their 4 single-position-gapped variants, over the trimmed CDR3.
    BenchRep-T Appendix B.2.2. '_' marks the gapped position."""
    out = []
    n = len(trimmed)
    for i in range(n - 3):
        km = trimmed[i:i + 4]
        out.append(km)
        for pos in range(4):
            out.append(km[:pos] + "_" + km[pos + 1:])
    return out


# Column-name aliases seen across TRUST4 versions/invocations on this cohort's real
# report.tsv files (confirmed on the VM 2026-09-23: `#count, frequency, CDR3nt, CDR3aa, V,
# D, J, C, cid, cid_full_length` -- NOT the `read_count/CDR3_amino_acids/...` names this
# script's docstring assumed from the protocol doc). Map whichever alias is present to the
# canonical names used throughout the rest of this file.
COLUMN_ALIASES = {
    "read_count": ["read_count", "#count", "count"],
    "CDR3_amino_acids": ["CDR3_amino_acids", "CDR3aa", "CDR3_aa"],
    "CDR3_dna": ["CDR3_dna", "CDR3nt", "CDR3_nt"],
}


def _resolve_column(df, canonical):
    for alias in COLUMN_ALIASES.get(canonical, [canonical]):
        if alias in df.columns:
            return alias
    return None


def load_report_tsv(path, chains=("TRB", "IGH")):
    """Load one person's TRUST4 report.tsv, productive-only, restricted to `chains`.
    Returns a DataFrame with columns: V, J, CDR3_amino_acids, chain, read_count."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    rc_col = _resolve_column(df, "read_count")
    aa_col = _resolve_column(df, "CDR3_amino_acids")
    missing = ({"V", "J"} - set(df.columns)) | ({"read_count"} if rc_col is None else set()) \
        | ({"CDR3_amino_acids"} if aa_col is None else set())
    if missing:
        die(f"{path}: report.tsv missing expected column(s) {missing}. "
            f"Actual columns: {list(df.columns)}. If TRUST4's schema changed, update "
            f"load_report_tsv()/COLUMN_ALIASES.")
    if aa_col != "CDR3_amino_acids":
        df["CDR3_amino_acids"] = df[aa_col]
    if rc_col != "read_count":
        df["read_count"] = df[rc_col]
    df = df[df["CDR3_amino_acids"].apply(is_productive)].copy()
    df["chain"] = df["V"].apply(chain_of)
    df = df[df["chain"].isin(chains)].copy()
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0).astype(int)
    return df[["V", "J", "CDR3_amino_acids", "chain", "read_count"]]


def rarefy_clones(df, target_reads, seed):
    """Downsample a person's clone table to `target_reads` total reads, without replacement,
    nested-subset draw (BenchRep-T Appendix B.4: fixed seed, D reads at depth D are a subset
    of any larger depth for the SAME person -- achieved here by drawing a single fixed-seed
    permutation of that person's reads once and taking prefixes, rather than redrawing per
    depth). Returns None if the person has fewer reads than target_reads (excluded, not
    padded)."""
    total = int(df["read_count"].sum())
    if total < target_reads:
        return None
    rng = np.random.default_rng(seed)
    # Expand to one row per read (memory-bounded: TRUST4 repertoires are thousands, not
    # millions, of clones/reads at the CDR3 level), permute once, take the first target_reads,
    # recount per-clone.
    read_idx = np.repeat(np.arange(len(df)), df["read_count"].to_numpy())
    perm = rng.permutation(read_idx)
    keep = perm[:target_reads]
    counts = np.bincount(keep, minlength=len(df))
    out = df.copy()
    out["read_count"] = counts
    return out[out["read_count"] > 0].reset_index(drop=True)


def build_vj_kmer_features(person_clone_tables, max_kmers=4000):
    """person_clone_tables: {person_id: DataFrame (V, J, CDR3_amino_acids, chain, read_count)}.

    Returns (X: DataFrame [person_id index, feature columns], feature_kinds: {col: 'V'|'J'|'kmer'}).
    V usage and J usage are each normalized to their own relative-frequency simplex, PER CHAIN
    (TRB_V_*, TRB_J_*, IGH_V_*, IGH_J_*, plus a combined TRB+IGH block) -- separate namespaces so
    a TRBV gene column is never confused with an IGHV gene column. k-mer vocabulary is built
    from the FULL cohort (label-blind -- gene/kmer identity does not depend on any phenotype),
    then capped to the top `max_kmers` by cross-sample variance, which only uses X itself, never
    y -- keeps the "features never see test-fold labels" guarantee true regardless of which CV
    fold a person later lands in.
    """
    per_person_rows = {}
    v_vocab, j_vocab, kmer_vocab = set(), set(), set()
    raw = {}
    for pid, df in person_clone_tables.items():
        if df.empty:
            raw[pid] = {"v": {}, "j": {}, "kmer": {}}
            continue
        v_counts = df.groupby("V").size()
        j_counts = df.groupby("J").size()
        v_freq = (v_counts / v_counts.sum()).to_dict() if v_counts.sum() else {}
        j_freq = (j_counts / j_counts.sum()).to_dict() if j_counts.sum() else {}
        kmer_counter = {}
        for cdr3 in df["CDR3_amino_acids"]:
            trimmed = trim_cdr3(cdr3)
            for km in gapped_4mers(trimmed):
                kmer_counter[km] = kmer_counter.get(km, 0) + 1
        total_kmers = sum(kmer_counter.values())
        kmer_freq = ({k: c / total_kmers for k, c in kmer_counter.items()}
                     if total_kmers else {})
        raw[pid] = {"v": v_freq, "j": j_freq, "kmer": kmer_freq}
        v_vocab.update(v_freq)
        j_vocab.update(j_freq)
        kmer_vocab.update(kmer_freq)

    v_vocab, j_vocab = sorted(v_vocab), sorted(j_vocab)
    kmer_vocab = sorted(kmer_vocab)

    pids = sorted(raw)
    v_mat = pd.DataFrame(
        [[raw[p]["v"].get(g, 0.0) for g in v_vocab] for p in pids],
        index=pids, columns=[f"V__{g}" for g in v_vocab])
    j_mat = pd.DataFrame(
        [[raw[p]["j"].get(g, 0.0) for g in j_vocab] for p in pids],
        index=pids, columns=[f"J__{g}" for g in j_vocab])
    kmer_mat = pd.DataFrame(
        [[raw[p]["kmer"].get(k, 0.0) for k in kmer_vocab] for p in pids],
        index=pids, columns=[f"KMER__{k}" for k in kmer_vocab])

    if kmer_mat.shape[1] > max_kmers:
        variances = kmer_mat.var(axis=0)
        keep_cols = variances.sort_values(ascending=False).index[:max_kmers]
        kmer_mat = kmer_mat[sorted(keep_cols)]

    X = pd.concat([v_mat, j_mat, kmer_mat], axis=1)
    feature_kinds = {c: ("V" if c.startswith("V__") else "J" if c.startswith("J__") else "kmer")
                      for c in X.columns}
    return X, feature_kinds


# ---------------------------------------------------------------------------
# Real-mode data acquisition (VM only)
# ---------------------------------------------------------------------------
def download_report(research_id, results_bucket, scratch_dir):
    """`gcloud storage cp` (falls back to gsutil) the person's report.tsv to scratch_dir;
    returns the local path or None if not found. Deletes nothing here -- caller deletes after
    feature extraction, same disk-bounded discipline as run_rnaseq_batch_local.sh.

    If results_bucket is NOT a gs:// URL, it's a local (e.g. gcsfuse-mounted) directory --
    `gcloud storage`/`gsutil` both refuse local-to-local copies ("meant for cloud operations"),
    so read the file in place instead of copying. The caller must not delete a path returned
    this way (it's the real mounted file, not a scratch copy) -- see the `is_local` check in
    main()'s per-person loop."""
    remote = f"{results_bucket.rstrip('/')}/{research_id}/{research_id}_report.tsv"
    if not results_bucket.startswith("gs://"):
        return remote if os.path.exists(remote) else None
    local = os.path.join(scratch_dir, f"{research_id}_report.tsv")
    for cmd in (["gcloud", "storage", "cp", remote, local], ["gsutil", "cp", remote, local]):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and os.path.exists(local):
                return local
        except FileNotFoundError:
            continue
    print(f"  WARN: could not fetch {remote}", file=sys.stderr)
    return None


def list_bucket_sample_ids(results_bucket):
    """`gcloud storage ls` (falls back to gsutil) the top-level per-sample directories under
    results_bucket -- Aleix's repertoire_results/<sample_id>/ layout. Read-only: never writes
    to results_bucket. Returns a sorted list of sample-directory basenames.

    If results_bucket is a local (e.g. gcsfuse-mounted) directory rather than a gs:// URL,
    lists it directly -- `gcloud storage ls`/`gsutil ls` both refuse local paths."""
    if not results_bucket.startswith("gs://"):
        if not os.path.isdir(results_bucket):
            die(f"{results_bucket}: not a gs:// URL and not a local directory.")
        return sorted(d for d in os.listdir(results_bucket)
                     if os.path.isdir(os.path.join(results_bucket, d)))
    prefix = results_bucket.rstrip("/") + "/"
    for cmd in (["gcloud", "storage", "ls", prefix], ["gsutil", "ls", prefix]):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0:
                ids = []
                for line in r.stdout.splitlines():
                    line = line.strip().rstrip("/")
                    if not line or line == prefix.rstrip("/"):
                        continue
                    ids.append(line.rsplit("/", 1)[-1])
                if ids:
                    return sorted(set(ids))
        except FileNotFoundError:
            continue
    die(f"Could not list {prefix} with gcloud storage or gsutil.")


def find_rnaseq_manifest():
    """Locate the RNA-seq manifest.tsv under the controlled-tier workspace mount. ENVIRONMENT.md
    quirk #35: never use ~/mnt/aou-controlled (stale, hangs). Searches under
    ~/workspace/vwb-aou-datasets-controlled-v9/*/multiomics/rnaseq/ for a file literally named
    manifest.tsv, since the exact subpath has moved before."""
    root = os.path.expanduser("~/workspace")
    if not os.path.isdir(root):
        die(f"{root} not found -- is the controlled-tier bucket mounted (RUNBOOK.md)?")
    hits = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if "aou-controlled" in dirpath and "mnt" in dirpath:
            continue  # quirk #35: never touch the stale manual mount
        if "manifest.tsv" in filenames and "rnaseq" in dirpath:
            hits.append(os.path.join(dirpath, "manifest.tsv"))
    if not hits:
        die(f"No rnaseq manifest.tsv found under {root}. Searched recursively for */multiomics/"
            f"rnaseq/manifest.tsv-shaped paths; none matched.")
    return sorted(hits)[0]


def build_rna_lr_cohort(results_bucket, manifest_path, cohort_membership_path,
                        relatedness_path, kin_min=KIN_MIN, strict_min=STRICT_MIN,
                        skip_unrelated=False):
    """Build the RNA-seq x long-read HLA-calling overlap cohort from scratch (this VM instance
    has no pre-built lr_rnaseq_overlap_cohort.tsv -- Aleix's is per-instance, not shared).

    1. List sample dirs under results_bucket (read-only).
    2. Join sample -> research_id via the RNA-seq manifest (flexible column names: tries
       sample_id/sample/aou_sample_id x research_id/person_id).
    3. Intersect with cohort_membership.tsv's in_lr==True people (long-read HLA calls exist).
    4. Restrict to the greedy-unrelated set (24's algorithm, same kinship threshold) unless
       skip_unrelated (local/dry-run testing only).
    5. Attach ancestry_pred + admixture proportions from cohort_membership.tsv.

    Returns a DataFrame: research_id, ancestry (uppercased ancestry_pred).
    """
    sample_ids = list_bucket_sample_ids(results_bucket)
    print(f"  {len(sample_ids)} sample dirs under {results_bucket}", file=sys.stderr)

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    sample_col = next((c for c in ("sampleid", "sample_id", "sample", "aou_sample_id", "sample_name")
                       if c in manifest.columns), None)
    rid_col = next((c for c in ("research_id", "person_id") if c in manifest.columns), None)
    if sample_col is None or rid_col is None:
        die(f"{manifest_path}: could not find a sample-id / research-id column pair. "
            f"Columns: {list(manifest.columns)}")
    manifest = manifest[[sample_col, rid_col]].dropna().drop_duplicates()
    sample_to_rid = dict(zip(manifest[sample_col].astype(str), manifest[rid_col].astype(str)))

    matched_ids = sorted({sample_to_rid[s] for s in sample_ids if s in sample_to_rid})
    unmatched = len(sample_ids) - sum(1 for s in sample_ids if s in sample_to_rid)
    print(f"  {len(matched_ids)}/{len(sample_ids)} sample dirs matched a research_id via "
          f"manifest ({unmatched} unmatched)", file=sys.stderr)

    if HAVE_VIZ_COMMON:
        cohort = load_cohort_membership(cohort_membership_path)
    else:
        cohort = pd.read_csv(cohort_membership_path, sep="\t", dtype=str)
    cohort["person_id"] = cohort["person_id"].astype(str)
    if "in_lr" in cohort.columns:
        lr_ids = set(cohort.loc[cohort["in_lr"].astype(str).str.lower().isin(
            ("true", "1", "1.0")), "person_id"])
    else:
        lr_ids = set(cohort["person_id"])  # degrade gracefully if column absent
    overlap_ids = sorted(set(matched_ids) & lr_ids)
    print(f"  {len(overlap_ids)}/{len(matched_ids)} repertoire-matched people also have "
          f"long-read HLA calls (in_lr)", file=sys.stderr)

    if not skip_unrelated:
        pairs = m24().load_relatedness_pairs(relatedness_path)
        kept, removed = m24().greedy_unrelated(overlap_ids, pairs, kin_min=kin_min)
        print(f"  greedy-unrelated: kept {len(kept)}, removed {len(removed)} "
              f"(kin >= {kin_min})", file=sys.stderr)
        overlap_ids = sorted(kept)

    cohort_ind = cohort.set_index("person_id")
    ancestry_pred = cohort_ind["ancestry_pred"].astype(str).str.upper()
    rows = []
    for pid in overlap_ids:
        anc = ancestry_pred.get(pid)
        rows.append({"research_id": pid,
                    "ancestry": anc if isinstance(anc, str) and anc in ANCESTRY_ORDER else None})
    return pd.DataFrame(rows)


def get_cdr(explicit):
    if explicit:
        return explicit
    for var in ("WORKSPACE_CDR", "CDR_STORAGE_PATH"):
        v = os.environ.get(var)
        if v:
            return v
    die("Could not resolve the CDR dataset. Pass --cdr explicitly (ENVIRONMENT.md quirk #5: "
        "not reliably exposed via shell env vars outside a notebook kernel).")


def run_bq(sql, project):
    try:
        import pandas_gbq
        return pandas_gbq.read_gbq(sql, dialect="standard", progress_bar_type=None,
                                   project_id=project)
    except ImportError:
        pass
    from google.cloud import bigquery
    client = bigquery.Client(project=project) if project else bigquery.Client()
    return client.query(sql).to_dataframe()


MIN_EHR_YEARS = 1.0  # controls must have >=1 year of EHR observation (Aleix's access-gradient note)


def fetch_phenotypes_and_covariates(ids, cdr, project):
    """BigQuery pull: demographics (year_of_birth, sex_at_birth) + EHR-depth covariates
    (n_visit_dates, ehr_years from observation_period) + all condition_occurrence rows for
    `ids`, matched against IMMUNE_PHENOTYPES. Mirrors 12_disease_phenotypes.py /
    aleix/RNA-seq/scripts/query_overlap_phenotypes.py's pattern exactly (same CDR-access
    idiom, same ICD3+SNOMED-substring matching).

    EHR-depth covariates matter here because Aleix documented an ancestry-correlated
    healthcare-access gradient in this cohort (more EHR contact -> more condition rows ->
    looks like "more disease" independent of biology) -- so n_visit_dates/ehr_years MUST be
    in the covariate-only baseline, not just available as an afterthought.

    Returns (covariates_df, labels_df) -- labels_df is one row per (research_id, phenotype,
    is_case: bool), case = >=2 matching condition_occurrence rows.
    """
    id_list = ",".join(ids)
    sql_person = f"""
    SELECT p.person_id AS research_id, p.year_of_birth,
           c_sex.concept_name AS sex_at_birth
    FROM `{cdr}.person` p
    LEFT JOIN `{cdr}.concept` c_sex ON c_sex.concept_id = p.sex_at_birth_concept_id
    WHERE p.person_id IN ({id_list})
    """
    sql_cond = f"""
    SELECT co.person_id AS research_id, co.condition_source_value, c.concept_name AS condition_name
    FROM `{cdr}.condition_occurrence` co
    LEFT JOIN `{cdr}.concept` c ON c.concept_id = co.condition_concept_id
    WHERE co.person_id IN ({id_list})
    """
    # EHR depth: distinct visit dates (visit_occurrence) and observation-window years
    # (observation_period end - start, summed across periods per person -- most people have one).
    sql_ehr_depth = f"""
    SELECT person_id AS research_id, COUNT(DISTINCT visit_start_date) AS n_visit_dates
    FROM `{cdr}.visit_occurrence`
    WHERE person_id IN ({id_list})
    GROUP BY person_id
    """
    sql_obs_window = f"""
    SELECT person_id AS research_id,
           SUM(DATE_DIFF(observation_period_end_date, observation_period_start_date, DAY)) / 365.25
             AS ehr_years
    FROM `{cdr}.observation_period`
    WHERE person_id IN ({id_list})
    GROUP BY person_id
    """
    covariates = run_bq(sql_person, project)
    cond = run_bq(sql_cond, project)
    ehr_depth = run_bq(sql_ehr_depth, project)
    obs_window = run_bq(sql_obs_window, project)

    # BigQuery returns person_id/research_id as int64; every id elsewhere in this script
    # (manifest, cohort_membership, --cohort-tsv, the `ids` list) is a string. Normalize to
    # str immediately so every downstream merge/lookup/isin keyed on research_id actually
    # matches instead of silently returning 0 rows (an int64-vs-str merge either raises, as
    # it did against `cohort`, or -- worse -- succeeds-but-empty against a plain dict .get()).
    for _df in (covariates, cond, ehr_depth, obs_window):
        _df["research_id"] = _df["research_id"].astype(str)

    covariates = covariates.merge(ehr_depth, on="research_id", how="left")
    covariates = covariates.merge(obs_window, on="research_id", how="left")
    covariates["n_visit_dates"] = covariates["n_visit_dates"].fillna(0)
    covariates["ehr_years"] = covariates["ehr_years"].fillna(0.0)

    icd3 = cond["condition_source_value"].astype(str).str.upper().str[:3]
    name = cond["condition_name"].astype(str).str.lower()
    rows = []
    for label, icd_codes, substrings in IMMUNE_PHENOTYPES:
        mask = icd3.isin(icd_codes) if icd_codes else pd.Series(False, index=cond.index)
        for s in substrings:
            mask = mask | name.str.contains(s, na=False, regex=False)
        match_counts = cond.loc[mask, "research_id"].value_counts()
        for pid in ids:
            rows.append({"research_id": pid, "phenotype": label,
                         "is_case": int(match_counts.get(pid, 0) >= 2)})
    labels = pd.DataFrame(rows)
    return covariates, labels


def apply_control_ehr_gate(covariates, y, min_years=MIN_EHR_YEARS):
    """Controls (y==0) must have >= min_years of EHR observation, else they're indistinguishable
    from "never showed up to be diagnosed" rather than genuinely disease-free -- an
    ancestry-correlated confound per Aleix's access-gradient finding. Cases are NOT filtered on
    this (a case already proved they have >=2 matching condition rows, i.e. enough EHR contact
    to be diagnosed). Returns a boolean keep-mask aligned to y's index."""
    is_case = pd.Series(np.asarray(y).astype(bool), index=covariates.index)
    if "ehr_years" not in covariates.columns:
        return pd.Series(True, index=covariates.index)  # synthetic/no-EHR mode: gate is a no-op
    ehr_years = covariates["ehr_years"].reindex(covariates.index).fillna(0.0)
    return is_case | (ehr_years >= min_years)


def fetch_hla_carriage(ids, hla_calls_path, alleles=("DRB1*15:01", "B*57:01")):
    """Biology-positive-control label: carriage (>=1 copy, 2-field) of a common HLA allele, from
    the local hla_calls_rich.tsv (never BigQuery -- HLA calls are pipeline output, not in the
    CDR). Returns {allele: {research_id: is_carrier}} for alleles with enough carriers to model;
    picks the first allele in `alleles` with >= MIN_CASES carriers among `ids`."""
    if not os.path.exists(hla_calls_path):
        print(f"  WARN: {hla_calls_path} not found -- skipping HLA positive control.",
              file=sys.stderr)
        return None
    df = pd.read_csv(hla_calls_path, sep="\t", dtype=str)
    pid_col = next((c for c in ("person_id", "research_id") if c in df.columns), None)
    allele_col = next((c for c in ("allele", "allele_2field", "call") if c in df.columns), None)
    if pid_col is None or allele_col is None:
        print(f"  WARN: {hla_calls_path} missing person/allele columns "
              f"({list(df.columns)}) -- skipping HLA positive control.", file=sys.stderr)
        return None
    df[pid_col] = df[pid_col].astype(str)
    ids_set = set(ids)
    two_field = df[allele_col].astype(str).str.extract(r"^([A-Z0-9]+\*\d+:\d+)")[0]
    for allele in alleles:
        carriers = set(df.loc[(two_field == allele) & df[pid_col].isin(ids_set), pid_col])
        if len(carriers) >= MIN_CASES:
            return {"allele": allele, "carriers": carriers}
    return None


# ---------------------------------------------------------------------------
# CV + models
# ---------------------------------------------------------------------------
def make_folds(y, n_splits, seed):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(skf.split(np.zeros(len(y)), y))


def fit_predict_lr(X_train, y_train, X_test, C_grid=(1.0, 0.2, 0.1, 0.05, 0.03), seed=0):
    """L1-LR with an internal 5-fold-stratified C sweep (BenchRep-T Appendix B.2.2), scored by
    AUROC on the internal folds, refit at the winning C on the full training fold."""
    best_c, best_score = C_grid[0], -np.inf
    inner = StratifiedKFold(n_splits=min(5, max(2, np.bincount(y_train).min())),
                            shuffle=True, random_state=seed)
    for C in C_grid:
        scores = []
        for tr_idx, va_idx in inner.split(X_train, y_train):
            if len(np.unique(y_train[va_idx])) < 2:
                continue
            clf = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=2000)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                clf.fit(X_train[tr_idx], y_train[tr_idx])
            p = clf.predict_proba(X_train[va_idx])[:, 1]
            scores.append(roc_auc_score(y_train[va_idx], p))
        if scores and np.mean(scores) > best_score:
            best_score, best_c = np.mean(scores), C
    clf = LogisticRegression(penalty="l1", solver="liblinear", C=best_c, max_iter=2000)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clf.fit(X_train, y_train)
    return clf.predict_proba(X_test)[:, 1], clf, best_c


def fit_predict_xgb(X_train, y_train, X_test, seed=0, n_jobs=2):
    if not HAVE_XGBOOST:
        return None, None
    # Faithful to BenchRep-T's two-stage grid, collapsed to a smaller sweep for tractability
    # on cohorts of our size (n_cases as low as 100 means an 43-combo full grid overfits the
    # inner-CV selection step itself) -- documented deviation, protocol doc notes the full grid.
    best_params, best_score, best_model = None, -np.inf, None
    inner = StratifiedKFold(n_splits=min(5, max(2, np.bincount(y_train).min())),
                            shuffle=True, random_state=seed)
    for max_depth in (3, 5):
        for lr in (0.03, 0.1):
            scores = []
            for tr_idx, va_idx in inner.split(X_train, y_train):
                if len(np.unique(y_train[va_idx])) < 2:
                    continue
                model = xgb.XGBClassifier(max_depth=max_depth, learning_rate=lr,
                                          n_estimators=200, subsample=0.8,
                                          colsample_bytree=0.8, min_child_weight=3,
                                          eval_metric="logloss", verbosity=0,
                                          random_state=seed, n_jobs=n_jobs)
                model.fit(X_train[tr_idx], y_train[tr_idx])
                p = model.predict_proba(X_train[va_idx])[:, 1]
                scores.append(roc_auc_score(y_train[va_idx], p))
            if scores and np.mean(scores) > best_score:
                best_score = np.mean(scores)
                best_params = dict(max_depth=max_depth, learning_rate=lr)
    model = xgb.XGBClassifier(n_estimators=200, subsample=0.8, colsample_bytree=0.8,
                              min_child_weight=3, eval_metric="logloss", verbosity=0,
                              random_state=seed, n_jobs=n_jobs, **best_params)
    model.fit(X_train, y_train)
    return model.predict_proba(X_test)[:, 1], model


def bootstrap_ci(y_true, y_score, metric_fn, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        vals.append(metric_fn(yt, ys))
    if not vals:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def run_cv_pooled(X, y, seeds, n_splits=3, model="lr", feature_names=None, n_jobs=2):
    """Repeated 3-fold stratified CV. SAME fold assignment reused across every model variant
    the caller passes the same (y, seeds) into -- callers must pass identical y/seed ordering
    for covariate-only / repertoire-only / combined so folds line up (checked by the unit
    test). Returns per-repeat pooled OOF (y_true, y_score) plus aggregated feature importances.
    """
    per_repeat = []
    importances = np.zeros(X.shape[1]) if feature_names is not None else None
    n_imp_folds = 0
    for seed in seeds:
        folds = make_folds(y, n_splits, seed)
        oof_score = np.zeros(len(y))
        for tr_idx, te_idx in folds:
            if model == "lr":
                p, clf, _ = fit_predict_lr(X[tr_idx], y[tr_idx], X[te_idx], seed=seed)
                if feature_names is not None:
                    importances += np.abs(clf.coef_.ravel())
                    n_imp_folds += 1
            elif model == "xgb":
                p, clf = fit_predict_xgb(X[tr_idx], y[tr_idx], X[te_idx], seed=seed,
                                         n_jobs=n_jobs)
                if p is None:
                    return None
                if feature_names is not None:
                    importances += clf.feature_importances_
                    n_imp_folds += 1
            else:
                raise ValueError(model)
            oof_score[te_idx] = p
        per_repeat.append((y.copy(), oof_score))
    if n_imp_folds:
        importances /= n_imp_folds
    return per_repeat, importances


def summarize_repeats(per_repeat, seed_for_ci=0):
    aurocs, auprcs = [], []
    pooled_y, pooled_s = [], []
    for y_true, y_score in per_repeat:
        aurocs.append(roc_auc_score(y_true, y_score))
        auprcs.append(average_precision_score(y_true, y_score))
        pooled_y.append(y_true)
        pooled_s.append(y_score)
    pooled_y = np.concatenate(pooled_y)
    pooled_s = np.concatenate(pooled_s)
    auroc_ci = bootstrap_ci(pooled_y, pooled_s, roc_auc_score, seed=seed_for_ci)
    auprc_ci = bootstrap_ci(pooled_y, pooled_s, average_precision_score, seed=seed_for_ci)
    return {
        "auroc_mean": float(np.mean(aurocs)), "auroc_ci_lo": auroc_ci[0], "auroc_ci_hi": auroc_ci[1],
        "auprc_mean": float(np.mean(auprcs)), "auprc_ci_lo": auprc_ci[0], "auprc_ci_hi": auprc_ci[1],
    }


# ---------------------------------------------------------------------------
# Synthetic data (for --synthetic / --dry-run)
# ---------------------------------------------------------------------------
def make_synthetic_cohort(n_people, seed=20260922):
    rng = np.random.default_rng(seed)
    trbv = [f"TRBV{i}" for i in range(1, 15)]
    trbj = [f"TRBJ{i}-{j}" for i in range(1, 3) for j in range(1, 4)]
    ighv = [f"IGHV{i}" for i in range(1, 8)]
    ighj = [f"IGHJ{i}" for i in range(1, 5)]
    person_clone_tables = {}
    covariates_rows = []
    signal_gene = "TRBV3"  # planted repertoire signal for cases
    phenotype_case_rate = {}
    labels_rows = []
    for i in range(n_people):
        pid = f"synthetic_{i:05d}"
        age = int(rng.integers(20, 80))
        sex = rng.choice(["Male", "Female"])
        ancestry = rng.choice(["AFR", "AMR", "EAS", "EUR", "MID", "SAS"])
        depth_reads = int(rng.integers(5000, 60000))
        is_case_signal = rng.random() < 0.25
        n_clones = int(rng.integers(150, 900))
        v_pool = trbv + ighv
        j_lookup = {**{v: trbj for v in trbv}, **{v: ighj for v in ighv}}
        rows = []
        for _ in range(n_clones):
            if is_case_signal and rng.random() < 0.35:
                v = signal_gene
            else:
                v = rng.choice(v_pool)
            j = rng.choice(j_lookup[v])
            cdr3_len = int(rng.integers(10, 18))
            cdr3 = "C" + "".join(rng.choice(AA_ALPHABET, cdr3_len - 2)) + "F"
            reads = int(rng.integers(1, 25))
            rows.append({"V": v, "J": j, "CDR3_amino_acids": cdr3, "read_count": reads})
        df = pd.DataFrame(rows)
        total_reads = int(df["read_count"].sum())
        # rescale to target depth_reads roughly, keep it simple/synthetic
        df["chain"] = df["V"].apply(chain_of)
        person_clone_tables[pid] = df[["V", "J", "CDR3_amino_acids", "chain", "read_count"]]
        covariates_rows.append({"research_id": pid, "age": age, "sex_at_birth": sex,
                                "ancestry": ancestry, "depth_reads": total_reads})
        for label, _, _ in IMMUNE_PHENOTYPES:
            base_rate = 0.08
            p_case = base_rate + (0.35 if (label == "Rheumatoid arthritis" and is_case_signal) else 0.0)
            is_case = int(rng.random() < p_case)
            labels_rows.append({"research_id": pid, "phenotype": label, "is_case": is_case})
    covariates = pd.DataFrame(covariates_rows)
    labels = pd.DataFrame(labels_rows)
    return person_clone_tables, covariates, labels


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
MODEL_ORDER = ["covariate_only", "repertoire_only_TRB", "repertoire_only_IGH",
              "repertoire_only", "repertoire_plus_covariate", "repertoire_rarefied"]
MODEL_COLORS = {"covariate_only": "#999999", "repertoire_only_TRB": "#66c2a5",
                "repertoire_only_IGH": "#8da0cb", "repertoire_only": "#1b9e77",
                "repertoire_plus_covariate": "#d95f02", "repertoire_rarefied": "#7570b3"}


def _bar_panel(ax, sub_df, phenotypes, models, title):
    x = np.arange(len(phenotypes))
    n_models = max(len(models), 1)
    bar_w = 0.8 / n_models
    for mi, m in enumerate(models):
        sub = sub_df[sub_df["model"] == m].set_index("phenotype")
        means = [sub.loc[p, "auroc_mean"] if p in sub.index else np.nan for p in phenotypes]
        lo = [max(sub.loc[p, "auroc_mean"] - sub.loc[p, "auroc_ci_lo"], 0) if p in sub.index else 0
             for p in phenotypes]
        hi = [max(sub.loc[p, "auroc_ci_hi"] - sub.loc[p, "auroc_mean"], 0) if p in sub.index else 0
             for p in phenotypes]
        xpos = x + (mi - n_models / 2) * bar_w + bar_w / 2
        ax.bar(xpos, means, width=bar_w, label=m, color=MODEL_COLORS.get(m, "#333333"),
              yerr=[lo, hi], capsize=2, error_kw={"linewidth": 0.5})
    ax.axhline(0.5, color="black", linewidth=0.5, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(phenotypes, rotation=45, ha="right")
    ax.set_ylabel("Pooled OOF AUROC")
    ax.set_ylim(0.3, 1.0)
    ax.set_title(title, fontsize=6)


def _draw_figure(results_df, out_path, plt, use_nature_style):
    is_pc = results_df["phenotype"].astype(str).str.startswith("positive_control_")
    disease_df, pc_df = results_df[~is_pc], results_df[is_pc]
    disease_phenotypes = sorted(disease_df["phenotype"].unique())
    pc_phenotypes = sorted(pc_df["phenotype"].unique())
    models = [m for m in MODEL_ORDER if m in set(results_df["model"])]

    ctx = nature_style() if use_nature_style else None
    if ctx is not None:
        ctx.__enter__()
    try:
        width = mm(NATURE_DOUBLE_COL_MM) if use_nature_style else 12
        n_panels = 1 + (1 if pc_phenotypes else 0)
        fig, axes = plt.subplots(1, n_panels, figsize=(width, 4),
                                 gridspec_kw={"width_ratios": [3, 1][:n_panels]})
        axes = [axes] if n_panels == 1 else list(axes)
        _bar_panel(axes[0], disease_df, disease_phenotypes, models,
                  "Disease phenotypes (EHR-derived)")
        if pc_phenotypes:
            _bar_panel(axes[1], pc_df, pc_phenotypes, models, "Positive controls")
        axes[0].legend(fontsize=4.5, loc="upper right", frameon=False)
        fig.tight_layout()
        if use_nature_style:
            save_fig(fig, out_path)
        else:
            fig.savefig(out_path + ".png", dpi=200, bbox_inches="tight")
            plt.close(fig)
    finally:
        if ctx is not None:
            ctx.__exit__(None, None, None)


def make_figure(results_df, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if HAVE_VIZ_COMMON:
        try:
            _draw_figure(results_df, out_path, plt, use_nature_style=True)
            return
        except Exception as exc:
            # Nature style requests Helvetica/Arial via matplotlib's font-substitution
            # machinery -- on some local dev machines (confirmed here: a fresh pip-installed
            # matplotlib in a scratch venv) glyph lookup for that family raises
            # `RuntimeError: failed to load glyph` at layout time, environment-specific and
            # unrelated to the statistics. Never let a font issue block getting aggregate
            # metrics out -- fall back to plain matplotlib defaults and say so loudly, so a
            # real run on the VM (where Arial has historically resolved fine for 36/38b's
            # nature_style figures) is not silently masked if this ever recurs there too.
            print(f"WARN: nature_style figure rendering failed ({exc!r}); retrying with "
                 f"plain matplotlib style.", file=sys.stderr)
            plt.close("all")
    _draw_figure(results_df, out_path, plt, use_nature_style=False)


# ---------------------------------------------------------------------------
# Main pipeline (shared by real + synthetic mode)
# ---------------------------------------------------------------------------
MIN_TRB_CLONOTYPES = 100  # depth floor: below this, a sample's TRB repertoire is too shallow
                          # to trust for VJ-kmer features (documented in the QC section of the
                          # README, not just hardcoded silently).


def filter_chain(person_clone_tables, chain):
    return {pid: df[df["chain"] == chain].reset_index(drop=True)
           for pid, df in person_clone_tables.items()}


def qc_clonotype_stats(person_clone_tables):
    """Per-sample productive clonotype counts (unique CDR3+V+J+C rows) per chain -- the
    sanity/QC step BEFORE any modeling. Returns (quantiles_df, frac_below_floor, qc_ids_kept)."""
    rows = []
    for pid, df in person_clone_tables.items():
        n_trb = int((df["chain"] == "TRB").sum())
        n_igh = int((df["chain"] == "IGH").sum())
        rows.append({"research_id": pid, "n_trb_clonotypes": n_trb, "n_igh_clonotypes": n_igh})
    qc = pd.DataFrame(rows)
    if qc.empty:
        return qc, np.nan, []
    quantiles = qc[["n_trb_clonotypes", "n_igh_clonotypes"]].quantile(
        [0.05, 0.25, 0.5, 0.75, 0.95]).reset_index().rename(columns={"index": "quantile"})
    frac_below_floor = float((qc["n_trb_clonotypes"] < MIN_TRB_CLONOTYPES).mean())
    kept_ids = qc.loc[qc["n_trb_clonotypes"] >= MIN_TRB_CLONOTYPES, "research_id"].tolist()
    return quantiles, frac_below_floor, kept_ids


def min_detectable_auroc(n_cases, n_total, alpha=0.05, power=0.80):
    """Minimum AUROC distinguishable from 0.5 at the given (n_cases, n_total), Hanley-McNeil
    normal approximation: Var(AUC | AUC=0.5) ~= (n1+n0+1) / (12*n1*n0). Two-sided alpha, given
    power. Used to caveat null results honestly in the README rather than just reporting
    "AUROC ~ 0.5, no signal" with no sense of how much signal could have been missed."""
    from scipy.stats import norm
    n1 = n_cases
    n0 = n_total - n_cases
    if n1 <= 1 or n0 <= 1:
        return np.nan
    var0 = (n1 + n0 + 1) / (12.0 * n1 * n0)
    z_a = norm.ppf(1 - alpha / 2)
    z_b = norm.ppf(power)
    return 0.5 + (z_a + z_b) * np.sqrt(var0)


def model_one_label(label_name, y_series, X_variants, seeds, args, metrics_rows,
                    top_features_rows, extra_cols=None):
    """Run the shared LR-then-XGBoost variant loop for one binary label (disease phenotype OR
    a positive-control label: sex, ancestry, HLA carriage). X_variants: {model_name: (X, fnames)}
    (fnames=None for covariate-only). All variants share the SAME y/seeds so folds line up
    across variants for this label (BenchRep-T Section 2.2). LR is run for every variant BEFORE
    any XGBoost call, across the whole variants dict, so a slow XGBoost sweep never blocks LR
    numbers for other model variants of this phenotype (VM is shared, <=2 cores)."""
    extra_cols = extra_cols or {}
    y_series = np.asarray(y_series).astype(int)
    lr_results = {}
    for model_name, (Xv, fnames) in X_variants.items():
        if Xv is None or len(np.unique(y_series)) < 2:
            continue
        result = run_cv_pooled(Xv, y_series, seeds, n_splits=args.n_splits, model="lr",
                               feature_names=fnames, n_jobs=args.n_jobs)
        if result is None:
            continue
        lr_results[model_name] = result
        per_repeat, importances = result
        summ = summarize_repeats(per_repeat, seed_for_ci=args.seed)
        n_cases = int(y_series.sum())
        summ.update({"phenotype": label_name, "model": model_name,
                    "n_cases": n_cases if n_cases >= MIN_CELL else f"<{MIN_CELL}",
                    "n_total": len(y_series) if len(y_series) >= MIN_CELL else f"<{MIN_CELL}",
                    "min_detectable_auroc": min_detectable_auroc(n_cases, len(y_series))})
        summ.update(extra_cols)
        metrics_rows.append(summ)
        if fnames is not None and importances is not None:
            order = np.argsort(-importances)[:20]
            for rank, idx in enumerate(order):
                top_features_rows.append({
                    "phenotype": label_name, "model": model_name, "rank": rank + 1,
                    "feature": fnames[idx], "importance": float(importances[idx]),
                })

    if args.with_xgboost and HAVE_XGBOOST:
        for model_name in ("repertoire_only", "repertoire_only_TRB", "repertoire_only_IGH",
                           "repertoire_plus_covariate"):
            if model_name not in X_variants or model_name not in lr_results:
                continue
            Xv, fnames = X_variants[model_name]
            result_xgb = run_cv_pooled(Xv, y_series, seeds, n_splits=args.n_splits, model="xgb",
                                       feature_names=fnames, n_jobs=args.n_jobs)
            if result_xgb is None:
                continue
            per_repeat_xgb, _ = result_xgb
            summ_xgb = summarize_repeats(per_repeat_xgb, seed_for_ci=args.seed)
            lr_summ = metrics_rows[-1] if metrics_rows else {}
            n_cases = int(y_series.sum())
            summ_xgb.update({"phenotype": label_name, "model": model_name + "_xgb",
                            "n_cases": n_cases if n_cases >= MIN_CELL else f"<{MIN_CELL}",
                            "n_total": len(y_series) if len(y_series) >= MIN_CELL else f"<{MIN_CELL}",
                            "min_detectable_auroc": min_detectable_auroc(n_cases, len(y_series))})
            summ_xgb.update(extra_cols)
            metrics_rows.append(summ_xgb)


def run_pipeline(person_clone_tables, covariates, labels, args):
    if not HAVE_SKLEARN:
        die("scikit-learn is required. pip install scikit-learn (xgboost optional).")

    # --- QC first: clonotype depth, before any modeling ---
    qc_quantiles, frac_below_floor, qc_kept_ids = qc_clonotype_stats(person_clone_tables)
    print(f"QC: fraction of samples below the {MIN_TRB_CLONOTYPES}-TRB-clonotype depth floor: "
          f"{frac_below_floor:.3f}", file=sys.stderr)
    person_clone_tables = {pid: df for pid, df in person_clone_tables.items()
                           if pid in set(qc_kept_ids)}
    print(f"QC: {len(person_clone_tables)} samples retained after the depth floor.",
          file=sys.stderr)

    X_full, feature_kinds = build_vj_kmer_features(person_clone_tables, max_kmers=args.max_kmers)
    X_trb, _ = build_vj_kmer_features(filter_chain(person_clone_tables, "TRB"),
                                      max_kmers=args.max_kmers)
    X_igh, _ = build_vj_kmer_features(filter_chain(person_clone_tables, "IGH"),
                                      max_kmers=args.max_kmers)
    covariates = covariates.set_index("research_id")

    # depth covariate: prefer an explicit column, else derive from clone tables
    if "depth_reads" not in covariates.columns:
        covariates["depth_reads"] = pd.Series(
            {pid: int(df["read_count"].sum()) for pid, df in person_clone_tables.items()})

    cov_cols_numeric = [c for c in ("age", "year_of_birth", "depth_reads", "n_visit_dates",
                                    "ehr_years") if c in covariates.columns]
    cov_cat_cols = [c for c in ("sex_at_birth", "ancestry") if c in covariates.columns]
    cov_encoded = covariates[cov_cols_numeric].copy()
    for c in cov_cat_cols:
        dummies = pd.get_dummies(covariates[c], prefix=c)
        cov_encoded = pd.concat([cov_encoded, dummies], axis=1)
    cov_encoded = cov_encoded.reindex(X_full.index).fillna(0.0)

    # rarefied repertoire: downsample every person to a common target read count
    target_reads = int(np.percentile(
        [int(df["read_count"].sum()) for df in person_clone_tables.values()], args.rarefy_pctile))
    rarefied_tables = {}
    for pid, df in person_clone_tables.items():
        r = rarefy_clones(df, target_reads, seed=args.seed)
        if r is not None:
            rarefied_tables[pid] = r
    X_rare, _ = build_vj_kmer_features(rarefied_tables, max_kmers=args.max_kmers) if rarefied_tables else (None, None)

    metrics_rows = []
    top_features_rows = []
    seeds = [args.seed + k for k in range(args.n_repeats)]

    phenotype_counts = labels.groupby("phenotype")["is_case"].sum()
    modeled_phenotypes = [p for p, n in phenotype_counts.items() if n >= MIN_CASES]
    print(f"Phenotypes with >= {MIN_CASES} cases: {modeled_phenotypes}", file=sys.stderr)
    for label, _, _ in IMMUNE_PHENOTYPES:
        n = int(phenotype_counts.get(label, 0))
        flag = n if n >= MIN_CELL else f"<{MIN_CELL}"
        print(f"  n_cases[{label}] = {flag} (modeled: {label in modeled_phenotypes})",
              file=sys.stderr)

    for phenotype in modeled_phenotypes:
        y_series_full = labels[labels["phenotype"] == phenotype].set_index("research_id")["is_case"]
        common_ids = sorted(set(X_full.index) & set(y_series_full.index) & set(cov_encoded.index))
        y_full = y_series_full.loc[common_ids]

        # Control EHR-depth gate: controls (y==0) need >= MIN_EHR_YEARS of observation, else
        # "control" may just mean "never showed up" (Aleix's ancestry-correlated access
        # gradient). Cases are exempt (already proved >=2 condition rows of EHR contact).
        keep_mask = apply_control_ehr_gate(covariates.loc[common_ids], y_full)
        common_ids = [i for i, k in zip(common_ids, keep_mask) if k]
        y_series = y_series_full.loc[common_ids]
        y = y_series.to_numpy().astype(int)
        if len(np.unique(y)) < 2 or y.sum() < MIN_CASES:
            continue

        X_rep = X_full.loc[common_ids].to_numpy(dtype=float)
        X_cov = StandardScaler().fit_transform(cov_encoded.loc[common_ids].to_numpy(dtype=float))
        X_comb = np.concatenate([X_rep, X_cov], axis=1)
        cov_fnames = list(X_full.columns) + [f"cov_{i}" for i in range(X_cov.shape[1])]

        variants = {"covariate_only": (X_cov, None),
                   "repertoire_only": (X_rep, X_full.columns),
                   "repertoire_plus_covariate": (X_comb, cov_fnames)}
        trb_ids = [i for i in common_ids if i in X_trb.index]
        if len(trb_ids) >= MIN_CASES and len(set(y_series.loc[trb_ids])) > 1:
            variants["repertoire_only_TRB"] = (
                X_trb.loc[trb_ids].to_numpy(dtype=float), X_trb.columns)
        igh_ids = [i for i in common_ids if i in X_igh.index]
        if len(igh_ids) >= MIN_CASES and len(set(y_series.loc[igh_ids])) > 1:
            variants["repertoire_only_IGH"] = (
                X_igh.loc[igh_ids].to_numpy(dtype=float), X_igh.columns)
        model_one_label(phenotype, y_series, variants, seeds, args, metrics_rows,
                        top_features_rows)

        if X_rare is not None:
            rare_ids = [i for i in common_ids if i in X_rare.index]
            if len(rare_ids) >= MIN_CASES and len(set(y_series.loc[rare_ids])) > 1:
                y_rare_series = y_series.loc[rare_ids]
                model_one_label(phenotype, y_rare_series,
                                {"repertoire_rarefied": (X_rare.loc[rare_ids].to_numpy(dtype=float),
                                                         X_rare.columns)},
                                seeds, args, metrics_rows, top_features_rows)

    # --- Positive controls: sex, ancestry, HLA carriage from TRB -- shows the pipeline finds
    # signal when it exists, and quantifies ancestry-leakage risk (Dewitt 2018 eLife / Emerson
    # 2017: HLA restriction imprints TCR repertoires -- biology-positive control for the HLA
    # half of the paper). Repertoire-only + covariate-only, same variant machinery. ---
    pc_common = sorted(set(X_full.index) & set(cov_encoded.index) & set(covariates.index))
    if "sex_at_birth" in covariates.columns:
        sex = covariates.loc[pc_common, "sex_at_birth"].astype(str)
        top2 = sex.value_counts().index[:2]
        sex_ids = [i for i in pc_common if sex.get(i) in top2]
        y_sex = (sex.loc[sex_ids] == top2[0]).astype(int)
        if y_sex.sum() >= MIN_CASES and (len(y_sex) - y_sex.sum()) >= MIN_CASES:
            X_cov_sex = StandardScaler().fit_transform(
                cov_encoded.loc[sex_ids].drop(
                    columns=[c for c in cov_encoded.columns if c.startswith("sex_at_birth")],
                    errors="ignore").to_numpy(dtype=float))
            model_one_label(f"positive_control_sex({top2[0]})", y_sex,
                            {"covariate_only": (X_cov_sex, None),
                             "repertoire_only": (X_full.loc[sex_ids].to_numpy(dtype=float), X_full.columns),
                             "repertoire_only_TRB": (X_trb.reindex(sex_ids).fillna(0).to_numpy(dtype=float),
                                                      X_trb.columns)},
                            seeds, args, metrics_rows, top_features_rows)
    if "ancestry" in covariates.columns:
        anc = covariates.loc[pc_common, "ancestry"].astype(str).str.upper()
        for target_anc in ("EUR", "AFR"):
            y_anc = (anc == target_anc).astype(int)
            if y_anc.sum() >= MIN_CASES and (len(y_anc) - y_anc.sum()) >= MIN_CASES:
                model_one_label(f"positive_control_ancestry({target_anc})", y_anc,
                                {"repertoire_only": (X_full.loc[pc_common].to_numpy(dtype=float), X_full.columns),
                                 "repertoire_only_TRB": (X_trb.reindex(pc_common).fillna(0).to_numpy(dtype=float),
                                                          X_trb.columns),
                                 "repertoire_only_IGH": (X_igh.reindex(pc_common).fillna(0).to_numpy(dtype=float),
                                                          X_igh.columns)},
                                seeds, args, metrics_rows, top_features_rows)

    if args.hla_calls and os.path.exists(args.hla_calls):
        hla = fetch_hla_carriage(pc_common, args.hla_calls)
        if hla is not None:
            y_hla = pd.Series({i: int(i in hla["carriers"]) for i in pc_common})
            if y_hla.sum() >= MIN_CASES and (len(y_hla) - y_hla.sum()) >= MIN_CASES:
                model_one_label(f"positive_control_HLA({hla['allele']})", y_hla,
                                {"repertoire_only_TRB": (X_trb.reindex(pc_common).fillna(0).to_numpy(dtype=float),
                                                          X_trb.columns)},
                                seeds, args, metrics_rows, top_features_rows)

    metrics_df = pd.DataFrame(metrics_rows)
    top_features_df = pd.DataFrame(top_features_rows)
    return metrics_df, top_features_df, qc_quantiles, frac_below_floor


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--synthetic", "--dry-run", dest="synthetic", action="store_true",
                    help="Run end to end on an in-memory synthetic cohort. No VM, no network, "
                         "no real data.")
    ap.add_argument("--n-people", type=int, default=300, help="Synthetic mode only.")
    ap.add_argument("--cohort-tsv", default=None,
                    help="Real mode: TSV with a research_id column (e.g. "
                         "lr_rnaseq_overlap_cohort.tsv).")
    ap.add_argument("--results-bucket", default="gs://aleix-rnaseq-wb-cordial-leechee-9743/repertoire_results",
                    help="Real mode: bucket prefix holding <research_id>/<research_id>_report.tsv.")
    ap.add_argument("--scratch", default=os.path.expanduser("~/pipeline_outputs/rnaseq/ws6_scratch"),
                    help="Real mode: VM-local staging dir for downloaded report.tsv files "
                         "(deleted per-person after feature extraction).")
    ap.add_argument("--cdr", default=None, help="Real mode: fully-qualified CDR dataset.")
    ap.add_argument("--project", default=os.environ.get("GOOGLE_PROJECT"))
    ap.add_argument("--out-dir", required=True, help="Where aggregate outputs are written.")
    ap.add_argument("--max-kmers", type=int, default=4000)
    ap.add_argument("--n-splits", type=int, default=3)
    ap.add_argument("--n-repeats", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20260922)
    ap.add_argument("--rarefy-pctile", type=float, default=20.0,
                    help="Common target depth = this percentile of the cohort's per-person "
                         "total read counts (default: 20th pctile, so most people qualify).")
    ap.add_argument("--with-xgboost", action="store_true",
                    help="Also fit XGBoost for repertoire-only / repertoire+covariate variants "
                         "(run only after L1-LR has completed for all phenotypes -- see "
                         "model_one_label).")
    ap.add_argument("--n-jobs", type=int, default=2,
                    help="XGBoost thread count. VM is shared with other agents -- keep <=2.")
    ap.add_argument("--manifest", default=None,
                    help="Real mode, cohort auto-build: RNA-seq manifest.tsv path. Auto-located "
                         "under ~/workspace/*/multiomics/rnaseq/ if omitted.")
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT_MEMBERSHIP,
                    help="Real mode, cohort auto-build: cohort_membership.tsv (in_lr + ancestry).")
    ap.add_argument("--relatedness", default=DEFAULT_RELATEDNESS,
                    help="Real mode, cohort auto-build: kinship pairs table for greedy_unrelated.")
    ap.add_argument("--hla-calls", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"),
                    help="Real mode: rich HLA calls table, for the HLA-carriage positive control.")
    ap.add_argument("--skip-unrelated", action="store_true",
                    help="Cohort auto-build: skip the greedy-unrelated filter (local/test only).")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.synthetic:
        person_clone_tables, covariates, labels = make_synthetic_cohort(args.n_people, seed=args.seed)
    else:
        if args.cohort_tsv:
            cohort = pd.read_csv(os.path.expanduser(args.cohort_tsv), sep="\t", dtype=str)
            if "research_id" not in cohort.columns:
                die(f"{args.cohort_tsv}: no research_id column. Columns: {list(cohort.columns)}")
        else:
            print("No --cohort-tsv given -- building the RNA-seq x long-read overlap cohort "
                  "from scratch (this VM instance).", file=sys.stderr)
            manifest_path = args.manifest or find_rnaseq_manifest()
            print(f"  manifest: {manifest_path}", file=sys.stderr)
            cohort = build_rna_lr_cohort(args.results_bucket, manifest_path,
                                         args.cohort_membership, args.relatedness,
                                         skip_unrelated=args.skip_unrelated)
        ids = sorted(set(cohort["research_id"].dropna()))
        cdr = get_cdr(args.cdr)
        print(f"Cohort: {len(ids)} people. CDR: {cdr}", file=sys.stderr)
        covariates, labels = fetch_phenotypes_and_covariates(ids, cdr, args.project)
        if "ancestry" in cohort.columns:
            covariates = covariates.merge(cohort[["research_id", "ancestry"]], on="research_id", how="left")

        os.makedirs(args.scratch, exist_ok=True)
        person_clone_tables = {}
        for pid in ids:
            local = download_report(pid, args.results_bucket, args.scratch)
            if local is None:
                continue
            person_clone_tables[pid] = load_report_tsv(local)
            # disk-bounded: never keep more than one downloaded report.tsv on disk at once --
            # but never delete a local passthrough path (download_report returns the real
            # mounted file unchanged when --results-bucket isn't gs://; only a genuine scratch
            # copy lives under args.scratch and is safe to remove).
            if os.path.commonpath([os.path.abspath(local), os.path.abspath(args.scratch)]) == \
                    os.path.abspath(args.scratch):
                os.remove(local)
        print(f"Loaded repertoires for {len(person_clone_tables)}/{len(ids)} people.",
              file=sys.stderr)

    metrics_df, top_features_df, qc_quantiles, frac_below_floor = run_pipeline(
        person_clone_tables, covariates, labels, args)

    metrics_path = os.path.join(args.out_dir, "42_repertoire_baseline_metrics.tsv")
    features_path = os.path.join(args.out_dir, "42_repertoire_baseline_top_features.tsv")
    qc_path = os.path.join(args.out_dir, "42_repertoire_baseline_qc.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    top_features_df.to_csv(features_path, sep="\t", index=False)
    if qc_quantiles is not None and not qc_quantiles.empty:
        qc_quantiles["frac_below_trb_depth_floor"] = frac_below_floor
        qc_quantiles["trb_depth_floor"] = MIN_TRB_CLONOTYPES
        qc_quantiles.to_csv(qc_path, sep="\t", index=False)
        print(f"Wrote {qc_path}", file=sys.stderr)
    print(f"Wrote {metrics_path}", file=sys.stderr)
    print(f"Wrote {features_path}", file=sys.stderr)

    if not metrics_df.empty:
        fig_path = os.path.join(args.out_dir, "42_repertoire_baseline_auroc")
        make_figure(metrics_df, fig_path)
        print(f"Wrote figure to {fig_path}(.pdf/.png)", file=sys.stderr)
    else:
        print("No phenotype cleared MIN_CASES -- no figure written.", file=sys.stderr)


if __name__ == "__main__":
    main()
