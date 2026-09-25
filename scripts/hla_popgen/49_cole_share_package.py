#!/usr/bin/env python3
"""49 -- Cole data-share package builder (S04 WS-E).

Builds, ON THE VM, a small, tidy, documented data package for Cole under
`~/s04/share_release_2026-09/` (never under `~/pipeline_outputs*` -- this script only reads
those paths, it never writes into them). This script never uploads anything; upload is a
separate, Marc-consented step (see "UPLOAD (TEXT ONLY -- DO NOT RUN FROM HERE)" below).

Inputs (VM paths, all read-only):
  - HLA Table 1: `~/pipeline_outputs/hla_calls_rich.tsv` (`--hla-table1`)
  - HLA cis pairs: `~/pipeline_outputs/hla_cis_pairs.tsv` (`--hla-cis-pairs`, optional)
  - Cohort membership (Table 4): `~/pipeline_outputs/cohort_membership.tsv` (`--cohort-membership`)
  - KIR per-person GTFs: `~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz`
    (`--kir-outroot`), optionally `hap{1,2}/cds.fa.gz` next to them for the CDS/protein hash join
  - Relatedness pairs: AoU v9 `samples_relatedness.tsv` (`--relatedness-table`)
  - DQ/DP phasing-confidence file (37's WS1 artifact-check output, VM-local, never pulled off the
    VM into git): `~/s03/results/37/phasing_confidence_per_person.tsv` (`--phasing-confidence`).
    Optional -- if missing, the corresponding output columns are all NA and the README says so;
    this file is per-person and stays inside the perimeter (VM -> in-perimeter bucket is not
    egress, matching ORCHESTRATOR_HANDOFF.md S4 WS-E and context/ENVIRONMENT.md's "Share bucket"
    section), but it must never leave the perimeter or land in git.

Reuses, rather than re-derives:
  - `24_novelty_by_field.py`: `load_table1`, `add_call_labels` (field_class/artifact_label),
    `normalize_allele_name`/`n_field_name` (two-field truncation), `greedy_unrelated`/
    `load_relatedness_pairs`, `sha8`/`protein_info`.
  - `34_novel_protein_recurrence.py`: novelty-tier vocabulary cross-check only (no re-derivation;
    34 operates on already-clustered `novel_alleles.tsv`, which this script does not consume).
  - `41_kir_pilot.py`: `parse_hap_gtf`/`classify_novelty_tier` (KIR novelty tiers), `KIR_GENES`.
  - `44_kir_recurrence_saturation.py`: `parse_hap_gtf_full` (contig-aware KIR parsing) and
    `resolve_kir_cds_sequences` (the `cds.fa.gz` join), used verbatim for the KIR CDS/protein hash.
  - `43_kir_full_aggregate.py`: `discover_people` (KIR person discovery).
  - `37_dq_g1g2_signed_ld.py`: only as documentation of `phasing_confidence_per_person.tsv`'s
    schema (`build_phasing_confidence_table`'s docstring) -- this script reads that file's output,
    it does not regenerate it.

Disclosure: this package is participant-level, in-perimeter data (ORCHESTRATOR_HANDOFF.md S4 WS-E
"already in-perimeter controlled data, not egress"). It must NEVER be committed to git, printed in
full to stdout, or pasted into a report/chat. This script's own stdout is aggregate-only (row
counts, byte counts, column names) -- no person_id and no allele string is ever printed. A
disclosure convention is still applied to novel-allele carrier counts if/when this script is asked
to print any (it currently prints none): counts 1-19 render as the literal string "<20".

Modes:
  --dry-run   Build every table in memory, print the planned file list with estimated gzip bytes
              (measured by actually gzip-compressing each table, not guessed) and row counts
              (aggregate only), then exit without writing anything to `--out-dir`.
  --verify    Skip building; re-read `MANIFEST.tsv` already in `--out-dir` and recompute md5/sha256
              for every listed file, reporting only filename + match/mismatch (never any file's
              contents).
  (default)   Build and write all five deliverables plus MANIFEST.tsv into `--out-dir`.

Output format: gzip-compressed TSV by default (`--format tsv`, most portable -- Cole's stated
tooling is base R / simple scripts). `--format parquet` is available if `pyarrow` is present on
the VM (smaller, typed, but a heavier dependency for a collaborator who may not have it) -- picked
per-run, not silently upgraded.

===============================================================================================
UPLOAD (TEXT ONLY -- DO NOT RUN FROM HERE; a human runs these after reviewing the built package
and Marc has given explicit consent naming the exact files and sizes, per ORCHESTRATOR_HANDOFF.md
S4 WS-E and CLAUDE.md's "explicit permission required" rule for shared-state VM/bucket writes):

  # 1. Dry-run listing of what would be uploaded (no data moves):
  gsutil -m rsync -n -r ~/s04/share_release_2026-09/ \\
      gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/

  # 2. Actual upload, once Marc has confirmed the file list/sizes printed by --dry-run above:
  gsutil -m rsync -r ~/s04/share_release_2026-09/ \\
      gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/

  # 3. Verify the upload landed intact (compares local vs remote checksums; this bucket is
  #    workspace-owned and NOT requester-pays, so no -u/--billing-project flag is needed, unlike
  #    the AoU source bucket):
  gsutil -m rsync -n -r -c ~/s04/share_release_2026-09/ \\
      gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/ 2>&1 | grep -v "^Building"

  # 4. Alternative single-file verification (md5 stored by GCS, compare to MANIFEST.tsv's md5):
  gcloud storage objects describe \\
      gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/MANIFEST.tsv \\
      --format="value(md5Hash)"
===============================================================================================

Usage (on the VM):
    python3 scripts/hla_popgen/49_cole_share_package.py --dry-run
    python3 scripts/hla_popgen/49_cole_share_package.py
    python3 scripts/hla_popgen/49_cole_share_package.py --verify

Tests (local, synthetic data, no VM): scripts/hla_popgen/tests/test_49_cole_share_package.py
"""
import argparse
import gzip
import hashlib
import importlib.util
import io
import os
import sys
import time

import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _viz_common as vc  # noqa: E402

VERSION = "2026-09"
GENERATED_BY = "scripts/hla_popgen/49_cole_share_package.py"
DQ_DP_PAIRS = (("DQA1", "DQB1"), ("DPA1", "DPB1"))
HLA_NOVELTY_LABELS = {
    "known": "known",
    "f4_noncoding": "non_coding",
    "f3_synonymous": "synonymous",
    "f2_protein": "protein",
    "f1_undetermined": "undetermined",
    "uncalled": "uncalled",
}


def log(msg):
    print(f"[49 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Lazy module loading (mirrors 43/44's own `_load_module`/`m41()`/`m03()` pattern -- each script
# stays independently importable/testable, no package restructuring needed).
# ---------------------------------------------------------------------------
def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_m24 = _m41 = _m43 = _m44 = None


def m24():
    global _m24
    if _m24 is None:
        _m24 = _load_module("24_novelty_by_field.py", "novelty_by_field_for_49")
    return _m24


def m41():
    global _m41
    if _m41 is None:
        _m41 = _load_module("41_kir_pilot.py", "kir_pilot_for_49")
    return _m41


def m43():
    global _m43
    if _m43 is None:
        _m43 = _load_module("43_kir_full_aggregate.py", "kir_full_aggregate_for_49")
    return _m43


def m44():
    global _m44
    if _m44 is None:
        mod = _load_module("44_kir_recurrence_saturation.py", "kir_recurrence_saturation_for_49")
        _m44 = mod
    return _m44


# ---------------------------------------------------------------------------
# Disclosure helpers (same convention as 43/44: counts 1-19 -> "<20", true 0 stays "0").
# ---------------------------------------------------------------------------
def suppressed(n):
    n = int(n)
    if 1 <= n < 20:
        return "<20"
    return str(n)


# ---------------------------------------------------------------------------
# 1. hla_calls.tsv.gz
# ---------------------------------------------------------------------------
def two_field_from_consensus(consensus, gene_bare, m24mod):
    """'HLA-A*01:01:01:01' -> 'A*01:01'. None for uncalled/undetermined/unparseable consensus."""
    if not isinstance(consensus, str) or consensus.strip().upper() in m24mod.vc.NULL_TOKENS:
        return None
    gene, full = m24mod.normalize_allele_name(consensus, fallback_gene_bare=gene_bare)
    if not full:
        return None
    two = m24mod.n_field_name(full, 2)
    return two or None


def load_phasing_confidence(path, pairs=DQ_DP_PAIRS):
    """Reads 37's per-person phasing-confidence file (VM-local, `~/s03/results/37/
    phasing_confidence_per_person.tsv`, schema per `build_phasing_confidence_table`'s docstring
    in 37_dq_g1g2_signed_ld.py) and returns a dict {(person_id, hap, gene): "PAIR:status"},
    restricted to the DQ/DP gene_a/gene_b pairs this package cares about. Returns {} (not a
    FATAL) if the file is absent -- it is explicitly optional ("where available") per the task,
    and 37's own README documents it as VM-local-only, so a fresh VM session without a completed
    37 artifact-check run must still be able to build the rest of the package."""
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"person_id", "hap", "gene_a", "gene_b", "status"}
    missing = required - set(df.columns)
    if missing:
        log(f"WARNING: phasing-confidence file lacks {missing}; skipping (no crash).")
        return {}
    wanted_pairs = {(a, b) for a, b in pairs}
    out = {}
    for row in df.itertuples(index=False):
        if (row.gene_a, row.gene_b) not in wanted_pairs:
            continue
        tag = f"{row.gene_a}~{row.gene_b}:{row.status}"
        for gene in (row.gene_a, row.gene_b):
            key = (row.person_id, row.hap, gene)
            out[key] = out.get(key, "") + (";" if key in out else "") + tag
    return out


def build_hla_calls(table1_path, phasing_confidence_path, limit=None):
    """One row per (person_id, hap, contig, gene, copy_index) -- Table 1's own grain, per
    SCHEMA.md. Returns (df, stats: dict aggregate-only)."""
    m24mod = m24()
    log("loading HLA Table 1 ...")
    t1 = m24mod.load_table1(table1_path, limit=limit)
    t1 = m24mod.add_call_labels(t1)

    gene_bare = t1["gene"].astype(str).str.replace("^HLA-", "", regex=True)
    two_field = [two_field_from_consensus(c, g, m24mod)
                 for c, g in zip(t1["consensus"].tolist(), gene_bare.tolist())]

    phasing = load_phasing_confidence(phasing_confidence_path)
    dq_dp_status = [phasing.get((pid, hap, g)) for pid, hap, g in
                     zip(t1["person_id"], t1["hap"], gene_bare)]

    novelty_tier = t1["field_class"].map(HLA_NOVELTY_LABELS).fillna("uncalled")
    already_shared = (novelty_tier == "known") & pd.Series(two_field).notna().values

    out = pd.DataFrame({
        "person_id": t1["person_id"].values,
        "hap": t1["hap"].values,
        "contig": t1["contig"].values,
        "gene": t1["gene"].values,
        "gene_class": t1["gene_class"].values,
        "copy_index": t1["copy_index"].astype("Int64").values,
        "allele": t1["consensus"].values,
        "two_field": two_field,
        "novelty_tier": novelty_tier.values,
        "artifact_label": t1["artifact_label"].values,
        "dq_dp_phasing_confidence": dq_dp_status,
        "already_shared_two_field_common": already_shared.values,
    })
    stats = {
        "n_rows": len(out),
        "n_persons": out["person_id"].nunique(),
        "n_novel": int((out["novelty_tier"] != "known").sum()
                       - (out["novelty_tier"] == "uncalled").sum()),
        "phasing_confidence_available": bool(phasing),
    }
    return out, stats


# ---------------------------------------------------------------------------
# 2. kir_calls.tsv.gz
# ---------------------------------------------------------------------------
def build_kir_calls(kir_outroot, limit=None):
    """One row per (person_id, hap, contig, gene, copy_index). Reuses 44's contig-aware GTF
    parser and cds.fa.gz join verbatim; falls back gracefully (cds_sha8/protein_sha8 = NA,
    cds_available = False) when cds.fa.gz is absent for a haplotype -- matches 44's own
    documented "count, never guess" contract for that join.

    Returns (df, stats: dict aggregate-only, kir_call_status: {person_id: bool} -- whether the
    person had >=1 valid, parseable hap with >=1 KIR gene call, used by build_persons() for
    kir_called / exclusion_reason without a second GTF read pass).
    """
    kir_mod = m41()
    m44mod = m44()
    m03mod = _load_module("03_novel_alleles.py", "novel_alleles_for_49")
    m24mod = m24()

    pids = m43().discover_people(kir_outroot, limit=limit)
    log(f"KIR: discovered {len(pids)} person directories under {kir_outroot!r}")

    rows = []
    kir_call_status = {}
    n_hap_seen = n_hap_cds_present = 0
    for pid in pids:
        person_had_call = False
        for hap in ("hap1", "hap2"):
            gtf_path = os.path.join(kir_outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
            hap_rows = m44mod.parse_hap_gtf_full(gtf_path, kir_mod)
            kir_rows = [r for r in hap_rows if r["gene"] in kir_mod.KIR_GENES]
            if not kir_rows:
                continue
            n_hap_seen += 1
            person_had_call = True
            cds_fa_path = os.path.join(kir_outroot, str(pid), "immuannot_output", hap,
                                        "cds.fa.gz")
            cds_available = os.path.exists(cds_fa_path)
            if cds_available:
                n_hap_cds_present += 1
                seq_by_key, _stats = m44mod.resolve_kir_cds_sequences(kir_rows, cds_fa_path,
                                                                       m03mod)
            else:
                seq_by_key = {}
            for r in kir_rows:
                # copy_index comes from parse_hap_gtf_full's own extraction (the gene_id ".N"
                # suffix, per reference/IMMUANNOT_GTF_SPEC.md) -- the correct semantic (a real
                # second mapping cluster), not an occurrence-order counter.
                seq = seq_by_key.get((r["gene"], r["consensus"]))
                cds_sha8 = protein_sha8 = None
                if seq:
                    cds_sha8 = m24mod.sha8(seq)
                    protein_sha8 = m24mod.sha8(m24mod.protein_info(seq)["protein"])
                rows.append({
                    "person_id": str(pid),
                    "hap": hap,
                    "contig": r["contig"],
                    "gene": r["gene"],
                    "copy_index": r["copy_index"],
                    "allele": r["consensus"],
                    "novelty_tier": r["novelty_tier"],
                    "cds_sha8": cds_sha8,
                    "protein_sha8": protein_sha8,
                    "cds_available": cds_available,
                })
        kir_call_status[str(pid)] = person_had_call

    df = pd.DataFrame(rows, columns=[
        "person_id", "hap", "contig", "gene", "copy_index", "allele", "novelty_tier",
        "cds_sha8", "protein_sha8", "cds_available"])
    stats = {
        "n_rows": len(df),
        "n_persons_attempted": len(pids),
        "n_persons_with_calls": sum(kir_call_status.values()),
        "n_hap_with_calls": n_hap_seen,
        "n_hap_cds_available": n_hap_cds_present,
        "cds_join_available": n_hap_cds_present > 0,
    }
    return df, stats, kir_call_status


# ---------------------------------------------------------------------------
# 3. persons.tsv
# ---------------------------------------------------------------------------
def build_persons(cohort_membership_path, relatedness_path, kir_call_status, hla_person_ids,
                   kin_min=0.0442, skip_relatedness=False):
    """One row per person in cohort_membership (the long-read/HLA-called cohort). `kir_call_status`
    and `hla_person_ids` come from the already-built hla_calls/kir_calls tables (no re-parsing).
    """
    m24mod = m24()
    cohort = vc.load_cohort_membership(cohort_membership_path)
    pcols = [f"p_{a.lower()}" for a in vc.ANCESTRY_ORDER]
    keep_cols = ["person_id", "ancestry_pred", "platform"] + [
        c for c in pcols if c in cohort.columns]
    people = cohort[[c for c in keep_cols if c in cohort.columns]].copy()
    people["person_id"] = people["person_id"].astype(str)

    if skip_relatedness:
        people["unrelated"] = True
    else:
        pairs = m24mod.load_relatedness_pairs(relatedness_path)
        kept, _removed = m24mod.greedy_unrelated(set(people["person_id"]), pairs, kin_min)
        people["unrelated"] = people["person_id"].isin(kept)

    people["hla_called"] = people["person_id"].isin(set(hla_person_ids))
    people["kir_called"] = people["person_id"].map(
        lambda p: kir_call_status.get(p, False))

    def exclusion_reason(row):
        if row["kir_called"] and row["hla_called"]:
            return ""
        reasons = []
        if not row["hla_called"]:
            reasons.append("hla_not_called")
        if not row["kir_called"]:
            plat = str(row.get("platform", "")).strip().lower()
            if plat == "sequel2":
                reasons.append("kir_excluded_sequel2_tier3_self_align")
            elif row["person_id"] not in kir_call_status:
                reasons.append("kir_dir_missing")
            else:
                reasons.append("kir_parse_empty")
        return ";".join(reasons)

    people["exclusion_reason"] = people.apply(exclusion_reason, axis=1)
    return people


# ---------------------------------------------------------------------------
# Writing + manifest
# ---------------------------------------------------------------------------
def gzip_bytes(df, fmt):
    """Serializes df to an in-memory gzip'd TSV (or parquet) buffer -- used both to write the
    real file and, in --dry-run, to MEASURE a real compressed size instead of guessing."""
    buf = io.BytesIO()
    if fmt == "parquet":
        df.to_parquet(buf, index=False, compression="gzip")
    else:
        with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
            df.to_csv(gz, sep="\t", index=False)
    return buf.getvalue()


def write_table(df, out_dir, stem, fmt):
    ext = "parquet" if fmt == "parquet" else "tsv.gz"
    path = os.path.join(out_dir, f"{stem}.{ext}")
    data = gzip_bytes(df, fmt)
    with open(path, "wb") as f:
        f.write(data)
    return path, len(data)


def hashes_of(path):
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            md5.update(chunk)
            sha256.update(chunk)
    return md5.hexdigest(), sha256.hexdigest()


def write_manifest(out_dir, files_meta):
    """files_meta: list of (relative_path, n_rows). Hashes are computed from the already-written
    files on disk. Never writes outside out_dir."""
    rows = []
    for rel_path, n_rows in files_meta:
        full = os.path.join(out_dir, rel_path)
        md5, sha256 = hashes_of(full)
        rows.append({
            "file": rel_path,
            "bytes": os.path.getsize(full),
            "rows": n_rows,
            "md5": md5,
            "sha256": sha256,
        })
    manifest = pd.DataFrame(rows, columns=["file", "bytes", "rows", "md5", "sha256"])
    manifest_path = os.path.join(out_dir, "MANIFEST.tsv")
    manifest.to_csv(manifest_path, sep="\t", index=False)
    return manifest_path, manifest


def verify_manifest(out_dir):
    manifest_path = os.path.join(out_dir, "MANIFEST.tsv")
    if not os.path.exists(manifest_path):
        sys.exit(f"FATAL: no MANIFEST.tsv found at {manifest_path!r}. Build the package first.")
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    ok = True
    for row in manifest.itertuples(index=False):
        full = os.path.join(out_dir, row.file)
        if not os.path.exists(full):
            print(f"MISSING  {row.file}")
            ok = False
            continue
        md5, sha256 = hashes_of(full)
        good = (md5 == row.md5) and (sha256 == row.sha256)
        print(f"{'OK      ' if good else 'MISMATCH'} {row.file}")
        ok = ok and good
    print("VERIFY: all files match manifest" if ok else "VERIFY: one or more mismatches found")
    return ok


# ---------------------------------------------------------------------------
# README.md / SCHEMA.md
# ---------------------------------------------------------------------------
def write_docs(out_dir, hla_stats, kir_stats, persons_df, fmt, date_str):
    ext = "parquet" if fmt == "parquet" else "tsv.gz"
    n_persons = len(persons_df)
    n_unrelated = int(persons_df["unrelated"].sum())
    n_hla_called = int(persons_df["hla_called"].sum())
    n_kir_called = int(persons_df["kir_called"].sum())
    n_kir_excluded_sequel2 = int(
        persons_df["exclusion_reason"].astype(str).str.contains(
            "kir_excluded_sequel2_tier3_self_align").sum())

    readme = f"""# Omni-HLA data package for Cole -- release {date_str}

Generated by `{GENERATED_BY}`, version `{VERSION}`.

**This is participant-level, in-perimeter controlled-tier data.** It stays inside the AoU
Verily Workbench VPC-SC perimeter (VM -> this workspace's own share bucket is not egress). It must
never be copied into git, a public report, a chat message, or anywhere outside the perimeter.

## What's new in this release vs. what Cole already has

Cole previously received (and confirmed using) **only common, two-field HLA alleles** -- see
`sprints/ALEIX_HANDOFF_2026-09-22.md` S3. **Novel HLA calls were never shared.** This release adds:
- novel HLA calls (all novelty tiers, not just the common two-field subset), with novelty tier and
  artifact label so they can be told apart from the previously-shared common calls
  (`already_shared_two_field_common` column in `hla_calls.{ext}`);
- KIR calls (chr19), never shared before now;
- ancestry probabilities/labels, sequencing platform, and unrelated-set membership
  (`persons.tsv`), never shared before now;
- the DQ/DP cis-pair phasing-confidence flag Cole asked for on the 2026-09-22 call
  (`dq_dp_phasing_confidence` column in `hla_calls.{ext}`), where available.

## Files

### `hla_calls.{ext}`
One row per (person_id, hap, contig, gene, copy_index) -- same grain as the internal pipeline's
Table 1 (`scripts/hla_popgen/SCHEMA.md`). {hla_stats['n_rows']} rows, {hla_stats['n_persons']}
people.

| Column | Meaning |
|---|---|
| `person_id` | AoU research id (in-perimeter controlled data, not a disclosive identifier once inside the perimeter) |
| `hap` | `hap1`/`hap2` -- which assembled haplotype |
| `contig` | assembly contig -- the cis key; two genes are only in physical cis if they share this |
| `gene` | HLA gene symbol |
| `gene_class` | classical_I / classical_II / nonclassical_I / pseudogene_I / etc. (SCHEMA.md) |
| `copy_index` | 1 unless a real second mapping cluster (segmental duplication / DRB copy number) |
| `allele` | full-resolution Immuannot consensus call (may carry a trailing novelty field) |
| `two_field` | two-field (protein-resolution) truncation, e.g. `A*01:01`; blank if uncalled/undetermined |
| `novelty_tier` | `known` / `non_coding` / `synonymous` / `protein` / `undetermined` / `uncalled` -- see caveat below |
| `artifact_label` | `clean` / `partial_cds` / `inframe_stop` / `homopolymer_indel` -- QC flag from the same logic as `24_novelty_by_field.py`/`34_novel_protein_recurrence.py` |
| `dq_dp_phasing_confidence` | semicolon-joined `PAIR:status` for DQA1~DQB1 / DPA1~DPB1 (status in physical/different_contig/missing_a/missing_b/missing_both); blank for other genes or when the phasing-confidence run wasn't available (see caveat) |
| `already_shared_two_field_common` | `True` if this exact two-field allele is part of what Cole already has from the earlier common-allele share |

### `kir_calls.{ext}`
One row per (person_id, hap, contig, gene, copy_index). {kir_stats['n_rows']} rows,
{kir_stats['n_persons_with_calls']} people with >=1 KIR call
({kir_stats['n_persons_attempted']} attempted).

| Column | Meaning |
|---|---|
| `person_id`, `hap`, `contig`, `gene`, `copy_index` | as above, KIR region (chr19) |
| `allele` | Immuannot consensus / closest IPD-KIR reference allele |
| `novelty_tier` | `known` / `novel_genomic_known_cds` / `novel_cds_synonymous` / `novel_protein` / `novel_unclassified` / `undetermined` (see `41_kir_pilot.classify_novelty_tier` docstring for the exact rule -- deliberately not force-mapped onto HLA's four-tier names, since the underlying evidence differs) |
| `cds_sha8` / `protein_sha8` | 8-char sha1 of the joined CDS / protein sequence (from `<hap>/cds.fa.gz`, when present) -- lets Cole tell "same underlying sequence" apart from "same consensus string" without us sharing raw sequence |
| `cds_available` | whether the `cds.fa.gz` join succeeded for this haplotype (cohort-wide: {"available" if kir_stats['cds_join_available'] else "NOT available -- cds_sha8/protein_sha8 are all NA this release"}) |

**"Use novel alleles with caution."** Per Cole's own words on the 2026-09-22 call: *"the novel
ones are interesting, but they're kind of dangerous, and we can't say much with them because we
don't have a lot of them"* (`sprints/ALEIX_HANDOFF_2026-09-22.md`). Any `novelty_tier` other than
`known` in either file is a single-cohort, unvalidated call -- do not use as a positive control or
as a stratifier without checking recurrence across unrelated people first
(`reports/hla_popgen/34_novel_protein_recurrence/`, `44`/`45`).

### `persons.tsv`
One row per person in the long-read (Immuannot) cohort. {n_persons} people; {n_unrelated} flagged
unrelated (KING kinship < 0.0442, third-degree-or-closer removed greedily); {n_hla_called} with
HLA calls; {n_kir_called} with >=1 KIR call ({n_kir_excluded_sequel2} of the excluded people are
the known sequel2/Tier-3 self-align gap, see caveats).

| Column | Meaning |
|---|---|
| `person_id` | AoU research id |
| `ancestry_pred` | predicted ancestry label (AFR/AMR/EAS/EUR/MID/SAS) |
| `p_afr`...`p_sas` | continuous admixture proportions |
| `platform` | sequencing platform (`revio` / `sequel2e` / `sequel2`) |
| `unrelated` | `True` if kept by the greedy third-degree-or-closer relatedness filter |
| `hla_called` | `True` if this person has >=1 row in `hla_calls.{ext}` |
| `kir_called` | `True` if this person has >=1 row in `kir_calls.{ext}` |
| `exclusion_reason` | semicolon-joined reason(s) when `hla_called`/`kir_called` is False (`hla_not_called`, `kir_excluded_sequel2_tier3_self_align`, `kir_dir_missing`, `kir_parse_empty`); blank when both are True |

### `README.md` / `SCHEMA.md`
This file, and a companion schema reference (versioned, dated `{date_str}`).

### `MANIFEST.tsv`
`file`, `bytes`, `rows`, `md5`, `sha256` for every deliverable above -- for Cole to verify the
copy he receives matches what was built.

## Provenance

- Pipeline: this repo's `hla_popgen` sub-project (Immuannot on AoU long-read WGS, v9 CDR).
- Reference databases: IPD-IMGT/HLA 3.55.0, IPD-KIR (version as pinned in
  `~/tools/Immuannot_refdata` at run time -- record the exact IPD-KIR version string here before
  upload if it differs from what a prior release used).
- Ancestry: AoU v9 `ancestry_preds.tsv` (six-population predicted ancestry + continuous
  probabilities).
- Relatedness: AoU v9 `samples_relatedness.tsv`, KING kinship, third-degree-or-closer (>=0.0442)
  removed greedily to build the `unrelated` flag.

## Caveats (read before analysis)

1. **Novel alleles are "interesting but dangerous"** (Cole, verbatim, 2026-09-22) -- single-cohort,
   not independently validated, often too rare to draw conclusions from. Treat `novelty_tier !=
   known` as provisional throughout both tables.
2. **KIR excludes 991 sequel2 people.** Tier-3/`self_align_needed` people on the `sequel2`
   platform were excluded upstream of the KIR production run by a known chr6-reference-cache bug
   in `run_immuannot_person.py` (`ensure_chr6_ref()`), not fixed as of this release. See
   `reports/hla_popgen/43_kir_full_cohort/README.md`. These people appear in `persons.tsv` with
   `kir_called=False`, `exclusion_reason` containing `kir_excluded_sequel2_tier3_self_align`.
3. **KIR2DL2/KIR2DL3 co-occurrence (~0.8% of haplotypes; biological expectation ~0)** is an open
   question -- miscall vs. real duplication vs. phase error, not yet resolved
   (`reports/hla_popgen/43_kir_full_cohort/README.md`, S4 open item). Do not treat carriers of
   both as confirmed duplications.
4. **HLA-A deletion.** Cole called this "a big result if real" but it is not yet independently
   confirmed; see `reports/hla_popgen/38_hla_a_deletion_validation/README.md` for the current
   validation status before using HLA-A absence as a genotype call rather than a QC flag.
5. **Relatedness.** `unrelated=True` uses a single genome-wide kinship threshold (0.0442, KING
   third-degree). It does not by itself guarantee independence for every possible downstream
   analysis (e.g. very close ancestry-matched pairs below this exact threshold are still possible)
   -- treat it as the pipeline's standard cohort filter, not a guarantee.
6. **`dq_dp_phasing_confidence` availability**: {"populated from `phasing_confidence_per_person.tsv` (37's WS1 artifact-check output)." if hla_stats['phasing_confidence_available'] else "**NOT populated this release** -- 37's per-person phasing-confidence file was not found on this VM at build time (it is produced by `37_dq_g1g2_signed_ld.py --artifact-checks` and never leaves the VM as a per-person file). Every row's `dq_dp_phasing_confidence` is blank. Re-run 37's artifact checks and rebuild this package to populate it."}
7. **Small-cell disclosure.** No small-cell counts (participant counts, not allele-cluster counts)
   appear in any file here -- these are per-person genotype rows, which is the standard grain for
   an in-perimeter collaborator share, not an aggregate publication. The still-open small-cell
   *publication* policy in `context/DECISIONS.md` applies only to material leaving the perimeter
   (reports, git, papers), not to this in-perimeter share.

## Upload

Draft `gsutil`/`gcloud` commands (text only, not run by this script or by any agent) are in this
script's own module docstring (`scripts/hla_popgen/49_cole_share_package.py`). Upload requires
Marc's explicit, in-chat consent naming the exact files and sizes -- see
`sprints/S04_kir_recurrence_style_share/ORCHESTRATOR_HANDOFF.md` S2 WS-E.
"""
    with open(os.path.join(out_dir, "README.md"), "w") as f:
        f.write(readme)

    schema = f"""# SCHEMA -- Cole data package, release {date_str}

Companion to `README.md` in this same directory. Column-by-column type/nullability reference.
See `README.md` for narrative meaning; this file is the quick-lookup contract.

## hla_calls.{ext}

| Column | Type | Nullable | Justification for inclusion |
|---|---|---|---|
| person_id | string | no | join key across all three tables |
| hap | string (`hap1`/`hap2`) | no | cis-phasing unit |
| contig | string | no | the only valid cis key (SCHEMA.md Table 1) |
| gene | string | no | which gene this row calls |
| gene_class | string (category) | no | lets Cole filter to classical genes without a lookup table of his own |
| copy_index | int | no | disambiguates segmental-duplication/DRB copy-number rows |
| allele | string | no | the actual typing call -- the point of the table |
| two_field | string | yes (uncalled/undetermined) | Cole's own prior use was two-field only; keeps that resolution available without redoing the truncation |
| novelty_tier | string (category) | no (uncalled is a valid value) | Cole's headline distinction: known vs. novel, and how novel |
| artifact_label | string (category) | no | lets Cole exclude QC-flagged rows from novel-allele analyses without re-deriving `template_warning` parsing himself |
| dq_dp_phasing_confidence | string | yes | directly answers Cole's 2026-09-22 request for a per-person DQ/DP phasing flag |
| already_shared_two_field_common | bool | no | tells Cole which rows duplicate what he already has vs. what's new |

Columns deliberately NOT included: `template_distance`/`cds_mut`/`n_aa_changes`/raw QC internals
-- these are pipeline-internal QC fields, not needed for Cole's stat-gen use, and every one of
them is one column that must be documented, kept schema-stable, and re-checked for disclosure
risk on every future release. `novelty_tier`+`artifact_label` already summarize what he needs.

## kir_calls.{ext}

| Column | Type | Nullable | Justification for inclusion |
|---|---|---|---|
| person_id, hap, contig, gene, copy_index | as above | no | same cis/grain contract as hla_calls |
| allele | string | no | the typing call |
| novelty_tier | string (category) | no | KIR's own five/six-tier novelty vocabulary (41_kir_pilot.classify_novelty_tier) -- kept distinct from HLA's because the underlying evidence (cds_distance/cds_mut from Immuannot's own KIR-CDS comparison) differs from HLA's RefIndex-based classification |
| cds_sha8 | string(8) | yes | sequence identity without sharing raw nucleotide sequence (same principle as SCHEMA.md Table 3's `cds_seq_sha1`) |
| protein_sha8 | string(8) | yes | same, at the protein level -- lets Cole group by protein identity across differing CDS/nucleotide calls |
| cds_available | bool | no | whether the cds.fa.gz join succeeded for this haplotype; False everywhere means cds_sha8/protein_sha8 are uninformative this release, not that sequences don't exist |

## persons.tsv

| Column | Type | Nullable | Justification for inclusion |
|---|---|---|---|
| person_id | string | no | join key |
| ancestry_pred | string (category) | yes (UNASSIGNED-style blanks possible upstream) | ancestry stratification, explicitly requested |
| p_afr...p_sas | float [0,1] | no | continuous admixture -- needed for the stricter ancestry-probability thresholds Cole has asked for elsewhere (e.g. 95-98%) |
| platform | string (category) | yes | needed to check for platform-driven batch effects, and to interpret `kir_excluded_sequel2_tier3_self_align` |
| unrelated | bool | no | avoids Cole re-deriving the same relatedness filter the rest of this pipeline uses |
| hla_called, kir_called | bool | no | lets Cole build denominators correctly without joining back to the per-row tables |
| exclusion_reason | string | no (blank = none) | makes non-calls legible instead of silently-missing rows |

## MANIFEST.tsv

| Column | Type | Justification |
|---|---|---|
| file | string | which deliverable |
| bytes | int | on-disk size, for a quick copy sanity check |
| rows | int | row count, for a quick copy sanity check |
| md5 | string(32 hex) | integrity check matching `gsutil`/`gcloud` reported hashes |
| sha256 | string(64 hex) | stronger integrity check for long-term archival |
"""
    with open(os.path.join(out_dir, "SCHEMA.md"), "w") as f:
        f.write(schema)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s04/share_release_2026-09"))
    ap.add_argument("--hla-table1", default=os.path.expanduser(
        "~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--hla-cis-pairs", default=os.path.expanduser(
        "~/pipeline_outputs/hla_cis_pairs.tsv"))
    ap.add_argument("--cohort-membership", default=os.path.expanduser(
        "~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--kir-outroot", default=os.path.expanduser("~/pipeline_outputs_kir"))
    ap.add_argument("--relatedness-table", default=os.path.expanduser(
        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
        "samples_relatedness.tsv"))
    ap.add_argument("--phasing-confidence", default=os.path.expanduser(
        "~/s03/results/37/phasing_confidence_per_person.tsv"))
    ap.add_argument("--format", choices=["tsv", "parquet"], default="tsv")
    ap.add_argument("--kin-min", type=float, default=0.0442)
    ap.add_argument("--skip-relatedness", action="store_true",
                     help="smoke-test only: skip the relatedness join")
    ap.add_argument("--limit", type=int, default=None,
                     help="smoke test: cap discovered people on each side")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verify", action="store_true")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.verify:
        ok = verify_manifest(args.out_dir)
        sys.exit(0 if ok else 1)

    forbidden = [os.path.expanduser("~/pipeline_outputs"), os.path.expanduser("~/pipeline_outputs_kir")]
    out_abs = os.path.abspath(args.out_dir)
    if any(out_abs == p or out_abs.startswith(p + os.sep) for p in forbidden):
        sys.exit(f"FATAL: --out-dir {args.out_dir!r} is under a production output path "
                 f"({forbidden}) -- refusing to write there. Use ~/s04/... instead.")

    date_str = time.strftime("%Y-%m-%d")

    hla_df, hla_stats = build_hla_calls(args.hla_table1, args.phasing_confidence,
                                         limit=args.limit)
    kir_df, kir_stats, kir_call_status = build_kir_calls(args.kir_outroot, limit=args.limit)
    persons_df = build_persons(args.cohort_membership, args.relatedness_table, kir_call_status,
                                set(hla_df["person_id"].unique()), kin_min=args.kin_min,
                                skip_relatedness=args.skip_relatedness)

    tables = [
        ("hla_calls", hla_df),
        ("kir_calls", kir_df),
    ]

    if args.dry_run:
        print("=== DRY RUN: planned files (aggregate-only; no IDs, no files written) ===")
        total_bytes = 0
        for stem, df in tables:
            data = gzip_bytes(df, args.format)
            ext = "parquet" if args.format == "parquet" else "tsv.gz"
            total_bytes += len(data)
            print(f"  {stem}.{ext}: {len(df)} rows, {len(df.columns)} columns, "
                  f"~{len(data)/1e6:.2f} MB compressed")
        persons_bytes = len(persons_df.to_csv(sep="\t", index=False).encode())
        total_bytes += persons_bytes
        print(f"  persons.tsv: {len(persons_df)} rows, {len(persons_df.columns)} columns, "
              f"~{persons_bytes/1e6:.2f} MB uncompressed (not gzip'd)")
        print("  README.md, SCHEMA.md: text, negligible size")
        print("  MANIFEST.tsv: negligible size")
        print(f"  --- total estimated package size: ~{total_bytes/1e6:.2f} MB ---")
        print(f"  HLA: {hla_stats['n_persons']} people, {hla_stats['n_rows']} rows, "
              f"phasing_confidence_available={hla_stats['phasing_confidence_available']}")
        print(f"  KIR: {kir_stats['n_persons_with_calls']}/{kir_stats['n_persons_attempted']} "
              f"people with calls, {kir_stats['n_rows']} rows, "
              f"cds_join_available={kir_stats['cds_join_available']}")
        print(f"  Persons: {len(persons_df)} total, "
              f"{int(persons_df['unrelated'].sum())} unrelated")
        return

    os.makedirs(args.out_dir, exist_ok=True)
    files_meta = []
    for stem, df in tables:
        path, nbytes = write_table(df, args.out_dir, stem, args.format)
        files_meta.append((os.path.basename(path), len(df)))
        log(f"wrote {path} ({nbytes} bytes, {len(df)} rows)")

    persons_path = os.path.join(args.out_dir, "persons.tsv")
    persons_df.to_csv(persons_path, sep="\t", index=False)
    files_meta.append(("persons.tsv", len(persons_df)))
    log(f"wrote {persons_path} ({os.path.getsize(persons_path)} bytes, {len(persons_df)} rows)")

    write_docs(args.out_dir, hla_stats, kir_stats, persons_df, args.format, date_str)
    files_meta.append(("README.md", None))
    files_meta.append(("SCHEMA.md", None))

    manifest_path, manifest = write_manifest(args.out_dir, files_meta)
    log(f"wrote {manifest_path}")
    print("=== BUILD COMPLETE (aggregate-only summary) ===")
    print(manifest.drop(columns=["md5", "sha256"]).to_string(index=False))


if __name__ == "__main__":
    main()
