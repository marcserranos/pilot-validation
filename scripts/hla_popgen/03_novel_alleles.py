#!/usr/bin/env python3
"""Produce Table 3 (`novel_alleles.tsv`) per scripts/hla_popgen/SCHEMA.md, plus a VM-only
`novel_alleles_seqs.fa`, by matching Table 1's novel calls (01_extract_rich.py) to the per-haplotype
observed CDS sequences already sitting on disk in `<outroot>/<person>/immuannot_output/hap{1,2}/
cds.fa.gz` (reference/IMMUANNOT_GTF_SPEC.md part D/E).

## Methodology (see scripts/hla_popgen/research/NOVEL_LIT.md for the full literature review --
## every QC/clustering choice below cites a specific section of it)

- **Clustering key = exact observed CDS sequence identity** (NOVEL_LIT.md section 2.3, "identical-
  sequence collapse ... is your primary clustering key"). `novel_id = <gene>_nov_<sha1(seq)[:8]>` --
  deterministic, recomputable, never a shifting counter (SCHEMA.md Table 3).
- **Recurrence across >=2 unrelated persons is the primary evidence of real biology, not assembly
  artifact** (NOVEL_LIT.md section 1.4 item 5, section 3.4's DRB/validation discussion, and Zhou
  et al. 2024's own cross-validation logic). `passes_qc` hard-requires `n_persons >= 2`.
- **Homopolymer-indel QC filter** reproduces the Huijse et al. 2023 pan-MHC graph paper's exact
  filter class (NOVEL_LIT.md section 1.3: 271/2580 raw candidates discarded because the only
  divergence was an indel inside a homopolymer run) -- see `is_homopolymer_indel_only()`'s docstring
  for the specific way this script approximates that check given what's actually recoverable from
  `cds_mut`'s minimap2 `cs`-string field (no full-contig sequence access at this stage -- a
  documented, deliberate deviation, see module-level DEVIATIONS section below).
- **Synonymous vs non-synonymous excess as a sanity signal**: NOVEL_LIT.md section 0 (Zhou et al.
  2024's own validation argument) and section 1.4 item 6 -- a novel set dominated by missense/
  nonsynonymous change (relative to a null expectation) is itself evidence of real biology under
  balancing selection, not sequencing noise. This script reports the synonymous:non-synonymous
  breakdown (from `novelty_class`/`n_aa_changes`) in its markdown report; it does not attempt a
  formal statistical test of "excess" (04_allele_saturation.py's remit is the population-genetic
  modeling; this script's remit is the candidate table itself).

## DEVIATIONS from NOVEL_LIT.md's letter (both flagged in NOVEL_LIT.md itself as adaptable)

1. NOVEL_LIT.md section 1.4 item 3 recommends checking "local sequence context for a >=7-10bp
   homopolymer run AT THE INDEL POSITION" in the full assembled contig. This script does not have
   contig-level sequence access (only `cds_mut`'s summary cs-string, not the raw alignment) at the
   point novel calls are being tabulated -- re-reading `hap{N}.trimmed.fa` by GTF coordinates for
   every novel call would be a second, much heavier disk pass, and reference/IMMUANNOT_GTF_SPEC.md
   part D explicitly flags that route as unconfirmed to agree byte-for-byte with `cds.fa.gz` in this
   codebase's own testing history. Instead: a novel call's cds-level diff is flagged
   `is_homopolymer_indel_only` when EVERY diff operation in its `cds_mut` cs-string is an indel
   (`+seq`/`-seq`, no `*` substitution operations) AND every such indel's inserted/deleted sequence
   is itself a homopolymer run (all one repeated base). This is a same-spirit, strictly weaker proxy
   (it can't see homopolymer context outside the indel itself) -- documented here so a future pass
   with contig access can tighten it, not silently assumed equivalent.
2. NOVEL_LIT.md section 2.2 recommends collapsing via official G-group/P-group logic (exon-2/2+3
   ARD-only identity) as the primary clustering axis over raw full-CDS identity, specifically to
   avoid overcounting alleles that differ only in a 3rd/4th-field (non-ARD) region. This script uses
   full observed-CDS identity instead (SCHEMA.md Table 3's own explicit instruction: "Grain: one row
   per distinct novel sequence cluster (identical observed CDS collapsed)") -- SCHEMA.md is the
   binding contract for this repo and is followed over the literature review's own secondary
   preference; the deviation is noted here because full-CDS clustering is a finer partition than
   G/P-group clustering would be (it can only ever split a G-group's members further, never merge
   two distinct G-groups), so Table 3's counts are, if anything, a slight overcount of "true" novel
   protein-level diversity relative to the G/P-group convention -- worth flagging in any downstream
   paper-writing pass.

## Ambiguous copy_index>1 matching (reference/IMMUANNOT_GTF_SPEC.md part D "Caveat")

`cds.fa.gz` headers are `{contig}_{gene}_{i}` where `i` is a DETECTION-order index across ALL genes
on that contig -- not necessarily `copy_index`. For a gene with exactly one row in Table 1 at
`(person_id, hap, contig, gene)`, matching by `(contig, gene)` prefix alone is unambiguous (only one
candidate can exist). For a gene with >1 row at that same key (copy_index>1 -- real segmental
duplication, SCHEMA.md), the join key that would disambiguate which numbered `cds.fa.gz` record
belongs to which copy (`tmp.gene.csv`) is deleted by Immuannot's own cleanup -- this script does NOT
guess. Those rows are counted separately (`n_ambiguous_copy`) and excluded from Table 3 entirely.

## Compliance (SCHEMA.md Table 3 / hard rule #5)

Actual nucleotide sequences go ONLY to `novel_alleles_seqs.fa`, written under `--outroot`
(per-person-data territory, VM-only, never `reports/`). `novel_alleles.tsv` and the markdown report
carry hashes/counts/lengths only -- never a sequence, never a bare person_id.

Usage:
    # Fixtures first, always (SCHEMA.md hard rule #6):
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 300
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures --sample
    python3 scripts/hla_popgen/02_build_cohorts.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-full /tmp/hla_fixtures/immuannot_cohort_full.tsv \\
        --ancestry-preds /tmp/hla_fixtures/ancestry_preds.tsv \\
        --hla-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --skip-mount-check --sample
    python3 scripts/hla_popgen/03_novel_alleles.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --out-dir /tmp/hla_fixtures/reports --sample

    # Real run on the VM (defaults match 01_extract_rich.py / 02_build_cohorts.py's real paths):
    python3 scripts/hla_popgen/03_novel_alleles.py
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

import pandas as pd

# Person-id directories live under people/, not directly at the top level -- ~12,000 top-level
# entries makes the Workbench Jupyter file browser unusably slow (real incident, 2026-09-04; see
# RUNBOOK.md for the one-time move). --outroot points at people/ (where cds.fa.gz lives);
# DEFAULT_DATA_ROOT is the unmoved top level, where hla_calls_rich.tsv/cohort_membership.tsv
# (written by 01/02) actually live -- do not derive those paths from --outroot any more.
DEFAULT_DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_OUTROOT = os.path.join(DEFAULT_DATA_ROOT, "people")
DEFAULT_REPORTS_DIR_NAME = os.path.join("reports", "hla_popgen")

# minimap2 short-form `cs` tag tokens: ':N' (N identical bases), '*xy' (substitution ref x, obs y),
# '+seq'/'-seq' (insertion/deletion of seq). See reference/IMMUANNOT_GTF_SPEC.md part A/C.
CS_TOKEN_RE = re.compile(r":\d+|\*[a-z]{2}|[+-][a-z]+", re.IGNORECASE)
HEADER_KEY_RE = re.compile(r"^(.*)_(\d+)$")

# SCHEMA.md "template_warning policy": 00_recon_vm.py measured template_warning present on 95.3% of
# transcript rows in the real 200-person production census (00b_warning_census.py, 2026-09-03) --
# but that 95% figure was itself a bug: it counted the attribute's mere PRESENCE, and Immuannot
# writes the literal string `template_warning "NA"` to mean *no warning* on 57.4% of rows (only
# 4.6% omit the attribute entirely). The TRUE warning rate is ~38%: `partial_CDS` 24.8%,
# `no-start_codon` 6.8%, `no-stop_codon` 6.3%, `inframe_stop` 0% (never observed in 200 people).
# `template_warning` describes whether the TEMPLATE's CDS could be cleanly reconstructed from the
# gene-level alignment (searchTemplate.py's checkCDScompleteness()), not whether the typing call
# itself is wrong. Pseudogenes legitimately have no valid start/stop codon and will always warn
# (HLA-N/HLA-S are 100% partial_CDS; HLA-P/HLA-T/HLA-W are dominated by paired no-start/no-stop) --
# `no-start_codon`/`no-stop_codon` are correct pseudogene biology, not disqualifying by default.
# `partial_CDS`, however, means the CDS reconstruction was truncated: a truncated CDS cannot
# support a novel-allele claim, and at only ~2.4% of classical-gene calls, excluding it is cheap.
# `inframe_stop` is retained despite never being observed in 200 people -- it is the token that
# would genuinely indicate a broken reconstruction if it appeared. Every filter must be
# TOKEN-AWARE and `NA`-aware (never gate on bare presence of the attribute). This default is a CLI
# flag (--disqualifying-warnings), never hardcoded past this one module-level constant.
DEFAULT_DISQUALIFYING_WARNINGS = frozenset(["partial_CDS", "inframe_stop"])


def warning_tokens(template_warning):
    """Split a (possibly comma-joined, possibly NaN/None) template_warning value into a set of
    tokens. Per reference/IMMUANNOT_GTF_SPEC.md part A: 'template_warning "no-start_codon,
    inframe_stop,..."' -- comma-joined, conditional (absent, not zero, when checkCDScompleteness()
    raised nothing)."""
    if template_warning is None or (isinstance(template_warning, float) and pd.isna(template_warning)):
        return frozenset()
    s = str(template_warning).strip()
    # "NA" (any case) is Immuannot's literal spelling of "no warning" -- it fires on 57.4% of
    # transcript rows (SCHEMA.md's "template_warning policy"), far more often than the attribute
    # is simply absent (4.6%). Excluding only "nan"/"none" (pandas/Python missing-value spellings)
    # and not "NA" was exactly the bug: it returned frozenset({"NA"}) for clean calls.
    if not s or s.lower() in ("nan", "none", "na"):
        return frozenset()
    return frozenset(t.strip() for t in s.split(",") if t.strip() and t.strip().lower() != "na")


def has_disqualifying_warning(tokens, disqualifying_set):
    return bool(tokens & disqualifying_set)


# ---------------------------------------------------------------------------
# cds.fa.gz parsing / matching
# ---------------------------------------------------------------------------
def parse_cds_fasta(path):
    """Parse one hap's cds.fa.gz. Returns {rest: [(i, seq), ...]} where header was
    '>{rest}_{i}' (rest = '{contig}_{gene}', i = the detection-order index, part D)."""
    index = defaultdict(list)
    if not os.path.exists(path):
        return index
    key = None
    seq_chunks = []

    def flush():
        if key is None:
            return
        seq = "".join(seq_chunks).upper()
        m = HEADER_KEY_RE.match(key)
        if m:
            index[m.group(1)].append((int(m.group(2)), seq))

    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                key = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    flush()
    return index


def is_homopolymer_indel_only(cds_mut):
    """QC flag reproducing (a weaker proxy of) the Huijse et al. 2023 pan-MHC graph paper's
    homopolymer-indel artifact filter (NOVEL_LIT.md section 1.3/1.4 item 3) -- see module docstring
    DEVIATIONS section 1 for exactly how and why this differs from the paper's full-context check.

    True iff cds_mut's cs-string contains >=1 indel token, ALL diff tokens are indels (no
    substitutions), and every indel's inserted/deleted run is a homopolymer (single repeated base).
    False (not flagged) for any substitution-containing or unparseable cds_mut -- fail toward
    NOT excluding a real candidate on missing/ambiguous evidence.
    """
    if not cds_mut or not isinstance(cds_mut, str):
        return False
    # cds_mut = "ref_allele|cs_string|aa_diff", comma-joined per tied candidate (part A). The cs
    # string is the middle pipe field of each comma-separated candidate entry.
    any_indel = False
    for entry in cds_mut.split(","):
        parts = entry.split("|")
        if len(parts) < 2:
            continue
        cs = parts[1]
        tokens = CS_TOKEN_RE.findall(cs)
        if not tokens:
            continue
        for tok in tokens:
            if tok.startswith("*"):
                return False  # a real substitution present -> not homopolymer-indel-only
            if tok[0] in "+-":
                any_indel = True
                run = tok[1:]
                if len(set(run.upper())) != 1:
                    return False  # indel present but not a homopolymer run
    return any_indel


def nearest_allele_of(row):
    """Best available 'closest documented IPD allele' for a novel row: the tied-best candidate list
    when a CDS-level search actually ran, else the gene-level template (SCHEMA.md Table 3)."""
    alleles = row.get("alleles")
    if isinstance(alleles, str) and alleles.strip():
        return alleles.split(",")[0].strip()
    ta = row.get("template_allele")
    return ta if isinstance(ta, str) else None


# ---------------------------------------------------------------------------
# Matching Table 1 novel rows to cds.fa.gz sequences
# ---------------------------------------------------------------------------
def match_novel_rows(table1, outroot):
    """Returns (matched_rows: list of dict, stats: dict).

    matched_rows carry exactly what's needed to cluster + report: person_id, hap, gene, gene_class,
    cds_seq_sha1, cds_len, nearest_allele, cds_distance, n_aa_changes, novelty_class,
    warning_tokens (frozenset -- see module docstring "template_warning policy": NEVER collapse this
    to a bare has-any-warning bool -- that gate misreads "NA" as a warning and rejects nearly
    every call), is_homopolymer_indel_only.
    NEVER the sequence itself (kept in a side dict, written only to the VM-only fasta).
    """
    table1 = table1.copy()
    table1["is_novel_bool"] = table1["is_novel"].astype(str).str.lower().isin(["true", "1"])
    novel = table1[table1["is_novel_bool"]]

    # Ambiguity: count Table-1 rows (ALL rows, novel or not) sharing (person_id, hap, contig, gene)
    # -- >1 means real copy_index>1 on this contig, i.e. the cds.fa.gz join is genuinely ambiguous
    # (reference/IMMUANNOT_GTF_SPEC.md part D "Caveat").
    group_sizes = table1.groupby(["person_id", "hap", "contig", "gene"]).size()

    matched_rows = []
    seqs_by_hash = {}
    n_ambiguous = 0
    n_missing_cds = 0
    n_matched = 0

    # Batch by (person_id, hap) so each cds.fa.gz is opened exactly once.
    for (person_id, hap), sub in novel.groupby(["person_id", "hap"]):
        cds_path = os.path.join(outroot, person_id, "immuannot_output", hap, "cds.fa.gz")
        fasta_index = None  # lazy -- don't touch disk for ambiguous-only haps
        for _, row in sub.iterrows():
            key = (row["person_id"], row["hap"], row["contig"], row["gene"])
            if group_sizes.loc[key] > 1:
                n_ambiguous += 1
                continue
            if fasta_index is None:
                fasta_index = parse_cds_fasta(cds_path)
            rest = f"{row['contig']}_{row['gene']}"
            candidates = fasta_index.get(rest, [])
            if len(candidates) == 0:
                n_missing_cds += 1
                continue
            if len(candidates) > 1:
                # Should not happen given group_sizes==1, but never guess (SCHEMA.md discipline).
                n_ambiguous += 1
                continue
            _i, seq = candidates[0]
            if not seq:
                n_missing_cds += 1
                continue
            sha1 = hashlib.sha1(seq.encode("ascii")).hexdigest()
            seqs_by_hash.setdefault(sha1, seq)
            matched_rows.append({
                "person_id": row["person_id"], "hap": row["hap"], "gene": row["gene"],
                "gene_class": row["gene_class"],
                "cds_seq_sha1": sha1, "cds_len": len(seq),
                "nearest_allele": nearest_allele_of(row),
                "cds_distance": row.get("cds_distance"),
                "n_aa_changes": row.get("n_aa_changes"),
                "novelty_class": row.get("novelty_class"),
                "warning_tokens": warning_tokens(row.get("template_warning")),
                "is_homopolymer_indel_only": is_homopolymer_indel_only(row.get("cds_mut")),
            })
            n_matched += 1

    stats = {
        "n_novel_table1_rows": int(len(novel)),
        "n_matched": n_matched,
        "n_ambiguous_copy": n_ambiguous,
        "n_missing_cds": n_missing_cds,
    }
    return matched_rows, seqs_by_hash, stats


# ---------------------------------------------------------------------------
# Clustering into Table 3
# ---------------------------------------------------------------------------
def build_table3(matched_rows, ancestry_by_person, disqualifying_warnings=DEFAULT_DISQUALIFYING_WARNINGS):
    """Group matched rows by (gene, cds_seq_sha1) -- the exact-sequence-identity clustering key
    (NOVEL_LIT.md 2.3). novel_id = <gene>_nov_<sha1[:8]> -- deterministic (SCHEMA.md Table 3).

    `passes_qc` is TOKEN-AWARE on template_warning (SCHEMA.md "template_warning policy") -- a
    cluster is disqualified on warnings only if some member carries a token in
    `disqualifying_warnings` (default: {'partial_CDS', 'inframe_stop'}), never on bare presence of
    any warning (and never on the literal token "NA", which means "no warning" -- see
    `warning_tokens()`).

    Clusters whose (majority) `novelty_class` is `undetermined` (novelty_depth==1 -- even the
    gene-level field unresolved, SCHEMA.md Table 1) are still emitted here (Table 3's grain is
    "every distinct novel sequence cluster") but the caller (`write_report`) must exclude them from
    headline novel-allele counts and report them as their own labelled category -- an allele whose
    gene-level identity is unresolved is not a defensible "novel allele" claim.
    """
    clusters = defaultdict(list)
    for r in matched_rows:
        clusters[(r["gene"], r["cds_seq_sha1"])].append(r)

    out_rows = []
    for (gene, sha1), members in clusters.items():
        novel_id = f"{gene}_nov_{sha1[:8]}"
        persons = {m["person_id"] for m in members}
        haps = {(m["person_id"], m["hap"]) for m in members}

        # Representative annotation fields: identical CDS should carry consistent CDS-level
        # annotation across all members; take the majority value and don't silently hide
        # disagreement -- report it as a warning field instead.
        def majority(field):
            vals = [m[field] for m in members if m[field] is not None and not pd.isna(m[field])]
            if not vals:
                return None
            return Counter(vals).most_common(1)[0][0]

        gene_class = majority("gene_class")
        nearest_allele = majority("nearest_allele")
        cds_distance = majority("cds_distance")
        n_aa_changes = majority("n_aa_changes")
        novelty_class = majority("novelty_class")
        cluster_warning_tokens = frozenset().union(*(m["warning_tokens"] for m in members)) \
            if members else frozenset()
        any_disqualifying_warning = has_disqualifying_warning(cluster_warning_tokens,
                                                                disqualifying_warnings)
        is_homopolymer = any(m["is_homopolymer_indel_only"] for m in members)

        ancestry_counts = Counter()
        for pid, _hap in haps:
            anc = ancestry_by_person.get(pid)
            if anc:
                ancestry_counts[anc] += 1

        n_persons = len(persons)
        passes_qc = (n_persons >= 2) and (not any_disqualifying_warning) and (not is_homopolymer)

        out_rows.append({
            "novel_id": novel_id,
            "gene": gene,
            "gene_class": gene_class,
            "nearest_allele": nearest_allele,
            "cds_distance": cds_distance,
            "n_aa_changes": n_aa_changes,
            "novelty_class": novelty_class,
            "cds_seq_sha1": sha1,
            "cds_len": members[0]["cds_len"],
            "n_haplotypes": len(haps),
            "n_persons": n_persons,
            "n_ancestries": len(ancestry_counts),
            "ancestry_counts": json.dumps(dict(sorted(ancestry_counts.items()))),
            "is_homopolymer_indel_only": is_homopolymer,
            "passes_qc": passes_qc,
        })

    return pd.DataFrame(out_rows)


TABLE3_COLUMNS = [
    "novel_id", "gene", "gene_class", "nearest_allele", "cds_distance", "n_aa_changes",
    "novelty_class", "cds_seq_sha1", "cds_len", "n_haplotypes", "n_persons", "n_ancestries",
    "ancestry_counts", "is_homopolymer_indel_only", "passes_qc",
]


# ---------------------------------------------------------------------------
# Ancestry-gradient integration check (the key assertion NOVEL_LIT.md / the task both call for)
# ---------------------------------------------------------------------------
def novel_rate_by_ancestry(table1, ancestry_by_person):
    """Fraction of (person, hap, gene) Table-1 rows that are novel, per ancestry. This is computed
    directly off Table 1 (not Table 3's clusters) so it measures the raw per-call novelty rate the
    fixtures deliberately encode (afr 0.22 vs eur 0.05) -- the primary integration assertion this
    whole pipeline must recover."""
    t = table1.copy()
    t["is_novel_bool"] = t["is_novel"].astype(str).str.lower().isin(["true", "1"])
    t["ancestry"] = t["person_id"].map(ancestry_by_person)
    t = t[t["ancestry"].notna()]
    g = t.groupby("ancestry")["is_novel_bool"]
    return (g.mean().rename("novel_rate").to_frame()
            .join(g.size().rename("n_calls"))
            .join(g.sum().rename("n_novel"))
            .sort_values("novel_rate", ascending=False))


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def load_table1(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 1 not found at {path!r}. Run 01_extract_rich.py first.")
    return pd.read_csv(path, sep="\t", dtype={"person_id": str})


def load_cohort_membership(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 4 (cohort_membership.tsv) not found at {path!r}. Run "
                 f"02_build_cohorts.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    if "ancestry_pred" not in df.columns:
        sys.exit(f"FATAL: {path!r} has no 'ancestry_pred' column. Actual: {list(df.columns)}")
    df["ancestry_pred"] = df["ancestry_pred"].astype(str).str.upper()
    df.loc[df["ancestry_pred"].isin(["NAN", "NONE", ""]), "ancestry_pred"] = None
    return dict(zip(df["person_id"], df["ancestry_pred"]))


def write_seqs_fasta(path, seqs_by_hash, gene_by_hash):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for sha1, seq in sorted(seqs_by_hash.items()):
            gene = gene_by_hash.get(sha1, "UNKNOWN")
            f.write(f">{gene}_nov_{sha1[:8]}\n")
            for k in range(0, len(seq), 60):
                f.write(seq[k:k + 60] + "\n")


def write_report(md_path, table3, match_stats, rate_df, out_paths, matched_rows,
                  disqualifying_warnings=DEFAULT_DISQUALIFYING_WARNINGS):
    # SCHEMA.md Fix 1 steer: a cluster whose gene-level identity is itself unresolved
    # (novelty_class == "undetermined", novelty_depth==1) is not a defensible "novel allele" claim.
    # It stays IN Table 3 (the TSV's grain is every distinct novel sequence cluster) but is split out
    # of every headline novel-allele count here and reported as its own labelled category, with the
    # count stated explicitly rather than left implicit in a combined total.
    is_undetermined = table3["novelty_class"] == "undetermined" if len(table3) else pd.Series(
        dtype=bool)
    resolved = table3[~is_undetermined] if len(table3) else table3
    undetermined = table3[is_undetermined] if len(table3) else table3

    total_clusters = len(table3)
    n_undetermined = len(undetermined)
    n_resolved = len(resolved)
    n_pass = int(resolved["passes_qc"].sum()) if n_resolved else 0
    n_undetermined_pass = int(undetermined["passes_qc"].sum()) if n_undetermined else 0
    n_homopolymer = int(table3["is_homopolymer_indel_only"].sum()) if total_clusters else 0
    n_singleton = int((table3["n_persons"] == 1).sum()) if total_clusters else 0

    md = []
    md.append("# Novel HLA allele discovery -- Table 3 (`novel_alleles.tsv`)\n")
    md.append(
        "Methodology: scripts/hla_popgen/research/NOVEL_LIT.md. Novel-allele candidates = "
        "Immuannot's `new`-tagged consensus calls (Table 1), clustered by exact observed-CDS "
        "sequence identity (sha1). Recurrence in >=2 unrelated persons is the primary evidence a "
        "candidate is real biology rather than an assembly artifact (NOVEL_LIT.md section 1.4 item "
        "5) -- this is the single hardest gate in `passes_qc`.\n"
    )

    md.append("\n## Matching Table 1 -> `cds.fa.gz` sequences\n")
    md.append("| Metric | Count |")
    md.append("|---|---|")
    for k, v in match_stats.items():
        md.append(f"| {k} | {v} |")
    md.append(
        f"\nAmbiguous-copy rows (`n_ambiguous_copy`) are genuinely unresolvable per "
        f"reference/IMMUANNOT_GTF_SPEC.md part D (the detection-order join key `tmp.gene.csv` is "
        f"deleted by Immuannot's own cleanup) -- excluded from Table 3 by design, not a bug.\n"
    )

    md.append("\n## `template_warning` QC policy for this run\n")
    md.append(
        f"Disqualifying warning token(s) used in this run's `passes_qc` gate: "
        f"**{', '.join(sorted(disqualifying_warnings)) or '(none -- every candidate passes on warnings)'}**. "
        "Per SCHEMA.md's `template_warning` policy, mere presence of ANY warning is NOT a gate -- "
        "Immuannot writes the literal string `NA` to mean *no warning* on 57.4% of real transcript "
        "rows, so the true warning rate is ~38%, not ~95% (a real warning describes whether the "
        "template's CDS could be cleanly reconstructed, not whether the typing call is wrong). Only "
        "the token(s) listed above disqualify a candidate; every other token is treated as benign by "
        "this run. Override with `--disqualifying-warnings`.\n"
    )
    warning_token_counts = Counter()
    n_candidates_any_warning = 0
    for r in matched_rows:
        toks = r.get("warning_tokens") or frozenset()
        if toks:
            n_candidates_any_warning += 1
        warning_token_counts.update(toks)
    md.append(f"\nWarning-token breakdown among the {len(matched_rows)} matched novel candidates "
               f"(a candidate may carry more than one token, so counts need not sum to the total; "
               f"{n_candidates_any_warning} carried at least one token):\n")
    md.append("| warning token | n candidates | disqualifying in this run? |")
    md.append("|---|---|---|")
    if warning_token_counts:
        for tok, n in warning_token_counts.most_common():
            flag = "YES" if tok in disqualifying_warnings else "no (benign by default)"
            md.append(f"| {tok} | {n} | {flag} |")
    else:
        md.append("| (none observed) | 0 | -- |")

    md.append("\n## Cluster-level summary\n")
    md.append(
        "Headline counts below cover only clusters with a RESOLVED gene-level identity "
        "(`novelty_class != \"undetermined\"`). Depth-1 (`undetermined`) clusters -- even the "
        "gene-level field unresolved -- are reported separately immediately after, per SCHEMA.md's "
        "Fix 1 steer: an allele whose gene identity is itself unresolved is not a defensible "
        "\"novel allele\" claim.\n"
    )
    md.append("| Metric | Count |")
    md.append("|---|---|")
    md.append(f"| Distinct novel-allele clusters, resolved identity (`novel_id`) | {n_resolved} |")
    md.append(f"| Passing all QC gates (`passes_qc`), resolved identity | {n_pass} |")
    md.append(f"| Singleton (1 person only) -- excluded by the recurrence gate | {n_singleton} |")
    md.append(f"| Flagged homopolymer-indel-only artifact | {n_homopolymer} |")
    md.append(
        f"\n**`undetermined` (gene-level-unresolved) clusters: {n_undetermined}** -- included in "
        f"the Table 3 TSV, EXCLUDED from every headline count above. Of those, {n_undetermined_pass} "
        f"also pass the recurrence/warning/homopolymer gates (i.e. would look like real novel "
        f"alleles by every gate except gene-level resolution) -- reported here explicitly rather "
        f"than folded into either total.\n"
    )
    md.append(f"\n(Total clusters across both categories: {total_clusters}.)\n")

    if n_resolved:
        md.append("\n## Synonymous vs non-synonymous breakdown (`novelty_class`, resolved only)\n")
        vc = resolved["novelty_class"].fillna("NA (no novelty_class)").value_counts()
        md.append("| novelty_class | n clusters | % |")
        md.append("|---|---|---|")
        for cls, n in vc.items():
            md.append(f"| {cls} | {n} | {100*n/n_resolved:.1f}% |")
        md.append(
            "\nA strong excess of `protein_altering` (non-synonymous) clusters over what a uniform "
            "random-error process would produce is itself evidence of real biology under balancing "
            "selection (Zhou et al. 2024's own validation argument, NOVEL_LIT.md section 0) -- this "
            "table reports the raw breakdown; 04_allele_saturation.py's neutral-model comparison is "
            "where that claim gets a formal statistical treatment.\n"
        )

        md.append("\n## By gene_class (resolved-identity clusters only)\n")
        gc = resolved.groupby("gene_class").agg(
            n_clusters=("novel_id", "count"), n_passing_qc=("passes_qc", "sum")).reset_index()
        md.append("| gene_class | n clusters | passing QC |")
        md.append("|---|---|---|")
        for _, r in gc.sort_values("n_clusters", ascending=False).iterrows():
            md.append(f"| {r['gene_class']} | {r['n_clusters']} | {r['n_passing_qc']} |")

    md.append(
        "\n## Ancestry gradient in raw novel-call rate (Table 1, key integration assertion)\n"
        "Fraction of `(person, hap, gene)` calls tagged novel, by ancestry. IPD-IMGT/HLA's "
        "European bias predicts non-European ancestries should show a materially higher novel "
        "rate -- against synthetic fixtures this must recover the encoded gradient "
        "(afr 0.22 vs eur 0.05); against real data this is the paper's central empirical claim.\n"
    )
    md.append("| ancestry | novel_rate | n_calls | n_novel |")
    md.append("|---|---|---|---|")
    for anc, r in rate_df.iterrows():
        md.append(f"| {anc} | {r['novel_rate']:.3f} | {int(r['n_calls'])} | {int(r['n_novel'])} |")

    md.append(f"\nOutputs: `{out_paths['table3']}`" +
              (f", `{out_paths['seqs']}` (VM-only, sequences never leave pipeline_outputs)"
               if out_paths.get("seqs") else "") + "\n")

    with open(md_path, "w") as f:
        f.write("\n".join(md))
    return "\n".join(md)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root holding <person_id>/immuannot_output/hap{1,2}/cds.fa.gz. Default: "
                         "~/pipeline_outputs/people (NOT ~/pipeline_outputs itself -- person_id "
                         "dirs live one level deeper, see RUNBOOK.md).")
    ap.add_argument("--table1", default=None,
                    help="Path to hla_calls_rich.tsv (Table 1). Default: "
                         "~/pipeline_outputs/hla_calls_rich.tsv")
    ap.add_argument("--cohort-membership", default=None,
                    help="Path to cohort_membership.tsv (Table 4, for ancestry). Default: "
                         "~/pipeline_outputs/cohort_membership.tsv")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write novel_alleles.tsv + the markdown report (aggregate-only "
                         "-- SCHEMA.md hard rule #5). Default: <repo_root>/reports/hla_popgen")
    ap.add_argument("--seqs-path", default=None,
                    help="Where to write the VM-only novel_alleles_seqs.fa (real nucleotide "
                         "sequences -- never belongs in a report or the repo). Default: "
                         "~/pipeline_outputs/novel_alleles_seqs[.sample].fa. Override this for any "
                         "local/fixture run so it doesn't write into a real ~/pipeline_outputs on "
                         "whatever machine happens to run the script.")
    ap.add_argument("--sample", action="store_true",
                    help="Write to sample-suffixed output paths (quirk #22b).")
    ap.add_argument("--disqualifying-warnings", default="partial_CDS,inframe_stop",
                    help="Comma-separated template_warning tokens that disqualify a novel-allele "
                         "cluster from passes_qc (SCHEMA.md 'template_warning policy' -- mere "
                         "presence of ANY warning is NOT a valid gate: Immuannot writes the "
                         "literal string 'NA' to mean 'no warning' on 57.4%% of rows, so the true "
                         "warning rate is ~38%%, not ~95%%). Default: 'partial_CDS,inframe_stop' "
                         "-- a truncated CDS can't support a novel-allele claim, and partial_CDS "
                         "is only ~2.4%% of classical-gene calls so excluding it is cheap. "
                         "'no-start_codon'/'no-stop_codon' are deliberately excluded from the "
                         "default: they are expected pseudogene biology (HLA-N/HLA-S are 100%% "
                         "partial_CDS; HLA-P/HLA-T/HLA-W are dominated by paired no-start/no-stop), "
                         "not evidence of a broken call. Pass '' for no warning-based "
                         "disqualification at all.")
    args = ap.parse_args()

    disqualifying_warnings = frozenset(
        t.strip() for t in args.disqualifying_warnings.split(",") if t.strip())

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    table1_path = args.table1 or os.path.join(DEFAULT_DATA_ROOT, "hla_calls_rich.tsv")
    cohort_path = args.cohort_membership or os.path.join(DEFAULT_DATA_ROOT, "cohort_membership.tsv")
    out_dir = args.out_dir or os.path.join(repo_root, DEFAULT_REPORTS_DIR_NAME)
    os.makedirs(out_dir, exist_ok=True)

    suffix = ".sample" if args.sample else ""
    table3_path = os.path.join(out_dir, f"novel_alleles{suffix}.tsv")
    md_path = os.path.join(out_dir, f"03_novel_alleles_report{suffix}.md")
    # Written under DEFAULT_DATA_ROOT, not --outroot -- keeps people/ containing ONLY
    # person-id-named directories, which is what keeps it fast to move/browse/rsync as a unit.
    seqs_path = args.seqs_path or os.path.join(DEFAULT_DATA_ROOT, f"novel_alleles_seqs{suffix}.fa")

    t0 = time.time()
    print(f"Loading Table 1 from {table1_path!r} ...", file=sys.stderr)
    table1 = load_table1(table1_path)
    print(f"Loading Table 4 (ancestry) from {cohort_path!r} ...", file=sys.stderr)
    ancestry_by_person = load_cohort_membership(cohort_path)

    print("Matching novel Table-1 rows to cds.fa.gz sequences ...", file=sys.stderr)
    matched_rows, seqs_by_hash, match_stats = match_novel_rows(table1, args.outroot)
    print(f"  {match_stats}", file=sys.stderr)

    print("Clustering by exact CDS sequence identity ...", file=sys.stderr)
    print(f"  Disqualifying template_warning token(s) for passes_qc: "
          f"{sorted(disqualifying_warnings) or '(none)'}", file=sys.stderr)
    table3 = build_table3(matched_rows, ancestry_by_person, disqualifying_warnings)
    gene_by_hash = {r["cds_seq_sha1"]: r["gene"] for r in matched_rows}

    for c in TABLE3_COLUMNS:
        if c not in table3.columns:
            table3[c] = pd.NA
    table3 = table3[TABLE3_COLUMNS]
    table3.to_csv(table3_path, sep="\t", index=False)
    write_seqs_fasta(seqs_path, seqs_by_hash, gene_by_hash)

    rate_df = novel_rate_by_ancestry(table1, ancestry_by_person)

    md_text = write_report(md_path, table3, match_stats, rate_df,
                            {"table3": table3_path, "seqs": seqs_path},
                            matched_rows, disqualifying_warnings)

    elapsed = time.time() - t0
    n_undetermined = int((table3["novelty_class"] == "undetermined").sum()) if len(table3) else 0
    n_resolved_pass = int(
        table3.loc[table3["novelty_class"] != "undetermined", "passes_qc"].sum()
    ) if len(table3) else 0
    print(f"\nWrote {len(table3)} novel-allele clusters to {table3_path!r} "
          f"({n_resolved_pass} resolved-identity clusters passing QC; "
          f"{n_undetermined} additional gene-level-undetermined clusters excluded from that "
          f"headline count -- see report).",
          file=sys.stderr)
    print(f"Wrote {len(seqs_by_hash)} sequences to {seqs_path!r} (VM-only).", file=sys.stderr)
    print(f"Wrote report to {md_path!r}. ({elapsed:.1f}s)", file=sys.stderr)
    print("\nAncestry novel-rate gradient (key integration check):", file=sys.stderr)
    print(rate_df, file=sys.stderr)


if __name__ == "__main__":
    main()
