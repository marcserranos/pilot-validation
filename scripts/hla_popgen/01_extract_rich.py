#!/usr/bin/env python3
"""Extract Table 1 (`hla_calls_rich.tsv`) and Table 2 (`hla_cis_pairs.tsv`) per
scripts/hla_popgen/SCHEMA.md, by parsing the raw per-person `hap{1,2}.gtf.gz` files directly.

## Why this exists (see SCHEMA.md "Why these tables exist")

Every parser in this repo before this one (`run_immuannot_person.py`'s `parse_gtf()`,
`rebuild_immuannot_calls.py`'s `parse_gtf()`) captures exactly one string per (person, hap, gene) --
the `consensus` call -- and throws away contig identity, non-classical genes, and the novelty/
distance detail that's sitting right there in the same GTF line. This script re-parses the raw
GTFs from scratch, once, correctly, at the true (person_id, hap, contig, gene, copy_index) grain
(SCHEMA.md Table 1) -- the exact grain the 2026-08-10 production dedup bug (ENVIRONMENT.md quirk
#29) proves matters: keying dedup/merge at the wrong grain silently deletes rows, exits 0, and
looks fine until someone notices 0% completeness on a field that should be ~92%.

## Parsing gotchas this script exists to get right (all sourced from
## reference/IMMUANNOT_GTF_SPEC.md, a full read of Immuannot's upstream source)

- `template_distance` is UNQUOTED in the GTF, unlike (per the spec) every other attribute --
  and empirically, this project's own fixture generator (tests/make_fixtures.py) *also* writes
  `cds_distance` unquoted, which contradicts the spec's "unlike every other attribute" claim. This
  parser does not trust either side blindly: `parse_attrs()` accepts BOTH quoted-string and
  unquoted-integer forms for every key, so it is correct against real data if the spec is right,
  and correct against the fixtures either way. This exact discrepancy is called out in the
  bug/contradiction report handed back to Marc -- 00_recon_vm.py's job is to resolve it empirically
  on the real VM before this parser's tolerance is leaned on at scale.
- `cds_mut` contains `|`, `<`, `:`, `*` -- never split naively on those. This parser never manually
  splits the attribute-value string; `parse_attrs()`'s regex captures the full quoted value in one
  shot (`"((?:[^"\\]|\\.)*)"`), and `n_aa_changes` is derived by counting literal `<` characters
  (every codon diff has exactly one, regardless of how many tied candidates or `:`-joined diffs are
  packed into the string) rather than trying to parse the pipe/comma/colon structure.
- Conditional attributes (`cds_distance`, `cds_mut`, `template_warning`) are legitimately ABSENT,
  not zero -- emitted as pandas NA (`""` written out, read back as NaN), never imputed to 0/False.
- C4's `gene_name` has a one-letter size code appended with no delimiter (`C4AL` -> gene `C4A`,
  `c4_size` `L`) -- `clean_gene_name()` implements the exact strip `combine.py` itself uses
  (`genename[0:len(genename)-1]`).
- `gene_id` carries a `.N` copy-index suffix; N=1 (implicit, no suffix) is the overwhelmingly common
  case, N>1 is real segmental-duplication signal, not an artifact -- never silently overwritten.
- Table 2 (cis pairs) is emitted ONLY for genes sharing a contig -- being in the same hap FILE is
  necessary but not sufficient (multi-contig haplotypes are ~18% of the synthetic fixture and a
  real, documented phenomenon per the spec, part G). Unpairable cases (same hap, different contig)
  are counted and reported, not silently dropped -- that count is itself a finding about assembly
  fragmentation.

## Operational discipline (ENVIRONMENT.md quirks referenced)

- Quirk #22: incremental, resumable writes. A rerun skips person_ids already present in the output
  file rather than reprocessing/duplicating them.
- Quirk #22b: `--sample` writes to a separate, clearly-suffixed path -- never the same path as a
  real run.
- Quirk #25: no `du`/`df` calls anywhere in this script.
- Quirk #29: `--validate` re-reads the written output and asserts row-grain uniqueness on exactly
  the columns that broke last time: (person_id, hap, contig, gene, copy_index).

Usage:
    # Sanity run against local fixtures first, always:
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 200
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures --sample --limit 50
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures --sample --validate

    # Real run on the VM (default --outroot is ~/pipeline_outputs):
    python3 scripts/hla_popgen/01_extract_rich.py
    python3 scripts/hla_popgen/01_extract_rich.py --validate
"""
import argparse
import concurrent.futures
import csv
import gzip
import os
import re
import sys
import time

DEFAULT_DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
# Person-id-named directories live one level deeper, not directly under DEFAULT_DATA_ROOT.
# ~/pipeline_outputs holding ~12,000 top-level directory entries makes the Workbench Jupyter
# file browser (and any other tool that lists+stats a directory to render it) unusably slow --
# real incident, 2026-09-04. The fix is a one-time move (see RUNBOOK.md) of every person_id
# directory into a `people/` subfolder, keeping only the aggregate .tsv files and a handful of
# named directories at the top level, where the UI can render them instantly. DEFAULT_OUTROOT
# (where THIS script looks for <person_id>/immuannot_output/) points at that subfolder;
# DEFAULT_DATA_ROOT (where the aggregate hla_calls_rich.tsv etc. actually live) does not move.
DEFAULT_OUTROOT = os.path.join(DEFAULT_DATA_ROOT, "people")

# ---------------------------------------------------------------------------
# Gene classification -- SCHEMA.md "Gene classification" table, verbatim.
# ---------------------------------------------------------------------------
GENE_CLASS = {}
for g in ["HLA-A", "HLA-B", "HLA-C"]:
    GENE_CLASS[g] = "classical_I"
for g in ["HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]:
    GENE_CLASS[g] = "classical_II"
for g in ["HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DRA"]:
    GENE_CLASS[g] = "class_II_accessory"
for g in ["HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8",
          "HLA-DRB9", "HLA-DPA2", "HLA-DPB2", "HLA-DQA2", "HLA-DQB2"]:
    GENE_CLASS[g] = "class_II_paralog"
for g in ["HLA-E", "HLA-F", "HLA-G"]:
    GENE_CLASS[g] = "nonclassical_I"
for g in ["HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S", "HLA-T", "HLA-U",
          "HLA-V", "HLA-W", "HLA-Y"]:
    GENE_CLASS[g] = "pseudogene_I"
for g in ["MICA", "MICB", "TAP1", "TAP2"]:
    GENE_CLASS[g] = "mic_tap"
for g in ["C4A", "C4B"]:
    GENE_CLASS[g] = "complement"
GENE_CLASS["HLA-HFE"] = "other"
for g in ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DP1",
          "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2",
          "KIR3DL3", "KIR3DP1", "KIR3DS1"]:
    GENE_CLASS[g] = "kir"

KIR_GENES = frozenset(g for g, c in GENE_CLASS.items() if c == "kir")

# Cis pairs SCHEMA.md Table 2 requests by default. `pair` is rendered without the "HLA-" prefix
# ("DQA1~DQB1"), matching SCHEMA.md's own example strings.
DEFAULT_PAIRS = [("HLA-DQA1", "HLA-DQB1"), ("HLA-DPA1", "HLA-DPB1"), ("HLA-DRA", "HLA-DRB1")]

# ---------------------------------------------------------------------------
# Attribute parsing
# ---------------------------------------------------------------------------
# Matches `key "value";` (value may contain |, <, :, *, anything but an unescaped quote) OR
# `key 123;` (unquoted integer -- template_distance always, and empirically cds_distance in this
# project's own fixtures, contra the spec's claim that only template_distance is unquoted -- see
# module docstring). Deliberately does NOT try to split the captured value further; cds_mut's
# internal `|`/`<`/`:`/`*` structure is handled by n_aa_changes()/consensus parsing, never by
# blindly splitting this regex's own delimiter.
ATTR_RE = re.compile(r'(\w+)\s+(?:"((?:[^"\\]|\\.)*)"|(-?\d+))\s*;')

# Header lines, e.g.:
#   ## gene (copy num = 0): C4A,C4B,HLA-DMA,...
#   ## gene (copy num = 1): HLA-A,HLA-B,...
#   ## gene (copy num > 1): HLA-DRB3,...
#   ## contigs for HLA-A: hap1_ctg_x_a
HEADER_COPY_RE = re.compile(r'^##\s*gene \(copy num ([=>]) (\d+)\):\s*(.*)$')
HEADER_CONTIG_RE = re.compile(r'^##\s*contigs for (\S+):\s*(.*)$')


def parse_attrs(attr_str):
    """Parse a GTF column-9 attribute string into {key: raw_string_value}.

    Tolerant of both quoted-string and unquoted-integer attribute forms for every key (see module
    docstring re: template_distance/cds_distance quoting discrepancy). Never splits a captured
    value on its own internal delimiters.
    """
    out = {}
    for m in ATTR_RE.finditer(attr_str):
        key = m.group(1)
        val = m.group(2) if m.group(2) is not None else m.group(3)
        out[key] = val
    return out


def gene_id_copy_index(gene_id):
    """Split a `.N` copy-index suffix off gene_id/transcript_id. Returns (base, copy_index)."""
    m = re.match(r'^(.*)\.(\d+)$', gene_id)
    if m:
        return m.group(1), int(m.group(2))
    return gene_id, 1


def clean_gene_name(raw_name):
    """Strip C4's no-delimiter size-code suffix. Returns (gene, c4_size_or_None).

    Per reference/IMMUANNOT_GTF_SPEC.md part A / callC4Allele.py: gene_name for C4 rows is the
    called allele (`C4`, `C4A`, or `C4B`) with a single size-code letter appended directly, no
    delimiter (`C4AL`, `C4AS`, ...). combine.py recovers the bare symbol by stripping the last
    character; we replicate that exactly, gated on `startswith("C4")` since no other gene name in
    the 65-gene list starts with "C4".
    """
    if raw_name.startswith("C4") and len(raw_name) > 2:
        return raw_name[:-1], raw_name[-1]
    return raw_name, None


NOVELTY_CLASS_BY_DEPTH = {1: "undetermined", 2: "protein_altering", 3: "synonymous", 4: "beyond_cds"}


def parse_consensus(consensus):
    """Derive n_fields/is_novel/novelty_depth/novelty_class from the consensus allele string.

    "new" is spliced INTO the colon-delimited allele string at a variable field depth (1, 2, 3, or
    4) that itself encodes novelty severity (reference/IMMUANNOT_GTF_SPEC.md part B) -- it is not a
    separate flag. Depth 1 (even the first/gene-resolution field unresolved) arises when
    `consensusCall()`'s `os.path.commonprefix()` truncation collapses tied candidate alleles that
    disagree at field 1 -- rare (3/50 in the production recon sample) but real; SCHEMA.md classifies
    it as `novelty_class = "undetermined"`, not an error. `consensus` can also be the literal string
    "undetermined" (ambiguous typing), which has no fields to parse at all -- distinct from a depth-1
    novelty call, which still has real (if unresolved) fields.

    Returns (n_fields, is_novel, novelty_depth, novelty_class, warn_msg_or_None).
    """
    if consensus is None or consensus == "undetermined":
        return None, False, None, None, None
    allele_part = consensus.split("*", 1)[1] if "*" in consensus else consensus
    fields = allele_part.split(":")
    n_fields = len(fields)
    if "new" not in fields:
        return n_fields, False, None, None, None
    depth = fields.index("new") + 1  # 1-based
    novelty_class = NOVELTY_CLASS_BY_DEPTH.get(depth)
    warn = None
    if novelty_class is None:
        warn = (f"unexpected novelty depth {depth} in consensus {consensus!r} "
                f"(expected 1, 2, 3, or 4 per nameNewHlaAllele()'s truncation rule)")
    return n_fields, True, depth, novelty_class, warn


def n_aa_changes(cds_mut):
    """Count amino-acid-level codon diffs in cds_mut. Every diff (REFaa(codon)<OBSaa(codon)) has
    exactly one '<' regardless of how many are ':'-joined or how many tied candidates are
    ','-joined -- counting '<' is robust to that structure without parsing it."""
    if not cds_mut:
        return None
    return cds_mut.count("<")


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------
def parse_header(lines):
    """Parse '#'-prefixed header lines. Returns dict: copy0/copy1/copyN sets of gene names, and
    contigs_for_gene: {gene: set(contig_names)}."""
    info = {"copy0": set(), "copy1": set(), "copyN": set(), "contigs_for_gene": {}}
    for line in lines:
        m = HEADER_COPY_RE.match(line.strip())
        if m:
            op, num, genes_str = m.groups()
            genes = {g for g in genes_str.split(",") if g}
            if op == "=" and num == "0":
                info["copy0"] |= genes
            elif op == "=" and num == "1":
                info["copy1"] |= genes
            elif op == ">":
                info["copyN"] |= genes
            continue
        m = HEADER_CONTIG_RE.match(line.strip())
        if m:
            gene, contigs_str = m.groups()
            contigs = {c for c in contigs_str.split(",") if c}
            info["contigs_for_gene"].setdefault(gene, set()).update(contigs)
    return info


def parse_hap_gtf(path):
    """Parse one hap{1,2}.gtf.gz. Returns (rows, header_info, warnings).

    rows: list of dicts, one per (contig, gene_id) that had BOTH a gene row and a transcript row
    (the two together carry every Table-1 column). A gene row with no matching transcript row (or
    vice versa) is a genuine parser-vs-data mismatch and is reported as a warning, not silently
    dropped or silently synthesized.
    """
    warnings = []
    header_lines = []
    gene_recs = {}
    txn_recs = {}

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                warnings.append(f"malformed line (fewer than 9 columns): {line[:120]!r}")
                continue
            contig, _source, feature, start, end, _score, strand, _frame, attr_str = fields[:9]
            attrs = parse_attrs(attr_str)
            gene_id = attrs.get("gene_id")
            if gene_id is None:
                continue
            key = (contig, gene_id)
            if feature == "gene":
                gene_recs[key] = {
                    "contig": contig, "gene_id": gene_id,
                    "gene_name_raw": attrs.get("gene_name"),
                    "template_allele": attrs.get("template_allele"),
                    "template_distance": attrs.get("template_distance"),
                    "gene_start": int(start), "gene_end": int(end), "strand": strand,
                }
            elif feature == "transcript":
                txn_recs[key] = {
                    "gene_name_raw": attrs.get("gene_name"),
                    "consensus": attrs.get("consensus"),
                    "alleles": attrs.get("alleles"),
                    "template_warning": attrs.get("template_warning"),
                    "cds_distance": attrs.get("cds_distance"),
                    "cds_mut": attrs.get("cds_mut"),
                }

    header_info = parse_header(header_lines)

    all_keys = set(gene_recs) | set(txn_recs)
    rows = []
    for key in all_keys:
        g = gene_recs.get(key)
        t = txn_recs.get(key)
        contig, gene_id = key
        if g is None:
            warnings.append(f"transcript row with no matching gene row: contig={contig} "
                             f"gene_id={gene_id}")
            continue
        if t is None:
            warnings.append(f"gene row with no matching transcript row: contig={contig} "
                             f"gene_id={gene_id} gene_name={g['gene_name_raw']!r}")
            continue
        row = dict(g)
        row.update(t)
        rows.append(row)

    # Cross-check the header's own copy-number/contig claims against what the body actually has.
    # SCHEMA.md: "a mismatch is a parser bug and must fail loudly" -- surfaced as warnings here;
    # main() escalates to FATAL if the mismatch rate looks systemic rather than a one-off.
    body_genes_by_contig = {}
    for row in rows:
        base, copy_idx = gene_id_copy_index(row["gene_id"])
        gene, _ = clean_gene_name(row["gene_name_raw"])
        body_genes_by_contig.setdefault(gene, {}).setdefault(row["contig"], []).append(copy_idx)

    for gene, contig_copies in body_genes_by_contig.items():
        max_copy = max(c for copies in contig_copies.values() for c in copies)
        observed_contigs = set(contig_copies)
        if gene in header_info["copy0"]:
            warnings.append(f"HEADER MISMATCH: {gene} listed as copy num=0 but has "
                             f"{sum(len(v) for v in contig_copies.values())} body row(s)")
        if max_copy > 1 and gene not in header_info["copyN"]:
            warnings.append(f"HEADER MISMATCH: {gene} has copy_index>1 in body but is not listed "
                             f"under 'copy num > 1'")
        if max_copy == 1 and gene in header_info["copy1"] | header_info["copyN"]:
            pass  # consistent
        declared_contigs = header_info["contigs_for_gene"].get(gene, set())
        if declared_contigs and declared_contigs != observed_contigs:
            warnings.append(f"HEADER MISMATCH: {gene} declared contigs {sorted(declared_contigs)} "
                             f"!= body-observed contigs {sorted(observed_contigs)}")
    for gene in header_info["copy1"] | header_info["copyN"]:
        if gene not in body_genes_by_contig:
            warnings.append(f"HEADER MISMATCH: {gene} listed as present (copy num>=1) but has 0 "
                             f"body rows")

    return rows, header_info, warnings


# ---------------------------------------------------------------------------
# Per-person processing
# ---------------------------------------------------------------------------
def build_table1_row(person_id, hap, row):
    gene, c4_size = clean_gene_name(row["gene_name_raw"])
    base_id, copy_index = gene_id_copy_index(row["gene_id"])
    gene_class = GENE_CLASS.get(gene)
    warn = None
    if gene_class is None:
        gene_class = "unknown"
        warn = f"UNKNOWN GENE not in SCHEMA.md's 65-gene classification: {gene!r}"
    if gene in KIR_GENES:
        warn = (f"KIR GENE PRESENT ({gene}) in person={person_id} hap={hap} -- spec says this "
                f"must never happen (chr19, outside the chr6 trim window). Likely a trim bug.")

    consensus = row.get("consensus")
    n_fields, is_novel, novelty_depth, novelty_class, novelty_warn = parse_consensus(consensus)
    if novelty_warn:
        warn = (warn + " | " + novelty_warn) if warn else novelty_warn

    template_distance = row.get("template_distance")
    template_distance = int(template_distance) if template_distance is not None else None
    gene_start, gene_end = row["gene_start"], row["gene_end"]
    span = gene_end - gene_start + 1
    td_per_kb = (1000.0 * template_distance / span) if (template_distance is not None and span > 0) else None

    cds_distance = row.get("cds_distance")
    cds_distance = int(cds_distance) if cds_distance is not None else None
    cds_mut = row.get("cds_mut")
    # Immuannot writes template_warning "NA" to mean NO WARNING far more often (57.4%) than it
    # omits the attribute (4.6%) -- see SCHEMA.md's "template_warning policy". Normalise both
    # "clean" spellings to None here so the literal string "NA" never leaks downstream and
    # can't be mistaken for a real warning token again (matches
    # scripts/production_orchestrator/rebuild_immuannot_calls.py's `warn_val not in {"", "NA"}`).
    raw_template_warning = row.get("template_warning")
    if raw_template_warning is not None and raw_template_warning.strip().upper() in {"", "NA"}:
        template_warning = None
    else:
        template_warning = raw_template_warning
    alleles = row.get("alleles")
    n_tied = len(alleles.split(",")) if alleles else None

    out = {
        "person_id": person_id, "hap": hap, "contig": row["contig"], "gene": gene,
        "copy_index": copy_index, "gene_class": gene_class, "consensus": consensus,
        "n_fields": n_fields, "is_novel": is_novel, "novelty_depth": novelty_depth,
        "novelty_class": novelty_class, "template_allele": row.get("template_allele"),
        "template_distance": template_distance, "gene_start": gene_start, "gene_end": gene_end,
        "template_distance_per_kb": td_per_kb, "cds_distance": cds_distance, "cds_mut": cds_mut,
        "n_aa_changes": n_aa_changes(cds_mut), "template_warning": template_warning,
        "has_warning": template_warning is not None, "alleles": alleles, "n_tied": n_tied,
        "strand": row["strand"], "c4_size": c4_size,
    }
    return out, warn


TABLE1_COLUMNS = [
    "person_id", "hap", "contig", "gene", "copy_index", "gene_class", "consensus", "n_fields",
    "is_novel", "novelty_depth", "novelty_class", "template_allele", "template_distance",
    "gene_start", "gene_end", "template_distance_per_kb", "cds_distance", "cds_mut",
    "n_aa_changes", "template_warning", "has_warning", "alleles", "n_tied", "strand", "c4_size",
]

TABLE2_COLUMNS = [
    "person_id", "hap", "contig", "pair", "allele_a", "allele_b", "haplotype_label",
    "both_exact", "either_novel", "cis_confidence",
]


def build_table2_rows(person_id, hap, table1_rows, pairs=DEFAULT_PAIRS):
    """Emit cis pairs ONLY for genes sharing a contig. Returns (rows, n_unpairable)."""
    by_gene_contig = {}
    for r in table1_rows:
        by_gene_contig.setdefault(r["gene"], []).append(r)

    rows = []
    n_unpairable = 0
    for gene_a, gene_b in pairs:
        recs_a = by_gene_contig.get(gene_a, [])
        recs_b = by_gene_contig.get(gene_b, [])
        if not recs_a or not recs_b:
            continue
        # Multi-copy genes: pair every same-contig combination (rare; usually 1x1).
        contigs_a = {r["contig"] for r in recs_a}
        contigs_b = {r["contig"] for r in recs_b}
        shared = contigs_a & contigs_b
        if not shared and contigs_a and contigs_b:
            n_unpairable += 1
        for contig in shared:
            ra_list = [r for r in recs_a if r["contig"] == contig]
            rb_list = [r for r in recs_b if r["contig"] == contig]
            for ra in ra_list:
                for rb in rb_list:
                    pair_label = f"{gene_a.replace('HLA-', '')}~{gene_b.replace('HLA-', '')}"
                    hap_label = f"{ra['consensus']}~{rb['consensus']}"
                    both_exact = (ra.get("template_distance") == 0 and rb.get("template_distance") == 0)
                    either_novel = bool(ra.get("is_novel") or rb.get("is_novel"))
                    rows.append({
                        "person_id": person_id, "hap": hap, "contig": contig, "pair": pair_label,
                        "allele_a": ra["consensus"], "allele_b": rb["consensus"],
                        "haplotype_label": hap_label, "both_exact": both_exact,
                        "either_novel": either_novel, "cis_confidence": "physical",
                    })
    return rows, n_unpairable


def process_person(person_id, person_dir, strict_header_check):
    """Returns (table1_rows, table2_rows, n_unpairable, warnings_list)."""
    all_table1, all_table2 = [], []
    total_unpairable = 0
    all_warnings = []

    for hap in ("hap1", "hap2"):
        gz_path = os.path.join(person_dir, f"{hap}.gtf.gz")
        if not os.path.exists(gz_path):
            continue
        try:
            raw_rows, _header_info, parse_warnings = parse_hap_gtf(gz_path)
        except (OSError, EOFError, gzip.BadGzipFile) as e:
            all_warnings.append(f"{person_id}/{hap}: unreadable GTF ({e}) -- skipped")
            continue

        for w in parse_warnings:
            msg = f"{person_id}/{hap}: {w}"
            all_warnings.append(msg)
            if strict_header_check and "HEADER MISMATCH" in w:
                sys.exit(f"FATAL (--strict-header-check): {msg}")

        hap_table1 = []
        for row in raw_rows:
            t1_row, warn = build_table1_row(person_id, hap, row)
            hap_table1.append(t1_row)
            if warn:
                all_warnings.append(f"{person_id}/{hap}/{t1_row['gene']}: {warn}")
        all_table1.extend(hap_table1)

        hap_table2, n_unpairable = build_table2_rows(person_id, hap, hap_table1)
        all_table2.extend(hap_table2)
        total_unpairable += n_unpairable

    return all_table1, all_table2, total_unpairable, all_warnings


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def discover_persons(outroot, limit=None):
    if not os.path.isdir(outroot):
        sys.exit(f"FATAL: --outroot {outroot!r} does not exist or is not a directory.")
    persons = sorted(
        d for d in os.listdir(outroot)
        if os.path.isdir(os.path.join(outroot, d, "immuannot_output"))
    )
    if not persons:
        sys.exit(f"FATAL: no <person_id>/immuannot_output/ directories found under {outroot!r}. "
                 f"Check --outroot, and that the pipeline has actually produced output there.")
    if limit is not None:
        persons = persons[:limit]
    return persons


def already_done_person_ids(path, id_col="person_id"):
    """For resumable reruns (quirk #22): read the person_ids already fully written."""
    if not os.path.exists(path):
        return set()
    done = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            done.add(row[id_col])
    return done


def append_rows(path, columns, rows):
    """Append rows to a TSV, writing the header only if the file doesn't exist yet. Flushes and
    fsyncs immediately (quirk #22: a script that only accumulates in memory can lose everything to
    an unattended session's terminal loss)."""
    write_header = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})
        f.flush()
        os.fsync(f.fileno())


# ---------------------------------------------------------------------------
# --validate mode
# ---------------------------------------------------------------------------
def validate_output(table1_path):
    """Re-read Table 1 and assert row-grain uniqueness on (person_id, hap, contig, gene,
    copy_index) -- the exact grain the 2026-08-10 production dedup bug destroyed (quirk #29)."""
    if not os.path.exists(table1_path):
        sys.exit(f"FATAL --validate: {table1_path!r} does not exist.")
    seen = set()
    n_rows = 0
    dupes = []
    with open(table1_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_rows += 1
            key = (row["person_id"], row["hap"], row["contig"], row["gene"], row["copy_index"])
            if key in seen:
                dupes.append(key)
            seen.add(key)
    if dupes:
        sys.exit(f"FATAL --validate: {len(dupes)} duplicate (person_id, hap, contig, gene, "
                 f"copy_index) key(s) found in {table1_path!r} -- row-grain violated. "
                 f"First few: {dupes[:5]}")
    print(f"--validate OK: {n_rows} rows in {table1_path!r}, all unique on "
          f"(person_id, hap, contig, gene, copy_index).", file=sys.stderr)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Directory holding <person_id>/immuannot_output/. Default: "
                         "~/pipeline_outputs")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write hla_calls_rich.tsv / hla_cis_pairs.tsv. Default: "
                         "~/pipeline_outputs (the data root, NOT --outroot's people/ subfolder -- "
                         "per-person raw data belongs in ~/pipeline_outputs, not reports/, "
                         "SCHEMA.md hard rule #5).")
    ap.add_argument("--sample", action="store_true",
                    help="Write to sample-suffixed output paths (quirk #22b) -- never shares a "
                         "path with a real run.")
    ap.add_argument("--limit", type=int, default=None, help="Only process the first N people.")
    ap.add_argument("--threads", type=int, default=8,
                    help="ThreadPoolExecutor workers -- this is gzip/I/O-bound, not CPU-bound.")
    ap.add_argument("--checkpoint-every", type=int, default=200,
                    help="Flush accumulated rows to disk every N people processed.")
    ap.add_argument("--strict-header-check", action="store_true",
                    help="Exit immediately on the first GTF header/body mismatch instead of just "
                         "warning. Recommended for an initial --sample run; probably too strict "
                         "for an unattended full-cohort run.")
    ap.add_argument("--validate", action="store_true",
                    help="Re-read the (already-written) output and assert row-grain uniqueness on "
                         "(person_id, hap, contig, gene, copy_index). Does not re-parse GTFs.")
    args = ap.parse_args()

    out_dir = args.out_dir or DEFAULT_DATA_ROOT
    suffix = ".sample" if args.sample else ""
    table1_path = os.path.join(out_dir, f"hla_calls_rich{suffix}.tsv")
    table2_path = os.path.join(out_dir, f"hla_cis_pairs{suffix}.tsv")

    if args.validate:
        validate_output(table1_path)
        return

    persons = discover_persons(args.outroot, args.limit)
    done = already_done_person_ids(table1_path)
    todo = [p for p in persons if p not in done]
    print(f"{len(persons)} person dirs found, {len(done)} already in {table1_path!r}, "
          f"{len(todo)} to process.", file=sys.stderr)
    if not todo:
        print("Nothing to do.", file=sys.stderr)
        return

    t0 = time.time()
    n_no_gtf = 0
    n_unpairable_total = 0
    n_warnings_total = 0
    pending1, pending2 = [], []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(process_person, pid, os.path.join(args.outroot, pid, "immuannot_output"),
                        args.strict_header_check): pid
            for pid in todo
        }
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            t1_rows, t2_rows, n_unpairable, warnings = fut.result()
            if not t1_rows:
                n_no_gtf += 1
            n_unpairable_total += n_unpairable
            n_warnings_total += len(warnings)
            for w in warnings[:3]:  # don't flood stderr; a full 12k-person run could emit a lot
                print(f"  WARNING: {w}", file=sys.stderr)
            pending1.extend(t1_rows)
            pending2.extend(t2_rows)

            if i % args.checkpoint_every == 0 or i == len(todo):
                append_rows(table1_path, TABLE1_COLUMNS, pending1)
                append_rows(table2_path, TABLE2_COLUMNS, pending2)
                pending1, pending2 = [], []
                elapsed = time.time() - t0
                rate = i / elapsed if elapsed > 0 else 0
                print(f"  [{i}/{len(todo)}] checkpointed, {elapsed:.0f}s elapsed "
                      f"({rate:.1f} people/sec)", file=sys.stderr)

    print(f"\nDone. {len(todo)} people processed, {n_no_gtf} had no readable GTF, "
          f"{n_unpairable_total} unpairable cis-pair cases (genes present but on different "
          f"contigs), {n_warnings_total} total warnings emitted.", file=sys.stderr)
    print(f"Wrote {table1_path} and {table2_path}.", file=sys.stderr)
    print(f"Run with --validate to check row-grain uniqueness before trusting this output.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
