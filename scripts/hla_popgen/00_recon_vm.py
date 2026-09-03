#!/usr/bin/env python3
"""CHEAP (target: under ~2 minutes on a sample of N=50 people) VM reconnaissance script. Run this
FIRST, on the real Verily Workbench VM, before anything else in this sub-project.

## Why this exists

`reference/IMMUANNOT_GTF_SPEC.md` is a full read of Immuannot's upstream source -- but it was
written with NO minimap2 available and NO access to this project's actual VM (ENVIRONMENT.md quirk
#28: VPC-SC blocks direct SSH entirely, so nothing could be run end-to-end while writing it). It
marks several load-bearing claims AMBIGUOUS and says so explicitly: whether the `hap{1,2}/`
intermediate folders (`cds.fa.gz`, `mm2.ipd.cds.paf.gz`, etc.) actually survived disk-space
pressure or manual cleanup on THIS project's real production run; whether the two independent CDS
sequence extraction routes agree; whether any GTF attribute keys exist in real output that weren't
seen in the source read. This script's whole job is to answer those questions empirically, cheaply,
on a small sample, so `01_extract_rich.py` and `02_build_cohorts.py` never have to build on an
assumption instead of a measurement.

**The single most important output**: presence rates and extrapolated total bytes for
`cds.fa.gz` / `mm2.ipd.cds.paf.gz` / `mm2.ipd.gen.paf.gz` / `gene.filtered.paf` per haplotype. The
entire future novel-allele workstream (Table 3, `novel_alleles.tsv` -- not built by this handoff,
but planned) depends on `cds.fa.gz` existing at scale; if it doesn't, that's the single fact this
whole sub-project most needs to know before anyone plans further work assuming it does.

## Design constraints (ENVIRONMENT.md quirks)

- Quirk #26: hard-fail immediately if `~/pipeline_outputs/` (or the sampled people's directories)
  don't exist, with the exact remediation. This script does NOT need the gcsfuse bucket mount --
  everything it inspects lives on the VM's local persistent disk under `~/pipeline_outputs/`.
- Quirk #25: never `du`/`df` a broad path. File sizes here are per-file `os.path.getsize()` calls
  on a small sample (N=50 people x ~2 haps x ~6 files = ~600 stat calls, not a tree walk).
- Cheap by construction: only reads N sampled people's files, and for the GTF attribute/gene
  inventory, reads each sampled GTF fully (they're small, a few hundred KB gzipped at most) rather
  than the whole cohort.

Usage:
    python3 scripts/hla_popgen/00_recon_vm.py                      # real run, N=50 people
    python3 scripts/hla_popgen/00_recon_vm.py --limit 20            # smaller/faster sample
    python3 scripts/hla_popgen/00_recon_vm.py --outroot /tmp/hla_fixtures  # against fixtures
"""
import argparse
import gzip
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_REPORT_DIR_NAME = "reports/hla_popgen"

# Files this script checks for, per haplotype -- exactly the set part E of
# reference/IMMUANNOT_GTF_SPEC.md says SHOULD survive Immuannot's own + this project's cleanup.
HAP_ROOT_FILES = ["{hap}.gtf.gz", "{hap}.trimmed.fa"]
HAP_SUBDIR_FILES = ["cds.fa.gz", "mm2.ipd.cds.paf.gz", "mm2.ipd.gen.paf.gz", "gene.filtered.paf"]

ATTR_RE = re.compile(r'(\w+)\s+(?:"((?:[^"\\]|\\.)*)"|(-?\d+))\s*;')

# Attribute keys documented in reference/IMMUANNOT_GTF_SPEC.md, per feature type -- anything seen
# in real output NOT in this set is a red flag worth surfacing loudly.
DOCUMENTED_ATTRS = {
    "gene": {"gene_id", "template_allele", "template_distance", "gene_name"},
    "transcript": {"gene_id", "transcript_id", "gene_name", "consensus", "alleles",
                   "template_warning", "cds_distance", "cds_mut", "pipetide_key"},
    "exon": {"gene_id", "transcript_id", "gene_name", "exon_number"},
    "CDS": {"gene_id", "transcript_id", "gene_name"},
    "UTR": {"gene_id", "transcript_id", "gene_name"},
    "start_codon": {"gene_id", "transcript_id", "gene_name", "codon"},
    "stop_codon": {"gene_id", "transcript_id", "gene_name", "codon"},
}

KIR_GENES = frozenset([
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DP1", "KIR2DS1",
    "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1",
    "KIR3DS1",
])

HEADER_COPY_RE = re.compile(r'^##\s*gene \(copy num ([=>]) (\d+)\):\s*(.*)$')
HEADER_CONTIG_RE = re.compile(r'^##\s*contigs for (\S+):\s*(.*)$')


def fail_if_no_outroot(outroot):
    if not os.path.isdir(outroot):
        sys.exit(
            f"FATAL: --outroot {outroot!r} does not exist.\n"
            f"Remediation: confirm you're on the right VM/session and that the Immuannot pipeline "
            f"has actually produced output there. Expected layout: "
            f"{outroot}/<person_id>/immuannot_output/hap{{1,2}}.gtf.gz"
        )


def parse_attrs(attr_str):
    out = {}
    for m in ATTR_RE.finditer(attr_str):
        key = m.group(1)
        val = m.group(2) if m.group(2) is not None else m.group(3)
        out[key] = val
    return out


def clean_gene_name(raw_name):
    if raw_name.startswith("C4") and len(raw_name) > 2:
        return raw_name[:-1], raw_name[-1]
    return raw_name, None


def gene_id_copy_index(gene_id):
    m = re.match(r'^(.*)\.(\d+)$', gene_id)
    if m:
        return m.group(1), int(m.group(2))
    return gene_id, 1


def scan_person(pid, person_dir, stats):
    """Populate `stats` (a big mutable dict of counters/accumulators) with everything about this
    one person. Kept as in-place mutation rather than returning a huge nested structure per person,
    since this needs to stay cheap across N people and we only care about the aggregate at the end.
    """
    stats["n_people_sampled"] += 1
    person_had_any_hap_output = False

    for hap in ("hap1", "hap2"):
        gtf_path = os.path.join(person_dir, f"{hap}.gtf.gz")
        trimmed_fa_path = os.path.join(person_dir, f"{hap}.trimmed.fa")
        hap_subdir = os.path.join(person_dir, hap)

        stats["file_presence"][f"{hap}.gtf.gz"]["checked"] += 1
        if os.path.exists(gtf_path):
            stats["file_presence"][f"{hap}.gtf.gz"]["present"] += 1
            stats["file_bytes"][f"{hap}.gtf.gz"].append(os.path.getsize(gtf_path))
            person_had_any_hap_output = True
        stats["file_presence"][f"{hap}.trimmed.fa"]["checked"] += 1
        if os.path.exists(trimmed_fa_path):
            stats["file_presence"][f"{hap}.trimmed.fa"]["present"] += 1
            stats["file_bytes"][f"{hap}.trimmed.fa"].append(os.path.getsize(trimmed_fa_path))

        for fname in HAP_SUBDIR_FILES:
            key = f"{hap}/{fname}"
            fpath = os.path.join(hap_subdir, fname)
            stats["file_presence"][key]["checked"] += 1
            if os.path.exists(fpath):
                stats["file_presence"][key]["present"] += 1
                stats["file_bytes"][key].append(os.path.getsize(fpath))

        if not os.path.exists(gtf_path):
            continue

        # Full parse of this hap's GTF -- attribute inventory, gene inventory, contig counts,
        # header cross-check, novelty-depth distribution.
        try:
            _scan_gtf_body(gtf_path, hap, stats)
        except (OSError, EOFError, gzip.BadGzipFile) as e:
            stats["unreadable_gtfs"].append(f"{pid}/{hap}: {e}")

    if not person_had_any_hap_output:
        stats["n_people_zero_output"] += 1


def _scan_gtf_body(gtf_path, hap, stats):
    header_lines = []
    contigs_seen = set()
    body_genes_by_contig = defaultdict(lambda: defaultdict(list))  # gene -> contig -> [copy_idx]

    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                stats["malformed_lines"] += 1
                continue
            contig, _source, feature, start, end, _score, strand, _frame, attr_str = fields[:9]
            contigs_seen.add(contig)
            attrs = parse_attrs(attr_str)

            # Attribute inventory, per feature type.
            fstats = stats["attrs_by_feature"][feature]
            for key, val in attrs.items():
                fstats["keys"][key] += 1
                examples = fstats["examples"][key]
                if len(examples) < 2 and val not in examples:
                    examples.append(val)
                if key not in DOCUMENTED_ATTRS.get(feature, set()):
                    stats["undocumented_attr_keys"].add((feature, key))

            gene_name_raw = attrs.get("gene_name")
            gene_id = attrs.get("gene_id")
            if gene_name_raw and gene_id:
                gene, c4_size = clean_gene_name(gene_name_raw)
                _, copy_idx = gene_id_copy_index(gene_id)
                if feature == "gene":
                    body_genes_by_contig[gene][contig].append(copy_idx)
                    stats["gene_freq"][gene] += 1
                    if gene in KIR_GENES:
                        stats["KIR_GENES_FOUND"].append((hap, contig, gene))

            if feature == "gene":
                td_raw = attrs.get("template_distance")
                if td_raw is not None:
                    try:
                        stats["template_distance_values"].append(int(td_raw))
                        stats["template_distance_by_gene"][gene_name_raw].append(int(td_raw))
                    except ValueError:
                        pass

            if feature == "transcript":
                consensus = attrs.get("consensus")
                if consensus:
                    depth = _novelty_depth(consensus)
                    if depth is not None:
                        stats["novelty_depth_counts"][depth] += 1
                stats["cds_distance_present"] += int("cds_distance" in attrs)
                stats["cds_mut_present"] += int("cds_mut" in attrs)
                # NOTE: the attribute's mere PRESENCE is NOT a warning indicator -- Immuannot
                # writes the literal string template_warning "NA" to mean *no warning* far more
                # often (57.4% of transcript rows) than it omits the attribute (4.6%). Counting
                # presence alone (the old `int("template_warning" in attrs)`) produced the
                # misleading "~95% of calls warned" headline this script used to report.
                # `template_warning_real_present` below counts only genuine warning tokens; the
                # full token breakdown (NA vs absent vs each real token) is also tracked so this
                # recon output can't be misread the same way again (SCHEMA.md's
                # "template_warning policy").
                warn_raw = attrs.get("template_warning")
                if warn_raw is None:
                    stats["template_warning_token_counts"]["<absent>"] += 1
                else:
                    warn_norm = warn_raw.strip().upper()
                    if warn_norm in ("", "NA"):
                        stats["template_warning_token_counts"]["NA"] += 1
                    else:
                        stats["template_warning_real_present"] += 1
                        for tok in warn_raw.split(","):
                            tok = tok.strip()
                            if tok:
                                stats["template_warning_token_counts"][tok] += 1
                stats["n_transcript_rows"] += 1

    stats["n_contigs_per_hapfile"].append(len(contigs_seen))

    # Header cross-check.
    header_info = {"copy0": set(), "copy1": set(), "copyN": set(), "contigs_for_gene": {}}
    for line in header_lines:
        m = HEADER_COPY_RE.match(line.strip())
        if m:
            op, num, genes_str = m.groups()
            genes = {g for g in genes_str.split(",") if g}
            if op == "=" and num == "0":
                header_info["copy0"] |= genes
            elif op == "=" and num == "1":
                header_info["copy1"] |= genes
            elif op == ">":
                header_info["copyN"] |= genes
            continue
        m = HEADER_CONTIG_RE.match(line.strip())
        if m:
            gene, contigs_str = m.groups()
            header_info["contigs_for_gene"][gene] = {c for c in contigs_str.split(",") if c}

    for gene, contig_copies in body_genes_by_contig.items():
        max_copy = max(c for copies in contig_copies.values() for c in copies)
        observed_contigs = set(contig_copies)
        if gene in header_info["copy0"]:
            stats["header_mismatches"].append(f"{gene}: listed copy0 but has body rows")
        if max_copy > 1 and gene not in header_info["copyN"]:
            stats["header_mismatches"].append(f"{gene}: copy_index>1 in body, not in copyN header")
        declared = header_info["contigs_for_gene"].get(gene)
        if declared and declared != observed_contigs:
            stats["header_mismatches"].append(
                f"{gene}: header contigs {sorted(declared)} != body contigs {sorted(observed_contigs)}")
    for gene in header_info["copy1"] | header_info["copyN"]:
        if gene not in body_genes_by_contig:
            stats["header_mismatches"].append(f"{gene}: header says present but 0 body rows")


def _novelty_depth(consensus):
    if consensus == "undetermined" or "*" not in consensus:
        return None
    allele_part = consensus.split("*", 1)[1]
    fields = allele_part.split(":")
    if "new" not in fields:
        return None
    return fields.index("new") + 1


def discover_persons(outroot, limit):
    persons = sorted(
        d for d in os.listdir(outroot)
        if os.path.isdir(os.path.join(outroot, d, "immuannot_output"))
    )
    if not persons:
        sys.exit(f"FATAL: no <person_id>/immuannot_output/ directories found under {outroot!r}.")
    if limit is not None and limit < len(persons):
        step = len(persons) / limit
        persons = [persons[int(i * step)] for i in range(limit)]
    return persons


def pct(n, d):
    return f"{100.0 * n / d:.1f}%" if d else "n/a"


def summarize_bytes(byte_list, n_people_total, n_haps_per_person=2):
    """Extrapolate observed average file size to a full-cohort byte estimate. n_people_total is
    the number the sample rate should be extrapolated to (passed in by the caller as
    --extrapolate-to, default 12000 per the task brief's '~12,000 people')."""
    if not byte_list:
        return {"n_observed": 0, "mean_bytes": None, "extrapolated_total_bytes": None}
    mean_b = sum(byte_list) / len(byte_list)
    return {
        "n_observed": len(byte_list),
        "mean_bytes": round(mean_b, 1),
        "extrapolated_total_bytes": round(mean_b * n_people_total * n_haps_per_person),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--limit", "-n", type=int, default=50,
                    help="Number of people to sample (default 50, spread evenly across the sorted "
                         "person-directory list rather than just the first N, so the sample isn't "
                         "biased toward whatever sorts first).")
    ap.add_argument("--extrapolate-to", type=int, default=12000,
                    help="Full-cohort size to extrapolate byte totals to (default 12000, per the "
                         "project's stated cohort size).")
    ap.add_argument("--report-dir", default=None,
                    help="Where to write the markdown+JSON report. Default: "
                         f"<repo_root>/{DEFAULT_REPORT_DIR_NAME}/ if running from inside the repo, "
                         f"else ./{DEFAULT_REPORT_DIR_NAME}/")
    args = ap.parse_args()

    t0 = time.time()
    fail_if_no_outroot(args.outroot)
    persons = discover_persons(args.outroot, args.limit)
    print(f"Sampling {len(persons)} people from {args.outroot!r} ...", file=sys.stderr)

    stats = {
        "n_people_sampled": 0,
        "n_people_zero_output": 0,
        "file_presence": defaultdict(lambda: Counter()),
        "file_bytes": defaultdict(list),
        "attrs_by_feature": defaultdict(lambda: {"keys": Counter(), "examples": defaultdict(list)}),
        "undocumented_attr_keys": set(),
        "gene_freq": Counter(),
        "KIR_GENES_FOUND": [],
        "template_distance_values": [],
        "template_distance_by_gene": defaultdict(list),
        "novelty_depth_counts": Counter(),
        "cds_distance_present": 0,
        "cds_mut_present": 0,
        "template_warning_real_present": 0,
        "template_warning_token_counts": Counter(),
        "n_transcript_rows": 0,
        "n_contigs_per_hapfile": [],
        "header_mismatches": [],
        "malformed_lines": 0,
        "unreadable_gtfs": [],
    }

    for pid in persons:
        scan_person(pid, os.path.join(args.outroot, pid, "immuannot_output"), stats)

    elapsed = time.time() - t0

    # ---- Build the report ----
    file_presence_report = {}
    for key, counts in stats["file_presence"].items():
        checked, present = counts["checked"], counts["present"]
        byte_summary = summarize_bytes(stats["file_bytes"].get(key, []), args.extrapolate_to)
        file_presence_report[key] = {
            "checked": checked, "present": present, "presence_rate": pct(present, checked),
            **byte_summary,
        }

    attrs_report = {}
    for feature, d in stats["attrs_by_feature"].items():
        attrs_report[feature] = {
            key: {"count": count, "examples": d["examples"][key][:2]}
            for key, count in d["keys"].most_common()
        }

    n_contigs = stats["n_contigs_per_hapfile"]
    n_multi_contig = sum(1 for c in n_contigs if c > 1)

    genes_multi_copy = {}
    # Recompute copy>1 frequency from gene_freq isn't directly tracked per-copy; approximate via
    # header_mismatches is wrong -- instead track directly: re-derive from template_distance_by_gene
    # keys presence is not copy-index-aware, so report what we *can* cheaply say: gene frequency and
    # whether any KIR appeared. Full copy-index census is 01_extract_rich.py's job at cohort scale;
    # recon only needs to sanity-check "does copy_index>1 exist at all in this sample."

    td_vals = stats["template_distance_values"]

    report = {
        "meta": {
            "outroot": args.outroot, "n_people_sampled": stats["n_people_sampled"],
            "elapsed_seconds": round(elapsed, 1),
            "extrapolated_to_n_people": args.extrapolate_to,
        },
        "MOST_IMPORTANT_file_presence_and_bytes": file_presence_report,
        "n_people_zero_output": stats["n_people_zero_output"],
        "gtf_attribute_inventory": attrs_report,
        "undocumented_attr_keys": sorted(f"{feat}.{key}" for feat, key in stats["undocumented_attr_keys"]),
        "gene_inventory": dict(stats["gene_freq"].most_common()),
        "KIR_GENES_FOUND": stats["KIR_GENES_FOUND"],
        "contigs_per_hapfile": {
            "n_hapfiles_seen": len(n_contigs),
            "n_multi_contig": n_multi_contig,
            "fraction_multi_contig": pct(n_multi_contig, len(n_contigs)),
            "max_contigs_seen": max(n_contigs) if n_contigs else None,
        },
        "template_distance": {
            "n": len(td_vals),
            "min": min(td_vals) if td_vals else None,
            "max": max(td_vals) if td_vals else None,
            "mean": round(sum(td_vals) / len(td_vals), 2) if td_vals else None,
            "n_zero_exact": sum(1 for v in td_vals if v == 0),
        },
        "conditional_field_presence": {
            "n_transcript_rows": stats["n_transcript_rows"],
            "cds_distance_present": stats["cds_distance_present"],
            "cds_distance_present_rate": pct(stats["cds_distance_present"], stats["n_transcript_rows"]),
            "cds_mut_present": stats["cds_mut_present"],
            "cds_mut_present_rate": pct(stats["cds_mut_present"], stats["n_transcript_rows"]),
            # "real" = a genuine warning token, i.e. NOT the literal string "NA" (Immuannot's
            # spelling of "no warning") and NOT attribute-absent. See
            # "template_warning_token_breakdown" below for the NA-vs-absent-vs-real split that
            # makes this unambiguous.
            "template_warning_real_present": stats["template_warning_real_present"],
            "template_warning_real_present_rate": pct(stats["template_warning_real_present"], stats["n_transcript_rows"]),
        },
        "template_warning_token_breakdown": {
            "n_transcript_rows": stats["n_transcript_rows"],
            "counts": dict(stats["template_warning_token_counts"].most_common()),
            "rates": {
                tok: pct(count, stats["n_transcript_rows"])
                for tok, count in stats["template_warning_token_counts"].most_common()
            },
        },
        "novelty_depth_distribution": dict(stats["novelty_depth_counts"]),
        "header_body_mismatches": {
            "n_mismatches": len(stats["header_mismatches"]),
            "examples": stats["header_mismatches"][:20],
        },
        "malformed_lines": stats["malformed_lines"],
        "unreadable_gtfs": stats["unreadable_gtfs"],
    }

    # ---- Write outputs ----
    report_dir = args.report_dir or DEFAULT_REPORT_DIR_NAME
    os.makedirs(report_dir, exist_ok=True)
    json_path = os.path.join(report_dir, "recon_report.json")
    md_path = os.path.join(report_dir, "recon_report.md")
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    with open(md_path, "w") as f:
        f.write(_render_markdown(report))

    print(f"\nWrote {json_path} and {md_path} in {elapsed:.1f}s.", file=sys.stderr)
    if elapsed > 120:
        print(f"WARNING: took {elapsed:.0f}s, over the ~2 minute target -- consider a smaller "
              f"--limit.", file=sys.stderr)

    _print_compact_summary(report)


def _render_markdown(report):
    lines = ["# hla_popgen VM recon report\n"]
    lines.append(f"Sampled **{report['meta']['n_people_sampled']}** people from "
                 f"`{report['meta']['outroot']}` in {report['meta']['elapsed_seconds']}s.\n")

    lines.append("## Most important: file presence + extrapolated bytes\n")
    lines.append("| file | present/checked | rate | mean bytes | extrapolated total "
                 f"(n={report['meta']['extrapolated_to_n_people']}) |")
    lines.append("|---|---|---|---|---|")
    for key, d in report["MOST_IMPORTANT_file_presence_and_bytes"].items():
        lines.append(f"| `{key}` | {d['present']}/{d['checked']} | {d['presence_rate']} | "
                     f"{d['mean_bytes']} | {d['extrapolated_total_bytes']} |")
    lines.append("")

    cds_row = report["MOST_IMPORTANT_file_presence_and_bytes"].get("hap1/cds.fa.gz", {})
    if cds_row.get("present", 0) == 0:
        lines.append("**cds.fa.gz WAS NOT FOUND IN THE SAMPLE.** The entire novel-allele "
                     "sequence-level workstream (Table 3, planned) depends on this file. "
                     "Investigate before planning further work that assumes it exists.\n")

    lines.append(f"## KIR gene check\n")
    if report["KIR_GENES_FOUND"]:
        lines.append(f"**LOUD WARNING: {len(report['KIR_GENES_FOUND'])} KIR gene occurrence(s) "
                     f"found** -- spec says this must never happen (chr19, outside the chr6 trim "
                     f"window). Indicates a trim bug. Examples: "
                     f"{report['KIR_GENES_FOUND'][:10]}\n")
    else:
        lines.append("No KIR genes found in sample -- consistent with spec expectation.\n")

    lines.append("## Undocumented GTF attribute keys\n")
    if report["undocumented_attr_keys"]:
        lines.append(f"**{len(report['undocumented_attr_keys'])} undocumented key(s) found:** "
                     + ", ".join(f"`{k}`" for k in report["undocumented_attr_keys"]) + "\n")
    else:
        lines.append("None -- every attribute key seen matches reference/IMMUANNOT_GTF_SPEC.md.\n")

    lines.append("## Contigs per hap file\n")
    cp = report["contigs_per_hapfile"]
    lines.append(f"{cp['n_multi_contig']}/{cp['n_hapfiles_seen']} hap files "
                 f"({cp['fraction_multi_contig']}) span more than one contig "
                 f"(max seen: {cp['max_contigs_seen']}). This is the ceiling on cis-pairing loss.\n")

    lines.append("## template_distance distribution\n")
    td = report["template_distance"]
    lines.append(f"n={td['n']}, min={td['min']}, max={td['max']}, mean={td['mean']}, "
                 f"exact (0) = {td['n_zero_exact']}\n")

    lines.append("## Conditional field presence\n")
    cf = report["conditional_field_presence"]
    lines.append(f"- cds_distance present: {cf['cds_distance_present']}/{cf['n_transcript_rows']} "
                 f"({cf['cds_distance_present_rate']})")
    lines.append(f"- cds_mut present: {cf['cds_mut_present']}/{cf['n_transcript_rows']} "
                 f"({cf['cds_mut_present_rate']})")
    lines.append(f"- template_warning REAL (excludes Immuannot's literal \"NA\" and absent "
                 f"attribute, both of which mean *clean*): {cf['template_warning_real_present']}/"
                 f"{cf['n_transcript_rows']} ({cf['template_warning_real_present_rate']})\n")

    lines.append("## template_warning token breakdown (NA and <absent> are both CLEAN)\n")
    twb = report["template_warning_token_breakdown"]
    for tok, count in twb["counts"].items():
        lines.append(f"- {tok}: {count}/{twb['n_transcript_rows']} ({twb['rates'][tok]})")
    lines.append("")

    lines.append("## Novelty depth distribution (which field 'new' lands in)\n")
    lines.append(str(report["novelty_depth_distribution"]) + "\n")

    lines.append("## Header vs body cross-check\n")
    hm = report["header_body_mismatches"]
    if hm["n_mismatches"]:
        lines.append(f"**{hm['n_mismatches']} mismatch(es) found -- indicates a parser bug, "
                     f"investigate before trusting the header as ground truth.**\n")
        for ex in hm["examples"]:
            lines.append(f"- {ex}")
    else:
        lines.append("No mismatches -- header claims agree with the parsed body.\n")

    lines.append("\n## Full gene inventory\n")
    for gene, n in report["gene_inventory"].items():
        lines.append(f"- {gene}: {n}")

    return "\n".join(lines) + "\n"


def _print_compact_summary(report):
    """Print a compact, paste-back-friendly summary to stdout (per the task brief: 'print a
    compact summary to stdout that Marc can paste back in one block')."""
    print("\n" + "=" * 70)
    print("COMPACT RECON SUMMARY (paste this back)")
    print("=" * 70)
    print(json.dumps({
        "n_people_sampled": report["meta"]["n_people_sampled"],
        "elapsed_seconds": report["meta"]["elapsed_seconds"],
        "n_people_zero_output": report["n_people_zero_output"],
        "file_presence_rates": {k: v["presence_rate"]
                                 for k, v in report["MOST_IMPORTANT_file_presence_and_bytes"].items()},
        "file_extrapolated_bytes": {k: v["extrapolated_total_bytes"]
                                     for k, v in report["MOST_IMPORTANT_file_presence_and_bytes"].items()},
        "undocumented_attr_keys": report["undocumented_attr_keys"],
        "KIR_genes_found_count": len(report["KIR_GENES_FOUND"]),
        "n_genes_in_inventory": len(report["gene_inventory"]),
        "fraction_multi_contig": report["contigs_per_hapfile"]["fraction_multi_contig"],
        "template_distance": report["template_distance"],
        "conditional_field_presence_rates": {
            k: v for k, v in report["conditional_field_presence"].items() if k.endswith("_rate")
        },
        "novelty_depth_distribution": report["novelty_depth_distribution"],
        "header_body_mismatch_count": report["header_body_mismatches"]["n_mismatches"],
        "malformed_lines": report["malformed_lines"],
        "unreadable_gtfs": report["unreadable_gtfs"],
    }, indent=2))
    print("=" * 70)


if __name__ == "__main__":
    main()
