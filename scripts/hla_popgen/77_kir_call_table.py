#!/usr/bin/env python3
"""S06 77: per-call KIR table for the call-quality replication (VM; reads only the existing Immuannot GTFs, no re-annotation).

One row per (person, hap, contig, gene, copy_index), read from ~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz with
44's own parser (parse_hap_gtf_full, which already applies 41's novelty tiers) plus one light pass for the fields 44 does not expose
(gene-row span / strand / template_distance, transcript `alleles` -> n_tied).

Labels (KIR alleles are digit strings, not colon fields; 61's protein_field is reused):
  label3  protein field  'KIR2DL1*0030204' -> 'KIR2DL1*003'   valid for tiers known / novel_genomic_known_cds / novel_cds_synonymous
  label5  coding field   'KIR2DL1*0030204' -> 'KIR2DL1*00302' valid for tiers known / novel_genomic_known_cds (CDS identical to a named allele)
  novel_protein and undetermined calls have no label (they are the separate "novel" stream in 77c).
Haplotype completeness uses 51's own rule (all calls on one contig and both 3DL3 and 3DL2 present): hap_class in
complete / one_contig_open_end / multi_contig / empty.  Contig names are extraction ranges `ctg:S-E`; the distance of a gene to the
two ends of that range (dist_start, dist_end) is kept as a proxy for "the gene sits at a contig end".
n_called = number of haplotypes (0-2) carrying the gene (copy_index 1 counted once per hap); person_n_haps = haps with any KIR call.

Output: INTERNAL ~/s06/internal/77/{kir_calls,kir_haps}.INTERNAL.tsv (+ parts/ checkpoint chunks), aggregate counts
(1-19 -> "<20") in ~/s06/results/77/.

    cd ~/repos/pilot-validation-main && pixi run -e spechla -- python3 ~/repos/pv-s06/scripts/hla_popgen/77_kir_call_table.py --workers 4 --resume
    ... --limit 40 --out-root /path/synthetic   # smoke
"""
import argparse
import glob
import gzip
import importlib.util
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

OUT = os.path.expanduser("~/s06/results/77/")
INT = os.path.expanduser("~/s06/internal/77/")
CHUNK = 200
LABEL3_TIERS = ("known", "novel_genomic_known_cds", "novel_cds_synonymous")
LABEL5_TIERS = ("known", "novel_genomic_known_cds")
RANGE_RE = re.compile(r"^(.*):(\d+)-(\d+)$")
TD_RE = re.compile(r"template_distance (\d+)")
ALLELES_RE = re.compile(r'alleles "([^"]*)"')
GID_RE = re.compile(r'gene_id "([^"]+)"')
_MODS = {}


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def mods():
    """Lazy (per process) import of the KIR scripts reused here, never re-derived."""
    if not _MODS:
        _MODS["m41"] = _load("kir41_for_77", "41_kir_pilot.py")
        _MODS["m44"] = _load("kir44_for_77", "44_kir_recurrence_saturation.py")
        _MODS["m51"] = _load("kir51_for_77", "51_kir_structure_qc.py")
    return _MODS


def protein_field(consensus):
    """61's protein_field verbatim: 'KIR2DL1*0030204' -> 'KIR2DL1*003'; 'KIR3DL2*00201new' -> 'KIR3DL2*002'."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    gene, rest = consensus.split("*", 1)
    digits = re.match(r"\d+", rest.replace("new", ""))
    if not digits:
        return None
    return f"{gene}*{digits.group(0)[:3]}"


def coding_field(consensus):
    """First five digits (protein + synonymous-coding field): 'KIR2DL1*0030204' -> 'KIR2DL1*00302'."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    gene, rest = consensus.split("*", 1)
    digits = re.match(r"\d+", rest.replace("new", ""))
    if not digits or len(digits.group(0)) < 5:
        return None
    return f"{gene}*{digits.group(0)[:5]}"


def labels(consensus, tier):
    l3 = protein_field(consensus) if tier in LABEL3_TIERS else None
    l5 = coding_field(consensus) if tier in LABEL5_TIERS else None
    return l3, l5


def parse_extra(gtf_path):
    """gene_id -> span / strand / template_distance (gene row) and `alleles` (transcript row)."""
    ex = {}
    if not gtf_path or not os.path.exists(gtf_path):
        return ex
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] not in ("gene", "transcript"):
                continue
            m = GID_RE.search(x[8])
            if not m:
                continue
            d = ex.setdefault(m.group(1), {})
            if x[2] == "gene":
                try:
                    d["start"], d["end"] = int(x[3]), int(x[4])
                except ValueError:
                    continue
                d["strand"] = x[6]
                t = TD_RE.search(x[8])
                if t:
                    d["template_distance"] = int(t.group(1))
            else:
                a = ALLELES_RE.search(x[8])
                if a:
                    d["n_tied"] = len([s for s in a.group(1).split(",") if s])
                if "template_distance" not in d:
                    t = TD_RE.search(x[8])
                    if t:
                        d["template_distance"] = int(t.group(1))
    return ex


def contig_len(contig):
    m = RANGE_RE.match(str(contig))
    return int(m.group(3)) - int(m.group(2)) + 1 if m else None


def person_tables(pid, root):
    """-> (call rows, hap rows) for one person."""
    M = mods()
    kir_genes = set(M["m41"].KIR_GENES)
    calls, haps = [], []
    for hap in ("hap1", "hap2"):
        gtf = os.path.join(root, str(pid), "immuannot_output", hap + ".gtf.gz")
        if not os.path.exists(gtf):
            continue
        rows = [r for r in M["m44"].parse_hap_gtf_full(gtf, M["m41"]) if r["gene"] in kir_genes]
        ex = parse_extra(gtf)
        tup = [(r["contig"], r["gene"], (ex.get(r["gene_id"]) or {}).get("start", 0),
                (ex.get(r["gene_id"]) or {}).get("end", 0)) for r in rows]
        hap_class = M["m51"].classify_hap(tup)
        contigs = sorted({r["contig"] for r in rows})
        haps.append(dict(person_id=str(pid), hap=hap, n_calls=len(rows), n_genes=len({r["gene"] for r in rows}),
                         n_contigs=len(contigs), hap_class=hap_class))
        for r in rows:
            e = ex.get(r["gene_id"]) or {}
            s, en = e.get("start"), e.get("end")
            clen = contig_len(r["contig"])
            l3, l5 = labels(r["consensus"], r["novelty_tier"])
            calls.append(dict(
                person_id=str(pid), hap=hap, contig=r["contig"], gene=r["gene"], copy_index=r["copy_index"],
                consensus=r["consensus"], novelty_tier=r["novelty_tier"], label3=l3, label5=l5,
                gene_start=s, gene_end=en, strand=e.get("strand"), template_distance=e.get("template_distance"),
                cds_distance=r["cds_distance"], cds_mut=r["cds_mut"], template_warning=r["template_warning"],
                n_tied=e.get("n_tied"), hap_class=hap_class, n_contigs_hap=len(contigs), contig_len=clen,
                dist_start=(s - 1) if s is not None else None,
                dist_end=(clen - en) if (clen is not None and en is not None) else None))
    return calls, haps


def add_person_level(calls, haps):
    """n_called per person-gene (haps carrying the gene) and person_n_haps."""
    if calls.empty:
        return calls
    c1 = calls[calls["copy_index"] == 1]
    nc = c1.groupby(["person_id", "gene"])["hap"].nunique().rename("n_called")
    calls = calls.join(nc, on=["person_id", "gene"])
    nh = haps[haps["n_calls"] > 0].groupby("person_id")["hap"].nunique().rename("person_n_haps")
    return calls.join(nh, on="person_id")


def run_chunk(args):
    pids, root = args
    calls, haps = [], []
    for p in pids:
        c, h = person_tables(p, root)
        calls += c
        haps += h
    return pd.DataFrame(calls), pd.DataFrame(haps)


def aggregates(calls, haps):
    t = calls.groupby(["gene", "novelty_tier"]).size().rename("n").reset_index()
    t["n"] = [C.mask_n(v) for v in t["n"]]
    C.write_tsv(t, OUT + "gene_by_novelty_tier.tsv")
    h = haps.groupby("hap_class").size().rename("n_haps").reset_index()
    h["n_haps"] = [C.mask_n(v) for v in h["n_haps"]]
    C.write_tsv(h, OUT + "hap_class_counts.tsv")
    c1 = calls[calls["copy_index"] == 1].drop_duplicates(["person_id", "gene"])
    g = c1.groupby(["gene", "n_called"]).size().rename("n_people").reset_index()
    g["n_people"] = [C.mask_n(v) for v in g["n_people"]]
    C.write_tsv(g, OUT + "gene_n_called.tsv")
    lab = calls.assign(has_label3=calls["label3"].notna()).groupby(["gene", "has_label3"]).size().rename("n").reset_index()
    lab["n"] = [C.mask_n(v) for v in lab["n"]]
    C.write_tsv(lab, OUT + "gene_label3_available.tsv")


def main(argv=None):
    global OUT, INT
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", default=C.DEFAULT_KIR_OUTROOT)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--resume", action="store_true", help="reuse finished chunks in the internal parts/ folder")
    ap.add_argument("--int-dir", default=INT)
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    OUT, INT = a.out_dir, a.int_dir
    for d in (OUT, INT, os.path.join(INT, "parts")):
        os.makedirs(d, exist_ok=True)
    pids = sorted(d for d in os.listdir(a.out_root) if os.path.isdir(os.path.join(a.out_root, d)))
    if a.limit:
        pids = pids[:a.limit]
    chunks = [pids[i:i + CHUNK] for i in range(0, len(pids), CHUNK)]
    C.log(f"{len(pids)} people, {len(chunks)} chunks")
    todo = []
    for k, ch in enumerate(chunks):
        part = os.path.join(INT, "parts", f"chunk_{k:05d}")
        if a.resume and os.path.exists(part + ".calls.tsv") and os.path.exists(part + ".haps.tsv"):
            continue
        todo.append((k, ch))
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [(k, ex.submit(run_chunk, (ch, a.out_root))) for k, ch in todo]
        for k, fu in futs:
            c, h = fu.result()
            part = os.path.join(INT, "parts", f"chunk_{k:05d}")
            c.to_csv(part + ".calls.tsv.tmp", sep="\t", index=False)
            h.to_csv(part + ".haps.tsv.tmp", sep="\t", index=False)
            os.replace(part + ".haps.tsv.tmp", part + ".haps.tsv")
            os.replace(part + ".calls.tsv.tmp", part + ".calls.tsv")
            C.log(f"chunk {k + 1}/{len(chunks)} calls {len(c)} haps {len(h)}")
    keep = {f"chunk_{k:05d}" for k in range(len(chunks))}
    cf = [f for f in sorted(glob.glob(os.path.join(INT, "parts", "chunk_*.calls.tsv"))) if os.path.basename(f)[:11] in keep]
    hf = [f.replace(".calls.", ".haps.") for f in cf]
    dt = dict(person_id=str, contig=str, hap=str)
    calls = pd.concat([pd.read_csv(f, sep="\t", dtype=dt) for f in cf], ignore_index=True)
    haps = pd.concat([pd.read_csv(f, sep="\t", dtype=dt) for f in hf], ignore_index=True)
    calls = add_person_level(calls, haps)
    for df, nm in ((calls, "kir_calls"), (haps, "kir_haps")):
        p = os.path.join(INT, nm + ".INTERNAL.tsv")
        with open(p, "w") as f:
            f.write("# INTERNAL: person-level, never leaves the VM\n")
        df.to_csv(p, sep="\t", index=False, mode="a")
        C.log(f"wrote {p} ({len(df)} rows)")
    aggregates(calls, haps)
    C.log("DONE 77")


if __name__ == "__main__":
    main()
