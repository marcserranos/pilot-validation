"""105 raw-read evidence package for candidate novel HLA alleles (S09-12; VM side, aggregate output only).

For every substitution signature that is recurrent and NOT a named allele in the current IPD release (91 --ipd-nuc output), plus every nonsense / start-loss
signature in HLA-E/F/G, find the exact GRCh38 position (GENCODE CDS model; the strand-aware base is checked against the assembly-derived ref/obs bases, so a wrong
mapping is detected), then pile up the raw LONG-READ alignments (the AoU grch38 BAMs, independent of the assembly the call came from) in up to --carriers carriers
and --controls random non-carriers:
  carriers must show the alternate base in ~half the reads (heterozygous) or ~all (homozygous); non-carriers must not.
Output (aggregates only): novel_allele_evidence_pub.tsv  (counts < 20 -> '<20'; per-signature verdict ALL / MOST / SOME / NONE of examined carriers supported,
pooled read fractions, control false-support) and ..._raw.private.tsv (mode 600, exact). Per-person pileups are checkpointed to evidence105.private.jsonl (resumable).
"""
import argparse
import csv
import gzip
import json
import os
import random
import re
import subprocess
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
GENES_OK = {"A", "B", "C", "E", "F", "G", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"}


def log(m):
    print("[105] " + m, file=sys.stderr, flush=True)


def supp(n, floor=20):
    return "<%d" % floor if 0 < n < floor else str(n)


def load_cds_models(gtf_gz, genes):
    """gene symbol -> dict(chrom, strand, exons=[(start,end)] in transcript order); transcript with the longest CDS wins."""
    tx = defaultdict(lambda: {"chrom": None, "strand": None, "cds": [], "gene": None})
    want = set(genes)
    with gzip.open(gtf_gz, "rt") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            x = ln.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] != "CDS":
                continue
            m = re.search(r'gene_name "([^"]+)"', x[8])
            if not m or m.group(1) not in want:
                continue
            t = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
            d = tx[t]
            d["chrom"], d["strand"], d["gene"] = x[0], x[6], m.group(1)
            d["cds"].append((int(x[3]), int(x[4])))
    best = {}
    for t, d in tx.items():
        L = sum(e - s + 1 for s, e in d["cds"])
        if d["gene"] not in best or L > best[d["gene"]][0]:
            ex = sorted(d["cds"], reverse=(d["strand"] == "-"))
            best[d["gene"]] = (L, {"chrom": d["chrom"], "strand": d["strand"], "exons": ex, "transcript": t, "cds_len": L})
    return {g: v[1] for g, v in best.items()}


def cds_to_genomic(model, pos1):
    """1-based CDS position -> 1-based genomic coordinate (strand-aware), or None if beyond the CDS."""
    off = pos1 - 1
    for s, e in model["exons"]:
        n = e - s + 1
        if off < n:
            return (s + off) if model["strand"] == "+" else (e - off)
        off -= n
    return None


def genomic_bases(model, ref, obs):
    """Bases as they appear on the + strand of the genome."""
    return (ref, obs) if model["strand"] == "+" else (ref.translate(COMP), obs.translate(COMP))


def parse_pileup_bases(bases):
    """Count A/C/G/T/ref/del from an samtools mpileup base string."""
    s = re.sub(r"\^.", "", bases).replace("$", "")
    out = {"A": 0, "C": 0, "G": 0, "T": 0, "ref": 0, "del": 0}
    i = 0
    while i < len(s):
        c = s[i]
        if c in "+-":
            j = i + 1
            k = j
            while k < len(s) and s[k].isdigit():
                k += 1
            n = int(s[j:k] or 0)
            i = k + n
            continue
        if c in ".,":
            out["ref"] += 1
        elif c == "*":
            out["del"] += 1
        elif c.upper() in "ACGT":
            out[c.upper()] += 1
        i += 1
    return out


def pileup(samtools, bam, fasta, chrom, pos, mapq):
    r = subprocess.run([samtools, "mpileup", "-r", "%s:%d-%d" % (chrom, pos, pos), "-q", str(mapq), "-Q", "0", "-B", "-d", "100000", "-f", fasta, bam],
                       capture_output=True, text=True, timeout=300)
    for ln in r.stdout.splitlines():
        x = ln.split("\t")
        if len(x) >= 5 and int(x[1]) == pos:
            c = parse_pileup_bases(x[4])
            refb = x[2].upper()
            c[refb] = c.get(refb, 0) + c["ref"]
            return c, refb
    return None, None


def read_manifest_bams(path, mount_from, mount_to):
    bam = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b = r.get("grch38_bam")
            if b in (None, "", "NA"):
                continue
            p = b.replace(mount_from, mount_to)
            bam.setdefault(r["research_id"], []).append((p, r.get("platform", "NA")))
    return bam


def resolve(bams, pid):
    for p, plat in bams.get(pid, []):
        if os.path.exists(p) and os.path.exists(p + ".bai"):
            return p, plat
    return None, None


def run(a):
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    src = os.path.expanduser(a.tables)
    raw = list(csv.DictReader(open(os.path.join(src, "recurrence_raw.private.tsv")), delimiter="\t"))
    carriers = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(src, "carriers.private.tsv")), delimiter="\t"):
        carriers[(r["gene"], r["ref_allele"], r["cds_events"])][r["person_id"]] = r["haps"].split(",")
    cands = []
    for r in raw:
        g = r["gene"]
        if g not in GENES_OK or r["ipd_match"] != "none":
            continue
        lof = g in ("E", "F", "G") and ("*" in r["protein"] or re.match(r"^M1", r["protein"]))
        suffix = bool(re.search(r"[A-Z]$", r["ref_allele"].split("*")[-1]))
        if int(r["persons"]) >= 2 and not suffix or lof:
            cands.append(r)
    log("candidate signatures: %d" % len(cands))
    models = load_cds_models(os.path.expanduser(a.gencode), ["HLA-" + g for g in GENES_OK])
    membership = {}
    with open(os.path.expanduser(a.membership)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("in_lr", "").lower() in ("true", "1"):
                membership[r["person_id"]] = (r.get("ancestry_pred") or "NA").upper()
    bams = read_manifest_bams(os.path.expanduser(a.manifest), "gs://vwb-aou-datasets-controlled/", os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/"))
    jp = os.path.join(out, "evidence105.private.jsonl")
    done = set()
    if os.path.exists(jp):
        done = {(d["sig"], d["person_id"]) for d in map(json.loads, open(jp))}
    rnd = random.Random(21)
    allp = sorted(membership)
    sites = []
    for r in cands:
        key = (r["gene"], r["ref_allele"], r["cds_events"])
        model = models.get("HLA-" + r["gene"])
        evs = [e for e in r["cds_events"].split(";") if e]
        if model is None or not evs:
            continue
        m = re.match(r"(\d+):([ACGT])>([ACGT])", evs[0])      # first substitution represents the signature
        if not m:
            continue
        pos1, ref, obs = int(m.group(1)), m.group(2), m.group(3)
        gpos = cds_to_genomic(model, pos1)
        if gpos is None:
            continue
        gref, gobs = genomic_bases(model, ref, obs)
        cs = sorted(carriers[key])
        rnd.shuffle(cs)
        cs = cs[:a.carriers]
        pool = [p for p in allp if p not in carriers[key]]
        ctl = rnd.sample(pool, min(a.controls, len(pool)))
        sid = "%s|%s|%s" % (r["gene"], r["protein"], r["cds_events"])
        sites.append({"sig": sid, "gene": r["gene"], "protein": r["protein"], "chrom": model["chrom"], "pos": gpos, "gref": gref, "gobs": gobs,
                      "carriers": cs, "controls": ctl, "haps": {p: len(carriers[key][p]) for p in cs}, "persons": int(r["persons"]),
                      "unrelated": int(r["unrelated"]), "anc": r["anc_unrelated"]})
    log("sites to pile up: %d" % len(sites))
    n = 0
    for s in sites:
        for role, plist in (("carrier", s["carriers"]), ("control", s["controls"])):
            for p in plist:
                if (s["sig"], p) in done:
                    continue
                bam, plat = resolve(bams, p)
                rec = {"sig": s["sig"], "person_id": p, "role": role, "platform": plat, "haps": s["haps"].get(p)}
                if bam is None:
                    rec["err"] = "no_bam"
                else:
                    try:
                        c, refb = pileup(os.path.expanduser(a.samtools), bam, os.path.expanduser(a.fasta), s["chrom"], s["pos"], a.mapq)
                        rec.update({"counts": c, "hg38_ref": refb} if c else {"err": "no_coverage"})
                    except Exception as e:  # noqa: BLE001
                        rec["err"] = type(e).__name__
                with open(jp, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                os.chmod(jp, 0o600)
                n += 1
                if n % 25 == 0:
                    log("pileups done: %d" % n)
    # summarise
    recs = defaultdict(list)
    for d in map(json.loads, open(jp)):
        recs[d["sig"]].append(d)
    rows = []
    for s in sites:
        rs = recs.get(s["sig"], [])
        row = {"gene": s["gene"], "protein": s["protein"], "chrom": s["chrom"], "pos": s["pos"], "persons": s["persons"], "unrelated": s["unrelated"], "anc_unrelated": s["anc"]}
        sup_c = sup_ctl = ex_c = ex_ctl = 0
        pooled_alt = pooled_tot = 0
        wrongbase = 0
        for d in rs:
            c = d.get("counts")
            if not c:
                continue
            alt = c.get(s["gobs"], 0)
            other = sum(c.get(b, 0) for b in "ACGT" if b not in (s["gobs"],))
            tot = alt + other
            if tot < a.min_depth:
                continue
            supported = alt >= 3 and alt / tot >= 0.2
            if d["role"] == "carrier":
                ex_c += 1
                sup_c += supported
                pooled_alt += alt
                pooled_tot += tot
                if d.get("hg38_ref") not in (s["gref"], s["gobs"]):
                    wrongbase += 1
            else:
                ex_ctl += 1
                sup_ctl += supported
        verdict = "NONE" if ex_c == 0 or sup_c == 0 else ("ALL" if sup_c == ex_c else ("MOST" if sup_c / ex_c >= 0.75 else "SOME"))
        row.update({"carriers_examined": ex_c, "carriers_supported": sup_c, "verdict": verdict, "pooled_alt_fraction": (pooled_alt / pooled_tot) if pooled_tot else "",
                    "controls_examined": ex_ctl, "controls_supported": sup_ctl, "coord_check_fail": wrongbase})
        rows.append(row)
    cols = ["gene", "protein", "chrom", "pos", "persons", "unrelated", "anc_unrelated", "carriers_examined", "carriers_supported", "verdict", "pooled_alt_fraction",
            "controls_examined", "controls_supported", "coord_check_fail"]
    for suffix, pub in (("pub", True), ("raw.private", False)):
        pth = os.path.join(out, "novel_allele_evidence_%s.tsv" % suffix)
        with open(pth, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                cells = []
                for c in cols:
                    v = r[c]
                    if pub and c in ("persons", "unrelated", "carriers_examined", "carriers_supported", "controls_examined", "controls_supported", "coord_check_fail"):
                        v = supp(int(v))
                    elif pub and c in ("anc_unrelated", "pos", "chrom"):
                        v = v if c != "anc_unrelated" else ";".join(p.split("=")[0] + "=" + supp(int(p.split("=")[1])) for p in str(v).split(";") if "=" in p)
                    elif isinstance(v, float):
                        v = "%.3f" % v
                    cells.append(str(v))
                f.write("\t".join(cells) + "\n")
        if not pub:
            os.chmod(pth, 0o600)
    log("wrote %d signature rows" % len(rows))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--tables", default="~/s08/r366")
    p.add_argument("--out-dir", default="~/s09/e105")
    p.add_argument("--gencode", default="~/ref/gencode.v46.basic.annotation.gtf.gz")
    p.add_argument("--fasta", default="~/ref/hg38_ucsc.fa")
    p.add_argument("--manifest", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--samtools", default="~/repos/pilot-validation-main/.pixi/envs/spechla/bin/samtools")
    p.add_argument("--carriers", type=int, default=8)
    p.add_argument("--controls", type=int, default=8)
    p.add_argument("--mapq", type=int, default=20)
    p.add_argument("--min-depth", type=int, default=8)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
