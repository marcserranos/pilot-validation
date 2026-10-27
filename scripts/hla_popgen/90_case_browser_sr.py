"""90 sr (VM side): do the short reads support the long-read variants in a case?

For each selected case (person, haplotype, gene, copy; from cases.json + the private id map) the step
  1. pulls the person's srWGS reads over the gene's hg38 locus (samtools on the mounted CRAM, one small region),
  2. realigns them (minimap2 -ax sr, single-end) to the CALLED haplotype's own gene sequence, so a read carrying the
     other haplotype's base shows up as a mismatch there,
  3. counts, at every difference of the case, reads carrying the called base / the other haplotype's base / anything else,
     with an indel-aware rule, and gives a verdict (sr_verdict in 90_case_browser.py),
  4. keeps a small pileup (window of the called CDS, other-haplotype row, up to ~20 reads).
Output: one JSON line per case in <out-dir>/sr.jsonl (checkpointed, resumable); `build --sr` merges it into the viewer.
Read-only on shared paths; temp files live in a private temp dir that is removed afterwards.
"""
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict

FLANK = 300          # bp of haplotype sequence kept around the gene for read placement
LOCUS_PAD = 1000     # bp around the hg38 gene locus when pulling reads
WINDOW = 31          # CDS bases shown in a pileup
MAX_SITES = 6
MAX_READS_SHOWN = 20
MOUNT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/")
COMP = str.maketrans("ACGTN", "TGCAN")


def log(m):
    print(m, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ coordinates
def cds_contig_positions(call, n):
    """1-based contig coordinate of every base of the observed CDS (strand-aware; the stop codon, which the GTF CDS rows
    leave out, is extended past the last segment)."""
    minus = call["strand"] == "-"
    pos = []
    for s, e in sorted(call["cds"], reverse=minus):
        pos.extend(range(e, s - 1, -1) if minus else range(s, e + 1))
    if pos and n > len(pos):
        step, last = (-1 if minus else 1), pos[-1]
        pos.extend(last + step * (k + 1) for k in range(n - len(pos)))
    return pos[:n]


# ------------------------------------------------------------------ reads
class Read(object):
    """One aligned read as three lookups over contig coordinates: base at a position, deleted positions, and
    insertions (keyed by the reference position they follow)."""
    __slots__ = ("name", "strand", "mapq", "start", "end", "base", "dele", "ins")

    def __init__(self, name, strand, mapq, start, end, base, dele, ins):
        self.name, self.strand, self.mapq, self.start, self.end = name, strand, mapq, start, end
        self.base, self.dele, self.ins = base, dele, ins


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_sam_line(line, ref_start):
    """SAM line -> Read in CONTIG coordinates (reference = contig slice starting at 1-based `ref_start`), or None for
    unmapped / secondary / supplementary records."""
    x = line.rstrip("\n").split("\t")
    if len(x) < 11 or line.startswith("@"):
        return None
    flag = int(x[1])
    if flag & (4 | 256 | 2048) or x[5] == "*":
        return None
    p = ref_start + int(x[3]) - 1
    seq, base, dele, ins, q = x[9], {}, set(), {}, 0
    for n, op in CIGAR_RE.findall(x[5]):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                base[p + k] = seq[q + k]
            p, q = p + n, q + n
        elif op == "I":
            ins[p - 1] = ins.get(p - 1, 0) + n
            q += n
        elif op in "DN":
            if op == "D":
                dele.update(range(p, p + n))
            p += n
        elif op == "S":
            q += n
    return Read(x[0], "-" if flag & 16 else "+", int(x[4]), ref_start + int(x[3]) - 1, p - 1, base, dele, ins)


# ------------------------------------------------------------------ one site
def _other_map(called, other, banded_align):
    """called CDS idx -> other-haplotype base ('-' = gap there), and the set of called idx with a gap run after them
    (the other haplotype has extra bases there)."""
    ra, rb = banded_align(other, called)       # ra = other row, rb = called row
    ob, extra, ci = {}, set(), -1
    for a, b in zip(ra, rb):
        if b != "-":
            ci += 1
            ob[ci] = a
        elif a != "-":
            extra.add(ci)                       # other has a base the called haplotype lacks, after called idx ci
    return ob, extra


def classify_site(ev, reads, pos, called, omap, strand):
    """(obs, alt, oth, het, other_base) read counts for one event. `pos` = contig coord per called CDS idx;
    `omap` = (other base per called idx, idx set after which the other hap has extra bases) or None."""
    ob, extra = omap if omap else ({}, set())
    minus = strand == "-"
    fix = (lambda b: b.translate(COMP)) if minus else (lambda b: b)
    obs = alt = oth = 0
    if ev["type"] == "snv":
        i = ev["qpos"] - 1
        p, want = pos[i], called[i]
        ob_base = ob.get(i)
        het = ob_base is not None and ob_base != want
        for r in reads:
            if not (r.start <= p <= r.end):
                continue
            b = None if p in r.dele else r.base.get(p)
            b = None if b is None else fix(b)
            if b == want:
                obs += 1
            elif het and b == ob_base:
                alt += 1
            else:
                oth += 1
        return obs, alt, oth, het, ob_base or ""
    if ev["type"] == "ins":       # the called haplotype carries extra bases: called idx qpos-1 .. qpos+len-2
        i, k = ev["qpos"] - 1, ev["len"]
        idx = list(range(i, i + k))
        het = any(ob.get(j, "-") == "-" for j in idx)
        l, r = shift_range(called, i, k)
        zone = list(range(l, r + k))                    # every position the inserted bases could equivalently occupy
        zp = [pos[j] for j in zone]
        for rd in reads:
            if not (rd.start <= min(zp) and max(zp) <= rd.end):
                continue
            if any(p in rd.dele for p in zp):
                alt += 1                                # a deletion anywhere in the zone: this read lacks the extra bases
            elif all(fix(rd.base.get(p, "?")) == called[j] for p, j in zip(zp, zone)):
                obs += 1
            else:
                oth += 1
        return obs, alt, oth, het, ""
    # del: the called haplotype LACKS bases that the reference has; the junction sits before called idx qpos-1
    i1 = ev["qpos"] - 1
    if i1 <= 0 or i1 >= len(pos):
        return 0, 0, 0, False, ""
    het = (i1 - 1) in extra
    l, r = shift_range(called[:i1] + ev["ref"] + called[i1:], i1, len(ev["ref"]))
    l, r = max(1, l), min(len(pos) - 1, r)
    juncs = [(min(pos[s - 1], pos[s]), max(pos[s - 1], pos[s])) for s in range(l, r + 1)]
    juncs = [(lo, hi) for lo, hi in juncs if hi - lo == 1]
    if not juncs:
        return 0, 0, 0, het, ""
    lo0, hi0 = min(j[0] for j in juncs), max(j[1] for j in juncs)
    for rd in reads:
        if not (rd.start <= lo0 and hi0 <= rd.end):
            continue
        if any(rd.ins.get(lo, 0) > 0 for lo, _ in juncs):
            alt += 1                                    # extra bases inserted somewhere in the zone: the reference-like haplotype
        elif any(p in rd.dele or p not in rd.base for p in range(lo0, hi0 + 1)):
            oth += 1
        else:
            obs += 1
    return obs, alt, oth, het, ""


def shift_range(seq, start, k):
    """Equivalent placements of the k-base indel segment seq[start:start+k]: returns (leftmost, rightmost) start index."""
    l = start
    while l > 0 and seq[l - 1] == seq[l + k - 1]:
        l -= 1
    r = start
    while r + k < len(seq) and seq[r] == seq[r + k]:
        r += 1
    return l, r


def pileup_window(site_idx, called, omap, reads, pos, strand, rng_n=MAX_READS_SHOWN):
    """Strings over WINDOW called-CDS bases centred on the site: called row, other-haplotype row, read rows."""
    minus = strand == "-"
    fix = (lambda b: b.translate(COMP)) if minus else (lambda b: b)
    lo = max(0, min(site_idx - WINDOW // 2, len(called) - WINDOW))
    idx = range(lo, min(len(called), lo + WINDOW))
    ob = omap[0] if omap else {}
    out = {"obs": "".join(called[i] for i in idx), "other": "".join(ob.get(i, " ") for i in idx), "focus": site_idx - lo, "reads": []}
    ps = [pos[i] for i in idx]
    for r in reads:
        if not any(r.start <= p <= r.end for p in ps):
            continue
        s = "".join(" " if not (r.start <= p <= r.end) else "-" if p in r.dele else fix(r.base.get(p, " ")) if p in r.base else " " for p in ps)
        if s.strip():
            out["reads"].append({"s": s, "st": r.strand, "mq": r.mapq, "_i": site_idx})
    return out


def pick_reads(window, want_idx, n=MAX_READS_SHOWN):
    """Balanced display subset: reads classed by the base at the focus column (called / other-haplotype / other)."""
    f, called, other = window["focus"], window["obs"], window["other"]
    groups = defaultdict(list)
    for r in window["reads"]:
        b = r["s"][f] if f < len(r["s"]) else " "
        r["a"] = "obs" if b == called[f] else "alt" if other[f:f + 1] not in (" ", called[f]) and b == other[f] else "oth"
        groups[r["a"]].append(r)
    keep = groups["obs"][: n // 2] + groups["alt"][: n // 2] + groups["oth"][: max(2, n // 5)]
    for r in keep:
        r.pop("_i", None)
    window["reads"] = keep
    return window


# ------------------------------------------------------------------ one case
def evaluate_case(core, case, call, called_cds, other_cds, hapA_contig_seq, sam_lines, ref_start, ref_end):
    """Pure core of the step (no subprocesses): returns the case's sr block."""
    pos = cds_contig_positions(call, len(called_cds))
    reads = [r for r in (parse_sam_line(l, ref_start) for l in sam_lines) if r]
    omap = _other_map(called_cds, other_cds, core.banded_align) if other_cds else None
    al = case["alignments"][0]
    order = sorted(range(len(al["events"])), key=lambda i: (al["events"][i]["type"] == "snv" and al["events"][i].get("effect") == "synonymous", i))
    sites = []
    for i in order[:MAX_SITES]:
        ev = al["events"][i]
        try:
            obs, alt, oth, het, ob = classify_site(ev, reads, pos, called_cds, omap, call["strand"])
            site_idx = ev["qpos"] - 1 if ev["type"] != "del" else max(0, ev["qpos"] - 1)
            win = pick_reads(pileup_window(min(site_idx, len(called_cds) - 1), called_cds, omap, reads, pos, call["strand"]), site_idx)
        except Exception as e:          # one odd site must not lose the case
            sites.append({"ev": i, "het": False, "other_base": "", "obs": 0, "alt": 0, "oth": 0, "verdict": "error: " + str(e)[:40], "window": None})
            continue
        sites.append({"ev": i, "het": bool(het), "other_base": ob, "obs": obs, "alt": alt, "oth": oth,
                      "verdict": core.sr_verdict(obs, alt, oth, het if omap else None), "window": win})
    span = max(1, ref_end - ref_start + 1)
    depth = sum(min(r.end, ref_end) - max(r.start, ref_start) + 1 for r in reads) / float(span)
    return {"sample": "srWGS", "sites": sites, "depth_sr": round(depth, 1), "n_reads": len(reads),
            "note": ("Reads (single-end) realigned to this person's called haplotype gene sequence. Counts: reads with the called base / "
                     "the other haplotype's base / any other base. Homozygous site: nearly all reads should carry the called base; "
                     "heterozygous: about half." if omap else
                     "Other haplotype has no call for this gene, so every site is treated as homozygous.")}


# ------------------------------------------------------------------ plumbing
def read_contig(path, name):
    seq, take = [], False
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if take:
                    break
                take = ln[1:].split()[0] == name
            elif take:
                seq.append(ln.strip())
    if not seq:
        raise KeyError("contig not in trimmed fasta")
    return "".join(seq).upper()


def load_gene_loci(path):
    loci = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or ln.startswith("chrom"):
                continue
            x = ln.rstrip("\n").split("\t")
            if len(x) > 5 and x[5] == "gene":
                loci[x[4]] = (x[0], int(x[1]), int(x[2]))
    return loci


def add_gencode_loci(loci, gtf_path, genes=("E", "F", "G", "DRB3", "DRB4", "DRB5")):
    """Fill hg38 gene loci that S06's table lacks (HLA-E/F/G...) from GENCODE; genes absent from the primary
    assembly (DRB3/DRB4 live on alt haplotypes) simply stay missing."""
    want = {"HLA-" + g: g for g in genes if g not in loci}
    if not want or not os.path.exists(gtf_path):
        return loci
    with gzip.open(gtf_path, "rt") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            x = ln.split("\t", 8)
            if len(x) < 9 or x[2] != "gene":
                continue
            m = re.search(r'gene_name "([^"]+)"', x[8])
            if m and m.group(1) in want and x[0].startswith("chr") and "_" not in x[0]:
                loci[want[m.group(1)]] = (x[0], int(x[3]), int(x[4]))
    return loci


def cram_paths(manifest, persons):
    out = {}
    with open(manifest) as f:
        for r in csv.DictReader(f):
            if r["person_id"] in persons:
                u = r["cram_uri"]
                out[r["person_id"]] = MOUNT + u.split("/", 3)[-1] if u.startswith("gs://") else u
    return out


def run(a, core):
    import importlib.util as iu
    spec = iu.spec_from_file_location("cb_extract_sr", os.path.join(os.path.dirname(os.path.abspath(__file__)), "90_case_browser_extract.py"))
    ext = iu.module_from_spec(spec)
    spec.loader.exec_module(ext)
    out = os.path.expanduser(a.out_dir)
    cases = json.load(open(os.path.join(out, "cases.json")))["cases"]
    priv = {r["case"]: r for r in csv.DictReader(open(os.path.join(out, "case_map.private.tsv")), delimiter="\t")}
    ck = os.path.join(out, "sr.jsonl")
    done = set()
    if os.path.exists(ck):
        done = {json.loads(l)["id"] for l in open(ck) if l.strip()}
    cand = [c for c in cases if c["tab"] == "ladder" and c["species"] == "HLA" and c["alignments"][0]["events"] and c["id"] in priv]
    cand.sort(key=lambda c: (c.get("set") != "all-novel-snv", c["alignments"][0]["stats"]["nm"] > 10, c["id"]))
    cand = [c for c in cand if c["id"] not in done][: max(0, a.max_cases - len(done))]
    store = ext.Store(os.path.expanduser(a.hla_root), "", os.path.expanduser(a.refdata))
    loci = add_gencode_loci(load_gene_loci(os.path.expanduser(a.features)), os.path.expanduser(a.gencode))
    crams = cram_paths(os.path.expanduser(a.manifest), {priv[c["id"]]["person"] for c in cand})
    tmp = tempfile.mkdtemp(prefix="cb_sr_")
    t0, ok, fail = time.time(), 0, defaultdict(int)
    try:
        for c in cand:
            if a.max_minutes and time.time() - t0 > 60 * a.max_minutes:
                log("time budget reached")
                break
            p = priv[c["id"]]
            try:
                person, hap, gene, copy = p["person"], p["hap"], p["gene"], int(p["copy"])
                call = store.calls("HLA", person, hap)[("HLA-" + gene, copy)]
                called = store.observed("HLA", person, hap, call)[0].upper()
                h2 = "2" if hap == "1" else "1"
                c2 = store.calls("HLA", person, h2).get(("HLA-" + gene, 1))
                other = store.observed("HLA", person, h2, c2)[0].upper() if c2 else None
                seq = read_contig(os.path.join(store.odir("HLA", person, hap), "hap%s.trimmed.fa" % hap), call["contig"])
                rs, re_ = max(1, call["start"] - FLANK), min(len(seq), call["end"] + FLANK)
                chrom, s, e = loci[gene]
                with open(os.path.join(tmp, "ref.fa"), "w") as f:
                    f.write(">ref\n%s\n" % seq[rs - 1:re_])
                region = "%s:%d-%d" % (chrom, max(1, s - LOCUS_PAD), e + LOCUS_PAD)
                cmd = "set -o pipefail; %s view -u -T %s %s %s | %s fastq -F 0x900 - > %s/r.fq" % (
                    a.samtools, a.hg38, crams[person], region, a.samtools, tmp)
                sh = subprocess.run(cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
                if sh.returncode:
                    raise RuntimeError("samtools: " + (sh.stderr.strip().splitlines()[0][:90] if sh.stderr.strip() else "rc %d" % sh.returncode))
                mm = subprocess.run([a.minimap2, "-ax", "sr", "-t", "1", os.path.join(tmp, "ref.fa"), os.path.join(tmp, "r.fq")],
                                    capture_output=True, text=True)
                if mm.returncode:
                    raise RuntimeError("minimap2: " + (mm.stderr.strip().splitlines() or ["rc %d" % mm.returncode])[-1][:90])
                sam = mm.stdout.splitlines()
                blk = evaluate_case(core, c, call, called, other, seq, sam, rs, re_)
                blk["id"] = c["id"]
                with open(ck, "a") as f:
                    f.write(json.dumps(blk, separators=(",", ":")) + "\n")
                    f.flush()
                    os.fsync(f.fileno())
                ok += 1
            except Exception as ex_:
                fail[type(ex_).__name__ + ":" + str(ex_)[:60]] += 1
            if (ok + sum(fail.values())) % 10 == 0:
                log("sr: %d ok, %d failed, %.0fs" % (ok, sum(fail.values()), time.time() - t0))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    log("sr done: %d ok; failures: %s; %.0fs" % (ok, dict(fail), time.time() - t0))


def register(sub):
    p = sub.add_parser("sr", help="VM: short-read support for the variants of the ladder cases -> <out-dir>/sr.jsonl")
    p.add_argument("--out-dir", default="~/case_browser")
    p.add_argument("--hla-root", default="~/pipeline_outputs/people")
    p.add_argument("--refdata", default="~/tools/Immuannot_refdata")
    p.add_argument("--features", default="~/s06/results/70/features.tsv", help="hg38 gene loci (S06 70)")
    p.add_argument("--gencode", default="~/ref/gencode.v46.basic.annotation.gtf.gz", help="fills loci S06 70 lacks (HLA-E/F/G)")
    p.add_argument("--manifest", default=MOUNT + "v9/wgs/cram/manifest.csv")
    p.add_argument("--hg38", default="~/ref/hg38_ucsc.fa")
    p.add_argument("--samtools", default="~/envs/yak/bin/samtools")
    p.add_argument("--minimap2", default="~/repos/pilot-validation/.pixi/envs/spechla/bin/minimap2")
    p.add_argument("--max-cases", type=int, default=250)
    p.add_argument("--max-minutes", type=float, default=0, help="stop cleanly after this many minutes (0 = no limit)")
    p.set_defaults(func=lambda a: (setattr(a, "hg38", os.path.expanduser(a.hg38)), setattr(a, "samtools", os.path.expanduser(a.samtools)),
                                   setattr(a, "minimap2", os.path.expanduser(a.minimap2)), run(a, _CORE)))


_CORE = None  # set by 90_case_browser.main()
