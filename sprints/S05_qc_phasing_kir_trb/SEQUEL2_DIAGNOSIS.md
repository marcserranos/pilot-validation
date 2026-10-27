# Why the 991 `sequel2` people have no HLA or KIR calls (S05 WS-K, script 63)

*2026-10-12. Read from `claude/project-thread-a3l14r` at `422d885`; the cited code files are
identical on `origin/main`. `origin/aleix/hla-resolve-phase1` was read without checking it out: it has nothing
on the Tier-3 trim or the chr19 bug. Its only `sequel2` notes concern read BAMs (the platform's
GRCh38 BAM is under `grch38_haplotagged_bam`) and the RNA-seq overlap (940/991). No participant
identifiers are in this note.*

## 1. Summary

- `sequel2` people have the phased assembly FASTAs but **no assembly-to-hg38 PAF or BAM**. The
  cohort builder therefore puts all of them in trim tier `self_align_needed` (Tier 3).
- **Proximate cause, for both HLA and KIR:** they were **never attempted**. The HLA production run
  and the S03 KIR run were both launched with `--single-phase --skip-trim-tier self_align_needed`.
  They were left out on purpose, because Tier 3 had failed its only real test.
- **Root cause (why Tier 3 could not be switched on), found in `scripts/run_immuannot_person.py` at HEAD:**
  1. **Fatal failure:** `trim_assembly()` runs `samtools faidx` on the source `.fa.gz` (l.251).
     For `sequel2` this file is (very likely) **plain gzip, not BGZF**, and faidx cannot index plain
     gzip. So every haplotype stops at l.445 (`SKIP -- samtools faidx trim failed`).
  2. **Underlying design flaw, hidden behind failure 1:** the self-align runs against **chr6
     only** (l.214, l.230). `regions_from_paf()` has no filters and uses a `min(qstart)..max(qend)`
     span (l.193–198). Together they selected 15–18k contigs and 1.3–1.9 GB per haplotype, which is
     the whole assembly.
  3. **KIR-specific:** `ensure_chr6_ref()` hardcodes `"chr6"` (l.214). The orchestrator's phase 2
     does not pass `--chr6-ref-cache` (`run_production_orchestrator.py` l.366–377). A chr19 KIR
     region would therefore be self-aligned to chr6 and could never overlap. The S03 KIR run
     excluded the tier for this reason.
  4. **Latent:** all concurrent workers share one `<cache>.tmp` file (l.212), which is a race.
- **None of the other suspects** in the brief is the cause:
  - The empty manifest columns are the PAF/BAM columns, and that is what defines the tier. The
    FASTA column is filled.
  - The assembly product (path root) is not known to differ in a way that matters.
  - Contig names do not matter, because no step depends on them.
  - Plain gzip is the one format difference with evidence behind it (an htslib error message), but
    it has **not been confirmed with `htsfile`**. See §5.

## 2. The tiers

Tier assignment is in `scripts/legacy/experiments/build_immuannot_cohort.py:91-100`. It is
computed per manifest row and then the best tier is taken per person (`best_tier()`, l.103–110).

| Tier | `trim_tier` | Needs (existence-checked on the mount) | How contigs are chosen | Trim |
|---|---|---|---|---|
| 1 | `paf_region` | hap FASTA + `assembly_hap{1,2}_aln2_hg38_paf` | blocks overlapping `--region` in the person's own PAF | `samtools faidx` sub-range ±100 kb |
| 2 | `bam_whole_contig` | hap FASTA + `assembly_hap{1,2}_aln2_hg38_bam` | `samtools view` of the BAM on `--region` | whole contigs |
| 3 | `self_align_needed` | hap FASTA only | self-align the FASTA with minimap2, opt-in via `--enable-self-align-fallback` | as Tier 1, from the synthetic PAF |

Where the `sequel2` people fall: `reports/pre_sprint/lr_data_census/README.md:34` lists `sequel2` as 991 rows
with "Partial (hap1/hap2/primary FASTA+GFA only; **no alignment-to-hg38 files**)". `row_trim_tier()`
returns `self_align_needed` when both FASTAs resolve but neither the PAFs nor the BAMs do (l.95–100).
So the whole platform lands in Tier 3. The other platforms (12,261 people in the production cohort,
EXPERIMENTS.md 2026-08-05) are `paf_region`.

## 3. The exact code path

All line numbers are at HEAD `422d885`, before the patch.

**Production (why there are no calls at all).**
- `scripts/production_orchestrator/PREFLIGHT.md:116`: "**So: launch with `--single-phase
  --skip-trim-tier self_align_needed`.** … leaves the 991 sequel2 people untouched in the cohort
  file, already tagged, for a fast-follow once Tier 3 actually works."
- `run_production_orchestrator.py:335-342` (`filter_people`) drops the tier before dispatch.
- KIR: `sprints/S03_call8_figures_kir_prediction/LOG.md:79`: "KIR run uses `--single-phase
  --skip-trim-tier self_align_needed`: the orchestrator's phase 2 (Tier-3 self-align) calls
  run_immuannot_person.py without --chr6-ref-cache, so it would align chr19 contigs against the
  cached chr6 slice → wrong/empty calls." The full run's target was the 12,261 `paf_region` people
  (LOG.md:83).

**Tier 3 itself, as tested on 2 real sequel2 people (EXPERIMENTS.md 2026-08-05, l.608–669; 0 genes
called).**

`process_person()` (l.484–513) resolves the FASTAs and finds no BAM or PAF. In
`process_haplotype()`, `targets is None and not aln_rel` is true (l.375), so the haplotype enters
Tier 3:

1. `ensure_chr6_ref()` (l.202–221) writes the output of `samtools faidx <hg38> chr6` to `~/ref/chr6.fasta`.
   - `"chr6"` is a literal (l.214).
   - The temp path is the shared `cache_path + ".tmp"` (l.212).
2. `self_align_paf()` (l.224–233) runs `minimap2 -x asm5 -t 4 chr6.fasta <hap.fa.gz>`. The whole
   ~3.1 Gb haplotype is aligned to chr6 alone.
3. `regions_from_paf()` (l.160–199):
   - It accepts every block on chr6 that overlaps the window, with no MAPQ, length, primary or
     identity filter (cols 10–12 are not read).
   - It spans `min(qstart)-pad .. max(qend)+pad` per contig (l.193–198).
   - Result: 15,070–18,447 contigs and 1,337–1,882 MB per haplotype (EXPERIMENTS table, l.616–621).
4. `trim_assembly()` (l.236–254) runs `samtools faidx --fai-idx … --gzi-idx … <hap.fa.gz> <targets>`
   (l.251). This fails with `Cannot index files compressed with gzip, please use bgzip`
   (EXPERIMENTS l.626–631). The haplotype is then SKIPped at l.444–446 and no `hapN.gtf.gz` is
   written.

EXPERIMENTS l.656–662 states why step 4 failing was the safer outcome: "fixing A alone would have been worse than failing … producing a 1.3-1.9 GB trimmed FASTA and feeding it to `immuannot.sh` at 24-way concurrency".

**The KIR variant.** Suppose phase 2 had run with `--region chr19:54600000-54920000`:
- Step 1 would still have built chr6. The orchestrator command (l.366–377) has no
  `--chr6-ref-cache` or `--hg38-ref`.
- Step 3 would then have seen only `chr6` target names. With `accept_chrom={chr19,19}` it returns
  `{}`, so every haplotype is SKIPped at l.408–411.

`41_kir_pilot.py:33-43` built a chr19-cache workaround, but it was **never exercised**: "the
candidate cohort had zero sequel2 people" (`reports/hla_popgen/41_kir_scoping/README.md:89`, l.364).

**Missing prerequisite (secondary).** `check_phase2_prereqs()` (orchestrator l.283–308) notes that
`bootstrap_vm.sh` does not fetch `~/ref/Homo_sapiens_assembly38.fasta`. Tier 3 cannot run without it.

## 4. What the patch changes

The fix is the list the project had already agreed. PREFLIGHT.md:120–124: "at minimum: align
against whole hg38 (not chr6 alone) so MAPQ is meaningful again, add length/identity filters to
`regions_from_paf()`, replace its `min(qstart)/max(qend)` span with per-contig block clustering, and
fix the `ensure_chr6_ref()` tmp-file race". To these four items the patch adds the plain-gzip
extraction and the region-derived chromosome.

**Changes in `scripts/run_immuannot_person.py` (Tier 3 only; Tiers 1 and 2 are byte-for-byte the
same code path):**
- **Extraction (cause 1).** `stream_extract_fasta()` reads the source once and writes the same
  `ctg:S-E` headers and 60-column lines that faidx writes. It accepts plain gzip, BGZF or plain
  FASTA and never writes an index. `is_bgzf()` is added for the VM check.
- **Two stages (cause 2).**
  - **Stage 1** aligns to a slice of the region's chromosome around the target window (±1 Mb). It
    keeps only contigs with a primary block of at least 10 kb (`candidate_contigs_from_paf`), and
    refuses to go on if the candidates total more than 600 MB.
  - **Stage 2** extracts only those contigs. It aligns them to the **whole primary GRCh38**
    (chr1–22, X, Y, M; a prebuilt `asm5` `.mmi`, passed as `--self-align-confirm-ref`), so MAPQ is
    meaningful again.
  - The trim comes from stage-2 blocks that pass all of: primary, ≥10 kb, MAPQ ≥20, col10/col11 ≥0.5,
    and on `--region`. Blocks are clustered per contig (gap ≤250 kb) and the largest cluster is kept.
  - Tier 3 now **refuses** to run without the stage-2 reference.
  - A size guard SKIPs any trimmed input above 25 MB per haplotype. The expected size is about 4 MB
    for HLA and about 1 MB for KIR.
- **Filters (cause 2).** `regions_from_paf()` gains optional filters and clustering. With the
  default arguments it returns exactly what it returned before (this is tested).
- **KIR (cause 3).** The stage-1 slice is derived from `--region` (`ensure_ref_slice`). Any stage-1
  reference whose sequence names do not include the region's chromosome is refused with a clear
  SKIP. A cached slice is reused only if its header matches.
- **Race (cause 4).** The temp file is unique per process and thread, and is renamed into place
  atomically.

The driver `scripts/hla_popgen/63_sequel2_recall.py` runs this for the `self_align_needed` people
only, into `~/s05/internal/63/{hla,kir}/pipeline_outputs`. It never writes into
`~/pipeline_outputs*`. Tests are in `scripts/hla_popgen/tests/test_63_sequel2_recall.py` (synthetic
plain-gzip haplotype with a segmental-duplication decoy, a repeat contig and a stray block).

**What the patch does not fix:**
- The production orchestrator's phase 2 is unchanged. Run as it was, it will now SKIP Tier 3
  cleanly ("needs --self-align-confirm-ref") instead of producing wrong calls.
- `41_kir_pilot.py`'s `--chr6-ref-cache` route also SKIPs for the same reason. It is historical.

## 5. Still unknown: check on the VM (read-only)

All of these are read-only. Run them in the specimmune env. Each prints counts or formats, not
identifiers. **Do not touch `~/mnt/aou-controlled`**: it hangs (ENVIRONMENT.md quirk #35).

```bash
PY=~/repos/pilot-validation/.pixi/envs/specimmune/bin/python3; BIN=$(dirname $PY)
M=~/workspace/vwb-aou-datasets-controlled-v9
# U1. Is the 991 gap really "never attempted"? Count tier and platform in the cohort, and how many
#     self_align_needed people have any output in the production trees (expect 0 calls; up to 2
#     dirs from the 2026-08-05 test).
$PY - <<'EOF'
import os, pandas as pd
c = pd.read_csv(os.path.expanduser("~/pipeline_outputs/immuannot_cohort_full.tsv"), sep="\t", dtype=str)
s = c[c.trim_tier == "self_align_needed"]
print("tiers", c.trim_tier.value_counts().to_dict()); print("platform of tier3", s.platform.value_counts().to_dict())
for root, f in [("~/pipeline_outputs", "immuannot_calls.tsv"), ("~/pipeline_outputs_kir", "immuannot_calls.kir.tsv")]:
    r = os.path.expanduser(root); ids = set(s.person_id)
    calls = set(pd.read_csv(os.path.join(r, f), sep="\t", dtype=str, usecols=["person_id"]).person_id)
    print(root, "tier3 with calls:", len(ids & calls), "tier3 with a person dir:", sum(os.path.isdir(os.path.join(r, p)) for p in ids))
EOF
# U2. Manifest: which assembly columns are filled for sequel2, and from which bucket and path root.
$PY - <<'EOF'
import pandas as pd
import os
m = pd.read_csv(os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv"), sep="\t", dtype=str)
for plat in ("sequel2", "revio"):
    d = m[m.platform == plat]
    cols = [c for c in m.columns if c.startswith("assembly_")]
    print(plat, len(d), {c: int(d[c].notna().sum()) for c in cols})
    fa = d["assembly_hap1_fa"].dropna()
    print("  bucket/root:", fa.str.split("/").str[:6].str.join("/").value_counts().head(5).to_dict())
    print("  suffix:", fa.str.extract(r"(\.[A-Za-z.]+)$")[0].value_counts().to_dict())
EOF
# U3. Compression: BGZF vs plain gzip, on 20 sequel2 and 20 revio hap1 FASTAs (counts only).
$PY - <<'EOF'
import os, sys, pandas as pd
sys.path.insert(0, os.path.expanduser("~/repos/pilot-validation-main/scripts"))
from run_immuannot_person import is_bgzf
M = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/")
m = pd.read_csv(M + "v9/wgs/long_read/manifest.tsv", sep="\t", dtype=str)
for plat in ("sequel2", "revio"):
    n = {"bgzf": 0, "plain_gzip": 0, "not_gzip": 0, "missing": 0}
    for u in m[m.platform == plat]["assembly_hap1_fa"].dropna().head(20):
        p = M + u.split("/", 3)[-1]
        if not os.path.exists(p): n["missing"] += 1; continue
        h = open(p, "rb").read(2)
        n["bgzf" if is_bgzf(p) else ("plain_gzip" if h == b"\x1f\x8b" else "not_gzip")] += 1
    print(plat, n)
EOF
#     and the htslib view of one file of each (prints only the format line):
$PY -c "import pandas as pd,os;M=os.path.expanduser('~/workspace/vwb-aou-datasets-controlled-v9/');m=pd.read_csv(M+'v9/wgs/long_read/manifest.tsv',sep='\t',dtype=str);[print(M+m[m.platform==p]['assembly_hap1_fa'].dropna().iloc[0].split('/',3)[-1]) for p in ('sequel2','revio')]" \
  | while read f; do $BIN/htsfile "$f" | sed 's#.*:\t##'; done
# U4. Contig naming and size of one sequel2 haplotype (names and lengths only; the first 3 names):
F=$($PY -c "import pandas as pd,os;M=os.path.expanduser('~/workspace/vwb-aou-datasets-controlled-v9/');m=pd.read_csv(M+'v9/wgs/long_read/manifest.tsv',sep='\t',dtype=str);print(M+m[m.platform=='sequel2']['assembly_hap1_fa'].dropna().iloc[0].split('/',3)[-1])")
zcat "$F" | awk '/^>/{if(n)print n"\t"l; n=$1; l=0; next}{l+=length($0)}END{print n"\t"l}' | sort -k2,2nr \
  | awk '{s+=$2; c++; if(c<=3) print "top:", $1, $2} END{print "contigs", c, "total_bp", s}'
# U5. Reference: present, indexed, and which non-primary sequences it carries (ALT/HLA/decoy).
ls -la ~/ref/Homo_sapiens_assembly38.fasta ~/ref/Homo_sapiens_assembly38.fasta.fai ~/ref/chr6.fasta 2>&1
cut -f1 ~/ref/Homo_sapiens_assembly38.fasta.fai | awk '{t="primary"} /_alt$/{t="alt"} /^HLA-/{t="hla"} /decoy|^chrEBV/{t="decoy"} /_random$|^chrUn/{t="unplaced"} {n[t]++} END{for(k in n) print k, n[k]}'
head -c 200 ~/ref/chr6.fasta | head -1
# U6. Box and tools.
nproc; free -g | head -2; df -h ~ | tail -1
$BIN/minimap2 --version; $BIN/samtools --version | head -1
# U7. The 2026-08-05 test logs, if still on disk (counts only).
ls ~/pipeline_outputs/*/immuannot_output/hap1.self_align.paf 2>/dev/null | wc -l
grep -l "Cannot index files compressed with gzip" ~/pipeline_outputs/*/immuannot_output/.orchestrator_last_run.log 2>/dev/null | wc -l
```

**What each answer changes:**
- **U1.** If any self_align_needed person already has calls, that person is skipped by the merge
  overlap check.
- **U2 and U3.** If `sequel2` FASTAs are BGZF after all, failure 1 is something else. The stream
  extractor still works either way. If the bucket differs from `gs://vwb-aou-datasets-controlled/`,
  `strip_bucket()` would treat every FASTA as missing. The smoke run would then show
  `no_assembly_fasta` in `failure_reasons.tsv`.
- **U4.** Names are not used by any step. Size gives the stage-2 input bound.
- **U5.** The `.fai` must exist, or `--prepare-refs` stops rather than writing next to the shared
  reference. ALT and HLA contigs are why stage 2 uses a primary-only copy.
- **U6.** Stage 2 holds a ~5–8 GB index per concurrent worker, so RAM limits `--workers`.

**What only the smoke run can show:**
- The real stage-1 and stage-2 timings.
- The number of candidate contigs per haplotype.
- That HLA trimmed sizes are a few MB (and KIR about 1 MB), and that genes are called on both
  haplotypes.

Check those in `~/s05/results/63_smoke/trim_summary.tsv` before the full run.
