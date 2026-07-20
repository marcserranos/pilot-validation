# Smoke test — HLA-Resolve vs Lai truth, HG002 demo

**2026-07-18 · install validation, NOT an accuracy benchmark**

Ran HLA-Resolve v0.8.5 on its own bundled `demo/HG002.hifi_reads.fastq.gz` (88,955 reads, PacBio HiFi, hybrid-capture), scored against `../reference/Lai_Supplementary-6.xlsx`.

```
--input_file demo/HG002.hifi_reads.fastq.gz --sample_name HG002_demo
--platform pacbio --scheme hybrid_capture --threads 12
```

Runtime **21m25s** (DeepVariant `make_examples`: 11m18s wall / 87m54s CPU across 8 shards).

## Result: 16/16 alleles exact at 4-field

| Gene | Truth (Lai) | HLA-Resolve | Match |
|---|---|---|---|
| A | `26:01:01:01` + `01:01:01:01` | `01:01:01:01` + `26:01:01:01` | exact |
| B | `38:01:01:01` + `35:08:01:01` | `35:08:01:01` + `38:01:01:01` | exact |
| C | `12:03:01:01` + `04:01:01:06` | `04:01:01:06` + `12:03:01:01` | exact |
| DPA1 | `01:03:01:04` + `01:03:01:02` | `01:03:01:02` + `01:03:01:04` | exact |
| DPB1 | `04:01:01:03` + `04:01:01:01` | `04:01:01:01` + `04:01:01:03` | exact |
| DQA1 | `03:01:01:01` + `01:05:01:01` | `01:05:01:01` + `03:01:01:01` | exact |
| DQB1 | `03:02:01:01` + `05:01:01:05` | `03:02:01:01` + `05:01:01:05` | exact |
| DRB1 | `04:02:01` + `10:01:01:03` | `04:02:01` + `10:01:01:03` | exact |

## Notes that matter for the real runs

- **Haplotype order differs on 6 of 8 genes.** Not a discordance — a genotype is an *unordered* pair, and hap1/hap2 labelling is arbitrary. Marc hit the identical thing with SpecHLA and logged it as cosmetic. **The scorer must compare unordered pairs** or it will manufacture 6 false errors.
- **`DRB1*04:02:01` is 3-field in both truth and call** — that allele is only defined to 3 fields in IMGT. Not a truncation; the scorer must not treat missing 4th fields as a miss.
- **Output covers the 8 classical genes only.** Truth also lists `DRB4*01:03:01:10` for HG002 (and no DRB3/DRB5), which the CSV has no column for — even though the run *did* classify 9,240 reads as DRB paralogs. DRB3/4/5 can't be scored from this output.

## What this does and doesn't establish

**Does:** the install is correct end-to-end — pixi env, DeepVariant container, augmented reference, IMGT matching all functional. That was the smoke test's purpose.

**Does not:** demonstrate real-world accuracy. This is the tool's own curated demo file, near-certainly its best case; beating its published 4-field figure (7/8→100%) on it means little. The real benchmark is **whole-genome** HiFi for HG002 / HG01258 / HG03579 pulled independently, run through **both** HLA-Resolve and SpecImmune, scored against this same truth table.
