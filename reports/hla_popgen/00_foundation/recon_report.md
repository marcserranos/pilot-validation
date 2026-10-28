# hla_popgen VM recon report

Sampled **50** people from `/home/jupyter/pipeline_outputs` in 324.2s.

## Most important: file presence + extrapolated bytes

| file | present/checked | rate | mean bytes | extrapolated total (n=12000) |
|---|---|---|---|---|
| `hap1.gtf.gz` | 50/50 | 100.0% | 9257.1 | 222171360 |
| `hap1.trimmed.fa` | 50/50 | 100.0% | 9696626.9 | 232719045600 |
| `hap1/cds.fa.gz` | 50/50 | 100.0% | 7833.7 | 188009280 |
| `hap1/mm2.ipd.cds.paf.gz` | 50/50 | 100.0% | 851785.7 | 20442856800 |
| `hap1/mm2.ipd.gen.paf.gz` | 50/50 | 100.0% | 793614.2 | 19046739840 |
| `hap1/gene.filtered.paf` | 50/50 | 100.0% | 7094784.5 | 170274828480 |
| `hap2.gtf.gz` | 50/50 | 100.0% | 8958.6 | 215005920 |
| `hap2.trimmed.fa` | 50/50 | 100.0% | 9699819.2 | 232795659840 |
| `hap2/cds.fa.gz` | 50/50 | 100.0% | 7681.1 | 184345920 |
| `hap2/mm2.ipd.cds.paf.gz` | 50/50 | 100.0% | 804414.0 | 19305936000 |
| `hap2/mm2.ipd.gen.paf.gz` | 50/50 | 100.0% | 756873.4 | 18164962560 |
| `hap2/gene.filtered.paf` | 50/50 | 100.0% | 6704709.1 | 160913018400 |

## KIR gene check

No KIR genes found in sample -- consistent with spec expectation.

## Undocumented GTF attribute keys

None -- every attribute key seen matches reference/IMMUANNOT_GTF_SPEC.md.

## Contigs per hap file

80/100 hap files (80.0%) span more than one contig (max seen: 16). This is the ceiling on cis-pairing loss.

## template_distance distribution

n=3421, min=0, max=123, mean=2.61, exact (0) = 2275

## Conditional field presence

- cds_distance present: 1155/3589 (32.2%)
- cds_mut present: 331/3589 (9.2%)
- template_warning present: 3421/3589 (95.3%)

## Novelty depth distribution (which field 'new' lands in)

{4: 734, 2: 314, 3: 104, 1: 3}

## Header vs body cross-check

No mismatches -- header claims agree with the parsed body.


## Full gene inventory

- HLA-DQA2: 101
- HLA-DQB2: 100
- TAP2: 99
- HLA-DOB: 99
- HLA-DPA2: 98
- HLA-DPB1: 98
- HLA-DPA1: 98
- HLA-DOA: 98
- HLA-DMB: 98
- HLA-E: 98
- HLA-N: 98
- HLA-DQA1: 97
- HLA-L: 97
- HLA-DQB1: 96
- HLA-DMA: 96
- HLA-S: 96
- HLA-W: 96
- HLA-A: 96
- HLA-DRA: 95
- HLA-B: 95
- HLA-J: 95
- HLA-G: 95
- HLA-V: 95
- HLA-DRB1: 94
- HLA-DPB2: 94
- TAP1: 94
- HLA-F: 94
- HLA-P: 94
- C4A: 94
- HLA-C: 93
- MICA: 92
- MICB: 88
- HLA-H: 87
- HLA-T: 86
- HLA-U: 85
- HLA-K: 84
- C4B: 74
- HLA-DRB3: 38
- HLA-DRB4: 28
- HLA-Y: 16
- HLA-HFE: 11
- HLA-DRB5: 9
