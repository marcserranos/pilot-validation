# 65 reads mode: TRB insertion (TRBV4-3 / TRBV3-2) genotype from GRCh38 HiFi reads

How it was made:
- VM run on 2026-10-12, 12:02 → 16:59 UTC, 2 workers, 17,726 s.
- Script at 3e0df40; aggregates regenerated at 5fa84b0 (`--aggregate-only`, PR #25).
- Probes and breakpoint from `--build-probes` (probe_classes.tsv): breakpoint chr7:142,346,200–142,352,104, insertion about 21,666 bp, insertion classes TRBV4-3 and TRBV3-2.
- Calls: reads spanning the breakpoint ±1 kb; min_reads 6; allele-fraction cuts 0.2 / 0.8. The median person has 14 informative reads.
- Revio and Sequel IIe GRCh38 BAMs carry no HP tags. Calls are per-person dosage (0/1/2 copies), and every haplotype is `no_hp_tag` (haplotype_completeness.tsv).

Folders:
- `reads_all/`: 12,233 people, 11,513 called.
- `reads_unrelated/`: 11,856 unrelated people, 11,151 called.

Headline (unrelated, TRBV4-3): people with 0 copies, i.e. germline absence:
- SAS 55.1%
- AFR 38.5%
- EAS 37.1%
- MID ~32% (counts masked in the public table)
- AMR 30.1%
- EUR 19.5%
- ALL 33.0%
- Insertion-allele frequency: 0.464 overall.

Caveat: there are too few heterozygotes in every ancestry (HWE p from 1e-31 to 0.002), consistent with heterozygotes being called homozygous at about 14 reads. Validation against the stratified-2,000 assembly run (65b) is pending.

Disclosure:
- No ids. Counts 1–19 are masked by the script, with cross-row complementary suppression.
- Checked at pull: no single-masked allele-fraction partition; no ALL-vs-platform or ALL-vs-ancestry single mask. That check missed the cross-table sums (genotype counts = sums of AF bins in [0,0.2), [0.2,0.8), [0.8,1]), people not called, and the 23 people with no listed ancestry.
- Re-masked 2026-10-13 by `scripts/hla_popgen/66_cross_table_suppression.py` (cross-table complementary suppression, `<20*`). Earlier versions of these tables remain in git history.
