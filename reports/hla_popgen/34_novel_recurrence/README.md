# 34 — novel protein alleles by recurrence and gene

*Marc's ask, 2026-09-22.* Built from committed aggregates only.

## Headline

- **1404** distinct novel protein alleles
- **231** seen in ≥2 unrelated people (the 'recurrent' set)
- **29** seen in ≥20 unrelated people — the threshold above which All of Us lets us publish a count, so these are the nameable ones

## Why the x-axis is not 1, 2, 3, … 19, 20+

All of Us suppresses every participant count between 1 and 19, so script 24's committed table writes them all as `<20`. The individual counts do not exist outside the VM, so a per-count histogram cannot be built here. The three-class split (seen once / 2–19 / ≥20) is the same information at the resolution the disclosure rule permits.

To get the true per-count histogram, one VM run is needed against the uncensored cluster table — the counts are computed there and only suppressed on the way out.

## Files

- `novel_protein_recurrence_by_gene.tsv` — the counts behind both panels.
- `fig_novel_recurrence.png` / `.pdf`.
