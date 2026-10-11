# Sprint S05 — QC and phasing, KIR, TRB germline, canonical counts (from the 2026-10-06 call)

*Opened 2026-10-11 on `main`. Owner: Marc. Source of the asks:
[`sprints/_calls/2026-10-06_CALL_SUMMARY.md`](../_calls/2026-10-06_CALL_SUMMARY.md) and the
`REPORT #11` deck agenda ("what are the counts exactly?", phasing QC, other branches). Plan
agreed with Marc in the project thread on 2026-10-11 ("complete everything").
Log: [`LOG.md`](LOG.md) (append-only). Report ledger: [`METHODS_LEDGER.md`](METHODS_LEDGER.md).*

## Reporting rule for this sprint

Marc will need to explain every method and finding himself, from a high signal-to-noise but
exhaustive final report. So every workstream records, as it goes, in `METHODS_LEDGER.md`:
the question, why it matters, exact inputs (file + version), every parameter and threshold with
the reason it was chosen, what was run (script + commit + command), the numbers that came out,
what they mean, what they do not mean, and caveats. Nothing is "obvious enough to skip". The
final report is assembled from the ledger, the per-script READMEs and `LOG.md`.

## Board

| WS | Script | Ask (source) | Needs VM | Status |
|---|---|---|---|---|
| A | `50_regional_depth.py` | Regional coverage QC: long-read depth over the MHC and the KIR region relative to each person's genome-wide depth (Cole + Marc, 10-06 call) | yes | code + tests written; VM run pending |
| B | `50_regional_depth.py` (TRB part), `51_trb_germline_usage.py` | Genotype TRBV gene presence/absence from depth (TRBV4-3, TRBV12-4 lead from Aleix's W1); join to Aleix's per-person V usage where readable on the VM | yes | 50 covers depth; 51 pending |
| C | `50b_depth_outcome_join.py` | Does low regional depth explain missing KIR framework genes, residual HLA "deletions", rare DQ G1/G2-incompatible calls? Propose a cutoff | yes | pending |
| D | — | Trio phasing QC with yak trioeval (REPORT #11 "Phasing QC"): find the run, finish it, write it up next to twin/relative QC | yes | who started it is unknown; find on VM |
| E | `52_canonical_counts.py` | One canonical counts table, the single source for every number in Figure 1 and the text (REPORT #11 "What are the counts exactly?") | yes | pending |
| F | `53_kir_by_ancestry_panel.py` | KIR gene content and allele frequencies by ancestry as a Figure 1 panel, after the coverage filter (10-06 call) | local from 43 aggregates + VM rerun | pending |
| G | `49_cole_share_package.py` (new release) | Phased haplotype units for Cole: DQA1~DQB1 (and DRB1~DQ) as one unit, with physical-phasing confidence and the new coverage flag (10-06 call: "treat DQ as a single gene unit") | yes | pending |
| H | — | Pull the 5 missing `44_public` tables; VM cleanup plan (Marc runs it) | yes | pending |
| I | — | Reading: Sharon et al. 2016 (V/J usage × MHC), AoU long-read flagship preprint (scoop check) | no | pending |
| J | `54_incompatibility_map.py` | Exploratory: generalize the DQ G1/G2 purge to every physically phased MHC gene pair | yes | after A–G |

## How the VM work runs

This cloud thread cannot drive Marc's browser. VM steps run either through a Claude Code session
on Marc's laptop (which has the Chrome tooling and the VM channel recipe in
`sprints/S03_call8_figures_kir_prediction/VM_CHANNEL.md`), or by Marc pasting the exact commands in
`VM_COMMANDS.md`. Every VM command writes PUBLIC aggregates under `~/s05/results/<NN>/` and
INTERNAL per-person tables under `~/s05/internal/<NN>/` (never pulled).
