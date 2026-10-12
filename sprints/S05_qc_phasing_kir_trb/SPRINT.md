# Sprint S05 — QC and phasing, KIR, TRB germline, canonical counts (supervisor asks, 2026-10-06)

*Opened 2026-10-11 on `main`. Owner: Marc. Source of the asks:
the 2026-10-06 supervisor meeting and the `REPORT #11` deck agenda ("what are the counts exactly?", phasing QC, other branches). Plan
agreed with Marc in the project thread on 2026-10-11 ("complete everything").
Log: [`LOG.md`](LOG.md) (append-only). Report ledger: [`METHODS_LEDGER.md`](METHODS_LEDGER.md).*

## Reporting rule for this sprint

Marc will need to explain every method and finding himself, from a high signal-to-noise but
exhaustive final report. So every workstream records, as it goes, in `METHODS_LEDGER.md`:
the question, why it matters, exact inputs (file + version), every parameter and threshold with
the reason it was chosen, what was run (script + commit + command), the numbers that came out,
what they mean, what they do not mean, and caveats. Nothing is "obvious enough to skip". The
final report is assembled from the ledger, the per-script READMEs and `LOG.md`.

## Board (re-planned 2026-10-11 after the 10-06 work, scripts 50–57c, landed on `main`)

The 10-06 work already answers several original asks. Its results are in `reports/hla_popgen/50_*`–`57_*` and
`context/EXPERIMENTS.md` (2026-10-06). New scripts start at 58.

| WS | Script | Ask | Needs VM | Status |
|---|---|---|---|---|
| A | `58_regional_depth.py` | Regional coverage QC: long-read depth over the MHC and the KIR region, relative to each person's genome-wide depth (supervisor ask) | yes | code + tests done; VM run pending |
| B | `58` (TRB part), `60_trb_germline_usage.py` | TRBV gene copy number from depth; join it to Aleix's per-person V usage. Prior art: the TRBV4-3/3-2/6-2 deletion is known (`READING.md` §B), so the contribution is scale, ancestry and the link to usage | yes | 58 running on VM; 60 code + tests done |
| B2 | `65_trb_assembly_genotype.py` | TRBV4-3 has no locus on GRCh38 primary chr7 (it sits on a ~21 kb insertion carried only on an ALT contig), so 58 cannot see it. Genotype the insertion per haplotype from each person's assemblies | yes | code + tests done; VM run pending |
| C | `59_depth_outcome_join.py` | Does low regional depth explain incomplete KIR haplotypes (51), structurally anomalous DQ contigs (50e) and residual absences? Propose the cutoff | yes | code + tests done; runs after 58 |
| D | `57`, `57b`, `57c` | Trio phasing QC with yak | — | **done 10-06** (179 trios). Optional: genome-wide yak; use trios to tune the confidence score |
| E | `50_qc_foundation.py` (a) | Canonical counts | — | **mostly done 10-06** (waterfall per ancestry). Remaining: haplotype / per-gene / novel-by-class rows in the same table, for the final report |
| F | `61_kir_by_ancestry_panel.py` | KIR gene content and allele frequencies by ancestry as a Figure 1 panel, on complete haplotypes only (51 rule), after the coverage filter | local draft + VM | code + tests done; local draft panel built; VM run after 59 |
| G | `62_cole_haplotype_units.py` | Phased units for fine-mapping: DQA1~DQB1 and DRB1~DQA1~DQB1 per haplotype, with the 57b physical-evidence rule and the coverage flag | yes | code + tests done; VM run after 59 |
| H | — | Pull the 5 missing `44_public` tables and the 53 ASE table | yes | in progress |
| I | — | Reading: Sharon 2016, TRBV deletion prior art, AoU long-read flagship scoop check | no | **done** (`READING.md`) |
| J | `52_mhc_haplotype_atlas.py` | Incompatibility map across MHC gene pairs | — | **done 10-06** (absent-haplotype atlas) |
| K | `63_*` | Fix the Tier-3 trim for the 991 Sequel II people (mostly AFR) and call HLA + KIR for them | yes, heavy | diagnosed + patched (`SEQUEL2_DIAGNOSIS.md`); VM checks + smoke pending; full run ~27 days on 4 vCPU, ~3 days on a resized app |
| L | `64_*` | Novel alleles: re-match against a current IPD-IMGT/HLA release and add read-level support (decides whether the protein-novel headline is demoted) | yes | code + tests done (21); VM run pending |
| M | `53` scale-up | ASE on all people with long reads and RNA | yes, needs a bigger app (Marc authorises the resize) | pending |

Open decisions for Marc (`context/DECISIONS.md`, S05): the exclusion list and ancestry-purity rule; demoting
classical protein-novel counts; KIR claims on complete haplotypes only; the physical-haplotype evidence rule from trios.

## How the VM work runs

This cloud thread cannot drive Marc's browser. VM steps run either through a Claude Code session
on Marc's laptop (which has the Chrome tooling and the VM channel recipe in
`sprints/S03_call8_figures_kir_prediction/VM_CHANNEL.md`), or by Marc pasting the exact commands in
`VM_COMMANDS.md`. Every VM command writes PUBLIC aggregates under `~/s05/results/<NN>/` and
INTERNAL per-person tables under `~/s05/internal/<NN>/` (never pulled). The 10-06 scripts use the same tree, so
numbers are never reused.
