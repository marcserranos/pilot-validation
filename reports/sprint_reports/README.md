# reports/sprint_reports/ — one finished report per sprint

The delivery-ready, self-contained document of each sprint lives here. Per-script result folders stay in
[`../hla_popgen/`](../hla_popgen/) (see its by-era table); older reports in [`../pre_sprint/`](../pre_sprint/).

| Sprint | Report | Built by | In git? |
|---|---|---|---|
| S05 QC / phasing / KIR / TRB | [`S05_REPORT.html`](S05_REPORT.html) | `sprints/S05_qc_phasing_kir_trb/report/build_report.py` (sources in `report/src/`) | yes |
| S06 call quality | `S06_FINAL_REPORT.html` (local) — built in `reports/hla_popgen/79_s06_final/`; copy it here when published | `scripts/hla_popgen/79_s06_final_report.py` | no (small-cell masking pass still open) |
| S08 case browser / novel recheck | `S08_REPORT.html` (local) | `scripts/hla_popgen/117_s08_full_report.py --out reports/sprint_reports/S08_REPORT.html` (branch `s08-case-browser`; shared figure library `116a_figs.py`) | no (participant-derived sequence) |
| S09 discovery loop | [`S09_REPORT.html`](S09_REPORT.html) | `scripts/hla_popgen/116_s09_full_report.py --out reports/sprint_reports/S09_REPORT.html` (branch `s09-discovery-loop`; shared figure library `116a_figs.py`) | yes (aggregate only, counts <20 masked) |

**Rule:** a sprint's report goes here as `S0N_REPORT.html` — never in a sprint folder or in `reports/hla_popgen/`. The builder's `--out` points here; per-script results stay in `reports/hla_popgen/NN_*/`.

S01–S04 have no single report: their plain-language findings are `FINDINGS_FOR_MARC.md` in each sprint folder
(`sprints/S0N_*/`), with the S04 figure catalogue in `FIGURES_INDEX.md`.
