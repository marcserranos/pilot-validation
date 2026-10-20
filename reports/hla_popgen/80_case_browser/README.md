# 80 — case browser (code only; outputs are never committed)

A small, dependency-free viewer for looking at the *bases* behind HLA/KIR call-quality numbers.

| Tab | What you see |
|---|---|
| Distance ladder | Immuannot's own CDS alignment (replayed from `mm2.ipd.cds.paf.gz`, not recomputed) for calls at ~1, 2, 5, 10, 50, 100, >=150 CDS edits from their reference allele: every differing base, codon and amino acid, exon / groove position, homopolymer or tandem-repeat indels, % of the CDS that differs. HLA and KIR. Tick cases to compare them on one axis. |
| Novel vs known | The same alignments split by Immuannot's novelty class. |
| Mendelian | A child haplotype against the parent 57b assigned it to: child vs its reference, each parent haplotype vs its reference (all Immuannot), plus a direct child-vs-parent CDS alignment. |
| Switches | Parent of origin along a child's haplotype contig at gene resolution (57b's existing k-mer assignment; no yak rerun). |

```
# VM (JupyterLab terminal), from a throwaway clone of this branch:
python3 scripts/hla_popgen/80_case_browser.py extract --out-dir ~/case_browser
python3 scripts/hla_popgen/80_case_browser.py build   --cases ~/case_browser/cases.json --out-dir ~/case_browser
# open ~/case_browser/index.html from JupyterLab (right-click -> Open in New Browser Tab)

# anywhere, synthetic data only:
python3 scripts/hla_popgen/80_case_browser.py demo --out-dir /tmp/cb_demo
python3 -m pytest scripts/hla_popgen/tests/test_80_case_browser.py scripts/hla_popgen/tests/test_80_extract.py -q
```

Local viewing of a downloaded `cases.json` (put it in the gitignored `case_browser/` folder, then):

```
python3 scripts/hla_popgen/80_case_browser.py build --cases case_browser/cases.json --out-dir case_browser
open case_browser/index.html
```

`build` also normalises an older `cases.json` (exon labels, Mendelian summary), so no re-extraction is needed.

Data rules: `cases.json` / `index.html` contain participant-derived sequence, so they stay in the Workbench; case ids are
`c001`-style, the id-to-person map is `case_map.private.tsv` (mode 600, gitignored). Parent A / B are the two parents in
sorted-ID order (57's convention), not sex.
