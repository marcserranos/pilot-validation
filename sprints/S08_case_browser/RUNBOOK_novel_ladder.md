# Runbook — novel-allele ladder (scripts 118/119), VM side

Written 2026-10-29 by a cloud session that had **no VM access**: 118/119 were tested only on synthetic fixtures (`scripts/hla_popgen/tests/test_118_ladder.py`, plus a browser render of 119 on synthetic data). Expect a first-run fix or two in 118's data-plumbing (reading `case_map.private.tsv`, `sr.jsonl`, the lrWGS manifest); the maths and the report are tested.

Newest IPD release: **3.66.0 (14 Oct 2026)**, verified from `release_version.txt` on the ANHIG/IMGTHLA `Latest` branch; the S08 log already used 3.66. Re-check `Allelelist.txt` header on the day.

Assume the VM restarted: remount the bucket and `ls` first (ENVIRONMENT quirk #14). Everything below is read-only on shared paths; outputs go to `~/s10` and `~/case_browser`.

1. Code: fresh `git clone` of `main` into `~/cb_code`.
2. IPD: `curl -O https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_nuc.fasta` (43 MB; needs GitHub raw reachable from the Workbench — otherwise upload it).
3. Short-read evidence for **every** substitution-only novel call (the S08 run capped F and G at 60):
   Use a NEW out-dir (`~/s10_cb`): a fresh extract renumbers case ids, so writing into `~/case_browser` would orphan its existing `sr.jsonl`.
   `python3 scripts/hla_popgen/90_case_browser.py extract --out-dir ~/s10_cb --exhaustive-snv 5000` (check the subcommand layout with `--help`; the 10-21 run used the same entry point), then
   `python3 scripts/hla_popgen/90_case_browser.py sr --out-dir ~/s10_cb --max-cases 5000` (~1.5 s/case; resumable; `--max-minutes` keeps it inside the one-day rule).
4. Ladder: `python3 scripts/hla_popgen/118_novel_ladder.py --ipd-nuc hla_nuc.fasta --ipd-release 3.66.0 --out-dir ~/s10 --sr-dir ~/s10_cb`
5. Exact report (stays in the Workbench): `python3 scripts/hla_popgen/119_novel_ladder_report.py --data ~/s10/ladder_data.private.json --out ~/s10/S08_NOVEL_LADDER.private.html` — open it from JupyterLab.
6. Cloud-safe report: `119 ... --data ~/s10/ladder_data.pub.json --out ~/s10/S08_REPORT.html`. `ladder_data.pub.json` is the only file meant to leave the VM, and **only after Marc OKs the allele-level class table** (DECISIONS.md, open item).
7. Delete the code clone and the IPD file afterwards; `novel_ladder.private.tsv` and `ladder_data.private.json` are mode 600 and stay.

Not done here (needs the VM): all distributions, counts, saturation curves. The report states this on its first page until data is loaded.
