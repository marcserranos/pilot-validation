# VM operator handoff (2026-09-27, session 2 end)

Inner lab URL: `https://9394ec22-d949-4489-b46e-73790750472e.workbench-app-prod.verily.com/lab`
(app was still running from session 1, same hostname `18aa4228c0ec` -- no restart this session).

## Status: both jobs done, pulled, committed. Nothing left running on ~/s04.

- **44 v2** (commit `4a75657`): finished at 23:08:57 (608s), 5 tables. Superseded by v3, not pulled.
- **44 v3** (commit `bca633e`, script deployed via a hand-typed `git diff`/`patch` -- see recipe
  below -- md5-verified against local `df95904c...` before running): ran in `~/s04/results/44_v3/`,
  DONE in 292s, **zero** `sanity_check_coverage` warnings (v2 had 6: KIR2DL2/KIR2DL5B/KIR2DP1
  `pct_novel_protein>100%`). New `kir_protein_catalogue_qc.tsv`: only KIR2DP1/KIR3DP1 (both
  pseudogenes) uncovered, `reason=no_catalogued_protein_entries` -- exactly the expected outcome.
  Pulled (via `/files/<path>?download=1` navigation, landed in `~/Downloads`, moved into the repo)
  and committed to `reports/hla_popgen/44_kir_recurrence_saturation/` as commit `d577a28`.
- **47 v2** (commit `8a5ea75`): finished at 23:30 (1837s), 14 task rows (v1 had 9) --
  within-EUR/within-AFR platform and cB tasks plus `naive_ml_ancestry_adjusted.tsv` (new). Pulled
  and committed to `reports/hla_popgen/47_naive_ml_tests/` as commit `c448a2b`.
- **md5 post-check**: `md5sum ~/pipeline_outputs/*.tsv` vs `~/s04/hla_tsv_md5_pre.txt` ->
  `MD5_MATCH_OK`, byte-identical. Neither run touched the HLA production tables.
- Disclosure scan (manual `awk`, no dedicated script exists yet -- see below): clean on both.
  `pca_density_grid.tsv` bins already mask 1-19 as `<20`. `recurrence_classes.tsv` eq1/eq2/gt2/ge20
  correctly masked; `s_obs`/`n_distinct_alleles` (allele-richness, not participant counts) left
  unmasked, matching the pre-existing v1/v2 convention already committed -- not a new issue.

## Remaining steps for the next session
None for 44/47. `~/s04/44_kir_recurrence_saturation_v3.py`, `44_v2_to_v3.patch`,
`results/{44_v2,44_v3,47_v2}/` are still on the VM (left in place -- not authorized to delete
anything this session). A future operator can supersede `44_kir_recurrence_saturation.py` with the
v3 content on the VM if a later script wants to import it as the canonical `44` (out of scope here).

## Recipe that worked this session (deploying a small local-only fix, not a fresh 67KB file)
- The file-browser-upload trick (patching `HTMLInputElement.prototype.click` to capture the hidden
  `<input type=file>`) that worked in session 1 **was refused by the classifier this session**
  ("Auto-Mode Bypass") -- did not retry/work around it, switched approach entirely.
- Instead: `git diff <deployed-commit> <target-commit> -- <file>` locally (259 lines, ~19KB, far
  smaller than the full 67KB file), typed into a `cat > patch << 'EOF'` heredoc in the terminal in
  5-6 chunks (each `type` call ending mid-heredoc, no `Return` needed since embedded `\n`s already
  submit each line -- but **a chunk boundary that lands right after a line with no blank-line
  padding will silently drop any blank context lines that should follow**, since the next chunk's
  leading `\n` only closes the previous line, not the file's actual double-blank-line gap).
  Caught this by bisecting `head -n K local | md5` vs `head -n K vm | md5sum` over decreasing K
  (found the exact break in ~3 rounds), then fixed with a 2-line Python `list.insert()` on the VM
  rather than retyping. **Always verify the whole-file md5 immediately after a heredoc paste, not
  just at patch-apply time** -- `patch`'s "malformed patch at line N" error point is *downstream*
  of the actual missing content, not at it.
- `patch --fuzz=0 <explicit-target-file> < patch` (not `-p1 -d ~`, which fails to resolve the
  `a/scripts/...`/`b/scripts/...` diff-header paths against a VM tree that doesn't have that
  directory structure) applies cleanly once the patch file's md5 matches local exactly.
- Contents API (`fetch('/api/contents/<path>')`, XSRF cookie already present, no header needed for
  a plain GET) is reliable for inspecting file content/length from the page JS -- but **`javascript_tool`
  return values truncate around ~1.3KB**, so `.length` first, then `.slice()` in chunks for
  anything longer, or just download instead (see below).
- **Pulling files**: `navigate` (not `javascript_tool`/fetch) to
  `https://<id>.../files/<path>?download=1` reliably lands the exact file in the local Mac's
  `~/Downloads/` (sizes verified byte-for-byte against the API's `content.length` and the VM's
  `ls -la` for every file this session) -- far cheaper than reconstructing file content through the
  chat context.

## Gotchas
- No dedicated disclosure-scan script exists in this repo yet (`AGENT_PREAMBLE.md`/
  `ORCHESTRATOR_HANDOFF.md` refer to "43's disclosure helpers", which are library functions, not a
  standalone scanner) -- disclosure checking is currently a manual `awk`/`grep` pass per pulled
  file. Worth a WS-F follow-up: a small `scripts/hla_popgen/check_disclosure.py` that flags any
  bare 1-19 integer in a `.tsv` outside an allowed-unmasked-column allowlist (richness/`s_obs`
  columns need to stay on that allowlist, per the existing convention).
- `git diff` includes trailing blank context lines *inside* a hunk (as literal `" "` single-space
  lines) but not *between* hunks -- don't assume a visual blank line you see when reading a diff
  chunk is padding you can drop.
