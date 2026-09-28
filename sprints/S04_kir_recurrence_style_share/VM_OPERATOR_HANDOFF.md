# VM operator handoff (2026-09-27, session 4 update)

## Session 4: 44 v4b full run -- already DONE when this session checked in, pulled + committed

Inner lab URL (unchanged, same hostname `18aa4228c0ec`, app restarted by orchestrator between
sessions but landed on the same instance id):
`https://9394ec22-d949-4489-b46e-73790750472e.workbench-app-prod.verily.com/lab`

Task handed to this session was "relaunch the dead v4b run" (per the prior session's ERROR log
entry: app stopped mid-run at ~33 min). On arrival: verified `~/s04/44_kir_recurrence_saturation_v4b.py`
md5 `0cdd97eef10d9692560e55fba9d41ed0` == local commit `e4211ad` byte-for-byte; HLA production
tsv md5 vs `~/s04/hla_tsv_md5_pre.txt` -> `MD5_MATCH_OK`. **Found the run had actually already
completed successfully** in `~/s04/results/44_v4b/` (`STATUS.txt`: started 12:50:24, `DONE in
2742s -- 9 tables written` at 13:36:06, ~8h before this session's clock time) -- no relaunch was
needed or performed. `ps aux` confirmed nothing running; `44_v4b_full.log` has zero
error/traceback/invariant-violation hits. `genomic_identity_qc.tsv`: both
`kir_genomic_available=True` and `hla_genomic_available=True` (the e4211ad HLA-prefix fix holds
on the full cohort). Post-run HLA tsv md5 re-checked: still `MD5_MATCH_OK` (nothing wrote to
`~/pipeline_outputs*`).

Pulled all 11 files in `~/s04/results/44_v4b/` via `/files/<path>?download=1` navigation
(md5-verified byte-identical VM vs downloaded vs committed for every file), replacing the v3
tables in `reports/hla_popgen/44_kir_recurrence_saturation/` and adding 4 new v4b-only QC tables
(`diagnostics_identity.tsv`, `artifact_qc.tsv`, `genomic_artifact_qc.tsv`,
`genomic_identity_qc.tsv`). Committed as `9ca7053`. **`README.md` and the 45/46
figures/`recurrence_stats.tsv` in that folder are now stale** (built against v3 data) -- next
session should re-run 45/46 against v4b before trusting the prose numbers in that README.

**Disclosure scan**: no person IDs found (`grep -E '[0-9]{6,}'` hits were only p-values in
scientific notation and large aggregate counts, not IDs). `recurrence_classes.tsv`'s
eq1/eq2/gt2/ge20 masked `<20` as usual (13 eq1 cells, 7 n_distinct_alleles cells masked even
within pooled ALL ancestry -- more than in v3, worth noting if anyone recomputes pooled
sums: coerce `<20` to NaN before summing, don't let pandas silently turn the whole column
to string). The **new** QC tables (`artifact_qc.tsv`, `genomic_artifact_qc.tsv`,
`diagnostics_identity.tsv`'s `n_artifact_*` columns) contain unmasked gene-level integer counts,
some 1-19 -- these are technical/QC tallies (calls flagged with an artifact label per
species+gene), the same category as the already-accepted unmasked `s_obs`
richness convention (LOG's prior ruling: richness/tally counts are not participant identifiers).
Flagged for Marc/Aleix awareness in the completion report; not blocked, but if this convention
is ever formalized in `context/DECISIONS.md`, these new tables should be reviewed against
whatever rule gets written down.

**Headline numbers** (pooled ALL ancestry, from `diagnostics_identity.tsv` +
`recurrence_classes.tsv` + `coverage_chao2.tsv`, masked cells coerced to NaN before summing):

| species | level | S_obs | %novel(singleton eq1) | %singleton | S_obs/Chao2 |
|---|---|---|---|---|---|
| kir | genomic | 38024 | 30178 | 79.4% | 19.2% |
| kir | any_novel | 37503 | 30090 | 80.2% | 18.9% |
| kir | protein | 1441 | 775 | 53.8% | 39.1% |
| kir | protein_novel | 1061 | 702 | 66.2% | 31.4% |
| hla | genomic | 17188 | 13491 | 78.5% | 17.6% |
| hla | any_novel | 14591 | 12468 | 85.4% | 12.3% |
| hla | protein | 1079 | 507 | 47.0% | 39.9% |
| hla | protein_novel | 163 | 119 | 73.0% | 6.5% |

Diagnostics pooled: KIR names=1121, genomic=38024, cds=2100, protein=1444, artifact
clean=151925/fs=709/homopolymer=5932/partial=64874. HLA names=4394, genomic=17188, cds=1363,
protein=1079, artifact clean=640305/fs=0/homopolymer=16036/partial=218838. All
`check_identity_invariants` orderings hold (protein <= cds <= genomic, novel <= baseline, for
both species).

Terminal used: `terminals/2` (opened this session, not reused). Closed at session end per scope.

## Note from orchestrator (received mid-session, addressed): checkpoint/resume contingency

Orchestrator sent local commit `22bf543` (adds `--resume`/`_checkpoints/` to 44) in case the v4b
run needed relaunching after dying again. **Not used** -- the run had already finished
successfully before this session started polling, so no relaunch/deploy of `22bf543` was needed.
Left as a note for whoever runs 44 next: if a future run dies mid-way, deploy `22bf543` (diff +
patch, md5-verify) into `~/s04/results/44_v4c/` instead of `44_v4b`, and never pull
`_checkpoints/` (per-person data, must stay on the VM).

---

# Prior handoff (2026-09-27, session 3 end) -- kept for the v3->v4 bug history below

Inner lab URL: `https://9394ec22-d949-4489-b46e-73790750472e.workbench-app-prod.verily.com/lab`
(app was still running from session 1/2, same hostname `18aa4228c0ec` -- no restart this session).

## Session 3: 44 v4 (true genomic identity) -- BLOCKED at smoke test, do not run full job yet

Deployed `44_kir_recurrence_saturation_v4.py` (local commit `55b3d83`) via the diff+patch-in-terminal
recipe below, md5-verified byte-for-byte against local (`08de53eb6f7ea08833dd3b13c05bc74c`). Input
checks: HLA `hap1.trimmed.fa` present for 12259 haplotypes under `~/pipeline_outputs/people`, KIR
same count under `~/pipeline_outputs_kir`, `~/tools/Immuannot_refdata/gen.fa.gz` present with 1530
KIR-matching headers.

**`--limit 300` smoke test ran to completion (132s, 9 tables, no exception, no
`check_identity_invariants` violation) but exported HLA `genomic`/`any_novel` as `NA` for every
row** (`genomic_identity_qc.tsv`: `hla_genomic_available=False, n_hap_seen_hla=0,
n_hap_trimmed_present_hla=0`, reported reason `no_trimmed_fasta` -- **misleading**, the files do
exist; `KIR2DP1`/`KIR3DP1` protein-catalogue warning is expected/pre-existing, not new).

**Root cause (confirmed, not yet fixed -- out of this session's scope):**
`build_person_true_genomic_identity()`/`build_true_genomic_identity_sets()` filter each GTF
transcript row by `r["gene"] in gene_names`, where `gene_names` for HLA is
`m39mod.CLASSICAL_GENES_BARE` (bare names: `A`, `B`, `C`, `DRB1`, ...). But the HLA GTF's own
`gene_name` attribute is **prefixed** (`HLA-A`, `HLA-B`, `HLA-DRB1`, ...) -- confirmed by direct
inspection of one person's `hap1.gtf.gz`. KIR's GTF uses bare gene names natively, so the same
filter happens to work for KIR (595/595 haplotypes seen in the smoke run) but silently zeroes out
every HLA row. This is a genuine bug in the new species-agnostic code path, not a smoke-sample
fluke -- it will reproduce on the full cohort. Needs a fix (e.g. strip an `HLA-` prefix before the
`gene_names` membership test, or normalize `gene_names` to include both bare and prefixed forms)
before the full v4 run is attempted.

**STOPPED here per task instruction** ("if any level is NA unexpectedly, STOP and report").
`~/s04/results/44_v4_smoke/` (9 tables) and `~/s04/44_v4_smoke.log` are left on the VM for the next
session to inspect. The full run was NOT started. `~/pipeline_outputs/*.tsv` were never touched
(read-only paths; the run only reads `hap{N}.gtf.gz`/`hap{N}.trimmed.fa`, never writes there) -- no
md5 pre/post check was needed since no write path was exercised, but a fresh check is cheap
insurance for the next session before any full run.

**Gotcha this session (new):** typing a large heredoc body across several separate `type` tool
calls into ONE open heredoc silently dropped ~17 lines somewhere in the run (root cause not fully
isolated, suspected a timing/echo race in the xterm.js websocket under sustained large paste).
Fix: write EACH chunk to its OWN small file with its OWN heredoc + immediate `md5sum`/`wc -l`
verification against the local chunk before moving to the next chunk, then `cat` all verified
chunks together and re-verify the concatenation's md5 against the local full diff before applying
`patch`. This caught two more transcription slips (a dropped line, a `+`/space typo) immediately,
each fixed with a targeted `sed` rather than retyping. Also: a large paste can take 10-20s to fully
land in the terminal before the next command is safe to send -- `wait` a beat before reading output
after any multi-KB `type` call, or the next command's output gets misread as still-echoing heredoc
body.

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
