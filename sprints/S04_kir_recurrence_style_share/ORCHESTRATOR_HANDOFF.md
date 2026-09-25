# S04 orchestrator handoff — read this first, then act

*Written 2026-09-26 by the S03 orchestrator (Opus 5.5) for the S04 orchestrator (Opus 5.5).
Branch `s04-kir-recurrence-style-share`, cut from `s03-call8-figures-kir-prediction` @ 5d5ef99 (it
contains S01–S03 and `fig1-drafts-and-research-map`). Owner: Marc Serrano. Nothing is pushed; Marc
pushes and merges.*

This file distils Marc's S04 asks, the operating rules, and a map of where everything lives. It is
the only file you need to read in full. Everything else is pointed to and should be read **by
sub-agents**, not by you.

---

## 0. Your first 15 minutes

1. Read this file. Then skim, do not study: `CLAUDE.md` (root), `sprints/README.md` (layered
   context), `sprints/S03_call8_figures_kir_prediction/SPRINT.md` (board + distilled findings). That
   is enough orientation. **Do not** read LOG.md files, VM_CHANNEL.md, ENVIRONMENT.md or report
   READMEs yourself. Point sub-agents at them.
2. Ask Marc to paste the **authorization block** in §7 (or his own wording). VM state changes need
   consent given to you directly, in chat. Sub-agents will not act on consent you relay, and that
   refusal is correct.
3. Ask Marc for the **figure-style guidelines document** (§2, WS-C) if it is not already at
   `reference/figure_style_guidelines.*`. It was mentioned but never arrived in the S03 session.
4. Fill in `SPRINT.md` (board) and start `LOG.md` (template in §6). Then launch work.

---

## 1. Why this sprint exists (Marc's words, distilled)

S03 called KIR on the full cohort (12,261 people; `reports/hla_popgen/43_kir_full_cohort/`). Before
pushing, Marc wants the KIR result brought up to the HLA standard, a figure-style pass, the data
shared with Cole, and the repo and VM cleaned up. His asks, in order:

- **A. KIR recurrence and saturation (the science).** Repeat for KIR the recurrence and saturation
  analyses already done for HLA, per ancestry and overall.
- **B. KIR vs HLA: which region the reference catalogue covers better.** Show KIR (IPD-KIR) against
  HLA (IPD-IMGT/HLA).
- **C. Figure polish with the supervisor's style guidelines.** A supervisor said the figures look
  "AI-like" and sent guidelines. Apply them to the key figures.
- **D. Explain and extend the prediction work, plus new experiment ideas.** Marc did not
  understand WS6. He wants in-scope experiments worth testing, and naive, easy-to-understand ML and
  simple tests that reveal the nature of the data.
- **E. Share data with Cole via the bucket.** Upload the smallest, most informative, structured
  data package, with documentation. The novel HLA calls were never shared.
- **F. Cleanup.** VM declutter after the upload, repo reorganisation (no loose files), and branch
  unification.
- **G. Reporting.** Everything is logged so that, after the sprint, one agent can build a full
  report: what was done, why, results, figures, trade-offs.

---

## 2. Workstreams (suggested board — you own the final plan)

Numbering convention: `scripts/hla_popgen/NN_description.py` → `reports/hla_popgen/NN_description/`.
The highest number used so far is **43** (43, 43b). Continue from **44**. Numbers are never reused;
gaps are fine.

### WS-A: KIR recurrence and saturation (VM + local), scripts 44/45
Mirror the HLA pipeline, reusing its code rather than re-deriving it:
- HLA recurrence: `03_novel_alleles.py` (novel allele clustering, recurrence across ≥2 unrelated
  people), `34_novel_protein_recurrence.py`, `24_novelty_by_field.py` (known / non-coding /
  synonymous / protein tiers).
- HLA saturation: `04_allele_saturation.py`, `18_discovery_rate.py`, `19_/27_allele_space_coverage.py`,
  and the corrected per-ancestry, equal-N method in `39_saturation_by_ancestry.py` (its README
  documents the S03 methodology fixes; read it before designing).
- KIR inputs: per-person GTFs on the VM at `~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz`
  (12,261 people). The parser and novelty tiers are in `41_kir_pilot.py` (`parse_hap_gtf`,
  `classify_novelty_tier`). The aggregate pattern is `43_kir_full_aggregate.py` (multiprocessing,
  disclosure helpers).

What Marc asked to see:
- **Discovery/saturation curves** overall and per ancestry (AFR AMR EAS EUR MID SAS; equal-N as in
  39), for KIR.
- **Saturation stratified by recurrence class**:
  - seen once (singletons);
  - seen twice;
  - seen more than twice;
  - seen at least N times. The threshold is your choice; ≥10 or ≥20 unrelated carriers is
    suggested. **≥20 has the advantage that it coincides with the disclosure floor.**

  Show this overall and for a few ancestries.
- **One more informative quantity** is welcome. Candidates:
  - expected undetected alleles (Chao1/Chao2);
  - Good–Turing sample coverage, i.e. "fraction of the allele space explored";
  - an iNEXT-style extrapolation to 2× N.

  Pick one or two and justify them in the README.
- **Keep both novelty levels visible.** Marc explicitly wants the nuance between general
  neo-alleles (any sequence difference) and new protein-level alleles. Show both together, e.g.
  paired panels or two line styles on one axis. Do not drop either.

### WS-B: KIR vs HLA catalogue coverage, script 46 (mostly local once A exports)
- Put KIR and HLA saturation, recurrence and novelty on one comparable footing: same unrelated set,
  same ancestry scheme, same subsampling, same recurrence classes.
- Answer: which region is better represented in its reference (IPD-KIR 1,530 genomic alleles vs
  IPD-IMGT/HLA)? Metrics:
  - % of calls novel (any level vs protein);
  - sample coverage;
  - slope at equal N;
  - Chao-estimated unseen.
- Report per gene where possible. Per S03: KIR novelty is 58.9% of calls, mostly non-coding; HLA
  classical-gene novelty is much lower. Quantify this properly.
- Candidate Figure 1/2 panel for the paper. Coordinate its style with WS-C.

### WS-C: figure style pass (local)
- **Input:** the supervisor's guidelines document (ask Marc; store it at
  `reference/figure_style_guidelines.*`; if it is a PDF, have a sub-agent read it and distil it
  into `reference/FIGURE_STYLE.md` as ≤60 lines of concrete, checkable rules).
- **Then:**
  - Encode the rules in `scripts/hla_popgen/_viz_common.py` (`nature_style()`, palette, fonts,
    `save_fig`) so every figure inherits them.
  - Re-render the key figures:
    - Figure 1 v5 (`40_figure1_v5.py compose`);
    - DQ G1/G2 (37/37c);
    - saturation (39);
    - KIR full cohort (43b);
    - the new WS-A/WS-B figures.
  - Put before/after PNGs in the sprint folder for Marc.
- **Known polish debt:**
  - 43b panel f: its y-axis starts at 86%, which exaggerates the bar differences; use points with
    CIs or a zero-based axis.
  - 43b panel c: the "3DP1 †" label collides with panel b's n= text.
  - Figure 1: panel b uses a 0.98 ancestry threshold, while panel c uses 0.9.
  - The figure-editor critique in `sprints/S03_call8_figures_kir_prediction/CRITIC_2.md`.
- **"Less AI-like"** usually means:
  - fewer panel titles that are sentences;
  - no redundant legends;
  - direct labelling instead of legends;
  - restrained colour (one accent);
  - no decorative gridlines or text boxes;
  - honest axes;
  - consistent number formats;
  - no emoji or "†" clutter.

  But **the guidelines document overrides this list**.

### WS-D: prediction baseline explained, plus experiment ideas (local; VM optional)
- **Explain WS6 in plain language** in the sprint report. Plain explanation:
  - We asked whether a person's immune repertoire (the T-cell receptor sequences TRUST4 pulled out
    of their blood RNA-seq, from Aleix's workstream) predicts their traits.
  - We used the BenchRep-T recipe: count V/J gene usage and short amino-acid "words" (4-mers) in
    each person's receptors, then train a simple L1-logistic regression with 3-fold
    cross-validation. AUROC 0.5 means chance; 1.0 means perfect.
  - Positive controls: ancestry is predictable (AFR AUROC 0.82). This is a warning that ancestry
    leaks into any disease signal. Sex is barely predictable (0.55).
  - No disease had enough cases in the 1,500-person subsample. The full 7,640-person run is still
    pending.

  Code and report: `42_repertoire_baseline.py`, `reports/hla_popgen/42_repertoire_baseline/`.
- **Naive ML and simple tests Marc asked for** (easy to understand, reveal data structure).
  Examples to evaluate; pick what is feasible with data on hand:
  - HLA/KIR genotype → ancestry with a one-feature-per-allele logistic regression or a small tree
    (which alleles are most ancestry-informative?);
  - PCA/UMAP of per-person HLA and KIR allele-carriage vectors, coloured by ancestry and by
    sequencing platform. Platform clustering would expose batch effects;
  - predict **sequencing platform or assembly quality** from the calls. If that is predictable,
    there is a technical artifact;
  - predict KIR gene content (cA/cB) from HLA genotype. KIR and HLA co-evolve through ligand
    pairs such as HLA-C1/C2 with KIR2DL1/2DL2/3, and HLA-Bw4 with KIR3DL1;
  - nearest-neighbour "does a person's closest genetic neighbour share their rare alleles?", a
    sanity check of relatedness filtering;
  - permutation baselines for every one of the above.
- **In-scope original experiments** (seed list; a sub-agent should rank them for feasibility,
  novelty and relevance to Cole's paper):
  - **KIR–HLA ligand co-occurrence by ancestry:** are functional receptor–ligand pairs (KIR3DL1 +
    HLA-Bw4; KIR2DL1 + C2; KIR2DL2/3 + C1) enriched above independence in AoU? Long reads give
    both regions from the same person;
  - **KIR novel-allele recurrence vs HLA:** is KIR novelty mostly private (singletons) or shared?
    This comes out of WS-A;
  - **KIR2DL2/KIR2DL3 co-occurrence (0.8%, expected ~0):** a per-haplotype look to separate
    miscall from duplication from phase error. This is an S03 open item;
  - **KIR framework-gene misses:** real deletion haplotypes vs miscalls (S03 open item);
  - HLA-C1/C2 and Bw4/Bw6 epitope frequencies per ancestry. Cheap from existing HLA calls, and
    feeds the ligand analysis;
  - WS6 full-cohort run (7,640 people plus XGBoost) and the HLA-carriage positive control. The
    control was skipped because this VM's `hla_calls_rich.tsv` lacked expected columns; check that.

### WS-E: data package for Cole in the share bucket (VM)
- **Bucket:** `gs://hla-calls-share-wb-cordial-leechee-9743`. It is inside the perimeter,
  workspace-owned and not requester-pays. Current layout: `aggregate/` holds small TSVs;
  `pipeline_outputs/` holds the per-person tree. Documentation: `context/ENVIRONMENT.md`, the "Share
  bucket" table around line 70.
- **Already shared:** two-field common HLA alleles (Cole used them).
- **Not yet shared:**
  - novel HLA calls (with novelty tier and artifact label from 24/34);
  - KIR calls;
  - ancestry proportions and labels;
  - unrelated-set membership;
  - the DQ/DP phasing-confidence flag (per person; it lives only at `~/s03/results/37/` on the VM).
- **Principle:** the smallest set of files carrying the most information for pop-gen/stat-gen.
  - One tidy long table per locus family, one row per person × haplotype × gene: allele, novelty
    tier, artifact flag, contig/phase info.
  - A person table: ancestry probabilities, predicted ancestry, platform, unrelated flag.
  - A `README.md` + `SCHEMA.md` at the bucket root (versioned, dated) describing every file and
    column and the known caveats (e.g. KIR excludes 991 sequel2 people; novel alleles are "dangerous"
    per Cole — flag them clearly).
  - Use a dated prefix, e.g. `release_2026-09/`, with checksums.
- **Disclosure:** this is controlled-tier data moving between in-perimeter locations, which is not
  egress. It is still participant-level: never copy it into git, a report or a chat reply. Only
  aggregates leave the VM for the repo. The small-cell publication question is still open in
  `context/DECISIONS.md`.
- **Needs Marc's consent** (shared-state write). State exactly which files and sizes before
  uploading.

### WS-F: cleanup and unification (last)
- **VM declutter**, only after WS-E is uploaded and verified.
  - Candidates:
    - `~/s03/` scratch (`_patch*.py`, old logs);
    - `~/pipeline_outputs_kir_test/`;
    - benchmark cohort files;
    - old WS6 scratch;
    - `*.corrupted-backup*` in `~/pipeline_outputs/`.
  - **Never** touch `~/pipeline_outputs/*.tsv` (the HLA production tables) or
    `~/pipeline_outputs_kir/`, unless Marc says so after the bucket copy is verified.
  - Write a dry-run list (paths + sizes) first, then get Marc's OK. Deletion is irreversible.
  - The VM repo `~/repos/pilot-validation` is on branch `needle-view-cds-diversity-density`, owned
    by other work. It has one untracked file from S03:
    `scripts/production_orchestrator/run_production_orchestrator_s03kir.py`. It can go once this
    branch is merged and the VM repo is updated.
- **Repo reorganisation.** Marc: "a lot of files are outside their folder".
  - Audit: loose scripts at `scripts/` top level (≈40 experiment_* / run_* / analyze_* files from
    before the numbered convention), `outputs/` (legacy), loose sprint-level files
    (`sprints/CALL_*`, `ALEIX_HANDOFF_*`), `scripts/hla_popgen/{NEEDLE_VIEW_BRIEF.md,RUNBOOK.md}`.
  - Propose a target layout in the LOG and apply it with `git mv`, so history is preserved.
  - Update every path reference: README.md, CLAUDE.md, context/*, script docstrings, tests.
  - Run the test suites after the moves.
  - Keep the four-tier doc system (`context/`, `reference/`, `scripts/`, `reports/`) intact.
- **Branch unification.** At the end, this branch holds everything. Leave one clean branch; Marc
  pushes and merges. Do not merge `aleix/*` branches; they are Aleix's. Check
  `context/DECISIONS.md` and STATUS for concurrent-editing notes first.

### WS-G: reporting (throughout, final)
- See §6. The final deliverable is `sprints/S04_*/FINDINGS_FOR_MARC.md`, written as a PI-style
  report: what ran, why, results with figures, costs, problems and how they were resolved,
  trade-offs, open items. Use an Artifact or document if Marc wants a shareable page.

---

## 3. Operating rules (hard-won in S01–S03; follow them)

### Context and agents
- **You orchestrate. Sub-agents implement.** Use Sonnet 5 sub-agents (`model: sonnet`) for code,
  analysis, figures, critique and doc reading. Each returns a compact result of ≤300 words: key
  numbers, file paths, commit hash, errors and fixes, open issues. **Never** read long logs,
  transcripts or big files yourself.
- **At most 3 agents at once.** More than 3 concurrent VM agents exhausted the account usage limit
  in about an hour in S03 (three outages).
- **Every agent brief is self-contained:** repo path, branch, the files to read (point to
  `AGENT_PREAMBLE.md` in this folder), exact outputs, commit rules, and return format.
- **One fresh-context critic agent per result folder before it is called done.** Format:
  `sprints/S03_*/CRITIC_2.md`. The critic checks disclosure, number traceability, figure
  legibility and overclaiming. The critic fixes text and figure issues and documents judgment
  calls. **Critics catch real errors** (S03: fabricated-zero cells, circular EM, an unreconciled
  "104").
- Review every agent's headline claims against its own curves and tables. S03 caught an agent
  claiming "near saturation" while its curves were still rising.

### Polling and cost
- Long VM jobs: launch with `setsid nohup`, write a small VM-local **status file**
  (aggregate counts only), and check it with background `sleep` waits of 1–3 h. Do not tail logs.
- VM-side **waiter scripts** chain steps without idle billing: wait for the PID to exit, then run
  the next step (S03: `launch_full_after_bench.sh`, `after_full.sh`).
- The app has a 1 h idle autostop. A running process counts as activity. Assume the VM restarted
  at the start of every session: remount gcsfuse and verify with `ls` before trusting reads.

### VM access (Chrome tools + JupyterLab REST/websocket)
- **Recipe:** `sprints/S03_call8_figures_kir_prediction/VM_CHANNEL.md` (the terminal websocket
  `bg`/`poll` helpers). ENVIRONMENT.md quirks #31–37 cover the rest.
- **Load the tools** with one ToolSearch: `select:mcp__claude-in-chrome__tabs_context_mcp,...navigate,
  computer,javascript_tool,get_page_text,find,tabs_create_mcp,tabs_close_mcp`. Use Marc's
  logged-in Chrome. You never sign in; if a Sign In page appears, ask Marc.
- **Getting into JupyterLab:** Workbench → workspace `full-cohort-hla-calling` → Apps → app
  `AoU_Jupyter_ComputeEngine_20260805_big_run` (the app-id has been `9394ec22-…` but can change).
  Navigate to the wrapper `workbench.verily.com/app/<id>/` once (it sets the cookie), then to the
  inner `https://<id>.workbench-app-prod.verily.com/lab`.
- **S03 quirks to know:**
  - the terminal websocket goes stale after ~1 h idle, so reopen a terminal on each check;
  - JS tool return values truncate at ~1.3 kB, and results that look like `KEY=value;` are
    blocked. Replace `=`/`;` before returning;
  - **pull files** via top-level navigation to `https://<id>.../files/<path>?download=1`. They land
    in `~/Downloads`; move them into the repo immediately and disclosure-check them;
  - **deploy files as plain text:** `PUT /api/contents/<path>` with JSON `{"type":"file","format":"text","content":...}`.
    **Never base64.** It got an agent's browser tools permanently refused. For large files, send a
    `diff` and `patch` it on the VM, then verify the md5 against the local copy;
  - `pixi run` buffers output, so use `PYTHONUNBUFFERED=1` / `python3 -u`;
  - the Workbench "Start" button sometimes needs a second click. The app card shows stale state;
    reload the page;
  - GCP capacity stockouts happen: n2-highcpu-96 was unavailable in us-central1-a for over 1.5 h.
    Marc chose n2-highcpu-80 ($2.98/h). Never switch machine family without Marc.
- **VM repo discipline:** the VM repo is on another branch owned by other work. Never checkout,
  commit or pull there. Deploy S04 files to `~/s04/` and run with
  `PYTHONPATH=~/s04:~/s03:~/repos/pilot-validation/scripts/hla_popgen`.
- **Classifier refusals** ("Auto-Mode Bypass", "Blind Apply"): do not work around them. One plain
  retry of a read-only action is fine. Otherwise stop and tell Marc.

### Data diligence (controlled-access All of Us data — non-negotiable)
- Only **aggregate** outputs leave the VM. No person IDs in files, logs you print, commit messages
  or chat. Mask IDs in any grep you run (`sed -E 's/[0-9]{6,}/<ID>/g'`). Listing
  `~/pipeline_outputs*/` prints IDs, so count instead.
- **Counts 1–19 are written as `<20`** and parsed as censored, **never 0**. A true 0 stays `0`.
  Blank any rate whose numerator or denominator is 1–19, because rate × denominator can
  back-reveal a count. Censored cells are hatched in figures, never drawn as 0. S02 shipped a
  fabricated-zero Figure 1 row; see memory `feedback_suppressed_counts_are_not_zero`.
- Run a disclosure check (a script that scans every pulled TSV for 1–19) before every commit.
  Include derived fractions: S03 caught a QC fraction equal to 1/1500.
- **Counts of alleles** (e.g. "number of novel alleles seen in exactly one person") are not counts
  of participants. The earlier convention (36's recurrence histogram) treated them as releasable,
  but that is a judgment call. An allele seen once identifies one person's rare variant. Keep it
  aggregate: never an allele name next to a count of 1–19 carriers. If unsure, flag it for Marc
  (open question in DECISIONS.md).
- **Never write to HLA output paths** (`~/pipeline_outputs/*.tsv`). Before any VM job that could
  touch them, record `md5sum` of `~/pipeline_outputs/*.tsv` and verify afterwards. S03 did this;
  the pre-run manifest is at `~/s03/kir_full/hla_tsv_md5_pre.txt`.
- The repo is **public on GitHub**. Participant-level data never goes in git.

### Git
- Commit locally with **explicit paths** (`git add <paths>`, never `-A`, never `.`). Sub-agents
  commit concurrently in a shared index.
- Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Never push.** Marc pushes. `git commit` over the VM terminal has been classifier-blocked
  before; hand the command to Marc.

---

## 4. Map: where everything lives

| What | Where |
|---|---|
| Rules of engagement | `CLAUDE.md` (root) |
| Doc tiers | `README.md`: `context/` (read fully), `reference/` (on demand), `scripts/` (use, don't read), `reports/` (on demand) |
| Live state / next commands | `context/STATUS.md` (rewrite at sprint end) |
| Environment and VM quirks | `context/ENVIRONMENT.md` (quirks #1–41+; share bucket around line 70) |
| Decisions + open questions | `context/DECISIONS.md` (small-cell publication question is open) |
| Run log (append-only) | `context/EXPERIMENTS.md` (one pointer entry per sprint) |
| Sprint convention | `sprints/README.md` (L0 STATUS → L1 SPRINT/LOG → L2 WS briefs → L3 scripts/reports) |
| Previous sprint | `sprints/S03_call8_figures_kir_prediction/`: SPRINT.md, LOG.md, FINDINGS_FOR_MARC.md, KIR_FULL_RUN.md, VM_CHANNEL.md, CRITIC_1/2.md |
| Call notes that drive the research | `sprints/CALL_SUMMARY_2026-09-22.md`, `sprints/ALEIX_HANDOFF_2026-09-22.md` (Aleix/Cole context: two-field alleles shared; novel alleles "dangerous"; KIR to be shared) |
| Research map | `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` (if present) |
| Pipeline scripts | `scripts/hla_popgen/NN_*.py` (README.md = pipeline table), shared helpers `_viz_common.py`, schema `SCHEMA.md`, tests in `scripts/hla_popgen/tests/` |
| Production orchestrator | `scripts/production_orchestrator/run_production_orchestrator.py` (has `--region/--pad/--out-suffix`; `person_done` needs both haps with ≥1 call) |
| Per-person caller | `scripts/run_immuannot_person.py` (Tier-3 self-align hardcodes a chr6 cache, which is the KIR sequel2 bug) |
| KIR results | `reports/hla_popgen/41_kir_scoping/` (pilots), `43_kir_full_cohort/` (full cohort, 6 aggregate TSVs + figure) |
| HLA novelty / recurrence / saturation | `reports/hla_popgen/{03,04,18,19,24,27,34,39}_*/` |
| Figure 1 | `reports/hla_popgen/40_figure1_v5/` (compose: `40_figure1_v5.py compose --panel-a-bins … --panel-c-totals …`) |
| VM data | `~/pipeline_outputs/` (HLA: `hla_calls_rich.tsv`, `hla_cis_pairs.tsv`, `cohort_membership.tsv`, `immuannot_cohort_full.tsv`), `~/pipeline_outputs_kir/` (KIR per-person GTFs, 75 GB), `~/s03/` (S03 scripts + results), relatedness at `~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv`, AoU bucket mount `~/mnt/aou-controlled` (gcsfuse recipe: ENVIRONMENT quirk #11) |
| App config | `AoU_Jupyter_ComputeEngine_20260805_big_run`: currently **n2-highmem-4, 2000 GB, autostop 1 h, Stopped**. Resizes are Marc-authorized per session (record the before-config, restore after) |
| Memory (auto) | `~/.claude/projects/-Users-marcserrano-WORK-STANFORD-pilot-validation/memory/` (index MEMORY.md; `project_s03_kir_full_run.md` has the S03 ops lessons) |

---

## 5. Suggested sequencing (adapt freely)

1. **Local, parallel:**
   - WS-C step 1: get the guidelines and distil them into `reference/FIGURE_STYLE.md`;
   - WS-D ideation: rank the experiments and naive-ML list; write the WS6 plain-language
     explainer;
   - WS-A design: read 03/04/24/34/39 and write the 44/45 spec + tests.
2. **VM session 1** (4 vCPU is probably enough; aggregation is I/O-bound. Ask Marc before any
   resize):
   - deploy 44/45 and run the KIR recurrence and saturation exports;
   - run matching HLA exports if 39's committed outputs aren't comparable enough;
   - run cheap experiments that need per-person data: ligand co-occurrence, naive ML on the VM
     (only aggregate metrics leave).
3. **Local:** WS-B comparison figures; WS-C re-render of the key figures in the new style; one
   critic per result.
4. **VM session 2:** WS-E bucket package (build, checksum, upload after Marc's OK, verify, write
   the bucket README). Then WS-F VM cleanup (dry-run list → Marc's OK → delete).
5. **Local:** WS-F repo reorganisation (a single agent, since `git mv` in a shared index needs one
   writer), tests, STATUS/EXPERIMENTS/DECISIONS updates, FINDINGS_FOR_MARC. Unify the branch.

---

## 6. Reporting protocol (so a later agent can build the final report)

- `sprints/S04_*/LOG.md`: append-only, `[time] TYPE — what / why / outcome`, with TYPE ∈ PLAN,
  RUN, RESULT, ERROR, FIX, DECISION, IDEA. Log every plan, launch, result (with numbers), error,
  fix and decision, **including trade-offs** and why. This is the source for the PI report.
- `SPRINT.md`: board (WS | script | ask | needs VM | status) plus a ≤10-line "Distilled" entry
  per WS. Keep it current.
- Each result folder `reports/hla_popgen/NN_*/README.md` covers:
  - question;
  - method, with the statistics named;
  - results, with numbers;
  - how to read each figure;
  - caveats;
  - a ≤10-line Distilled section.
- Figures are PNG (600 dpi) + PDF, in the result folder. **List every figure** with a one-line
  caption in `sprints/S04_*/FIGURES_INDEX.md`, so the final report agent can assemble them.
- Costs and VM time go in a `VM_RUNS.md` (like S03's `KIR_FULL_RUN.md`): config before/after,
  start/end times, $.
- At sprint end:
  - `FINDINGS_FOR_MARC.md` (PI-style, plain language);
  - `context/STATUS.md` rewritten;
  - one pointer entry in `context/EXPERIMENTS.md`;
  - durable calls in `context/DECISIONS.md`;
  - new quirks in `context/ENVIRONMENT.md`.
- **Final message to Marc:** what ran, cost, time, results, problems and how they were resolved,
  and what he must do (push/merge).

---

## 7. Authorization block for Marc to paste into the new session (draft; edit freely)

> I, Marc (data owner), authorize you, the orchestrator, directly for sprint S04:
> 1. start/stop the app `AoU_Jupyter_ComputeEngine_20260805_big_run`, and resize it if a job
>    needs it (record the config before; restore it after);
> 2. deploy scripts to `~/s04/` and run aggregation and analysis jobs on the VM;
> 3. build a data package for Cole and, after telling me the exact files and sizes, upload it to
>    `gs://hla-calls-share-wb-cordial-leechee-9743/` under a dated prefix;
> 4. after I confirm the upload, delete the VM scratch paths you list for me (never the HLA
>    production tables or `~/pipeline_outputs_kir/` unless I say so).
> Same rules as S03: aggregates only leave the VM, `<20` censoring, explicit-path local commits,
> never push.

---

## 8. Open items inherited from S03 (fold into the board)

- KIR2DL2 + KIR2DL3 co-occur on 0.8% of haplotypes; expected ~0.
- KIR novelty is not yet recurrence-validated (this is WS-A).
- 991 sequel2 people are missing from KIR. Fix `ensure_chr6_ref()` / `--chr6-ref-cache` for
  non-chr6 regions in `run_immuannot_person.py`, then run them. This needs a VM run of ~1 h at
  80 cores. It is optional for S04; ask Marc.
- WS6 full cohort (7,640 people + XGBoost) and the HLA-carriage positive control.
- Figure 1: the ancestry threshold differs between panel b (0.98) and panel c (0.9).
- The small-cell publication policy is still open (DECISIONS.md).
