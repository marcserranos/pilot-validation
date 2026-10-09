# Decisions & Open Questions — RNA-seq / repertoire workstream

> **Role:** the "why" register. Open questions at the top; resolved decisions with rationale
> below. Move an item down when settled — append, don't delete history.
> **Read:** to understand why the pipeline is shaped this way, or before reopening a call.

## Open questions

- **Real research_ids committed to a public GitHub repo, bound to individual-level findings**
  — in at least 6 files (incl. `results/1000291_trust4_smoke_test.md`,
  `../../context/ENVIRONMENT.md`). Flagged 2026-08-17, escalated to Marc, still undecided.
  **Do NOT add git push credentials to the VM** until resolved — paste results back so a
  human reviews every byte leaving the controlled environment. Same class of issue as the
  root repo's "bare person_id in a public repo" open question — don't resolve independently,
  needs a joint call with Marc + supervisors.
- **Ancestry recovery gap: reframed by the full cohort (2026-09-27), not yet tested.**
  At n=100: AFR high / EAS low, 1.63x, "consistent across all 5 major chain types," KW p=0.18.
  At n=7,922 (group means only, `results/rnaseq_cohort_ancestry_summary.csv`) **both halves of
  that framing fail**: (1) EUR, not EAS, is now lowest (AFR/EUR 1.44x; EAS mid-pack); (2) the
  gap is **carried by the B-cell chains** (IGH/IGK/IGL AFR/EUR 1.68-1.75x) while **TRB is
  near-flat** (AFR/EUR 1.06x, all-group spread 1.19x, AFR not even highest). So "reference/
  database bias against non-European receptor sequences" is now the *less* likely reading —
  that would hit TRB too. Leading hypothesis to test, not assert: a biological Ig-expression
  difference (plasmablast/plasma-cell Ig mRNA load; higher serum Ig in African-ancestry
  populations is documented) rather than a TRUST4 artifact. **Not yet a result:** means only,
  no dispersion, no per-person test, no depth normalization at this n. Needs the per-person
  per-chain test (VM-local `batch_summary_detail.tsv` exists) before any claim.
- **Group min/max in committed summaries are single-person values.** `cdr3_min`/`cdr3_max` in
  `results/rnaseq_cohort_ancestry_summary.csv` are each one individual's count. Not an
  identifier, but an n=1 statistic in a public repo — same small-cell disclosure class as the
  root repo's open item. Flag for Marc/supervisors; don't decide here. Leaning: drop min/max
  (or report p5/p95) in future committed summaries. **2026-10-04: the leaning is implemented
  going forward** (`aggregate_rnaseq_results.py` now writes `cdr3_p05`/`cdr3_p95`); the old
  values remain in git history, which is the part still for Marc/supervisors.
- **Rarefaction vs. post-hoc depth normalization for the recovery metric.** Flagged in
  `reference/TRUST4_DEEP_DIVE.md` §8.5 as probably the most important methodological fix
  outstanding — post-hoc normalization only closed 40% of the ancestry gap, which itself
  argues for switching the metric rather than continuing to patch around it. Not yet done.
- **Whether disease burden predicts CDR3 recovery.** New lead (2026-08-17) — both variables
  now exist per-person for the same people. Cheap to check, not yet run.
- **HLA-label attachment timing** — attach Marc's HLA calls (AoU-native + long-read) once the
  ancestry-recovery confound question is far enough along that it won't contaminate the join.
  Not yet started; gated on the full-cohort ancestry-gap result above, not on anything
  external.

## Resolved decisions

- **W0 analysis table (report 09), design decisions (2026-10-08).** One row per person for
  every disease/ancestry analysis; exclusions are flags, never deletions.
  - *Relatedness resolved within our cohort* from AoU's pairwise `samples_relatedness.tsv`
    (kin >= 0.0884, second degree), not AoU's global flagged list, which also removes people
    whose relative is not in our 7,922. Greedy maximal unrelated set keeps, in order of priority,
    the relative who passes the other main filters, then the one with more TRB reads (so a
    pair never loses its only usable member). `family_id` is kept too, so models can use grouped CV and keep relatives.
  - *Read floor = 823 TRB reads*, the same fixed depth as report 02's rarefaction, so every
    depth-normalised metric uses one number. The cost of 250-2,000 is reported by ancestry.
  - *Diversity at fixed depth*: exact rarefied richness, and Shannon / e^H / inverse
    Simpson / public fraction averaged over 10 seeded subsamples of 823 reads. Raw versions
    kept to compare with Cole's PheWAS (which uses raw entropy + total reads).
  - *Sex is a covariate (Female / Male / Unknown), not a filter*: requiring a recorded sex
    removed 13% of AFR but 1.5% of EAS. `has_sex` remains for sex-specific diseases.
  - *Ancestry purity (0.90 / 0.95) only for ancestry-stratified figures*; models use AoU
    PC1-10 instead, so admixed people are not dropped.
  - *Ancestry probability order* (AFR, AMR, EAS, EUR, MID, SAS) is asserted, then checked
    (argmax must reproduce ancestry_pred for >= 99%).
  - *rnaseq_metadata.tsv is profiled, not interpreted*: the draw-date column is chosen after
    seeing its fill rate (Cole: ~3,700 people have an RNA-seq date).

- **Identifiability headline = z-scored (decided 2026-10-06).** Quote report 06 as: SCEPTR
  mean + TRBV usage, **27.0% (95% CI 25.6-28.4%)** of 4,000 people matched to their own
  split-half repertoire, chance 0.025%; TRBV usage alone 16.5% (15.4-17.7%); SCEPTR mean
  alone 11.5% (10.6-12.6%). Reasons: (1) raw cosine on non-negative profiles is dominated
  by the mean direction every person shares, so it rewards whichever block has the larger
  magnitude (raw TRBV usage 23.8% drops to 16.5% once standardised, while the combination
  barely moves); (2) z-scoring makes the comparison scale-free, so raw and block-balanced
  concatenations agree exactly; (3) it is the conservative choice (27.0% < 29.4% raw);
  (4) removing depth leaves 26.9%, so the headline is not a depth effect. Raw and
  depth-removed values stay in the table as sensitivity analyses.
- **Per-person cap stays at 500 (decided 2026-10-06); no rerun.** After the tie-break fix
  the 500 are: every clonotype above the cutoff read count, then a seeded random sample of
  the clonotypes at the cutoff (1 read for 78% of people, 2 for 20%). A random sample of
  the low-count tail is unbiased, so the old objection was to the alphabetical fill, not the
  size. Lowering to ~300 would not remove ties (the median person's tied block is 625) and
  would leave each person's vector noisier. Stability supports 500: top-250 vs top-500
  distance correlation 0.76. The cap now only affects the descriptive pool reports (01, 03
  robustness panel, 05); 02, 04 and 06 use whole repertoires. Methods wording: "up to 500
  TRB clonotypes per person, ranked by read support; ties at the cutoff sampled at random
  (seeded)". This replaces the 2026-09-28 "500 is inherited, robustness check owed" item.

- **Methods review fixes (2026-10-04), from reading the code for the methods artifact.**
  All implemented, tested on synthetic data with planted effects, run on the VM via
  `scripts/run_polish.sh`:
  - *Top-500 tie-break.* Ties at the per-person cap were broken alphabetically by (chain,
    V, CDR3): `groupby` sorts its keys and the read sort was stable. Most clonotypes have
    1-2 reads, so the cutoff is almost always inside a tie and the pool favoured TRBV10-12
    over TRBV7/9. Now a seeded per-person hash (`embed_cdr3s.tie_key`). The old pool and
    embedding are kept in `embeddings/v1_alphabetical_ties/`; report 08 measures the shift.
  - *Report 06 covers whole repertoires.* It looked vectors up from the top-500 pool, so only
    sequences in somebody's top 500 had one. `embed_full_cache.py` embeds every distinct
    (TRBV, CDR3) of the full cache once; 06 uses it and writes the coverage.
  - *Concatenations are block-balanced* (`repfig.balanced_concat`). A raw `np.hstack`
    lets the larger-magnitude block dominate the cosine: in a synthetic check a perfectly
    identifying 48-d block went from 100% alone to 0.1% when hstacked with a 64-d
    large-magnitude noise block, and back to 50% balanced. Raw versions kept for comparison.
  - *Report 07 person axis.* It used k=100 CLR cluster profiles (0.5% in 06), not 06's
    winner as its docstring said. Now mean pooling (embedding only) and mean + TRBV usage,
    both z-scored.
  - *Confidence intervals.* Epitope AUROC / 5-NN in 04 and 07 get an epitope-level bootstrap
    (epitopes are the independent units) with paired differences vs the production model;
    06 identifiability gets Wilson CIs. The earlier "statistically indistinguishable" for
    synthetic_data in EXPERIMENTS.md had no test behind it.
  - *Report 05* adds a Mantel-Haenszel test stratified by TRBV gene, separating CDR3 motif
    from V-gene bias. *Report 04* writes how often the decoy Pgen tolerance had to widen.
  - *Aggregation* assigns chain from V, then J, then C, as the clonotype step does.
- **Correction (2026-10-04): the top-N cap-robustness check does not gate report 02.** 02
  uses the full clonotype cache, not the top-500 pool, so the cap cannot affect the ageing
  result. The cap matters for the embedding pool (01, 03, 05) and is covered by 03's
  robustness panel.

- **TRB is the primary chain for every HLA-facing analysis (Cole Shanks, 2026-09-27).**
  Cole: what matters most are chains that undergo V(D)J recombination, and he's "most
  interested specifically in the TCRBs." TRG/TRD don't really interact with HLA (γδ T cells
  are largely not MHC-restricted) — ignore for now; possibly an interesting disease signal
  later, not an immediate focus. Consequences: `--chain TRB` (already `embed_cdr3s.py`'s
  default) is the headline for embeddings and HLA joins; every ancestry/recovery result is
  reported per chain with TRB first, never only as an all-chain total (the all-chain total
  is dominated by Ig chains and gave the wrong story at n=7,922 — see open question above).
  TRA stays a secondary (it pairs with TRB, but bulk RNA-seq can't pair them). IG chains are
  kept for the recovery-confound analysis, not the HLA aim.
- **All figures use `cnsplots`, journal (Cell/Nature/Science) style (David Bonet,
  2026-09-27).** Current figures read as AI-generated; David asked that every figure from now
  on follow `github.com/faridrashidi/cnsplots` (matplotlib-based, BSD-3, `pip install
  cnsplots`, bundled Claude skill via `cnsplots skill install --agent claude --scope
  project`). Conventions to follow: sizes in points (`cns.figure(width, height)`),
  `cns.multipanel(max_width=...)` with lettered panels, built-in stats annotation
  (`test=`, `p_adjust=`) rather than hand-drawn asterisks, export SVG/PDF with editable text
  (`cns.savefig`), PNG 300 dpi only as a preview. **Where figures get rendered:** figures
  drawn purely from committed group-level CSVs can be built locally; any figure that needs
  per-person values (box/violin/ECDF/scatter) is rendered **on the VM** with cnsplots
  installed there, and only the finished figure leaves — no individual points drawn for
  groups/bins under n=20, no per-person data file ever egresses. So cnsplots goes into
  `aleix/RNA-seq/pixi.toml` (or the VM's base python) as well as locally.

- **SCEPTR adopted as the embedding model (2026-09-10, EXPERIMENTS.md).** Compared against
  ESMC-300M on the same 500-person CDR3 pool: SCEPTR gap (same- vs diff-V-gene cosine
  similarity) = 0.019 in 23s; ESMC-300M gap = 0.0006 (30x worse, embedding space nearly
  collapsed) in 3,666s (156x slower). Not a close call — SCEPTR wins on both signal and cost.
  ESMC is not planned to be re-run at full-cohort scale (would cost ~16 hours vs. SCEPTR's
  ~6 minutes for no discrimination benefit shown so far).
- **Embedding input for HLA/disease work = `b_sceptr` on TRBV + CDR3B, top-500 clonotypes
  per person by read support (2026-09-27).** SCEPTR's docs list a beta-chain-only variant
  (`variant.b_sceptr`); we had been running the paired-chain `default()` on CDR3B alone,
  which also discards CDR1/CDR2 (V-encoded, the loops that contact HLA). `cdr3_only` is kept
  for the V-gene benchmark only, where V as input would make the check circular. Clonotype =
  unique (chain, V, CDR3aa), canonical junction (C…F/W), CDR3_score ≥ 0.02 from `cdr3.out`;
  people without `cdr3.out` are skipped, not silently run unfiltered via `report.tsv` (that
  mix was ~70% of the first full-cohort pool — the shard VMs' `cdr3.out` were never pulled
  back). Chain is assigned from V, then J, then C, so un-V-called TRB CDR3s aren't lost from
  the counts; they are excluded from `b_sceptr` input (it needs a V) and the excluded fraction
  is reported per ancestry.
- **CDR3 length cap ≤ 30 aa, applied before the per-person top-500 (2026-09-28).** Found
  via Fig 1B of report 01: one ≥20-person cell averaged ~100 aa. Full-cohort pool lengths:
  median 14, 99.9% ≤ 22, only 97 clonotypes in 26–40 aa, then a separate mode of 2,036 at
  > 40 aa (0.053% of clonotypes, 1,175 people) — assembly artifacts that cluster together
  in embedding space. 30 aa sits in the empty gap, so it removes the artifact mode without
  trimming the natural tail. Applied pre-cap so affected people get a real clonotype in the
  freed slot; hence a full re-embed rather than dropping rows post hoc.
- **Why top-500 clonotypes per person, by read support (recorded 2026-09-28).** (1)
  Comparability: TRUST4 recovers hundreds to thousands of TRB clonotypes per person with
  depth/T-cell content; a fixed-size sample keeps the person vector from tracking depth
  (median 494 after V exclusion, 96% of people ≥ 400). (2) Top-by-reads keeps the most
  reliably assembled and most expanded (antigen-driven) clones; 1–2-read clonotypes are the
  noisiest part of bulk RNA-seq. (3) The value 500 itself is inherited, not tuned — a
  robustness check across caps (100/250/500/all) is owed before any paper claim.
- **Copy-local BAM staging (not gcsfuse-direct reads) for the TRUST4 batch runner
  (2026-09-10).** gcsfuse streams a BAM at ~11 MB/s and doesn't parallelize (one shared fuse
  pipe); native `gcloud storage cp` pulls the same BAM at ~670 MB/s and bypasses the fuse
  daemon entirely. Local-copy-then-delete makes extraction CPU-bound, so `--jobs` actually
  scales. Full rationale in `scripts/run_rnaseq_batch_local.sh` header and
  `context/ENVIRONMENT.md` quirk 5.
- **Sharded 4-VM batch (n1-highmem-16, `--jobs 8` each) for the full 7,922-person cohort
  (2026-09-19), not a cheaper single-VM run.** Explicit hard constraint from the 3-4 day
  deadline: only accept a cost reduction with zero performance/timeline risk. Sharding turns
  a multi-day sequential job into a same-order-of-magnitude parallel one; collision-safety is
  structural (globally-unique research_id + disjoint `split -l` shards), not best-effort, so
  the parallelism carries no correctness risk.
- **Per-person incremental sync to GCS (not end-of-batch only), added 2026-09-19.** Protects
  against VM preemption/shutdown/loss during a multi-day unattended run — a killed VM loses
  at most the one person it was mid-run on, not the whole shard. Retrofitted after the
  original 500-person run had already completed without it, which is why that batch needed a
  manual backfill (see ENVIRONMENT.md quirk 7) — the gap this decision closes going forward.
- **Results bucket is `gs://aleix-rnaseq-wb-cordial-leechee-9743`, never the disease-counts
  bucket (corrected 2026-09-19, commit `b9f7e65`).** The two buckets have different owners:
  disease-counts = Cole's phenotype/burden deliverable; rnaseq = LR calls + embeddings +
  repertoire results. Caught mid-batch-launch by a direct question ("why in disease counts,
  aren't we storing the rnaseq repertoires in the rnaseq bucket") — the default had been set
  wrong in the script; fixed before most of the batch ran.

Related: [[../../context/DECISIONS.md]] (root, Marc's — the shared public-repo/research_id
disclosure question is tracked jointly there and here).
