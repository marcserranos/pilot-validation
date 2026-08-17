# After TRUST4 — what we can do with a bag of CDR3s

**Written 2026-08-17. Deliberately shallow** — a map of the options and their entry costs,
not a recommendation to build any of them yet. Depth comes after we choose.

The input to everything below is what TRUST4 gives us per person: **a list of receptor
sequences with abundances, unpaired, from blood, at one timepoint.** Every option is
constrained by that sentence.

---

## 0. The constraint that shapes every option

**Bulk RNA-seq gives unpaired chains.** A real T-cell receptor is an **αβ heterodimer** —
one alpha chain and one beta chain, physically bound, jointly determining what the cell
recognizes. Single-cell methods keep them together via cell barcodes. Bulk data does not:
we get a bag of alphas and a bag of betas from thousands of cells, mixed.

Consequences, stated once and then assumed throughout:
- Anything requiring paired α/β is out unless it has a β-only mode.
- **TRB (beta chain) is the workhorse.** It carries more of the recognition specificity than
  alpha, and essentially all β-only tooling exists because bulk data is common.
- **IGH (heavy chain)** is the B-cell equivalent.
- Practical default: **build everything on TRB, treat IGH as the parallel B-cell track, use
  TRA/IGK/IGL as supporting signal only.**

---

## 1. Tier 1 — repertoire summary statistics (no model, immediate)

Collapse each person's whole repertoire into a handful of numbers. Every one is a candidate
feature for a downstream risk model, and they're computable today from files we already have.

| Metric | What it means intuitively |
|---|---|
| **Richness** | How many distinct sequences. Depth-dependent — needs rarefaction. |
| **Shannon entropy** | Diversity accounting for evenness: many clones at similar size = high. |
| **Clonality** (1 − normalized Shannon) | The inverse: is the repertoire dominated by a few big clones? **High clonality = an active or past strong immune response** (or age/immunosenescence). Probably the single most interpretable number here. |
| **Simpson / inverse Simpson** | Probability two random sequences are the same clone. Weighted to big clones. |
| **Chao1** | Estimate of the *true* richness including unseen clones. |
| **D50 / top-N fraction** | How many clones make up 50% of the repertoire. Blunt but very readable. |
| **T:B ratio, chain composition** | Crude cell-composition readout. **Better replaced by RSEM TPM of `CD3E`/`MS4A1`** (see data report §3.2 / item N5). |
| **Mean `CDR3_germline_similarity`** (IGH) | Somatic hypermutation load = how much affinity maturation. Real biology, currently unused. |
| **V/J gene usage frequencies** | ~50 TRBV × ~13 TRBJ → a fixed-length usage vector. **Depth-robust** (it's a proportion), unlike richness. Quietly one of the best feature sets available. |
| **CDR3 length distribution** | Shifts with thymic selection and with HLA. |

**Cost: near zero.** Pure pandas over files that already exist. **Standard tooling:**
`immunarch` (R), `VDJtools` (Java, and `report.tsv` is explicitly VDJtools-compatible),
`scirpy`/`pyrepseq` (Python).

**The catch, and it is the same catch as everywhere:** richness-family metrics are strongly
depth-dependent. **Rarefaction to a common depth is mandatory**, not optional. Usage
frequencies and length distributions are much safer.

---

## 2. Tier 2 — sequence clustering / motif discovery (no ML training)

Group CDR3s that look alike, on the theory that **similar CDR3s recognize similar
antigens**. Turns a bag of unique strings into a much smaller set of interpretable clusters,
which can then be counted per person.

| Tool | Approach | β-only OK? |
|---|---|---|
| **TCRdist3** | Biophysically-weighted distance between CDR3s; well-established | Yes |
| **GLIPH2** | Finds shared local motifs; explicitly models HLA association | Yes |
| **GIANA / clusTCR** | Fast approximate clustering, built for large repertoires | Yes |
| **VDJdb / McPAS-TCR / IEDB lookup** | Match our CDR3s against databases of sequences with **known antigen specificity** (CMV, EBV, influenza, SARS-CoV-2, autoantigens) | Yes |

**The database-lookup option deserves a flag: it is the cheapest possible source of
interpretable, biologically meaningful features.** "This person has 47 CDR3s matching known
CMV-reactive sequences" is a feature a clinician can read, and it directly encodes
*exposure history* — the axis that genotype structurally cannot reach. Downside: coverage
is thin and skewed toward well-studied common viruses in mostly-European cohorts, so
**match rates will likely be lower in non-EUR ancestry groups** — which, given our open
ancestry-recovery question, is a bias we would be importing on top of a bias we haven't
resolved. Worth knowing before leaning on it.

**GLIPH2 is notable** for explicitly integrating HLA — it is the closest existing tool to
the repertoire↔HLA link that our workstream is built on.

---

## 3. Tier 3 — embedding models (the direction the supervisor papers point at)

Turn each CDR3 (with V gene) into a fixed-length numeric vector, so that
**sequence similarity becomes geometric proximity** and standard ML can operate on it.

### SCEPTR — the model from the paper you were given
*Simple Contrastive Embedding of the Primary sequence of T cell Receptors*, Nagano et al.,
**Cell Systems, Jan 2025**.

Concrete facts, from its own documentation:
- BERT-style transformer, **only 153,108 parameters** — tiny. Runs fine on CPU.
- Trained with **autocontrastive learning + masked language modelling**. "Autocontrastive"
  = the model is shown two corrupted views of the *same* TCR (residues deleted, a chain
  dropped, dropout noise) and taught to place them close together, while pushing different
  TCRs apart. **Critically, dropping a chain is one of the training augmentations** — which
  is exactly why the model tolerates β-only input.
- Input: a DataFrame with `TRAV`, `CDR3A`, `TRBV`, `CDR3B`. **"Incomplete rows are
  allowed (e.g. only beta chain data available)"** — this is the sentence that makes SCEPTR
  usable on bulk data at all.
- Output: **64-dimensional** vectors (default model); 16-dim for the `tiny` variant.
- API: `calc_vector_representations()`, `calc_cdist_matrix()`, `calc_pdist_vector()`,
  `calc_residue_representations()`. `pip install sceptr`.

**Fit to our data: good, with one caveat.** We have `CDR3B` and `TRBV` from TRUST4 — exactly
the β-only case SCEPTR supports. The caveat is that SCEPTR was benchmarked on *paired* data
and its headline results are for the paired setting; β-only performance is a documented
capability, not the showcase. We would be using it in its degraded mode. That's fine and
honest — just not the same thing as the paper's numbers.

### Alternatives, briefly
- **TCR-BERT, catELMo, TCR2vec** — earlier protein-language-model approaches to the same
  problem. Generally larger and not obviously better; a 2025 review compares them.
- **ESM-2** — general protein language model, not TCR-specific. Bigger, slower, and CDR3s
  are short and hypervariable in ways general PLMs handle poorly. Not recommended as a first
  choice.
- **antiBERTy / AbLang / BALM** — the **B-cell** equivalents, for IGH. If we pursue the BCR
  track (and IGH is our highest-count chain) these are the analogues to SCEPTR.

### The aggregation problem — the genuinely unsolved bit

SCEPTR embeds **one receptor**. We need **one vector per person**, from thousands of
receptors. That mapping is a real modelling decision, not a detail:

| Approach | Note |
|---|---|
| Mean / abundance-weighted mean of embeddings | Trivial. Probably washes out exactly the rare specific clones that carry the signal. |
| Cluster embeddings, use cluster-occupancy histogram | Fixed-length, interpretable, robust. **Likely the sane default.** |
| Multiple-instance learning (**DeepRC**, **DeepLION2**) | Treats a repertoire as a bag of instances and learns which instances matter, with attention. Purpose-built for exactly this. Needs labels and real training data. |
| Attention/set-transformer pooling | Most flexible, most data-hungry. |

**This is where a nonlinear model actually earns its keep** — the person-level signal is
plausibly "presence of a few specific clone families," which is a genuinely nonlinear
function of the repertoire and is invisible to any linear summary.

---

## 4. Tier 4 — direct supervised repertoire classification

Skip hand-built features; train end-to-end from repertoire → label.

- **DeepRC** — the canonical multiple-instance-learning architecture for repertoire
  classification.
- **DeepLION2** — MIL + contrastive learning + motif attention; cancer-associated TCRs.
- **The precedent that matters:** Emerson et al.'s CMV-status classifier from TRB repertoires
  is the field's proof that a **repertoire alone can predict an exposure/serostatus label**,
  which is the exact shape of what we'd be attempting.

**Requires labels and n.** With 8,980 RNA-seq people total, and far fewer in any specific
disease, this tier is only reachable for **common** phenotypes. **Quantifying exactly which
phenotypes clear that bar is Task 2's job** — which is why Task 2 is not a side quest, it is
the feasibility gate for this entire tier.

---

## 5. The natural sequencing of work

```
  Tier 1 summary stats ──> validate against HLA labels ──> Tier 3 embeddings ──> Tier 4
   (today, free)            (the honest check)             (if Tier 1 holds up)
        │
        └── Tier 2 database lookup (cheap, interpretable, do in parallel)
```

**The validation step is the load-bearing one.** Before trusting any embedding for anything
harder, check whether the repertoire predicts **HLA type** — which we can do, because AoU
ships HLA calls for free, and because the published result (AUC 0.95 from β chain) tells us
what success looks like. If our data can't reproduce a known-true relationship, no amount of
model sophistication will rescue a harder one.

**And the standing precondition, unchanged:** if repertoire recovery is itself technically
uneven by ancestry — still open at p = 0.18 — then any model built on it risks learning a
measurement artifact. Rarefaction (TRUST4 deep-dive §8.5) plus the untested QC confounds
(data report §8, N4) are the way to close that, and both are cheap.

---

## Sources
- Nagano, Y. et al. *Contrastive learning of T cell receptor representations.* Cell Systems
  (2025). <https://www.cell.com/cell-systems/fulltext/S2405-4712(24)00369-7> ·
  docs <https://sceptr.readthedocs.io/> · code <https://github.com/yutanagano/sceptr>
- *TCR representation learning with protein language models: a comprehensive review* —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC12802949/>
- *Machine Learning Approaches to TCR Repertoire Analysis* (Frontiers in Immunology, 2022) —
  <https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.858057/full>
- DeepLION2 — <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10956474/>
- *A comprehensive evaluation of diversity measures for TCR repertoire profiling*,
  BMC Biology (2025) — <https://link.springer.com/article/10.1186/s12915-025-02236-5>
- *A multi-bin rarefying method for evaluating alpha diversities in TCR sequencing data* —
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11246167/>
- Tool catalogue: <https://github.com/slowkow/awesome-vdj>
