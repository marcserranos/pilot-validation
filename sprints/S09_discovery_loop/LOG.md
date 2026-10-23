# S09 log (append-only)
- [10-24] START — Marc asked for an iterative loop (think tank -> MVP -> critics/literature -> think tank) with mixed modes (moonshot / secure / naive / heuristic) and "Nature-cover" ambition. Split from S08 (closed: novel-allele story).
  Round 2 think tank run (9 Haiku agents); critic + literature check launched; HLA-F ASE (S08-B, script 96) still running on the VM at hand-off.
- [10-24] RESULT S09-01 (script 100, VM, 12,208 people; aggregates only) — native vs long-read: 2-field genotype concordance is low for DQA1 (0.52), DQB1 (0.78), DPB1 (0.87), C (0.90; AFR 0.86) but 0.98-0.99 at G-group level for every gene and ancestry;
  DQA1*05:01/05:03/05:05 and DQB1*02:01/02:02 share a G group (identical over exon 2), so the native caller is not wrong, it is exon-limited. Lesson: do NOT claim native errors from 2-field names; always score at G/P-group level. S09-02 (misclassification-corrected ORs) killed by this.
  Remaining value of long reads: variation outside the peptide-binding exons (S09-01b) and haplotype/SV/KIR context.
