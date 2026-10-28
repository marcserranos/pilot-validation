"""116 S09 full report (HTML, inline SVG): overview first, then one section per experiment with evidence, method and caveat drop-downs.

  python3 scripts/hla_popgen/116_s09_full_report.py --out reports/hla_popgen/S09_discovery_loop/S09_report.html [--review reviewer_notes.html]

Every number printed here was read from a run output recorded in sprints/S09_discovery_loop/{LOG,LEDGER}.md or from the aggregate tables on the VM
(counts of 1-19 are never shown). Person-level files stay on the VM. Figures come from the lists defined in this file; each figure carries a 'numbers behind
this figure' table. Items that were proposed but NOT run are labelled as such.
"""
import argparse
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


F = _load("f116", "116a_figs.py")
BG = _load("bg104", "104a_background.py")
fig, ev, badge, cards, table, section = F.fig, F.ev, F.badge, F.cards, F.table, F.section

# ===================================================================== data
GENES = ["HLA-A", "HLA-C", "DPB1", "DQA1", "DQB1", "DRB1"]
TWO_FIELD = [0.968, 0.899, 0.873, 0.522, 0.778, 0.965]
G_LEVEL = [0.987, 0.988, 0.993, 0.990, 0.989, 0.980]

RNA_FOREST = [
    ("#", "Positive controls (must be positive: a gene's own copy number drives its own RNA)"),
    ("KIR2DL1 dosage → KIR2DL1 RNA", 0.1719, 0.0135, "pos"), ("KIR2DL3 dosage → KIR2DL3 RNA", 0.3008, 0.0106, "pos"),
    ("C4A copies → C4A RNA", 0.0291, 0.0036, "pos"), ("C4B copies → C4B RNA", 0.0692, 0.0039, "pos"),
    ("#", "Reversed controls (KIR3DL1 / KIR2DS4 share one locus; reads mis-map between them)"),
    ("KIR3DL1 dosage → KIR3DL1 RNA", -0.1521, 0.0153, "rev"), ("KIR2DS4 dosage → KIR2DS4 RNA", -0.3774, 0.0169, "rev"),
    ("#", "Negative controls (must be zero)"),
    ("KIR2DL1 dosage → CD3E (T-cell gene)", -0.0174, 0.0217, "neg"), ("C4 copies → GAPDH (housekeeping)", -0.0000, 0.0036, "neg"),
    ("#", "The actual questions"),
    ("Activating-KIR count → NK-cell score", -0.0012, 0.0044, "test"), ("Total KIR genes → NK-cell score", -0.0045, 0.0036, "test"),
    ("C4 copies → interferon score", -0.0208, 0.0106, "test"), ("1 vs 2 HLA-A calls → HLA-A RNA", 0.0245, 0.0147, "test"),
]
WITHIN_G = [("HLA-DQA1", 80, 0.58, 0.002), ("HLA-DPA1", 6, 0.35, 0.002), ("HLA-DQB1", 75, 0.17, 0.024), ("HLA-C", 28, 0.04, 0.48)]

RAW_READS = [("Supported in every examined carrier, 0/8 controls", 18, "pos"), ("Uninformative (alt base is in all controls)", 3, "neg"), ("No raw-read support = assembly artefact", 5, "rev")]

HLAE_ANC = [("EUR (n=1,038)", -0.0860, 0.0265, "main"), ("AFR (n=466)", -0.1013, 0.0413, "main"), ("AMR (n=595)", -0.0657, 0.0366, "main"), ("EAS (n=428)", -0.0982, 0.0474, "main"),
            ("SAS (n=448)", -0.1212, 0.0447, "main"), ("Meta-analysis, 5 ancestries (n=2,975)", -0.0906, 0.0164, "test")]
HLAE_SPEC = [("HLA-E RNA itself (positive control)", -0.0659, 0.0068, "pos"), ("KLRC1 = NKG2A (the finding)", -0.0818, 0.0180, "test"), ("KLRD1 = CD94 (NKG2A's partner)", 0.0007, 0.0132, "neg"),
             ("KLRC2 = NKG2C (activating paralogue)", 0.0052, 0.0197, "neg"), ("KLRK1 = NKG2D", -0.0100, 0.0108, "neg"), ("NCR1", 0.0102, 0.0111, "neg"), ("KLRB1", 0.0022, 0.0153, "neg"),
             ("NCAM1 = CD56", 0.0158, 0.0090, "neg"), ("IFNG", -0.0062, 0.0119, "neg")]
HLAE_ROBUST = [("Full model (NK, interferon, CD8/γδ T, KLRC2, leader, age, sex)", -0.0904, 0.0160, "main"), ("+ 10 ancestry PCs, ancestry meta (n=2,975)", -0.0847, 0.0164, "test"),
               ("+ PCs + 3 KLRC1-region cis variants, pooled (n=3,140)", -0.0753, 0.0143, "main"), ("+ PCs + cis variants, ancestry meta (n=2,975)", -0.0729, 0.0147, "test")]
HLAE_PLASMA = [("KLRC1 (NKG2A) — all ancestries pooled", 0.0226, 0.0416, "test"), ("KLRC1 — meta over 5 ancestries", 0.0381, 0.0432, "test"), ("KLRD1 (CD94) — pooled", 0.0368, 0.0178, "neg"),
               ("KLRK1 (NKG2D) — pooled", 0.0329, 0.0142, "neg"), ("KLRK1 — meta", 0.0401, 0.0142, "neg"), ("KLRB1 — pooled", 0.0347, 0.0166, "neg")]

PHEWAS_SETS = [("HLA alleles + ligand groups", 1.123, 1.130, 6.98, "0", "calibrated"), ("C4A / C4B copy number", 1.105, 0.869, 4.07, "1 (C4A → I77)", "calibrated"),
               ("KIR gene content", 1.153, 1.788, 4.29, "6", "NOT calibrated: do not trust")]
PHEWAS_CTRL = [("DQB1*02:01 → type 1 diabetes (E10)", 5.46, "pos"), ("DRB1*03:01 → type 1 diabetes", 3.20, "pos"), ("DQB1*03:02 → type 1 diabetes", 1.75, "neg"), ("DRB1*15:01 → multiple sclerosis (G35)", 1.76, "neg"),
               ("A*31:01 → drug-reaction code (M05)", 1.47, "neg"), ("C*06:02 → psoriasis (L40)", 1.07, "neg")]
PHEWAS_TOP = [("DQA1*05:03", "E87 fluid/electrolyte disorders", 304, 1122, 5.88), ("DRB1*01:02", "L63 alopecia areata", 567, 35, 5.74), ("DQA1*05:03", "K72 hepatic failure", 304, 62, 5.63),
              ("DQA1*05:03", "K70 alcoholic liver disease", 304, 67, 5.56), ("DQB1*02:01", "E10 type 1 diabetes", 1414, 117, 5.46), ("DQA1*05:03", "A41 sepsis", 304, 290, 5.41),
              ("B*45:01", "N87 cervical dysplasia", 330, 57, 5.33), ("DQA1*05:01", "E10 type 1 diabetes", 1495, 117, 5.23), ("B*27:05", "F15 stimulant use", 309, 128, 4.90),
              ("DQA1*03:03", "E10 type 1 diabetes", 1210, 117, 4.75), ("B*27:05", "H20 iridocyclitis", 309, 58, 4.74), ("B*08:01", "O02 abnormal conception products", 918, 45, 4.68)]
ANC5 = ["AFR", "AMR", "EAS", "EUR", "SAS"]
DQA_BEFORE = {"E87 fluid/electrolyte": [0.91, 6.45, -0.77, 0.74, -1.43], "A41 sepsis": [3.33, 5.27, -0.13, -0.93, -0.36]}
DQA_AMR = [("E87 fluid/electrolyte", 6.45, 2.00), ("K72 hepatic failure", 5.41, 1.00), ("K70 alcoholic liver", 4.84, 1.26), ("A41 sepsis", 5.27, 1.82)]
C4A_BEFORE = [0.53, 0.76, 0.65, 4.38, 0.74, 0.23]
C4A_AFTER = [0.37, 0.63, 0.90, 4.32, 1.90, 0.40]
B27_BEFORE = [2.01, 3.23, -0.11, -0.83, 9.90]
B27_AFTER = [2.30, 4.07, -0.20, -0.58, 6.58]

PQTL_CTRL = [("KIR3DL1 gene content → KIR3DL1 protein", 50.99, "pos"), ("KIR2DS4 gene content → KIR3DL1 protein (same haplotype)", 50.96, "pos"), ("KIR2DL2 content → KIR2DL2/3 protein", 29.34, "pos"),
             ("KIR2DL3 content → KIR2DL2/3 protein (opposite sign: the other allele)", -28.52, "pos"), ("C4B copy number → C4B protein", 29.46, "pos"), ("A*01:01 → HLA-A protein (assay epitope effect)", 46.53, "ref"),
             ("A*02:01 → HLA-A protein (assay epitope effect)", -25.02, "ref"), ("E*01:03 (107G) → HLA-E protein (no effect)", 1.82, "neg")]
PQTL_ANC = [("B*14:02 → LRPAP1", [32.8, 39.6, 9.7, 55.5, 10.6]), ("C*08:02 → LRPAP1", [25.4, 39.0, 8.3, 53.0, 9.4]), ("KIR3DL1 content → KIR3DL1", [13.5, 19.2, 26.4, 26.1, 25.9]),
            ("A*01:01 → HLA-A", [17.4, 20.4, 8.5, 29.4, 19.9]), ("KIR2DL5B content → KIR2DL2/3", [5.6, 3.9, 4.8, 5.3, 8.8]), ("B*57:01 → MICA/MICB", [-3.5, -3.1, -0.9, -10.8, -9.0])]
PQTL_TOP = [("B*14:02", "LRPAP1", 78.8, 66.3), ("C*08:02", "LRPAP1", 73.7, 60.4), ("KIR3DL1 dosage", "KIR3DL1", 51.0, 49.7), ("KIR2DS4 dosage", "KIR3DL1", 51.0, 49.7), ("KIR2DS1 dosage", "KIR3DL1", -50.8, -49.4),
            ("A*01:01", "HLA-A", 46.5, 42.7), ("KIR2DL5A dosage", "KIR3DL1", -44.0, -41.9), ("B*14:02", "LRP1", 43.7, 34.7), ("KIR3DS1 dosage", "KIR3DL1", -42.9, -40.7), ("C*08:02", "LRP1", 40.8, 32.0),
            ("KIR2DS5 dosage", "KIR3DL1", -39.7, -38.1), ("B*14:02", "SGSH", 35.7, 29.1), ("C*08:02", "SGSH", 31.8, 25.1), ("G*01:06", "HLA-A", 31.3, 28.0), ("KIR2DS2 dosage", "KIR2DL2/3", 31.0, 29.8),
            ("DRB1*01:02", "LRPAP1", 29.7, 23.6), ("C4B copies", "C4B", 29.5, 28.9), ("KIR2DL2 dosage", "KIR2DL2/3", 29.3, 28.0), ("KIR2DL3 dosage", "KIR2DL2/3", -28.5, -27.2), ("C*07:02", "KIR2DL2/3", 27.5, 27.4),
            ("DQB1*03:02", "HLA-DRA", 26.4, 25.1), ("A*02:01", "HLA-A", -25.0, -22.9), ("DRB1*01:02", "LRP1", 19.0, None), ("B*07:02", "KIR2DL2/3", 18.9, None), ("C*06:02", "CD8A", 19.1, None),
            ("C*06:02", "CCHCR1", -16.3, None), ("B*14:02", "MICA/MICB", -16.9, None), ("B*40:01", "SFTPD", 15.8, None), ("B*27:05", "CLEC4E", 15.1, None), ("B*08:01", "MICA/MICB", 14.8, None)]

B14_SHIFT = [("LRPAP1", 2.59, -0.15), ("LRP1", 1.70, -0.10), ("SGSH", 1.37, -0.08)]
B14_ANC = [("AFR (66 carriers)", 2.74, -0.13), ("AMR (111)", 2.91, -0.18), ("EUR (253)", 2.49, -0.27)]
B14_SUB = [("#", "HLA-B alleles (SD shift in LRPAP1 per copy)"), ("B*14:02 (474 carriers)", 2.72, "warn"), ("B*14:01 (79)", 0.93, "test"), ("B*15:01 (398)", -0.08, "neg"), ("B*15:03 (216)", -0.10, "neg"),
           ("#", "HLA-C alleles"), ("C*08:02 (559) — tags B*14:02", 2.42, "warn"), ("C*08:01 (327)", -0.02, "neg"), ("C*08:03 (31)", -0.13, "neg"), ("C*07:02 (1,548)", -0.09, "neg"), ("C*07:01 (1,179)", -0.14, "neg")]
B14_JOINT = [("LRPAP1", 1.74, 1.20), ("LRP1", 1.15, 0.71), ("SGSH", 1.13, 0.48)]
B14_AA = [("HLA-B residue 110 W", 67.2, 9.5), ("HLA-B residue 129 F", 32.9, 3.3), ("HLA-B residue 72 C", 32.0, 3.2), ("HLA-B residue 126 Y", 22.5, 2.5), ("HLA-B residue 86 D", 20.6, 1.1),
          ("HLA-B residue 12 V", 18.2, 1.1), ("HLA-C residue 174 K", 45.6, 14.5), ("HLA-C residue 64 Q", 38.2, 11.9), ("HLA-C residue 174 T", -25.8, -9.0), ("HLA-C residue 64 R", -24.6, -8.6)]
B14_LIGANDS = [("#", "The three proteins that move"), ("LRPAP1", 78.8, "warn"), ("LRP1", 43.7, "warn"), ("SGSH", 35.7, "warn"),
               ("#", "Other LRP1 receptors / ligands / lysosomal enzymes on the same panel (all |z| < 2.1)"), ("LRP2", 0.3, "neg"), ("LRP11", -1.2, "neg"), ("LDLR", 0.9, "neg"), ("SORT1", -0.1, "neg"), ("SERPINE1", 0.7, "neg"),
               ("PLAU", 0.1, "neg"), ("PLAUR", -0.9, "neg"), ("MMP9", -0.3, "neg"), ("MMP13", -0.9, "neg"), ("THBS2", -0.4, "neg"), ("PCSK9", 0.6, "neg"), ("APP", 0.4, "neg"), ("PSAP", 0.7, "neg"),
               ("IDUA", -0.3, "neg"), ("IDS", -0.2, "neg"), ("CTSD", 0.1, "neg"), ("CTSB", -1.5, "neg"), ("GUSB", -2.0, "neg"), ("LAMP1", 0.7, "neg"), ("TPP1", -0.6, "neg")]
B14_RNA = [("#", "Genes of the three proteins and neighbours (log2 TPM shift in carriers; n = 8,000, 447 carriers)"), ("LRPAP1", 0.009, "neg", "p 0.68"), ("LRP1", 0.055, "neg", "p 0.13"), ("SGSH", 0.018, "neg", "p 0.50"),
           ("LRP2", 0.000, "neg", "p 0.35"), ("LRP10", -0.004, "neg", "p 0.85"), ("LDLR", 0.008, "neg", "p 0.67"), ("SORT1", 0.025, "neg", "p 0.40"), ("A2M", 0.021, "neg", "p 0.53"),
           ("#", "Positive control: the same haplotype does change its own neighbours' RNA"), ("C4A", 0.157, "pos", "p 3e-32"), ("C4B", 0.161, "pos", "p 7e-31"), ("HLA-C", -0.145, "pos", "p 2e-9"), ("HSPA1B", -0.111, "pos", "p 2e-5")]

IDEAS = [
    ("Native vs long-read HLA calls", "secure", "done", "S09-01/02 → section 1", "run"),
    ("Do KIR gene content / C4 copy number move blood RNA?", "naive", "done: null", "S09-10 → section 2", "run"),
    ("Class II expression outside the peptide groove", "heuristic", "weak signal", "S09-11 → section 2", "run"),
    ("Raw-read support for unnamed alleles; HLA-F loss-of-function", "secure", "done", "S09-12/13 → section 3", "run"),
    ("B*27 misfolding → unfolded-protein response in blood", "moonshot", "null", "S09-15 → section 4", "run"),
    ("HLA-E R107G → NKG2A per NK cell", "moonshot", "lead", "S09-16 → section 4", "run"),
    ("HLA / KIR / C4 × ICD-10 disease scan", "heuristic", "no discovery", "S09-14 → section 5", "run"),
    ("HLA / KIR / C4 → plasma proteome", "heuristic", "validation + atlas", "S09-17 → section 6", "run"),
    ("B*14:02 → LRPAP1 / LRP1 / SGSH", "serendipity", "open", "S09-18 → section 7", "run"),
    ("HLA-F T237I/P205T × KIR3DS1 vs viral diagnoses", "moonshot", "proposed", "carriers <50: needs power check first", "notrun"),
    ("HLA-A deletion shifts KIR content", "moonshot", "proposed", "carriers likely <20", "notrun"),
    ("Parent-of-origin effect on HLA expression", "moonshot", "proposed", "needs a trio × RNA overlap check", "notrun"),
    ("Balancing-selection halos around HLA genes", "moonshot", "proposed", "not started", "notrun"),
    ("MICA intron-5 repeat as an expression switch", "moonshot", "proposed", "not started", "notrun"),
    ("Phased MHC haplotype atlas", "secure", "proposed", "not started", "notrun"),
    ("Trio Mendelian-error rate per gene", "secure", "proposed", "partly covered by S06 trio tables", "notrun"),
    ("KIR atlas with a measured short-read false-negative rate", "secure", "proposed", "not started", "notrun"),
    ("Ten quick 'naive looks' (template-distance histograms, KIR gene-count outliers, DRB3/4/5 vs DRB1 consistency, genes-called by platform, haplotype block lengths)", "naive", "proposed", "NOT RUN: script 115_naive_looks.py is written and tested on synthetic data; the VM was unreachable when this report was rebuilt", "notrun"),
    ("Residue-level association scan across all HLA positions", "heuristic", "partly run", "done only for B*14:02 → LRPAP1 (section 7)", "notrun"),
    ("Haplotype-window scan; anomaly detection on per-person divergence", "heuristic", "proposed", "not started", "notrun"),
    ("Drug-hypersensitivity positive controls (abacavir/B*57:01 …), blood counts and labs, COVID severity from epitope load", "power-aware", "proposed", "needs the native 535k cohort", "notrun"),
]

LEDGER = [
    ("S09-01", "Native short-read HLA calls are wrong where long reads are right", "Killed", "2-field concordance 0.52–0.97, but 0.98–0.99 at G-group level: a resolution effect, not error", 1),
    ("S09-02", "Native errors distort HLA–disease odds ratios", "Killed", "follows from S09-01", 1),
    ("S09-10a", "KIR gene content shifts NK-cell transcripts", "Null", "n=5,120; p=0.79 / 0.22; negative control null", 2),
    ("S09-10b", "C4 copy number tunes interferon tone", "Null", "nominal p=0.049 only; controls behave", 2),
    ("S09-11", "Variation outside the peptide groove changes class II expression", "Weak", "0.2–0.6% of variance (permutation p 0.002–0.024); mapping bias not excluded", 2),
    ("S09-12", "Unnamed alleles are supported by raw long reads", "Strong", "18/26 supported in every carrier, 0/8 controls; 5 HLA-F loss-of-function alleles real; HLA-E/G stop-gains are artefacts", 3),
    ("S09-13", "HLA-E/F/G diversity landscape", "Done", "E balanced; F nearly monomorphic; G null haplotypes 4–5× more frequent in AFR", 3),
    ("S09-13b", "HLA-F loss-of-function lowers its own transcript", "Weak", "R294* haplotype 26% of reads; one carrier each", 3),
    ("S09-14", "HLA/KIR/C4 × ICD-10 phenome scan", "Lead", "calibrated for HLA/C4; no HLA family-wise hit; C4A → I77 survives PCs but EUR-driven; KIR not calibrated", 5),
    ("S09-15", "B*27 → unfolded-protein response in blood", "Null", "β −0.016 ± 0.024, n=5,642; controls flat", 4),
    ("S09-16", "HLA-E 107G lowers NKG2A per NK cell", "Strong", "meta p 3.6e-8; 5/5 ancestries; specific to KLRC1", 4),
    ("S09-16c", "…survives ancestry structure (10 PCs)", "Strong", "pooled −0.090 (p 1.7e-8); meta −0.085 (p 2.6e-7)", 4),
    ("S09-16d", "…and KLRC1 cis-eQTLs", "Strong", "3 cis leads: −0.075 (p 1.5e-7), meta −0.073 (p 7e-7)", 4),
    ("S09-16e", "…and plasma NKG2A (Olink)", "Open", "n=2,709: +0.023 NPX, p 0.59; plasma = soluble fraction, inconclusive", 4),
    ("S09-17", "HLA/KIR/C4 → plasma proteome atlas", "Strong", "5,415 assays, n=8,884, λ 1.10; KIR3DL1 z 51; C4B z 29; top-60 hits keep their sign in ≥4/5 ancestries", 6),
    ("S09-18", "B*14:02 / C*08:02 → LRPAP1, LRP1, SGSH", "Open", "+2.6/+1.7/+1.4 SD; allele-specific, same size in AFR/AMR/EUR, blood RNA null; protein link known from UK Biobank, allele unnamed there", 7),
]


# ===================================================================== helpers
def dt(headers, rows):
    return table(headers, rows)


def sec_link(n):
    return '<a href="#s%d">section %d</a>' % (n, n)


def sd(v):
    return "%+.2f" % v


# ===================================================================== sections
def overview(review_html):
    imp = """<div class="box imp"><h3>Most important or most novel findings</h3><ol>
<li><b>HLA-E R107G lowers the NKG2A (KLRC1) RNA in blood, about 6% per copy.</b> {} {}
<span class="why">The best new mechanistic lead. It is specific (8 neighbouring NK genes do not move), appears in all five ancestries, and survives ancestry and local-genetics controls (p ≈ 7e-7). The plasma protein test is null, but plasma only holds the soluble fraction, so the protein-level question is open. See {}.</span></li>
<li><b>One HLA allele, HLA-B*14:02, moves three unrelated plasma proteins by 1.4–2.6 SD.</b> {} {}
<span class="why">The biggest effect we saw (456 carriers among 8,884). It is allele-specific, the same size in AFR, AMR and EUR, and absent from blood RNA. UK Biobank already reported that HLA explains these proteins but named no allele; the cause (real biology or assay interference) is unresolved. See {}.</span></li>
<li><b>The long-read KIR and C4 calls are independently validated by plasma protein levels.</b> {}
<span class="why">KIR3DL1 gene content predicts the KIR3DL1 protein with z = 51 and C4B copy number predicts C4B protein (z = 29), in all five ancestries. This is a method-validation result, useful for the paper. See {}.</span></li>
<li><b>Native short-read HLA calls are accurate at the resolution that matters.</b> {}
<span class="why">DQA1 looks 52% concordant by allele name but 99% by peptide-binding groove (G group). This closed off a tempting idea ("native errors distort disease studies"). See {}.</span></li>
<li><b>Five real HLA-F loss-of-function alleles; the HLA-E and HLA-G ones were artefacts.</b> {}
<span class="why">Raw long reads support 18 of 26 recurrent unnamed variants in every carrier. See {}, and the S08 report for how the unnamed HLA-F alleles were found.</span></li></ol></div>""".format(
        badge("lead"), badge("open", "Protein not yet tested"), sec_link(4), badge("open"), badge("known", "Partly known"), sec_link(7), badge("known", "Validation"), sec_link(6),
        badge("known", "Validation"), sec_link(1), badge("strong"), sec_link(3))
    prom = """<div class="box prom"><h3>Most promising lines of research (in order)</h3><ol>
<li><b>Test HLA-E R107G → NKG2A at the protein level.</b> Flow cytometry for surface NKG2A on NK cells from donors chosen by HLA-E genotype. It is the only decisive test left and the effect is specific enough to justify it. <span class="why">Needs a wet-lab partner or a recontact cohort; cannot be done in the Workbench.</span></li>
<li><b>Settle HLA-B*14:02 → LRPAP1 / LRP1 / SGSH.</b> Measure LRPAP1 by an independent method (ELISA or mass spectrometry) in carriers, check the UK Biobank supplement and liver eQTLs. <span class="why">If real it explains the "unclear and intriguing" HLA–LRPAP1 link that the UK Biobank authors left open; if an assay artefact it is a useful warning for every Olink HLA study.</span></li>
<li><b>Write up the HLA/KIR/C4 → proteome atlas.</b> Five ancestries, long-read-resolved alleles, KIR content and C4 copy number on a 5,415-protein panel. <span class="why">Only the top 60 of 780,000 associations have been looked at; mining the rest (for example which KIR haplotypes move soluble immune proteins) is cheap.</span></li>
<li><b>Prepare the IPD submission for HLA-F T237I and P205T</b> (see the S08 report): needs raw-read confirmation in at least 3 unrelated carriers and an independent trio.</li>
<li><b>Replicate C4A copy number → arterial disorders (I77).</b> Weak but survived every check we ran. <span class="why">Needs C4 copy number in a second cohort.</span></li></ol></div>"""
    under = """<div class="box under"><h3>Underwhelming or negative (still informative)</h3><ul>
<li><b>No disease-association discovery.</b> The scan of HLA, KIR and C4 against ~600 diagnosis categories produced no family-wise HLA hit. Its strongest signal (DQA1*05:03 → liver/sepsis codes) was ancestry confounding inside one group and is retired. {}</li>
<li><b>KIR gene content does not change NK-cell RNA</b> (n=5,120) and <b>C4 copy number does not tune interferon</b>: clean nulls with working controls. {}</li>
<li><b>B*27 does not trigger an unfolded-protein-response signature in blood</b> (n=5,642). {}</li>
<li><b>Class II expression differences inside identical-groove alleles are tiny</b> (0.2–0.6% of variance) and mapping bias is not excluded. {}</li>
<li><b>Plasma NKG2A is null</b> (n=2,709), but plasma holds only the soluble fraction, so it is inconclusive rather than a refutation.</li>
<li><b>The ten "naive looks" were proposed but never run</b> (the VM was unreachable when this report was rebuilt). They are listed honestly in {} with a plan.</li></ul></div>""".format(
        badge("retired"), badge("null"), badge("null"), badge("lead", "Weak signal"), sec_link(8))
    cheat = """<details class="bg"><summary><span class="k">Cheat-sheet</span> How to read the numbers in this report (30 seconds)</summary><div class="in"><dl class="cheat">
<dt>z-score</dt><dd>How many "standard errors" an effect is from zero. About |z| &gt; 2 is conventionally "interesting" for one test; with thousands of tests we need |z| ≳ 6 (see the threshold printed in each section). |z| of 50 is absurdly strong.</dd>
<dt>p-value</dt><dd>Chance of a result at least this extreme if there were no effect. 3e-8 means 3 in 100 million.</dd>
<dt>SD</dt><dd>Standard deviation, the typical person-to-person spread. "+2.6 SD" means carriers sit about two and a half spreads above everyone else.</dd>
<dt>log2 TPM</dt><dd>RNA amount on a log2 scale. −0.09 ≈ 6% less RNA; −1 = half.</dd>
<dt>NPX</dt><dd>The Olink plasma-protein unit, also log2: +1 NPX = twice the signal.</dd>
<dt>Control</dt><dd>A test whose answer we know (positive = must be non-zero, negative = must be zero). A result is only believed if the controls behave.</dd>
<dt>Ancestry / PCs</dt><dd>People are grouped by genetic ancestry (AFR, AMR, EAS, EUR, SAS, MID). "PCs" are numbers summarising finer ancestry within a group; adding them removes ancestry confounding.</dd>
<dt>Replicates in 5/5 ancestries</dt><dd>The effect has the same sign in each ancestry separately, so it is unlikely to be one group's quirk.</dd></dl></div></details>"""
    return """<h2 id="glance">At a glance</h2>
<p class="sub">What this sprint asked: do HLA, KIR and complement C4 variants resolved by long-read sequencing explain anything in blood RNA, health records or plasma proteins that short-read data could not? Cohort: up to 12,233 All of Us long-read genomes, blood RNA-seq (3,141–8,000 people depending on the test), Olink plasma proteomics (up to 8,884) and health records.</p>
<div class="glance three">%s%s%s</div>%s%s""" % (imp, prom, under, cheat, review_html)


def s1():
    f1 = fig(F.grouped(GENES, [("2-field allele name", TWO_FIELD, "b"), ("G group (peptide-binding exons)", G_LEVEL, "c")], 1.1, "Genotype concordance, native short-read caller vs long-read truth (12,208 people with both)", height=300, fmt="{:.2f}"),
             "Left bar: do the two methods agree on the allele name? Right bar: do they agree on the peptide-binding groove? DQA1, DQB1 and DPB1 look poor by name and excellent by groove.",
             "script 100 → ~/s09/native_vs_lr_pub.tsv (VM); ledger S09-01",
             dt(["Gene", "2-field name", "G group"], [(g, "%.3f" % a, "%.3f" % b) for g, a, b in zip(GENES, TWO_FIELD, G_LEVEL)]))
    return section("s1", "1. Are the native short-read HLA calls wrong?", badge("known", "Idea killed"),
                   "<b>No, not where it matters.</b> The All of Us native calls agree with long-read truth for 98–99% of genotypes when judged at the level of the peptide-binding groove. The apparent errors at the allele-name level (DQA1 only 52%) are a resolution effect: short reads cannot tell apart alleles that differ outside the groove. The hypothesis that native errors distort HLA–disease studies was dropped.",
                   cards([("0.52 → 0.99", "DQA1 agreement, allele name vs groove (G group)"), ("12,208", "people with both native and long-read calls"), ("6 genes", "HLA-A, C, DPB1, DQA1, DQB1, DRB1")]),
                   f1 + ev("What 'resolution' and 'G group' mean", "<p>An allele name has up to four colon-separated fields; only the first two change the protein. A <b>G group</b> collects alleles that are identical over exons 2 and 3, the peptide-binding groove that T cells and NK cells read. A caller that sees only those exons can be right about the groove and unable to name the allele.</p><figure>%s</figure>" % BG.fig_resolution(), "bg") +
                   ev("Scripts, models and thresholds", "<p>For each person and gene, the genotype from the AoU native short-read callers was compared with the long-read assembly call (Immuannot, template distance 0 = ground-truth-like). Names were compared as 2-field strings; groove-level comparison maps each allele to its G group using IPD-IMGT/HLA <code>hla_nom_g.txt</code>.</p>", "how") +
                   ev("What could still be wrong", "<ul><li>G-group concordance is an <i>upper bound</i>: it ignores differences outside exons 2 and 3.</li><li>Long-read calls are the reference here but are not perfect; the 'template distance 0' subset is the most trustworthy.</li><li>This closes the 'native errors' idea for peptide-binding questions; it does not say native calls capture structural variation or phase, which long reads add.</li></ul>", "cav"))


def s2():
    f2 = fig(F.forest(RNA_FOREST, (-0.5, 0.5), "Effect of long-read-only gene features on blood RNA (95% CI)", label_w=380, xlabel="change in log2 TPM per unit of the feature"),
             "Green: a gene's own copy number drives its own RNA, so the pipeline sees real biology. Red: KIR3DL1 and KIR2DS4 are one locus with two names and reads mis-map between them, so KIR gene-level RNA cannot validate every gene. Grey: unrelated genes do not move. Orange: the three real questions are all flat.",
             "script 102 → ~/s09/expression_tests_pub.tsv (VM); ledger S09-10",
             dt(["Test", "beta", "se"], [(r[0], "%+.4f" % r[1], "%.4f" % r[2]) for r in RNA_FOREST if r[0] != "#"]))
    f3 = fig(F.hbars([("%s" % g, e, "b", "p = %g (groove alone: %d%%)" % (p, gg)) for g, gg, e, p in WITHIN_G], vmax=0.9, title="Expression variance explained beyond the peptide-binding groove (% of total)", label_w=150, fmt="{:.2f}%", right=330),
             "Among people homozygous for alleles that share a groove, which 3-field allele they carry explains a tiny extra share of RNA level. Real for HLA-DQA1 and DPA1, borderline for DQB1, absent for HLA-C. HLA-A, B, DPB1 and DRB1 could not be tested (no groove group with enough homozygotes).",
             "script 103 → ~/s09/within_g_expression_pub.tsv (VM); ledger S09-11",
             dt(["Gene", "G-group identity explains (%)", "extra explained within G (%)", "permutation p"], [(g, gg, e, p) for g, gg, e, p in WITHIN_G]))
    return section("s2", "2. Do KIR content, C4 copies and allele families change blood RNA?", badge("null", "Mostly null"),
                   "<b>The controls work and the questions come back empty.</b> Having a KIR gene, or more C4A/C4B copies, raises that gene's own RNA as it should. But KIR gene content does not shift NK-cell RNA (n=5,120) and C4 copy number does not shift interferon-response RNA. A small extra effect of non-groove variation on class II expression exists (0.2–0.6% of variance) but mapping bias is not excluded.",
                   cards([("n = 5,120", "people for the KIR → NK-cell test (p = 0.79 and 0.22)"), ("p = 0.049", "C4 → interferon, nominal only; both controls behave"), ("0.2–0.6%", "extra expression variance explained within a groove group")]),
                   f2 + ev("Why a positive control matters here", "<p>Short-read RNA maps poorly across KIR genes, which are near-identical. Before believing 'KIR content does not change NK RNA' we first check that KIR2DL1 and KIR2DL3 copy number <i>does</i> change KIR2DL1 and KIR2DL3 RNA. It does (z ≈ 13–28). KIR3DL1 and KIR2DS4 come out reversed (a mapping artefact between loci), so gene-level KIR RNA is only trusted for some genes.</p>", "bg") + f3 +
                   ev("Scripts, models and thresholds", "<p>Script 101 streams the RNA-SeQC2 TPM file for the blood RNA-seq cohort and keeps a panel of KIR, HLA, complement, lineage-marker and interferon genes. Script 102 regresses log2(TPM+1) on the long-read feature with ancestry dummies (groups under 20 dropped); NK and interferon scores are the mean of z-scored marker genes. Script 103 compares nested models (G group vs G group + 3-field allele) among homozygotes and builds a permutation null.</p>", "how") +
                   ev("What could still be wrong", "<ul><li>RNA is aligned to GRCh38, so haplotypes far from the reference can look lower for mapping reasons. Next step: re-quantify on personal two-haplotype references (S05/S08 machinery).</li><li>Bulk blood RNA mixes cell types; NK scores are composition proxies.</li><li>Nulls here exclude effects larger than roughly the confidence intervals shown, not tiny ones.</li></ul>", "cav"))


def s3():
    f_raw = fig(F.hbars([(a, b, c) for a, b, c in RAW_READS], vmax=24, title="Raw-read check of 26 recurrent unnamed signatures", label_w=380, fmt="{:.0f} signatures"),
                "Each signature was located exactly on GRCh38 from the GENCODE gene model (checked against assembly bases) and piled up in raw long reads of up to 8 carriers and 8 matched non-carriers. Real inherited variants appear in every carrier read and no control read.",
                "script 105 → ~/s09/e105 (VM, person-level; aggregate in LOG S09-12)",
                dt(["Outcome", "signatures"], [(a, b) for a, b, c in RAW_READS]))
    f_lof = fig(F.hbars([("R294* stop-gain (HLA-F)", 26, "rev", "expected ≈50% if normal"), ("M1V start-loss (HLA-F)", 40, "test", "transcript persists, protein lost"), ("missense variants, lowest carrier", 35, "neg"), ("missense variants, highest carrier", 50, "neg")], vmax=60,
                        title="Share of HLA-F RNA reads coming from the variant haplotype (single carriers)", label_w=380, fmt="{:.0f}%"),
                 "In a heterozygote half the RNA should come from each copy. The stop-gain copy contributes only ~26%, the signature of nonsense-mediated decay. One carrier per allele, so this is suggestive, not proven.",
                 "script 96 on personal two-haplotype references (VM); LOG S09-13b",
                 dt(["Variant", "variant-haplotype read fraction"], [("R294* stop-gain", "26%"), ("M1V start-loss", "40%"), ("missense variants", "35–50%")]))
    f_land = fig(F.hbars([("HLA-E expected heterozygosity (%), all 6 ancestries", 51, "c", "range 50–52%"), ("HLA-F top allele frequency, EAS", 96, "a"), ("HLA-F top allele frequency, AFR", 80, "a"),
                          ("HLA-G null/questionable haplotypes, AFR", 6.7, "warn"), ("HLA-G null/questionable haplotypes, SAS and EUR", 1.5, "mut", "range 1.4–1.6%")], vmax=110, title="HLA-E, F and G diversity at equal sample size (400 unrelated people per ancestry)", label_w=380, fmt="{:.1f}%"),
                "HLA-E is balanced and the same in every ancestry (its stable R107G polymorphism, studied in section 4). HLA-F is nearly monomorphic and ancestry-structured. HLA-G loses function 4–5 times more often in AFR than in SAS or EUR.",
                "script 106 → ~/s09/landscape (VM); LOG S09-13",
                dt(["Measure", "value"], [("HLA-E expected heterozygosity", "0.50–0.52 in all six ancestries"), ("HLA-F top allele", "96% EAS … 80% AFR"), ("HLA-G null/questionable haplotypes", "AFR 6.7%; SAS and EUR 1.4–1.6%")]))
    return section("s3", "3. HLA-E, F and G: which unnamed variants are real?", badge("strong"),
                   "<b>HLA-F, not HLA-E or HLA-G, is where real protein-truncating variation sits.</b> Raw long reads support 18 of 26 recurrent unnamed signatures in every carrier and in none of 8 controls, including five HLA-F loss-of-function alleles (stop-gains E70*, Q75*, W265*, R294*; start-loss M1V; each in fewer than 20 unrelated people). The HLA-E and HLA-G stop-gain calls and a DPA1 call had no read support: they were assembly artefacts.",
                   cards([("18 / 26", "unnamed signatures supported by raw reads in every carrier"), ("0 / 8", "matched control people with the variant"), ("26%", "of HLA-F RNA reads from the stop-gain copy (R294*)"), ("4–5×", "more HLA-G null haplotypes in AFR than SAS/EUR")]),
                   f_raw + ev("Background: raw reads, assemblies and nonsense-mediated decay", "<p>Long reads are stitched into an assembly, and alleles are called from the assembly. If the assembly is wrong at one base, a 'novel allele' appears that no read supports, so we return to the raw reads. A premature stop signal makes the cell destroy the message (nonsense-mediated decay); the carrier then shows fewer RNA reads from the mutated copy.</p><figure>%s</figure>" % BG.fig_reads(), "bg") +
                   f_lof + f_land +
                   ev("Scripts, models and thresholds", "<p>Script 105: GENCODE v46 CDS → GRCh38 coordinates (strand-aware, verified against assembly bases), samtools pileup on raw AoU long-read BAMs of carriers and matched controls. Script 96: RNA reads re-aligned to a personal two-haplotype transcript reference built from the carrier's own assembly, with a simulated-read bias control. Script 106: allele frequencies and heterozygosity from 400 unrelated people per ancestry.</p>", "how") +
                   ev("What could still be wrong", "<ul><li>Three 'supported' signatures are uninformative because the alternate base is present in every control.</li><li>'Supported' means real DNA variation, not 'new': names are checked against IPD (see the S08 report).</li><li>Loss-of-function frequencies per allele are below 20 people, so no counts are shown.</li></ul>", "cav"))


def s4():
    f_anc = fig(F.forest([("#", "Per ancestry")] + HLAE_ANC, (-0.25, 0.1), "HLA-E 107G dose → KLRC1 (NKG2A) RNA, fully adjusted (95% CI)", label_w=330, xlabel="log2 TPM per 107G allele (−0.09 ≈ 6% less)"),
                "Each row is a separate group of people. All five estimates are negative and overlap; the meta-analysis combines them.",
                "script 108 (S4) → ~/s09/mechanistic_tests_pub.tsv (VM); LOG S09-16",
                dt(["Group", "beta", "se"], [(a, "%+.4f" % b, "%.4f" % c) for a, b, c, d in HLAE_ANC]))
    f_spec = fig(F.forest(HLAE_SPEC, (-0.25, 0.1), "Specificity: which NK genes move with 107G? (95% CI)", label_w=330, xlabel="log2 TPM per 107G allele"),
                 "HLA-E's own RNA falls (positive control). Of the other NK genes only KLRC1 moves; its partner KLRD1, the activating paralogue KLRC2 and five other NK markers do not. A cell-count effect would have moved all of them.",
                 "script 108 (S3) ; LOG S09-16",
                 dt(["Gene", "beta", "se"], [(a, "%+.4f" % b, "%.4f" % c) for a, b, c, d in HLAE_SPEC]))
    f_rob = fig(F.forest(HLAE_ROBUST, (-0.2, 0.05), "Does the effect survive ancestry and local-genetics controls? (95% CI)", label_w=420, xlabel="log2 TPM per 107G allele"),
                "Adding ten ancestry PCs shrinks the estimate by about 7%. Adding the three strongest KLRC1-region variants shrinks it by about 17% in total. It stays at p ≈ 1e-7.",
                "script 108 (S5, S6); LOG S09-16c/d",
                dt(["Model", "beta", "se"], [(a, "%+.4f" % b, "%.4f" % c) for a, b, c, d in HLAE_ROBUST]))
    f_pl = fig(F.forest(HLAE_PLASMA, (-0.15, 0.2), "Plasma protein test (Olink): does 107G change soluble NK proteins? (95% CI)", label_w=330, xlabel="NPX (log2 protein units) per 107G allele"),
               "Plasma NKG2A (KLRC1) is flat, and the other soluble NK proteins move slightly the opposite way (nominal, not corrected for ten assays). Plasma holds only the shed, soluble protein, so this neither confirms nor refutes the RNA result.",
               "script 111 → ~/s09/olink/olink_tests_pub.tsv (VM, n = 2,709); LOG S09-16e",
               dt(["Protein / analysis", "beta (NPX)", "se"], [(a, "%+.4f" % b, "%.4f" % c) for a, b, c, d in HLAE_PLASMA]))
    return section("s4", "4. HLA-E R107G and NKG2A: the best mechanistic lead", badge("strong", "Strong at RNA") + badge("open", "Protein untested"),
                   "<b>People who carry the stable HLA-E variant 107G make about 6% less NKG2A (KLRC1) RNA per copy.</b> The effect is specific to KLRC1, has the same sign in all five ancestries, is independent of NK-cell abundance and interferon tone, and survives ancestry (ten PCs) and local-genetics (three KLRC1-region variants) controls. The only protein-level test available (plasma) is null, but plasma is the soluble fraction. Surface NKG2A on NK cells has not been measured.",
                   cards([("−0.09", "log2 TPM per 107G copy (≈ 6% less KLRC1 RNA)"), ("3.6e-8 → 7e-7", "meta p-value, before and after the extra controls"), ("1 of 9", "NK genes that move (KLRC1)"), ("n = 3,141", "people with blood RNA, HLA-E and ancestry data")]),
                   ev("Background: what HLA-E, NKG2A and 107G are", "<p>NK cells kill stressed cells but spare healthy ones. They check HLA class I on the target through HLA-E, which displays a fragment of other class I molecules' leader sequence; the NK receptor CD94/NKG2A (genes KLRD1/KLRC1) reads it and sends a brake signal. Everyone carries R107 (E*01:01) or G107 (E*01:03, more stable at the cell surface). We asked whether the version you carry changes how much NKG2A your NK cells make.</p><figure>%s</figure>" % BG.fig_nkg2a(), "bg", True) + f_anc + f_spec +
                   ev("Background: why ancestry PCs and a cis control", "<p>A skeptic can say: people with 107G may also tend to carry a KLRC1-region variant that lowers KLRC1. Two different chromosomes are only correlated through ancestry, so we test it directly: add ten ancestry PCs, then the strongest local variants (extracted from short-read genotypes; the lead variant alone has p = 3.7e-143, so the control is powerful). The measured correlation between 107G and those variants is ≤ 0.05.</p><figure>%s</figure>" % BG.fig_cis(), "bg") + f_rob + f_pl +
                   ev("A failed mechanism for contrast: HLA-B*27 and the unfolded-protein response", "<p>The 1999 misfolding hypothesis says B*27 misfolds and triggers an ER stress response that fuels spondyloarthritis. In blood RNA of 5,642 people, B*27 dose does not change a UPR module score (β = −0.016 ± 0.024, p = 0.49); ERAP1/2, IL23A, BiP and XBP1 are individually null, and B*07, B*57, housekeeping and T-cell controls are flat. A well-powered human null in whole blood (not in the affected joint tissue).</p>", "ev") +
                   ev("More detail on the HLA-E effect", "<ul><li><b>Dose pattern:</b> heterozygotes −0.144, homozygotes −0.159 (near-dominant).</li><li><b>Not mediated by HLA-E RNA</b> (adding HLA-E TPM as a covariate leaves the estimate at −0.091).</li><li><b>Leader-peptide hypothesis</b> (HLA class I signal-peptide methionine at −21 supplying HLA-E): only nominal (meta p ≈ 0.02), attenuated by NK abundance.</li><li><b>Closest prior work:</b> a study of about 100–200 donors found no HLA-E allele effect on NKG2A frequency. We found no report of the RNA effect; the literature search was not exhaustive.</li></ul>", "ev") +
                   ev("Scripts, models and thresholds", "<p>Script 108 regresses log2(TPM+1) of KLRC1 and control genes on E*01:03 dose with age, age², sex, ancestry dummies and a composition block (T/B/myeloid markers); the fully adjusted model adds NK abundance, interferon score, CD8/γδ T content, KLRC2 RNA and the HLA leader-peptide load. Per-ancestry results are combined by inverse-variance meta-analysis. S5 adds ten ancestry PCs (All of Us ancestry_preds). S6 adds the strongest independent KLRC1-region variants (1,756 variants from chr12:10.2–10.7 Mb in plink2). Script 111 repeats the test on Olink plasma proteins.</p>", "how") +
                   ev("What could still be wrong, and what would change our mind", "<ul><li><b>Not excluded:</b> CMV serostatus (adaptive NK cells; KLRC2 RNA is only a partial proxy), other chr6 variants in linkage with E*01:03, bulk-RNA cell-type attribution.</li><li><b>Protein level unknown.</b> A flow-cytometry result showing no difference in surface NKG2A per NK cell would make this a transcript-only curiosity.</li><li>The effect is small (6%). It matters mainly as a specific, replicated, mechanistic signal.</li></ul>", "cav"))


def s5():
    f_cal = dt(["Test set", "inflation λ (all)", "λ (negative-control chapters)", "family-wise |z| threshold", "hits past threshold", "verdict"], [(a, "%.3f" % b, "%.3f" % c, "%.2f" % d, e, f) for a, b, c, d, e, f in PHEWAS_SETS])
    f_ctrl = fig(F.hbars([(a, b, c) for a, b, c in PHEWAS_CTRL], vmax=8, title="Positive controls: known associations the scan should rediscover (z-score)", label_w=360, fmt="{:.2f}"),
                 "Green: detected. Grey: the scan has too few cases to detect them (62–167 cases), so z of 1–3 is expected. The scan finds the main type 1 diabetes association (DQB1*02:01, z 5.5).",
                 "script 107 → ~/s09/phewas/phewas_controls_hla.txt (VM)",
                 dt(["Control", "z"], [(a, "%.2f" % b) for a, b, c in PHEWAS_CTRL]))
    f_top = dt(["HLA allele", "Diagnosis category", "carriers", "cases", "z"], [(a, b, c, d, "%.2f" % e) for a, b, c, d, e in PHEWAS_TOP])
    f_dq = fig(F.dumbbell([(a, b, c) for a, b, c in DQA_AMR], (-2, 8), "DQA1*05:03 hits inside the AMR group, before and after adding ancestry PCs (z-score)", ("before PCs", "with 10 PCs"), ("warn", "c"), label_w=300, vlines=[(6.98, "family-wise threshold 6.98")], xlabel="z-score within AMR (234 carriers)"),
               "All four 'liver and sepsis' hits live only in AMR and collapse when finer ancestry is added: DQA1*05:03 marks admixture in that group, not an HLA effect.",
               "script 109 → ~/s09/phewas/followup_pcs_pub.txt (VM); LOG S09-14b/c",
               dt(["Outcome in AMR", "z before PCs", "z with 10 PCs"], [(a, "%.2f" % b, "%.2f" % c) for a, b, c in DQA_AMR]))
    f_c4 = fig(F.heat(["before adding PCs", "with 10 ancestry PCs"], ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"], [C4A_BEFORE, C4A_AFTER], vmax=5, title="C4A copy number → arterial disorders (I77): z-score by ancestry", label_w=190, cell_w=110,
                      note="Same sign in every ancestry; only EUR reaches significance (n of cases: EUR 114, others fewer)."),
               "The signal sits in EUR and keeps its size after PCs (4.38 → 4.32). Other ancestries are small and positive. Combined with PCs the meta z is 4.14 against a family-wise threshold of 4.07: borderline.",
               "script 109; LOG S09-14c",
               dt(["Ancestry", "z before PCs", "z with PCs"], [(a, "%.2f" % b, "%.2f" % c) for a, b, c in zip(["AFR", "AMR", "EAS", "EUR", "MID", "SAS"], C4A_BEFORE, C4A_AFTER)]))
    f_b27 = fig(F.heat(["before adding PCs", "with 10 ancestry PCs"], ANC5, [B27_BEFORE, B27_AFTER], vmax=8, title="B*27:05 → iridocyclitis (uveitis, H20): z-score by ancestry", label_w=190, cell_w=110,
                       note="Carriers: AFR 40, AMR 108, EAS <20, EUR 125 (no carrier cases), SAS 21 (z driven by <20 carrier cases)."),
                "The classic B*27 uveitis association looked rediscovered (pooled z 4.7). Inside ancestries it is inconsistent: nothing in EUR (zero carrier cases among 125 carriers) and the biggest z comes from a 21-carrier stratum. Inconclusive; not counted as a rediscovery.",
                "script 109; LOG S09-14c",
                dt(["Ancestry", "z before PCs", "z with PCs"], [(a, "%.2f" % b, "%.2f" % c) for a, b, c in zip(ANC5, B27_BEFORE, B27_AFTER)]))
    return section("s5", "5. Hypothesis-free scan of HLA, KIR and C4 against diagnoses", badge("lead", "Leads only") + badge("retired", "Top hit retired"),
                   "<b>The scan found no disease discovery.</b> It recovers the main type 1 diabetes association and its HLA-calibration is acceptable, but no HLA allele reaches family-wise significance and its strongest signal (DQA1*05:03 → liver/sepsis codes) was ancestry confounding within AMR. Only C4A copy number → arterial disorders (I77) survived every check, and it rests on one ancestry and one diagnosis code.",
                   cards([("11,856", "unrelated long-read people with health records"), ("594", "diagnosis categories with ≥30 cases"), ("λ = 1.12", "HLA inflation (negative-control chapters 1.13): calibrated"), ("0", "HLA hits past the family-wise threshold (|z| 6.98)")]),
                   ev("Background: PheWAS, calibration and why strong hits are often false", "<p>A phenome-wide scan tests every allele against hundreds of diagnoses. With ~100,000 tests some hits are chance, so we use <b>negative-control outcomes</b> (injury, external-cause and symptom chapters that should not depend on HLA), <b>label permutations</b> for a family-wise threshold, and the <b>inflation factor λ</b> (1.0 is ideal) as a calibration check. The biggest trap is <b>confounding</b>: within a broad group like AMR people differ in ancestral mix, which affects both which alleles they carry and which diagnoses they get; ancestry PCs are the standard fix.</p><figure>%s</figure>" % BG.fig_confound(), "bg") +
                   "<h3>Calibration</h3>" + ev("Calibration table: is the scan trustworthy?", f_cal + "<p>HLA and C4 are calibrated. KIR is not (negative-control λ 1.79), so the six KIR hits past its threshold are treated as ancestry-linked gene-content structure.</p>", "ev", True) +
                   f_ctrl + "<h3>The strongest hits and what happened to them</h3>" + ev("Top 12 HLA hits (none passes the family-wise threshold)", f_top + "<p>Odds ratios in sparse cells (tens to hundreds) are score-test artefacts and are not shown.</p>", "ev") + f_dq + f_c4 + f_b27 +
                   ev("Scripts, models and thresholds", "<p>Script 107: logistic score tests (ancestry, sex, age and age² as covariates) for each allele × 3-character ICD-10 category in unrelated long-read people; BH false-discovery control; family-wise threshold from 20 label permutations within ancestry×sex strata. The first version was discarded: negative-control chapters showed 269 false hits because near-constant common alleles broke the score-test arithmetic. Fixing it (centred exposures, minimum variance, conditioned permutations, regression test) restored calibration. <b>Lesson: a miscalibrated negative control is a stop sign; do not read the top hits.</b> Script 109 re-tests leads within each ancestry with a healthcare-utilisation covariate and ten ancestry PCs.</p>", "how") +
                   ev("What could still be wrong", "<ul><li>Binary diagnoses have little power below ~100 cases (the multiple sclerosis and psoriasis controls are underpowered).</li><li>Only the top hits per exposure set were followed up.</li><li>Native 535k cohort not used; it would give the power this scan lacks.</li></ul>", "cav"))


def s6():
    f_c = fig(F.sbars([(a, b, c) for a, b, c in PQTL_CTRL], (-60, 60), "Gene-to-protein controls: does gene content predict its own plasma protein? (z-score)", label_w=420, fmt="{:+.1f}",
                      marks=[(5.91, "+5.9"), (-5.91, "−5.9")]),
              "Long-read KIR gene content and C4 copy number predict the matching plasma proteins with z up to 51. HLA-allele effects on the HLA proteins are expected epitope (antibody-binding) effects. HLA-E 107G has no effect on soluble HLA-E.",
              "script 112 → ~/s09/olink/pqtl_hits_pub.tsv (VM)",
              dt(["Exposure → protein", "z"], [(a, "%+.2f" % b) for a, b, c in PQTL_CTRL]))
    f_a = fig(F.heat([a for a, _ in PQTL_ANC], ANC5, [v for _, v in PQTL_ANC], vmax=40, title="The same associations in each ancestry separately (z-score)", label_w=300, cell_w=100, fmt="{:+.0f}"),
              "Every association shown keeps its sign in all five ancestries (the 'replicates in 5/5 ancestries' criterion; across the top 60, two hits were 4/5). EAS values are smaller because that group has fewer carriers.",
              "script 112; per-ancestry re-estimation of the top hits",
              dt(["Association"] + ANC5, [(a,) + tuple("%+.1f" % x for x in v) for a, v in PQTL_ANC]))
    f_t = dt(["Allele or gene content", "Protein", "z (pooled)", "z (meta over ancestries)"], [(a, b, "%+.1f" % c, ("%+.1f" % d) if d is not None else "") for a, b, c, d in PQTL_TOP])
    return section("s6", "6. HLA, KIR and C4 against 5,415 plasma proteins", badge("known", "Validation") + badge("strong", "Atlas"),
                   "<b>The long-read KIR and C4 calls are independently validated, and a proteome-wide atlas is feasible.</b> In 8,884 people, KIR3DL1 gene content predicts plasma KIR3DL1 (z = 51), KIR2DL2/2DL3 content predicts the KIR2DL2/3 assay (z = ±29) and C4B copy number predicts C4B protein (z = 29). The scan is well calibrated (λ = 1.10) and every top-60 hit keeps its sign in at least four of the five ancestries (all five for nearly all). A UK Biobank atlas of HLA effects on plasma proteins (Europeans only) already exists; our additions are five ancestries, long-read-resolved alleles, KIR content and C4 copy number.",
                   cards([("144 × 5,415", "exposures × protein assays tested (780,000 associations)"), ("λ = 1.10", "genomic inflation: calibrated"), ("z = 51", "KIR3DL1 content → KIR3DL1 protein"), ("≥ 4 / 5", "ancestries with the same sign for every top-60 hit (5/5 for nearly all)")]),
                   ev("Background: how Olink measures a protein, what NPX and a pQTL are", "<p>Olink Explore HT measures ~5,400 plasma proteins with a proximity extension assay: two antibodies, each with a DNA tag, must bind the same protein before the tags pair and are counted. The result is NPX (log2 scale). A gene variant that shifts a protein's level is a <b>pQTL</b>. Plasma holds the <i>soluble</i> fraction, and because the readout depends on antibodies a variant can change the signal without changing the amount (an <b>epitope effect</b>).</p><figure>%s</figure>" % BG.fig_pea(), "bg", True) +
                   f_c + f_a + ev("Top 30 associations", f_t + "<p>All pass the family-wise threshold (|z| ≥ 5.9). Pooled z adjusts for age, sex, ancestry dummies and ten PCs; meta z combines per-ancestry estimates.</p>", "ev") +
                   ev("Scripts, models and thresholds", "<p>Script 112 builds a person × assay matrix (54 million Olink records, passing sample and assay QC, replicates averaged; assays with <70% coverage dropped), regresses every protein on every exposure (144: HLA alleles with ≥40 carriers, C4A/C4B copy number, KIR gene content) with age, age², sex, ancestry dummies and ten ancestry PCs by residual-projection OLS, and sets the family-wise threshold (|z| ≥ 5.91) by 10 permutations within ancestry×sex strata. Top hits are re-estimated within each ancestry.</p>", "how") +
                   ev("What could still be wrong", "<ul><li>Only the top 60 of 780,000 associations have been examined.</li><li>Plasma = soluble fraction, and antibodies can be sensitive to amino-acid changes: strong HLA→HLA-protein hits are probably epitope effects.</li><li>Cohort enrichment: 8,884 of the 9,969 Olink participants are in the long-read cohort, so power is high but the cohort is not population-representative.</li></ul>", "cav"))


def s7():
    f_shift = fig(F.grouped([a for a, _, _ in B14_SHIFT], [("B*14:02 carriers (456)", [b for _, b, _ in B14_SHIFT], "warn"), ("non-carriers (8,410)", [c for _, _, c in B14_SHIFT], "mut")], 3.2, "Plasma level, in standard deviations from the cohort mean", fmt="{:+.2f}", height=290, ymin=-0.8),
                  "Carriers sit 1.4–2.6 SD above everyone else for all three proteins (n = 456 heterozygotes; fewer than 20 homozygotes).",
                  "script 113 → ~/s09/olink/b1402_pub.txt (VM)",
                  dt(["Protein", "carriers (SD)", "non-carriers (SD)"], [(a, sd(b), sd(c)) for a, b, c in B14_SHIFT]))
    f_anc = fig(F.grouped([a for a, _, _ in B14_ANC], [("B*14:02 carriers", [b for _, b, _ in B14_ANC], "warn"), ("non-carriers", [c for _, _, c in B14_ANC], "mut")], 3.4, "LRPAP1 in each ancestry (SD from the cohort mean)", fmt="{:+.2f}", height=280, ymin=-0.8),
                "The same ~2.5–2.9 SD shift appears in AFR, AMR and EUR, whose HLA haplotype backgrounds differ. EAS has fewer than 20 carriers and SAS too few to show.",
                "script 113; per-ancestry carrier means",
                dt(["Ancestry", "carriers (SD)", "non-carriers (SD)"], [(a, sd(b), sd(c)) for a, b, c in B14_ANC]))
    f_sub = fig(F.sbars(B14_SUB, (-0.5, 3.2), "Which allele carries the effect? Shift in LRPAP1 per copy (SD)", label_w=330, fmt="{:+.2f}"),
                "The effect is allele-specific: B*14:02 +2.7 SD, its close relative B*14:01 about a third of that, B*14:02's HLA-C neighbour C*08:02 +2.4 SD only because it is inherited with B*14:02, while C*08:01 and other B and C alleles do nothing.",
                "script 113 (section B)",
                dt(["Allele (carriers)", "shift in LRPAP1 (SD)"], [(r[0], sd(r[1])) for r in B14_SUB if r[0] != "#"]))
    f_joint = fig(F.grouped([a for a, _, _ in B14_JOINT], [("B*14:02 (joint model)", [b for _, b, _ in B14_JOINT], "warn"), ("C*08:02 (joint model)", [c for _, _, c in B14_JOINT], "b")], 2.2, "Fitted together, how much of each effect belongs to each allele? (SD)", fmt="{:.2f}", height=280),
                  "B*14:02 carries the larger share for all three proteins; C*08:02 keeps a smaller share, so the true driver is B*14:02 or something tightly linked to it on that haplotype. DRB1*01:02 and DQB1*05:01 add nothing once B*14:02 and C*08:02 are modelled.",
                  "script 113 (section A)",
                  dt(["Protein", "B*14:02 (SD)", "C*08:02 (SD)"], [(a, "%.2f" % b, "%.2f" % c) for a, b, c in B14_JOINT]))
    f_aa = fig(F.dumbbell([(a, b, c) for a, b, c in B14_AA], (-30, 70), "Is it a single amino-acid position? z for each residue, alone vs after accounting for B*14:02", ("residue alone", "after accounting for B*14:02"), ("warn", "c"), label_w=260, xlabel="z-score for LRPAP1"),
               "Every single residue that tags B*14:02 (e.g. HLA-B position 110) collapses once B*14:02 itself is in the model. No shared residue explains the signal beyond the allele.",
               "script 113 (section C; IPD alignments)",
               dt(["Residue", "z alone", "z after B*14:02"], [(a, "%+.1f" % b, "%+.1f" % c) for a, b, c in B14_AA]))
    f_lig = fig(F.sbars(B14_LIGANDS, (-10, 10), "Specificity: B*14:02 against other LRP1 receptors, ligands and lysosomal enzymes (z-score)", label_w=300, bar_h=14, gap=5, fmt="{:+.1f}"),
                "Only LRPAP1, LRP1 and SGSH move. Every other protein on the panel in this pathway has |z| below 2.1; the first three are clipped at ±10 (their values are printed).",
                "script 113 (section D)",
                dt(["Protein", "z"], [(r[0], "%+.1f" % r[1]) for r in B14_LIGANDS if r[0] != "#"]))
    f_plr = fig(F.grouped(["LRPAP1", "LRP1", "SGSH"], [("plasma protein z", [78.8, 43.7, 35.7], "warn"), ("blood RNA z", [0.41, 1.51, 0.67], "a")], 85, "B*14:02 carriers vs non-carriers: z-score", fmt="{:.1f}", height=270),
                "The plasma shift is huge; the blood RNA of the same genes does not move at all (RNA z computed from the reported p-values).",
                "scripts 113 and 114 (b1402_rna_pub.txt on the VM)",
                dt(["Protein / gene", "plasma z", "blood-RNA z"], [("LRPAP1", 78.8, 0.41), ("LRP1", 43.7, 1.51), ("SGSH", 35.7, 0.67)]))
    f_rna = fig(F.sbars(B14_RNA, (-0.25, 0.25), "Blood RNA in B*14:02 carriers: change in log2 TPM", label_w=300, fmt="{:+.3f}"),
                "The haplotype has the expected cis effects on its own neighbours (C4A and C4B up, HLA-C down), proving the carriers are detected. The three proteins' genes show no change.",
                "script 114 → ~/s09/olink/b1402_rna_pub.txt (VM; 8,000 people, 447 carriers)",
                dt(["Gene", "log2 TPM change", "p"], [(r[0], "%+.3f" % r[1], r[3]) for r in B14_RNA if r[0] != "#"]))
    return section("s7", "7. HLA-B*14:02 and three unrelated plasma proteins", badge("open") + badge("known", "Partly known"),
                   "<b>One HLA allele shifts three proteins encoded on three other chromosomes, but the cause is unknown.</b> HLA-B*14:02 carriers (456 people) sit +2.6 SD for LRPAP1, +1.7 SD for LRP1 and +1.4 SD for SGSH. The effect is allele-specific, the same size in AFR, AMR and EUR, specific to these three of 5,415 proteins and not visible in blood RNA. The UK Biobank paper (Europeans) reports HLA explaining ~43% of LRPAP1, ~29% of SGSH and ~25% of LRP1 variance and calls it 'unclear and intriguing' without naming an allele; we attribute it to B*14:02 and show it replicates outside Europe.",
                   cards([("+2.6 SD", "LRPAP1 in B*14:02 carriers (LRP1 +1.7, SGSH +1.4)"), ("456", "carriers among 8,884 people"), ("8 of 5,415", "proteins past the family-wise threshold for B*14:02"), ("p ≥ 0.13", "blood RNA of all three genes")]),
                   ev("Background: haplotypes, trans effects and what the proteins are", "<p>HLA-B and HLA-C sit side by side on chromosome 6 and are inherited as a block (a <b>haplotype</b>): B*14:02 almost always comes with C*08:02, so either looks causal unless you model them together. The three affected proteins come from other chromosomes (4, 12, 17), so this is a <b>trans</b> effect. LRPAP1 is a chaperone of the LRP1 receptor family; SGSH is a lysosomal enzyme; they share no known HLA connection.</p><figure>%s</figure>" % BG.fig_haplotype(), "bg", True) +
                   f_shift + f_anc + ev("Which allele, and is it one amino acid?", f_sub + f_joint + f_aa, "ev", True) + ev("Specificity and the RNA test", f_lig + f_plr + f_rna, "ev") +
                   ev("More numbers", "<ul><li><b>Dose-response:</b> 0 copies −0.15 SD (8,410 people); 1 copy +2.59 SD (456); 2 copies: fewer than 20 people.</li><li><b>Co-variation:</b> residuals of LRPAP1/LRP1/SGSH correlate in carriers (0.61/0.46/0.24 for LRPAP1~LRP1, LRPAP1~SGSH, LRP1~SGSH) and in non-carriers (0.40/0.32/0.15): the three proteins move together in everyone, suggesting a shared factor.</li><li><b>Forward selection among 237 HLA alleles for LRPAP1:</b> B*14:02 (z 79) → C*08:02 (z 25.5) → C*07:18 (−3.8) → DRB1*01:02 (3.5).</li><li><b>Also affected (past the threshold):</b> MICA/MICB (−17), C4B (+10), HSPA1A (−7), GABARAP (+7), CDSN (−6), neighbouring MHC genes.</li></ul>", "ev") +
                   ev("Scripts, models and thresholds", "<p>Script 113 reuses the cached person × assay matrix of script 112. It models LRPAP1/LRP1/SGSH on B*14:02, C*08:02, DRB1*01:02 and DQB1*05:01 singly and jointly (ancestry, sex, age and ten PCs), forward-selects among all alleles, scans per-residue dosages built from IPD-IMGT/HLA protein alignments for HLA-B and HLA-C before and after conditioning on B*14:02, and lists the B*14:02 z-score for every assay. Script 114 extracts RNA for LRPAP1, LRP1, SGSH and related genes from the blood RNA-seq matrix and tests B*14:02, C*08:02 and B*14:01 dose.</p>", "how") +
                   ev("What it could be, and what would settle it", "<ul><li><b>Assay interference or antibody cross-reactivity:</b> a plasma factor specific to B*14:02 carriers (for example soluble HLA-B*14:02) could change the readout of several assays at once. The proteins also co-vary in non-carriers, which fits a shared technical factor too.</li><li><b>A real effect outside cells:</b> the haplotype could change how these proteins are secreted, cleared or shed (LRP1 is an endocytic receptor); blood RNA cannot see that.</li><li><b>Decisive:</b> measure LRPAP1 by ELISA or mass spectrometry in carriers; check the UK Biobank supplement (does it name B*14:02?); look for B*14:02-haplotype effects on LRPAP1 in liver or other tissue eQTL resources.</li></ul>", "cav") +
                   ev("Literature", "<p>Krishna et al., <i>Nature Communications</i> 2024 (UK Biobank, 45,330 Europeans, 2,940 Olink proteins): 504 proteins show HLA effects, >80% mapping to peptide-binding-groove positions; LRPAP1 (42.9% of variance explained by HLA), SGSH (29.2%) and LRP1 (24.9%) are among the top; no driving allele or residue is given for these in the text, and the authors note SGSH has no known immune role. Our result differs in attributing it to B*14:02 (with C*08:02 as tag) and in showing it outside Europeans. No assay-artefact analysis is reported there.</p>", "ev"))


def s8():
    rows = [(a, b, "<span class='badge %s'>%s</span> %s" % ("notrun" if k == "notrun" else "strong", "Not run" if k == "notrun" else "Run", c), d) for a, b, c, d, k in IDEAS]
    return section("s8", "8. Ideas from the think tanks: what was run and what was not", badge("notrun", "Partly not run"),
                   "<b>Of about 20 ideas from the think-tank rounds, 9 were run and the rest were not.</b> Several 'naive look' plots that the portfolio promised (template-distance histograms, KIR gene-count outliers, DRB3/4/5 consistency, genes called by platform, block lengths) were <b>never executed</b>: they need person-level data on the VM, which was unreachable when this report was rebuilt. They are the cheapest next step.",
                   "",
                   dt(["Idea", "Mode", "Status", "Where / why not"], rows) +
                   ev("What the 'naive looks' would be and how they would be made", "<ol><li><b>Template-distance histograms</b> (how far each person's HLA haplotypes are from the reference) by platform and ancestry: one column of the cohort membership table; a platform effect would show as shifted peaks.</li><li><b>KIR gene-count per person</b> by ancestry and platform: a missing-gene artefact would show as a low-count tail only on one platform.</li><li><b>DRB1 vs DRB3/4/5 consistency</b>: DRB1*03 should come with DRB3, DRB1*04 with DRB4, DRB1*15 with DRB5; exceptions are either real structural variants or calling gaps.</li><li><b>Genes called per person</b> by sequencing platform and trim tier.</li><li><b>Haplotype block lengths</b> between phase breaks, by ancestry.</li></ol><p>Panels 1–4 are implemented in <code>115_naive_looks.py</code> (written and tested on synthetic data; <b>not yet run on the real data</b>). Block lengths are not in the call tables. All outputs are aggregate histograms with counts under 20 suppressed, and take well under an hour of VM time.</p>", "how"))


def ledger_section():
    rows = [(a, b, "<span class='badge %s'>%s</span>" % ({"Strong": "strong", "Null": "null", "Killed": "known", "Weak": "lead", "Lead": "lead", "Open": "open", "Done": "known"}[c], c), d, '<a href="#s%d">%d</a>' % (n, n)) for a, b, c, d, n in LEDGER]
    return '<h2 id="ledger">Ledger of every hypothesis</h2>' + dt(["ID", "Hypothesis", "Status", "Evidence", "Section"], rows)


def methods():
    return '<h2 id="methods">Methods, data and reproduction</h2>' + \
        ev("Data used", "<ul><li><b>Long-read genomes:</b> up to 12,233 All of Us participants (PacBio Revio and Sequel II/IIe, ONT); HLA and KIR alleles called from assemblies with Immuannot against IPD-IMGT/HLA 3.55 (re-checked against 3.66 in S08).</li><li><b>Native HLA calls:</b> AoU short-read callers, 535k participants (12,208 overlap with long reads).</li><li><b>Blood RNA-seq:</b> RNA-SeQC2 gene TPM (3,141 people for the HLA-E tests; 8,000 for the B*14:02 RNA test).</li><li><b>Plasma proteomics:</b> Olink Explore HT, 9,969 participants, 5,415 assays; 8,884 unrelated long-read people used.</li><li><b>Health records:</b> OMOP ICD-10 3-character categories with at least two visits.</li><li><b>Ancestry:</b> AoU predicted ancestry and its first ten PCs.</li></ul>", "how") + \
        ev("Scripts (branch s09-discovery-loop)", "<p><code>100_native_vs_lr.py</code>, <code>101_rna_panel_extract.py</code>, <code>102_expression_dosage_tests.py</code>, <code>103_expression_within_g_group.py</code>, <code>105_novel_allele_evidence.py</code>, <code>106_nonclassical_landscape.py</code>, <code>107_phewas_engine.py</code>, <code>108_mechanistic_rna_tests.py</code>, <code>109_phewas_followup.py</code>, <code>111_olink_nkg2a_test.py</code>, <code>112_hla_kir_c4_pqtl_scan.py</code>, <code>113_b1402_pqtl_dissection.py</code>, <code>114_b1402_rna_test.py</code>; report: <code>116_s09_full_report.py</code> and <code>116a_figs.py</code>. Run log: <code>sprints/S09_discovery_loop/LOG.md</code>; hypothesis ledger: <code>LEDGER.md</code>.</p><p>Person-level files stay on the Workbench VM (mode 600). Everything shown here is aggregated; counts of 1–19 are written &lt;20.</p>", "how") + \
        ev("The loop used to generate and filter ideas", "<p>Haiku think-tank personas (moonshot, secure claim, naive, heuristic) proposed ideas, a stronger model with web search critiqued them, each surviving idea got a ledger row with a kill criterion before it ran, and every result passed a red-team checklist (positive and negative controls, chance baseline, relatedness, ancestry, platform, read support, independent truth, suppression).</p>", "bg")


def build(review_html=""):
    toc = [("glance", "At a glance"), ("s1", "1 Native vs long-read"), ("s2", "2 RNA nulls"), ("s3", "3 HLA-E/F/G"), ("s4", "4 HLA-E → NKG2A"), ("s5", "5 Disease scan"), ("s6", "6 Proteome atlas"),
           ("s7", "7 B*14:02"), ("s8", "8 Ideas run/not run"), ("ledger", "Ledger"), ("methods", "Methods")]
    lede = '<p class="sub">Omni-HLA · sprint S09 · rebuilt 28 Oct 2026 · aggregate results only · companion to the S08 report (are the "novel" alleles real?)</p>'
    body = overview(review_html) + s1() + s2() + s3() + s4() + s5() + s6() + s7() + s8() + ledger_section() + methods()
    return F.page("S09 discovery loop", "Omni-HLA · S09", "HLA, KIR and C4 in RNA, health records and plasma proteins: what held up", lede, toc, body)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--review", default="", help="HTML fragment with the independent reader's notes (placed in the overview)")
    a = ap.parse_args(argv)
    rv = open(a.review).read() if a.review and os.path.exists(a.review) else ""
    html = build(rv)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w").write(html)
    print("wrote", a.out, len(html) // 1024, "KB")


if __name__ == "__main__":
    main()
