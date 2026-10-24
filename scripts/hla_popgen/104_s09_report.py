"""104 S09 discovery-loop report (HTML, inline SVG; aggregate numbers only, every count shown is >= 20 or '<20').

  python3 scripts/hla_popgen/104_s09_report.py --out reports/hla_popgen/S09_discovery_loop/S09_report.html
Reuses the figure helpers of 94_s08_report.py.
"""
import argparse
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("r94", os.path.join(HERE, "94_s08_report.py"))
r94 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r94)
esc, grouped, hbars, CSS = r94.esc, r94.grouped, r94.hbars, r94.CSS

# native (AoU HLA-HD/Polysolver/OptiType) vs long-read truth, 12,208 people: genotype concordance at 2-field name vs at G-group level (upper bound)
GENES = ["HLA-A", "HLA-C", "DPB1", "DQA1", "DQB1", "DRB1"]
TWO_FIELD = [0.968, 0.899, 0.873, 0.522, 0.778, 0.965]
G_LEVEL = [0.987, 0.988, 0.993, 0.990, 0.989, 0.980]
# effect of long-read-only feature on RNA expression (log2 TPM per unit), beta +/- se
FOREST = [
    ("POS  KIR2DL1 dosage -> KIR2DL1", 0.1719, 0.0135, "pos"),
    ("POS  KIR2DL3 dosage -> KIR2DL3", 0.3008, 0.0106, "pos"),
    ("POS  C4A copies -> C4A", 0.0291, 0.0036, "pos"),
    ("POS  C4B copies -> C4B", 0.0692, 0.0039, "pos"),
    ("REV  KIR3DL1 dosage -> KIR3DL1", -0.1521, 0.0153, "rev"),
    ("REV  KIR2DS4 dosage -> KIR2DS4", -0.3774, 0.0169, "rev"),
    ("NEG  KIR2DL1 dosage -> CD3E", -0.0174, 0.0217, "neg"),
    ("NEG  C4 copies -> GAPDH", -0.0000, 0.0036, "neg"),
    ("TEST activating-KIR count -> NK score", -0.0012, 0.0044, "test"),
    ("TEST total KIR genes -> NK score", -0.0045, 0.0036, "test"),
    ("TEST C4 copies -> interferon score", -0.0208, 0.0106, "test"),
    ("TEST 1 vs 2 HLA-A calls -> HLA-A", 0.0245, 0.0147, "test"),
]
HLAE_ANC = [("EUR (n=1,038)", -0.0860, 0.0265, "pos"), ("AFR (n=466)", -0.1013, 0.0413, "pos"), ("AMR (n=595)", -0.0657, 0.0366, "pos"), ("EAS (n=428)", -0.0982, 0.0474, "pos"),
            ("SAS (n=448)", -0.1212, 0.0447, "pos"), ("META, 5 ancestries (n=2,975)", -0.0906, 0.0164, "test")]
HLAE_SPEC = [("HLA-E TPM (positive control)", -0.0659, 0.0068, "pos"), ("KLRC1 (NKG2A)", -0.0818, 0.0180, "test"), ("KLRD1 (CD94)", 0.0007, 0.0132, "neg"), ("KLRC2 (NKG2C)", 0.0052, 0.0197, "neg"),
             ("KLRK1 (NKG2D)", -0.0100, 0.0108, "neg"), ("NCR1", 0.0102, 0.0111, "neg"), ("KLRB1", 0.0022, 0.0153, "neg"), ("NCAM1 (CD56)", 0.0158, 0.0090, "neg"), ("IFNG", -0.0062, 0.0119, "neg")]
COL = {"pos": "--c", "rev": "--warn", "neg": "--mut", "test": "--b"}


def forest(rows, width=900, label_w=300, row_h=26, lim=0.45):
    h = len(rows) * row_h + 50
    x0 = label_w + (width - label_w - 20) / 2
    sc = (width - label_w - 40) / (2 * lim)
    out = ['<svg viewBox="0 0 %d %d" role="img">' % (width, h)]
    out.append('<line x1="%.1f" x2="%.1f" y1="8" y2="%d" stroke="currentColor" opacity=".35"/>' % (x0, x0, h - 30))
    for t in (-0.4, -0.2, 0, 0.2, 0.4):
        x = x0 + t * sc
        out.append('<text class="m" x="%.1f" y="%d" font-size="11" text-anchor="middle">%s</text>' % (x, h - 12, "%.1f" % t))
    for i, (lab, b, se, k) in enumerate(rows):
        y = 18 + i * row_h
        c = COL[k]
        out.append('<text x="0" y="%d" font-size="13">%s</text>' % (y + 4, esc(lab)))
        out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" style="stroke:var(%s)" stroke-width="3"/>' % (x0 + (b - 1.96 * se) * sc, x0 + (b + 1.96 * se) * sc, y, y, c))
        out.append('<circle cx="%.1f" cy="%d" r="5" style="fill:var(%s)"/>' % (x0 + b * sc, y, c))
    out.append("</svg>")
    return "".join(out)


def loop():
    boxes = [("Think tank", "Haiku personas, 4 modes", "--a"), ("Critic + literature", "stronger model, web search", "--b"), ("Ledger row", "hypothesis + kill criterion first", "--c"),
             ("MVP on the VM", "aggregates, controls, <20 rule", "--c"), ("Red-team", "artefact checklist", "--warn"), ("Next round", "orthogonal directions", "--a")]
    w, h, bw, bh = 960, 130, 142, 70
    out = ['<svg viewBox="0 0 %d %d" role="img">' % (w, h), '<defs><marker id="ar2" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6z" fill="currentColor" opacity=".6"/></marker></defs>']
    for i, (t, s, c) in enumerate(boxes):
        x = 8 + i * (bw + 20)
        out.append('<rect x="%d" y="20" width="%d" height="%d" rx="10" style="fill:none;stroke:var(%s);stroke-width:2"/>' % (x, bw, bh, c))
        out.append('<text x="%d" y="52" font-size="13" font-weight="600" text-anchor="middle">%s</text>' % (x + bw / 2, esc(t)))
        out.append('<text class="m" x="%d" y="72" font-size="10.5" text-anchor="middle">%s</text>' % (x + bw / 2, esc(s)))
        if i < len(boxes) - 1:
            out.append('<path d="M%d 55 h18" stroke="currentColor" opacity=".5" marker-end="url(#ar2)"/>' % (x + bw + 1))
    out.append("</svg>")
    return "".join(out)


LEDGER = [
    ("S08-A", "HLA-F T237I/P205T are real inherited alleles", "passed", "single flanking haplotype; no platform/centre skew; concordant relatives"),
    ("S08-B", "...and they are transcribed", "passed", "variant haplotype carries ~54% / ~65% of informative HLA-F reads (control F*01:25 ~58%)"),
    ("S09-01", "native short-read calls are wrong where long reads are right", "killed", "2-field concordance 0.52-0.97, but 0.98-0.99 at G-group level: resolution, not error"),
    ("S09-02", "native errors distort HLA-disease odds ratios", "killed", "follows from S09-01 (peptide-binding-domain calls are right)"),
    ("S09-10a", "KIR content shifts NK transcripts", "null", "n=5,120 people; p=0.79 / 0.22; negative control null"),
    ("S09-10b", "C4 copy number tunes interferon tone", "null", "nominal p=0.049 only; positive/negative controls behave"),
    ("S09-11", "non-peptide-binding variation changes class II expression", "weak signal", "0.2-0.6% of variance (perm p 0.002-0.024); mapping bias not excluded"),
    ("S09-12", "unnamed alleles are supported by raw long reads", "passed", "18/26 supported in all carriers, 0/8 controls; five HLA-F loss-of-function alleles real; E/G stop-gains are artefacts"),
    ("S09-13", "HLA-E/F/G diversity landscape", "done", "E balanced, F nearly monomorphic, G null haplotypes 4-5x higher in AFR"),
    ("S09-15", "B*27 -> UPR signature in blood", "null", "beta -0.016 +/- 0.024, n=5,642; controls flat"),
    ("S09-16", "HLA-E 107G lowers NKG2A per NK cell", "passed", "meta beta -0.091, p 3.6e-8, 5/5 ancestries, specific to KLRC1; survives 10 ancestry PCs (meta -0.085, p 2.6e-7) and 3 independent KLRC1-region cis variants (meta -0.073, p 7e-7)"),
    ("S09-17", "HLA-F LoF lowers its own transcript", "weak signal", "R294* haplotype 26% of reads; one carrier each"),
    ("S09-03", "HLA-F T237I/P205T x KIR3DS1 vs viral EHR codes", "proposed", "needs a power check first (carriers <50)"),
    ("S09-04", "HLA-A deletion shifts KIR content", "proposed", "carriers likely <20"),
    ("S09-05", "allele x phecode PheWAS with >= 50 negative controls", "proposed", "native 535k, HLA-prior phecodes first"),
    ("S09-07", "parent-of-origin effect on HLA ASE", "proposed", "check trio x RNA overlap"),
]
PHEWAS_TEXT = ("The HLA / KIR / C4 x ICD-10 phenome scan (594 outcome categories with at least 30 cases in 11,856 unrelated long-read people; score tests with ancestry, sex and age; negative-control outcome chapters; 20-permutation max|z| thresholds) first failed its own calibration check (numerical problem with near-constant common alleles) and was rerun. Corrected calibration is acceptable for HLA (lambda 1.12; negative chapters 1.13) and C4, but not for KIR (negative chapters lambda 1.79), so KIR hits are treated as ancestry-linked gene-content structure until shown otherwise. No HLA allele passes the family-wise threshold (max|z| 6.98). Known biology is recovered at the nominal level for type 1 diabetes (DQB1*02:01 -> E10, z 5.5; DQA1*05:01/*03:03 -> E10). The B*27:05 -> iridocyclitis hit (H20, z 4.7) looked like the classical B*27 anterior-uveitis association but does not hold within ancestry: it comes from AMR and one sparse SAS stratum (21 carriers), and shows nothing in EUR, so it is NOT counted as a recovery. The strongest HLA signal, DQA1*05:03 -> E87/K72/K70/A41 (z 5.6-5.9), hits four unrelated outcomes that are all severe-liver/critical-illness phenotypes, Within ancestry it is carried by AMR alone (234 of the 304 carriers; E87 z 6.45 in AMR, unchanged by healthcare-utilisation adjustment) and it collapses once ten ancestry PCs are added within AMR (E87 z 2.0; K72 1.0, K70 1.3, A41 1.8; meta z 1.8): DQA1*05:03 is a marker of ancestry structure inside the AMR group (a Native-American-ancestry share that tracks liver-disease and critical-illness risk), not an HLA effect, and is retired. C4A copy number -> I77 (other arterial disorders) is concordant in sign in all six ancestries and survives PC adjustment (EUR z 4.32 with PCs; PC- and utilisation-adjusted meta z 4.14 against the 4.07 threshold, negative-control chapters max|z| 2.5), but it rests on EUR (leave-EUR-out meta z 1.3) and a single outcome code, so it is a lead requiring replication (C4 copy number from short reads in the native cohort). Odds ratios in sparse cells (tens to hundreds) are score-test artefacts and are not interpreted. Honest summary: the scan validates the pipeline and yields leads, not a discovery.")
TAG = {"done": "known", "passed": "known", "killed": "art", "null": "art", "weak signal": "new", "proposed": "new"}


def build():
    ledger = "".join("<tr><td>%s</td><td>%s</td><td><span class='tag %s'>%s</span></td><td>%s</td></tr>" % (esc(a), esc(b), TAG[c], esc(c), esc(d)) for a, b, c, d in LEDGER)
    tpl = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>S09 discovery loop report</title><style>@@CSS@@</style></head><body><main>
<p class="sub">Omni-HLA · sprint S09 · 24 Oct 2026 · aggregate results only (counts of 1–19 are written &lt;20) · builds on S08</p>
<h1>Discovery loop: what survived, what died, what is next</h1>
<p>Marc asked for a loop — think tank → falsifiable MVP → critics and literature → think tank — with a mix of modes: moonshot, secure claims, naive exploration, heuristic search. Three rounds of Haiku think tanks (≈20 agents, ≈100 ideas) and a web-searching critic produced a portfolio; the first seven experiments have now run on the VM. <b>Four hypotheses died or came back null, two held up, and the method itself got sturdier.</b></p>
<div class="cards">
<div class="card"><b>0.52 → 0.99</b><span>DQA1 concordance of the AoU-native caller vs long reads, at 2-field name vs peptide-binding-exon (G-group) resolution</span></div>
<div class="card"><b>transcribed</b><span>both unnamed HLA-F alleles are expressed in blood (≈54% and ≈65% of informative reads)</span></div>
<div class="card"><b>null, n=5,120</b><span>KIR gene content does not shift NK transcripts; positive and negative controls behave</span></div>
<div class="card"><b>0.2–0.6%</b><span>of class II expression variance is explained by variation <i>outside</i> the peptide-binding exons within an allele family</span></div></div>

<h2>1. The loop</h2><figure>@@LOOP@@<figcaption>Every experiment gets a ledger row with its kill criterion before it runs; every result goes through a red-team checklist (positive and negative control, chance baseline, relatedness, ancestry, platform, read support, independent truth, pre-registered family, suppression).</figcaption></figure>

<h2>2. Native short-read calls are right where it matters, and blind where it doesn't</h2>
<figure>@@BARS@@<figcaption>Up to 12,208 people with both calls (fewer per gene). At the usual 2-field <i>name</i> level the native caller seems to fail at DQA1, DQB1 and DPB1; scored at G-group level (identity over the peptide-binding exons, the resolution exon-based typers can deliver) concordance is 98–99% for every gene and ancestry (a permissive upper bound). DQA1*05:01, *05:03 and *05:05 share one G group, as do DQB1*02:01/*02:02. So the "native errors distort disease associations" idea is dead — which is also what long reads <i>add</i>: variation outside exons 2/3, structural variation and phase.</figcaption></figure>

<h2>3. Do long-read-only features move blood RNA? Controls first</h2>
<figure>@@FOREST@@<figcaption>Effect on log2 TPM per unit of the predictor, with 95% CIs (green = positive controls that must be positive; red = sign-reversed controls; grey = negative controls; orange = tests). KIR2DL1/2DL3 and C4A/C4B copy number drive their own transcripts and the negative controls are flat, so the pipeline is valid. KIR2DS4 and KIR3DL1 come out reversed: they share one locus and reads mis-map between them, so KIR gene-level TPM cannot validate every KIR gene. On the actual questions: no effect of activating-KIR count or total KIR genes on an NK transcript score (n = 5,120), C4 copy number barely touches interferon tone (nominal p = 0.049), and one-vs-two HLA-A calls is a call-quality, not a deletion, signal. DRB5 presence drives DRB5 expression by +2.75 log2 (off-scale; positive control).</figcaption></figure>

<h2>4. Alleles identical in the groove can still differ in expression</h2>
<table><tr><th>Gene</th><th>G-group identity explains</th><th>Extra explained by 3-field allele within G group</th><th>Permutation p</th></tr>
<tr><td>HLA-DQA1</td><td>80%</td><td>0.58%</td><td>0.002</td></tr><tr><td>HLA-DPA1</td><td>6%</td><td>0.35%</td><td>0.002</td></tr><tr><td>HLA-DQB1</td><td>75%</td><td>0.17%</td><td>0.024</td></tr><tr><td>HLA-C</td><td>28%</td><td>0.04%</td><td>0.48</td></tr>
<tr><td>HLA-A, B, DPB1, DRB1</td><td>23%, 6%, 12%, —</td><td colspan="2">untestable: no G group with ≥2 allele classes of ≥20 homozygous people</td></tr></table>
<div class="note"><b>Guard.</b> RNA is mapped to GRCh38, so haplotypes that diverge from the reference may look lower for mapping reasons, and ancestry is only coarsely adjusted. Next step: re-quantify these class II groups on personal two-haplotype references (the S05/S08 machinery) before calling it regulatory.</div>

<h2>5. Second wave: raw reads, landscape, mechanisms</h2>
<h3>5a. Raw long-read evidence for the unnamed alleles</h3>
<p>For 26 candidate unnamed signatures we located the exact GRCh38 position from the GENCODE CDS model (checked against the assembly bases) and piled up the raw AoU long-read alignments, independent of the assembly the call came from, in up to 8 carriers and 8 matched non-carriers. <b>18 of 26 are supported in every examined carrier and in 0 of 8 controls</b>: HLA-F T237I and P205T, T279I, G128R, G37S, V309I, A322T, M119I, Q247R, three recurrent synonymous signatures, and <b>five HLA-F loss-of-function alleles (stop-gains E70*, Q75*, W265*, R294*; start-loss M1V; each in fewer than 20 unrelated people)</b>. The HLA-E (Y28*, G321*) and HLA-G (R181*, Y142*) stop-gain calls and DPA1 E210G have no raw-read support: assembly artefacts. Three “supported” signatures (HLA-C and HLA-G synonymous, DRB1 S29A) are uninformative because the alt base is in every control. So <b>HLA-F, not E or G, is where real protein-truncating variation sits</b>.</p>
<h3>5b. Do the HLA-F loss-of-function alleles lower their own transcript?</h3>
<p>One heterozygous carrier per allele on a personal two-haplotype reference: the stop-gain <b>R294*</b> haplotype carries only 26% of informative HLA-F reads (bias-corrected log2 −1.23, ≈2.3× lower; consistent with nonsense-mediated decay), the start-loss <b>M1V</b> 40% (transcript persists, translation lost), and the missense variants 35–50%. Single-carrier evidence only.</p>
<h3>5c. HLA-E, F and G: diversity landscape by ancestry</h3>
<p>At equal sample size (400 unrelated people per ancestry): <b>HLA-E is balanced</b> (expected heterozygosity 0.50–0.52 in all six ancestries — the stable R107G polymorphism); <b>HLA-F is nearly monomorphic and ancestry-structured</b> (top allele 96% in EAS, 80% in AFR); <b>HLA-G carries 4–5× more null/questionable-suffix haplotypes in AFR (6.7%) than in SAS or EUR (1.4–1.6%)</b>. For scale, HLA-B has 56–86 distinct 2-field alleles.</p>
<h3>5d. A mechanism that failed, and one that did not</h3>
<p><b>HLA-B*27 does not induce an unfolded-protein-response signature in blood RNA</b> (UPR module β = −0.016 ± 0.024, p = 0.49, n = 5,642; ERAP1/2, IL23A, BiP, XBP1 all null; B*07, B*57, housekeeping and T-cell controls flat): a well-powered human test of the 1999 misfolding hypothesis in whole blood, informative as a null.</p>
<p><b>HLA-E R107G (E*01:03) lowers NKG2A (KLRC1) expression per NK cell.</b> The effect is specific to KLRC1, independent of NK abundance, interferon tone, CD8/γδ T-cell content, KLRC2 (adaptive-NK proxy) and the HLA leader peptides, has the same sign in all five ancestries, is not mediated by HLA-E mRNA, and is near-dominant (heterozygote β −0.144, homozygote −0.159). No prior report was found (closest: ~100–200 donors, no effect on NKG2A frequency); novelty is unverified.</p>
<figure>@@A1@@<figcaption>Fully adjusted effect of HLA-E 107G dose on KLRC1 (log2 TPM per allele, 95% CI) by ancestry and fixed-effect meta-analysis (Q = 1.05 for 4 d.f., leave-one-ancestry-out meta p between 8e-6 and 1e-7).</figcaption></figure>
<figure>@@A2@@<figcaption>Specificity (age, sex, ancestry, composition adjusted): the positive control (HLA-E mRNA, p = 2e-20) and KLRC1 move; the CD94 partner (KLRD1), the activating paralogue (KLRC2), and other NK receptors do not.</figcaption></figure>
<div class="note"><b>Not yet excluded.</b> KLRC1 cis SNPs and KLRC2 deletion (need genotypes), CMV serostatus (adaptive NK), and bulk-RNA cell-type attribution. The leader −21M peptide-supply hypothesis is only nominal (meta p ≈ 0.02, attenuated by NK abundance).</div>
<h3>5e. Hypothesis-free association scan</h3>
<p>@@PHEWAS@@</p>

<h2>6. Ledger</h2><table><tr><th>ID</th><th>Hypothesis</th><th>Status</th><th>Evidence / note</th></tr>@@LEDGER@@</table>

<h2>7. What the think tanks produced (themes, not claims)</h2>
<ul><li><b>Moonshots:</b> HLA-F × KIR3DS1 epistasis; HLA-A deletion as NK-ligand compensation; parent-of-origin HLA expression; balancing-selection halos; MICA intron-5 STR as an expression switch; C4 as an interferon dial (tested: null).</li>
<li><b>Secure claims:</b> native-vs-long-read benchmark (done: it is a resolution paper); phased MHC haplotype atlas; trio Mendelian error per gene; KIR gene-content atlas with a measured short-read false-negative rate; IPD submission of the validated HLA-F alleles.</li>
<li><b>Naive looks:</b> ten 10-minute tables and plots (template-distance histograms, KIR gene-count outliers, DRB3/4/5 co-occurrence gaps, warning rate by sequencing date, block lengths).</li>
<li><b>Automated search:</b> allele × phecode PheWAS with permutation nulls and ≥50 negative-control phenotypes; residue-level scan; haplotype-window scan; anomaly detection on per-person divergence.</li>
<li><b>Power-aware disease designs:</b> residue/pocket burden two-stage (535k → 12k), blood counts and labs before disease, drug-hypersensitivity positive controls (abacavir/B*57:01, carbamazepine/B*15:02), COVID severity from predicted epitope load with a negative-control proteome.</li></ul>

<h2>8. Lessons that now run the project</h2>
<ul><li>A concordance number without a resolution definition is not an error rate (DQA1 0.52 vs 0.99).</li><li>"Novel" without a reference-release date is database age (S08).</li>
<li>Gene-level TPM validates some KIR genes and not others; controls decide which tests are allowed.</li><li>Quantitative traits (RNA) have power that binary EHR labels (&lt;20 carriers) do not.</li></ul>

<h2>9. Decisions for Marc</h2>
<ol><li>Commit the S08/S09 report HTMLs and suppressed tables? (open small-cell disclosure item.)</li>
<li>Next experiment: personal-reference quantification of the class II within-G effect (mapping-bias test), or the parent-of-origin ASE test (needs trio × RNA overlap check), or the PheWAS engine with negative controls.</li>
<li>IPD submission path for HLA-F T237I/P205T (needs ≥3 unrelated carriers confirmed in raw reads plus an independent trio).</li></ol>

<h2>10. Reproduce</h2><p><code>100_native_vs_lr.py</code>, <code>101_rna_panel_extract.py</code>, <code>102_expression_dosage_tests.py</code>, <code>103_expression_within_g_group.py</code>, <code>104_s09_report.py</code> on branch <code>s09-discovery-loop</code>; ledger in <code>sprints/S09_discovery_loop/</code>. VM runs: one core, nice 10, throwaway clones, private tables mode 600 on the VM only.</p>
</main></body></html>"""
    bars = grouped([TWO_FIELD, G_LEVEL], GENES, [("2-field name", "--b"), ("G-group level", "--c")], "Genotype concordance, native vs long-read truth (up to 12,208 people per gene)", ymax=1.1, height=290, fmt="{:.1f}")
    tpl = tpl.replace("@@A1@@", forest(HLAE_ANC, lim=0.2, label_w=260)).replace("@@A2@@", forest(HLAE_SPEC, lim=0.2, label_w=260))
    for k, v in (("@@CSS@@", CSS), ("@@LOOP@@", loop()), ("@@BARS@@", bars), ("@@FOREST@@", forest(FOREST)), ("@@LEDGER@@", ledger), ("@@PHEWAS@@", PHEWAS_TEXT)):
        tpl = tpl.replace(k, v)
    return tpl


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w").write(build())
    print("wrote", a.out)


if __name__ == "__main__":
    main()
