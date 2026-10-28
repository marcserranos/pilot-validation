"""104a Plain-language background drop-downs (inline SVG drawings) for the S09 report. Imported by 104_s09_report.py.

Each entry of BACKGROUND is a <details> block: a short explanation for a reader who is not an immunogeneticist plus one schematic drawing.
Drawings are schematic, not data. Colours use the report's CSS variables (--a, --b, --c, --warn, --mut)."""
from html import escape as esc

EXTRA_CSS = """
details.bg{border:1px solid rgba(128,128,128,.35);border-radius:10px;padding:2px 14px;margin:14px 0;background:rgba(128,128,128,.05)}
details.bg>summary{cursor:pointer;font-weight:600;padding:8px 0}
details.bg[open]>summary{border-bottom:1px solid rgba(128,128,128,.25);margin-bottom:8px}
details.bg p{margin:.5em 0}
details.bg figure{margin:10px 0}
details.bg svg{width:100%;height:auto}
"""


class S:
    """Tiny SVG builder."""
    _n = 0

    def __init__(self, w, h):
        S._n += 1
        self.w, self.h, self.id = w, h, "bgm%d" % S._n
        self.parts = ['<svg viewBox="0 0 %d %d" role="img"><defs><marker id="%s" markerWidth="9" markerHeight="9" refX="7" refY="3.5" orient="auto"><path d="M0 0 L7 3.5 L0 7z" fill="currentColor" opacity=".7"/></marker></defs>' % (w, h, self.id)]

    def rect(self, x, y, w, h, col="--a", fill=0.0, rx=8, sw=2, dash=False):
        self.parts.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%d" style="fill:var(%s);fill-opacity:%.2f;stroke:var(%s);stroke-width:%d%s"/>' % (x, y, w, h, rx, col, fill, col, sw, ";stroke-dasharray:5 4" if dash else ""))

    def circle(self, x, y, r, col="--a", fill=0.15, sw=2):
        self.parts.append('<circle cx="%.1f" cy="%.1f" r="%.1f" style="fill:var(%s);fill-opacity:%.2f;stroke:var(%s);stroke-width:%d"/>' % (x, y, r, col, fill, col, sw))

    def text(self, x, y, s, size=12, anchor="middle", bold=False, muted=False):
        self.parts.append('<text%s x="%.1f" y="%.1f" font-size="%s" text-anchor="%s"%s>%s</text>' % (' class="m"' if muted else "", x, y, size, anchor, ' font-weight="600"' if bold else "", esc(s)))

    def line(self, x1, y1, x2, y2, col="--mut", dash=False, arrow=False, sw=2):
        self.parts.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" style="stroke:var(%s);stroke-width:%d%s"%s/>' % (x1, y1, x2, y2, col, sw, ";stroke-dasharray:5 4" if dash else "", ' marker-end="url(#%s)"' % self.id if arrow else ""))

    def done(self):
        return "".join(self.parts) + "</svg>"


def block(title, paras, svg, caption):
    return '<details class="bg"><summary>%s</summary>%s<figure>%s<figcaption>%s</figcaption></figure></details>' % (esc(title), "".join("<p>%s</p>" % p for p in paras), svg, esc(caption))


def fig_resolution():
    s = S(900, 250)
    names = [("02", "allele group", "roughly the old blood-test family", "--a"), ("01", "protein", "different amino acids", "--b"), ("01", "silent DNA change", "same protein", "--c"), ("01", "outside the protein", "introns / UTRs", "--mut")]
    s.text(450, 22, "HLA-A*02:01:01:01 : an allele name has up to four fields", 14, bold=True)
    for i, (n, t, c, col) in enumerate(names):
        x = 40 + i * 215
        s.rect(x, 36, 195, 44, col, 0.12)
        s.text(x + 97, 65, "%s : %s" % (n, t), 13, bold=True)
        s.text(x + 97, 98, c, 11, muted=True)
    s.text(450, 128, "The gene is cut into 8 exons; the protein's peptide-binding groove is built from exons 2 and 3", 12.5, bold=True)
    for e in range(8):
        x = 90 + e * 100
        hot = e in (1, 2)
        s.rect(x, 140, 80, 34, "--warn" if hot else "--mut", 0.2 if hot else 0.05)
        s.text(x + 40, 162, "exon %d" % (e + 1), 12)
    s.parts.append('<path d="M190 182 v8 h180 v-8" fill="none" style="stroke:var(--warn);stroke-width:2"/>')
    s.text(280, 208, "G group = alleles identical over exons 2+3", 12, bold=True)
    s.text(280, 224, "(what T cells and NK cells actually 'see')", 11, muted=True)
    s.text(700, 208, "Two callers can disagree on the NAME", 12, bold=True)
    s.text(700, 224, "and still agree on the G group (our DQA1: 0.52 vs 0.99)", 11, muted=True)
    return s.done()


def fig_reads():
    s = S(900, 300)
    s.text(225, 20, "A. Is the variant in the raw reads?", 14, bold=True)
    s.text(225, 42, "carrier", 12, bold=True)
    for r in range(6):
        y = 52 + r * 13
        x0 = 40 + (r % 3) * 25
        s.rect(x0, y, 280, 8, "--b", 0.25, rx=3, sw=1)
        s.rect(215, y - 1, 8, 10, "--warn", 0.9, rx=1, sw=1)
    s.text(225, 140, "every read carries the variant (orange) at the same spot", 11, muted=True)
    s.text(225, 168, "non-carrier control", 12, bold=True)
    for r in range(6):
        y = 178 + r * 13
        x0 = 40 + (r % 3) * 25
        s.rect(x0, y, 280, 8, "--b", 0.25, rx=3, sw=1)
    s.text(225, 266, "no read does: the variant is real inherited DNA, not an assembly mistake", 11, muted=True)
    s.text(225, 286, "(an assembly artefact would show the change in the assembly but NOT in the reads)", 11, muted=True)
    s.text(660, 20, "B. A 'stop-gain' and nonsense-mediated decay", 14, bold=True)
    s.rect(500, 40, 320, 34, "--c", 0.12)
    s.text(660, 62, "normal mRNA  →  full-length protein", 12)
    s.rect(500, 100, 320, 34, "--warn", 0.12)
    s.text(660, 122, "premature STOP in mRNA  →  cell destroys it", 12)
    s.parts.append('<text x="810" y="64" font-size="18" style="fill:var(--c)">✓</text><text x="808" y="124" font-size="16" style="fill:var(--warn)">✗</text>')
    s.text(660, 164, "In a carrier, half the HLA-F mRNA comes from each copy.", 12)
    s.text(660, 182, "If the stop copy is destroyed we see only ~26% of reads from it", 12)
    s.text(660, 200, "(HLA-F R294*), instead of ~50%.", 12)
    s.text(660, 226, "That is the signature of nonsense-mediated decay,", 11, muted=True)
    s.text(660, 242, "i.e. the allele really is a loss-of-function allele.", 11, muted=True)
    return s.done()


def fig_nkg2a():
    s = S(900, 300)
    s.text(450, 22, "How a natural killer (NK) cell decides not to attack a healthy cell", 14, bold=True)
    s.rect(30, 50, 300, 170, "--c", 0.08, rx=20)
    s.text(180, 74, "healthy cell in your body", 13, bold=True)
    s.rect(125, 118, 110, 58, "--c", 0.25)
    s.text(180, 142, "HLA-E", 14, bold=True)
    s.text(180, 160, "(holds a peptide)", 11)
    s.circle(180, 196, 9, "--warn", 0.6)
    s.text(180, 242, "peptide = a piece cut from the leader (signal)", 11, muted=True)
    s.text(180, 258, "sequence of HLA-A / B / C / G: 'my class I is working'", 11, muted=True)
    s.rect(570, 50, 300, 170, "--b", 0.08, rx=20)
    s.text(720, 74, "NK cell", 13, bold=True)
    s.rect(650, 118, 140, 58, "--b", 0.25)
    s.text(720, 140, "CD94 / NKG2A", 12.5, bold=True)
    s.text(720, 158, "(genes KLRD1 / KLRC1)", 11)
    s.line(240, 148, 645, 148, "--a", arrow=True, sw=3)
    s.text(440, 138, "binds", 12, bold=True)
    s.text(440, 172, "'self is intact: do NOT kill'", 11, muted=True)
    s.text(720, 242, "an inhibitory brake; fewer receptors = weaker brake", 11, muted=True)
    s.rect(120, 272, 660, 24, "--warn", 0.1)
    s.text(450, 289, "Our result: the stable HLA-E variant 107G (R107G) goes with ~6% less KLRC1 (NKG2A) mRNA per copy", 12, bold=True)
    return s.done()

def fig_confound():
    s = S(900, 270)
    s.text(450, 22, "Why a strong association can be false: confounding by ancestry within a group", 14, bold=True)
    s.rect(300, 40, 300, 44, "--warn", 0.12)
    s.text(450, 67, "admixture proportion (varies inside 'AMR')", 12.5, bold=True)
    s.rect(60, 170, 270, 44, "--a", 0.12)
    s.text(195, 197, "allele frequency (DQA1*05:03)", 12.5, bold=True)
    s.rect(570, 170, 270, 44, "--b", 0.12)
    s.text(705, 197, "outcome (liver disease, sepsis)", 12.5, bold=True)
    s.line(400, 86, 230, 166, "--mut", arrow=True)
    s.line(500, 86, 670, 166, "--mut", arrow=True)
    s.line(335, 192, 565, 192, "--warn", dash=True)
    s.text(450, 184, "looks like an allele effect", 12, bold=True)
    s.text(450, 208, "(spurious)", 11, muted=True)
    s.text(450, 244, "Checks we run: ancestry principal components, tests inside each ancestry, negative-control outcomes (injuries, etc.), label permutations.", 11, muted=True)
    s.text(450, 260, "DQA1*05:03 failed the PC check (z 6.4 → 2.0); HLA-E 107G passed it.", 11, muted=True)
    return s.done()


def fig_cis():
    s = S(900, 290)
    s.text(450, 22, "Could a nearby variant at KLRC1 explain the HLA-E result?", 14, bold=True)
    s.text(60, 52, "chromosome 6", 12, anchor="start", bold=True)
    s.rect(60, 60, 360, 26, "--mut", 0.12, rx=13)
    s.rect(180, 56, 100, 34, "--a", 0.3)
    s.text(230, 78, "HLA-E 107G", 12, bold=True)
    s.text(60, 142, "chromosome 12", 12, anchor="start", bold=True)
    s.rect(60, 150, 780, 26, "--mut", 0.12, rx=13)
    s.rect(330, 146, 100, 34, "--b", 0.3)
    s.text(380, 168, "KLRC1 gene", 12, bold=True)
    s.rect(620, 146, 110, 34, "--c", 0.3)
    s.text(675, 168, "cis variant", 12, bold=True)
    s.line(250, 92, 360, 142, "--warn", arrow=True, sw=3)
    s.text(330, 104, "our effect: −0.09 log2 TPM", 11.5, bold=True, anchor="start")
    s.line(616, 163, 436, 163, "--c", arrow=True, sw=3)
    s.text(526, 140, "known strong effect (p = 3.7e-143)", 11.5, bold=True)
    s.text(640, 70, "different chromosomes: HLA-E 107G and the cis", 11, muted=True, anchor="start")
    s.text(640, 85, "variant are linked only through population", 11, muted=True, anchor="start")
    s.text(640, 100, "structure (measured |r| ≤ 0.05)", 11, muted=True, anchor="start")
    s.text(450, 212, "So we (1) included ten ancestry PCs and (2) added the three strongest cis variants as covariates.", 12, bold=True)
    s.text(450, 234, "If the HLA-E effect were really the cis variant's, it would vanish. It shrank by ~17% and stayed at p ≈ 1e-7.", 12)
    s.text(450, 262, "Because the cis effect is so strong, this control had plenty of power to catch a confound.", 11, muted=True)
    return s.done()

def fig_pea():
    s = S(900, 320)
    s.text(450, 22, "How Olink measures a protein in plasma (proximity extension assay)", 14, bold=True)
    s.circle(160, 150, 32, "--b", 0.3)
    s.text(160, 154, "protein", 12, bold=True)
    for (x, y, dx) in ((110, 105, -1), (210, 105, 1)):
        s.line(x, y, 160 + dx * 18, 130, "--a", sw=3)
        s.line(x - 8, y - 18, x + 8, y - 18, "--a", sw=4)
        s.line(x, y - 18, x, y, "--a", sw=3)
    s.text(112, 62, "antibody A", 11)
    s.text(210, 62, "antibody B", 11)
    s.text(160, 200, "each antibody carries a DNA tag; both only", 11, muted=True)
    s.text(160, 214, "when bound close together → tags pair up", 11, muted=True)
    s.line(250, 150, 330, 150, "--mut", arrow=True, sw=3)
    s.rect(335, 128, 150, 44, "--c", 0.12)
    s.text(410, 148, "copy + count the DNA", 12, bold=True)
    s.text(410, 164, "(sequencing)", 11)
    s.line(488, 150, 560, 150, "--mut", arrow=True, sw=3)
    s.rect(565, 128, 140, 44, "--a", 0.12)
    s.text(635, 148, "NPX value", 12, bold=True)
    s.text(635, 164, "(log2 scale, per person)", 11)
    s.text(730, 134, "1 SD  =  the typical", 11.5, bold=True, anchor="start")
    s.text(730, 150, "spread between people", 11.5, bold=True, anchor="start")
    s.text(450, 250, "A higher NPX normally means more protein. A genetic variant can also change it WITHOUT changing the amount:", 12, bold=True)
    s.text(450, 270, "an amino-acid change under an antibody's footprint (epitope effect), or something in plasma that interferes with the pair.", 12)
    s.text(450, 296, "HLA alleles measured on HLA proteins (A*01:01 → HLA-A, z 47) are the textbook epitope case.", 11, muted=True)
    return s.done()


def fig_haplotype():
    s = S(900, 340)
    s.text(450, 22, "A haplotype travels together; a 'trans' effect acts on genes elsewhere in the genome", 14, bold=True)
    s.text(40, 52, "chromosome 6 (one copy)", 12, anchor="start", bold=True)
    s.rect(40, 62, 460, 24, "--mut", 0.12, rx=12)
    for x, lab, col in ((70, "HLA-C*08:02", "--b"), (200, "HLA-B*14:02", "--warn"), (350, "DRB1*01:02", "--mut")):
        s.rect(x, 58, 100, 32, col, 0.3)
        s.text(x + 50, 79, lab, 11.5, bold=True)
    s.parts.append('<path d="M70 98 v8 h230 v-8" fill="none" style="stroke:var(--warn);stroke-width:2"/>')
    s.text(40, 128, "inherited as a block:", 12, anchor="start", bold=True)
    s.text(40, 146, "B*14:02 and C*08:02 almost always come together", 11.5, anchor="start")
    s.text(40, 164, "(C*08:02 'looks' causal only because it tags B*14:02)", 11, muted=True, anchor="start")
    chroms = [("chr 4", "LRPAP1", "+2.6 SD", 560, "--a"), ("chr 12", "LRP1", "+1.7 SD", 680, "--b"), ("chr 17", "SGSH", "+1.4 SD", 800, "--c")]
    for lab, g, v, x, col in chroms:
        s.rect(x - 45, 210, 90, 22, "--mut", 0.12, rx=11)
        s.rect(x - 32, 206, 64, 30, col, 0.3)
        s.text(x, 227, g, 11.5, bold=True)
        s.line(300, 92, x, 202, "--warn", dash=True, arrow=True)
        s.text(x, 256, lab, 11, muted=True)
        s.text(x, 274, v, 12, bold=True)
    s.text(680, 292, "plasma level in B*14:02 carriers", 11, muted=True)
    s.text(450, 316, "Blood RNA of these three genes does not change, so carriers do not simply make more of the message.", 12)
    s.text(450, 334, "Open question: a real effect on how the proteins are handled outside cells, or an assay/interferent effect.", 12)
    return s.done()

BACKGROUND = {
    "resolution": block(
        "Background: what is an HLA allele name, and what does 'resolution' mean?",
        ["HLA genes (HLA-A, -B, -C, DRB1, …) are the most variable genes in humans. Each version is an <b>allele</b>, named with up to four colon-separated fields. The first two fields change the protein; the last two usually do not. At the usual 2-field level a typing tool has to be right about every amino-acid-changing difference.",
         "What the immune system cares about is the <b>peptide-binding groove</b>, built from exons 2 and 3. A <b>G group</b> is a set of alleles that are identical over exactly those exons. Exon-based typing tools (like the short-read callers in All of Us) can only be judged fairly at G-group level, which is why the 'error' we first saw (DQA1 0.52) vanished at G-group level (0.99)."],
        fig_resolution(), "Schematic: the four fields of an allele name (top) and the exons that define a G group (bottom)."),
    "reads": block(
        "Background: raw reads, assembly artefacts, stop-gains and nonsense-mediated decay",
        ["A sequencing run produces millions of <b>reads</b>; software stitches them into an <b>assembly</b> and calls alleles from it. If the assembly is wrong at one position, a 'novel allele' can appear that no read supports. So we go back to the raw reads: a real inherited variant is present in <i>every</i> read at that position in carriers and in none in non-carriers.",
         "A <b>stop-gain</b> puts an early STOP signal in the gene. The cell often destroys such messages (<b>nonsense-mediated decay</b>), so the mutated copy contributes fewer RNA reads than the normal copy. Seeing only ~26% of RNA reads from the stop-gain copy is therefore evidence that the allele is truly loss-of-function."],
        fig_reads(), "Schematic. A: read pile-up in a carrier and a control. B: why a premature stop shows up as a deficit of RNA reads."),
    "nkg2a": block(
        "Background: HLA-E, NKG2A and why 107G matters",
        ["NK cells kill infected or stressed cells but must spare healthy ones. They check for <b>HLA class I</b> on the target. <b>HLA-E</b> displays a fragment of the leader (signal) sequence of the other class I molecules; the NK receptor <b>CD94/NKG2A</b> (genes KLRD1/KLRC1) reads it and sends a brake signal: 'self is intact'.",
         "Almost everyone carries one of two HLA-E versions, differing at amino acid 107: <b>R107</b> (E*01:01) or <b>G107</b> (E*01:03, the more stable one at the cell surface). We asked whether the version you carry changes how much NKG2A your blood's NK cells make. The <b>TPM</b> numbers in this report are RNA counts per million, shown on a log2 scale (−0.09 ≈ 6% less)."],
        fig_nkg2a(), "Schematic of the HLA-E / NKG2A 'don't kill me' handshake and our finding."),
    "admix": block(
        "Background: what is a phenome-wide scan, and why do 'strong' hits so often turn out false?",
        ["A <b>phenome-wide association scan (PheWAS)</b> tests every allele against hundreds of diseases from health records. With ~100,000 tests, some hits are chance, so we compare against <b>negative-control outcomes</b> (diagnosis chapters that should not depend on HLA, such as injuries) and against <b>label permutations</b>. The <b>genomic inflation λ</b> summarises whether the whole test cloud is miscalibrated (1.0 is ideal).",
         "The biggest trap is <b>confounding</b>: within a broad group like 'AMR' people differ in their ancestral mix, which affects both which alleles they carry and (through environment and access to care) which diagnoses they receive. <b>Ancestry principal components (PCs)</b> are numbers that summarise this fine ancestry; adding them is the standard fix."],
        fig_confound(), "Schematic of ancestry confounding and the checks that expose it."),
    "cis": block(
        "Background: ancestry principal components and 'cis' controls",
        ["A genetic effect on a gene's activity is called <b>cis</b> when the controlling variant sits next to the gene (here, variants within a few hundred thousand letters of KLRC1) and <b>trans</b> when it acts from far away (HLA-E is on a different chromosome). Cis variants usually have big, easy-to-detect effects.",
         "A skeptic's question: maybe people who carry HLA-E 107G also tend to carry a particular KLRC1-region variant, and <i>that</i> lowers KLRC1. Two chromosomes are only associated through shared ancestry, so we test it directly: add ten ancestry PCs, then add the strongest cis variants, and see whether the HLA-E effect survives."],
        fig_cis(), "Schematic of the cis control for the HLA-E → KLRC1 result."),
    "pea": block(
        "Background: plasma proteomics (Olink), NPX and pQTLs",
        ["<b>Olink Explore HT</b> measures ~5,400 proteins in a drop of blood plasma using a <b>proximity extension assay</b>: two antibodies, each carrying a DNA tag, must bind the same protein before the tags can pair; the paired tags are then counted. The result is <b>NPX</b> on a log2 scale. A person with a gene variant that changes a protein's plasma level is said to have a <b>pQTL</b> (protein quantitative trait locus).",
         "Plasma holds the <i>soluble</i> fraction of a protein, not what sits on a cell's surface. And because the readout depends on antibodies, a variant can change the signal without changing the amount (an <b>epitope effect</b>). That is why our gene-to-protein controls matter: KIR3DL1 gene content predicts its own protein with z ≈ 51 and C4B copy number predicts C4B protein, so the pipeline sees real biology."],
        fig_pea(), "Schematic of the Olink readout and the two ways a genetic variant can change it."),
    "haplotype": block(
        "Background: haplotypes, linkage, trans effects, and what 'z' and 'SD' mean here",
        ["A <b>haplotype</b> is a stretch of one chromosome inherited as a block. In the MHC region alleles of HLA-B, HLA-C and DRB1 are tightly linked, so when one looks associated, its neighbours do too; <b>conditional analysis</b> (asking whether an allele still matters after accounting for another) separates the driver from its passengers.",
         "<b>SD</b> (standard deviation) is a person-to-person spread; '+2.6 SD' means carriers sit far above nearly everyone else. <b>z</b> is the effect divided by its uncertainty (|z| above ~6 survives the correction for testing thousands of proteins). The B*14:02 effect is <b>trans</b>: the allele is on chromosome 6, the three proteins are encoded on chromosomes 4, 12 and 17, and in blood the three genes are not more active, which is what makes it puzzling."],
        fig_haplotype(), "Schematic of the B*14:02 haplotype and its trans associations with three plasma proteins."),
}
