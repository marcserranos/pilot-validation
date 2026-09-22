#!/usr/bin/env python3
"""Build the call report: every figure embedded as a data URI, with the explanations in blocks."""
import base64, os, sys

REPO = "/Users/marcserrano/WORK/STANFORD/pilot-validation"
R = os.path.join(REPO, "reports", "hla_popgen")
OUT = sys.argv[1]


def img(rel, alt, maxw=None):
    p = os.path.join(R, rel)
    if not os.path.exists(p):
        return '<p class="missing">missing: %s</p>' % rel
    b = base64.b64encode(open(p, "rb").read()).decode()
    style = ' style="max-width:%s"' % maxw if maxw else ""
    return ('<figure><img src="data:image/png;base64,%s" alt="%s"%s>'
            '<figcaption>%s<span class="path">%s</span></figcaption></figure>'
            % (b, alt, style, alt, rel))


HEAD = """<title>HLA Call Pack</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{--ground:#F6F7F8;--surface:#FFF;--ink:#15181C;--ink-soft:#4A525C;--rule:#D9DEE3;
  --rule-soft:#E8ECEF;--accent:#0F6E6E;--accent-soft:#E2F0EF;--warn:#9C5A00;--warn-soft:#FBEEDC;
  --good:#1B6B45;--good-soft:#E3F1E9;--mono-bg:#EFF2F4;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#121417;--surface:#191D21;
  --ink:#E8ECEF;--ink-soft:#A3ADB7;--rule:#2C333A;--rule-soft:#242A30;--accent:#4FC0BC;
  --accent-soft:#12302F;--warn:#E0A75C;--warn-soft:#2E2417;--good:#69C296;--good-soft:#142A20;
  --mono-bg:#20262B;}}
:root[data-theme="dark"]{--ground:#121417;--surface:#191D21;--ink:#E8ECEF;--ink-soft:#A3ADB7;
  --rule:#2C333A;--rule-soft:#242A30;--accent:#4FC0BC;--accent-soft:#12302F;--warn:#E0A75C;
  --warn-soft:#2E2417;--good:#69C296;--good-soft:#142A20;--mono-bg:#20262B;}
body{background:var(--ground);color:var(--ink);font-family:"IBM Plex Sans",-apple-system,sans-serif;
  font-size:16px;line-height:1.6;}
.wrap{max-width:82ch;margin:0 auto;padding-inline:20px;padding-block:44px 80px;}
h1,h2,h3{font-family:"IBM Plex Serif",Georgia,serif;text-wrap:balance;line-height:1.24;}
h1{font-size:2.1rem;font-weight:600;margin:0 0 .3em;}
h2{font-size:1.45rem;font-weight:600;margin:2.8em 0 .5em;padding-top:.8em;border-top:2px solid var(--accent);}
h3{font-size:1.06rem;font-weight:600;margin:2em 0 .35em;}
p{margin:0 0 1em;max-width:72ch;} ul,ol{margin:0 0 1.1em;padding-left:1.2em;max-width:72ch;}
li{margin-bottom:.4em;} strong{font-weight:600;}
code{font-family:"IBM Plex Mono",monospace;font-size:.85em;background:var(--mono-bg);
  padding:.1em .35em;border-radius:3px;}
.kicker{font-size:.72rem;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);
  font-weight:600;margin:0 0 .7em;}
.standfirst{font-size:1.05rem;color:var(--ink-soft);max-width:66ch;margin-bottom:1.4em;}
.meta{font-size:.8rem;color:var(--ink-soft);border-top:1px solid var(--rule);padding-top:.8em;
  margin-bottom:2em;font-family:"IBM Plex Mono",monospace;}
.block-tag{display:inline-block;font-family:"IBM Plex Mono",monospace;font-size:.7rem;
  font-weight:500;color:var(--accent);background:var(--accent-soft);border-radius:3px;
  padding:.12em .5em;vertical-align:.2em;margin-left:.5em;}
figure{margin:1.5em 0;background:var(--surface);border:1px solid var(--rule);border-radius:5px;
  padding:12px;}
figure img{display:block;width:100%;height:auto;border-radius:3px;background:#fff;}
figcaption{font-size:.83rem;color:var(--ink-soft);margin-top:.7em;line-height:1.45;}
figcaption .path{display:block;font-family:"IBM Plex Mono",monospace;font-size:.72rem;
  color:var(--ink-soft);opacity:.75;margin-top:.3em;word-break:break-all;}
.tablewrap{overflow-x:auto;margin:0 0 1.3em;border:1px solid var(--rule);border-radius:4px;
  background:var(--surface);}
table{border-collapse:collapse;width:100%;font-size:.87rem;}
th,td{text-align:left;padding:.48em .75em;border-bottom:1px solid var(--rule-soft);white-space:nowrap;}
th{font-size:.7rem;letter-spacing:.06em;text-transform:uppercase;color:var(--ink-soft);
  font-weight:600;background:var(--mono-bg);}
tbody tr:last-child td{border-bottom:none;}
td.num{font-family:"IBM Plex Mono",monospace;font-variant-numeric:tabular-nums;}
td.w{white-space:normal;min-width:20ch;}
.flag{padding:.85em 1.05em;border-radius:4px;margin:0 0 1.3em;font-size:.93rem;}
.flag.warn{background:var(--warn-soft);border:1px solid var(--warn);}
.flag.good{background:var(--good-soft);border:1px solid var(--good);}
.flag.say{background:var(--accent-soft);border:1px solid var(--accent);}
.flag p:last-child{margin-bottom:0;}
.flag .label{font-size:.68rem;letter-spacing:.11em;text-transform:uppercase;font-weight:600;
  display:block;margin-bottom:.3em;}
.flag.warn .label{color:var(--warn);} .flag.good .label{color:var(--good);}
.flag.say .label{color:var(--accent);}
.eq{font-family:"IBM Plex Mono",monospace;background:var(--mono-bg);padding:.7em .9em;
  border-radius:4px;margin:0 0 1em;font-size:.9rem;overflow-x:auto;}
footer{margin-top:3em;padding-top:1.1em;border-top:1px solid var(--rule);font-size:.8rem;
  color:var(--ink-soft);}
@media(max-width:480px){h1{font-size:1.65rem;}.wrap{padding-block:30px 52px;}}
</style>
"""


def build():
    h = [HEAD, '<div class="wrap">']
    A = h.append

    A('<p class="kicker">Omni-HLA · supervisor call pack</p>')
    A("<h1>Everything for this call, in blocks</h1>")
    A('<p class="standfirst">Every figure, what it shows, and the sentence to say when it is on '
      'screen. Ordered the way the call runs. Each figure names its file so you can open the '
      'original at full resolution.</p>')
    A('<p class="meta">11,856 unrelated participants · 22,000+ assembled haplotypes · '
      'scripts 24–35 · branch <code>fig1-drafts-and-research-map</code></p>')

    # ---------------------------------------------------------------- Figure 1
    A('<h2>Figure 1<span class="block-tag">ask A1</span></h2>')
    A("<p>The four panels Cole specified, assembled into one sheet. Panels (a) and (b) are the "
      "figures he already saw and approved; (c) and (d) are rebuilt on the corrected novelty "
      "definition.</p>")
    A(img("35_figure1/figure1.png",
          "Figure 1 — (a) cohort ancestry, (b) HLA-B allele ancestry simplex, "
          "(c) novelty by gene/ancestry/field with artifact control, (d) allele-frequency spectrum"))
    A('<div class="flag warn"><span class="label">Say this before they ask</span>'
      "<p>This is a <strong>layout</strong>, not a re-analysis — panels are composited from the "
      "rendered figures. Panels (a) and (b) are still at the <em>original</em> admixture "
      "threshold, not the stricter ~98% you asked for: scripts 06 and 10 have no such parameter "
      "today, so that is a small code change, not just a rerun. And panel (d) predates the "
      "corrected novelty definition — its shape is right, its absolute novel counts are "
      "over-estimates.</p></div>")

    # ---------------------------------------------------------------- Block A
    A('<h2>Block A · Novel allele discovery<span class="block-tag">asks A1c, A6</span></h2>')
    A("<h3>The chain of numbers — get these exactly right</h3>")
    A('<div class="tablewrap"><table><thead><tr><th>Step</th><th>Number</th><th>What it is</th>'
      "</tr></thead><tbody>"
      '<tr><td class="w">Calls flagged <code>new</code></td><td class="num">280,695</td>'
      '<td class="w">haplotype×gene <em>calls</em>, not alleles</td></tr>'
      '<tr><td class="w">Old "novel allele" clusters</td><td class="num">~3,000</td>'
      '<td class="w">clustered on CDS only — <strong>the artifact</strong></td></tr>'
      '<tr><td class="w">Distinct novel proteins</td><td class="num">1,404</td>'
      '<td class="w">re-clustered at protein level, sequence truth-checked</td></tr>'
      '<tr><td class="w">Seen in ≥2 unrelated people</td><td class="num">231</td>'
      '<td class="w">the recurrent set</td></tr>'
      '<tr><td class="w">Seen in ≥20 unrelated people</td><td class="num">29 proteins</td>'
      '<td class="w">+ 9 synonymous = <strong>38</strong> in the callout list</td></tr>'
      "</tbody></table></div>")
    A('<div class="flag warn"><span class="label">Two things not to garble</span>'
      "<p>It is <strong>1,404</strong>, not 1,000. And the 38 is <strong>29 proteins + 9 "
      "synonymous CDS</strong> — do not let it get called “38 new proteins”.</p></div>")

    A("<h3>The funnel: how 280,695 becomes 231</h3>")
    A(img("24_novelty_by_field/fig_novelty_funnel.png",
          "Every filtering step from flagged calls down to recurrent novel proteins"))
    A("<p>Two steps do most of the work: ~25,000 “novel” calls turn out to have a coding sequence "
      "that is <em>already catalogued</em> under another name, and the bulk of what is left is "
      "frameshift/stop calls dominated by homopolymer artifacts — the classic long-read error.</p>")

    A("<h3>Novel proteins by recurrence and gene — the figure you asked for</h3>")
    A(img("34_novel_recurrence/fig_novel_recurrence.png",
          "(a) novel protein alleles per gene, stacked by how many unrelated people carry them; "
          "(b) general vs strict — everything found, against what clears the reporting threshold"))
    A('<div class="flag say"><span class="label">The sentence</span>'
      "<p>“TAP2 has 201 novel proteins and 10 that clear the reporting threshold. "
      "<strong>HLA-B has 60 novel proteins, 1 recurrent, and 0 reportable.</strong> "
      "The classical genes are essentially catalogued; the discovery is everywhere else.”</p></div>")
    A('<div class="flag warn"><span class="label">A constraint to state, not improvise around</span>'
      "<p>All of Us suppresses every count from 1 to 19, so bars at 1, 2, 3 … 19 "
      "<strong>cannot</strong> be drawn from committed data — those numbers exist only on the VM. "
      "The three classes shown (seen once / 2–19 / ≥20) are the same information at the "
      "resolution the rule permits. A true per-count histogram is one VM run away.</p></div>")

    A("<h3>The callout list — what can be named in main text</h3>")
    A(img("32_novel_callouts/fig_novel_callouts.png",
          "The 38 novel alleles carried by ≥20 unrelated people, with their ancestry composition"))
    A("<p>The largest is a <strong>TAP1</strong> protein one residue from TAP1*01:01, carried by "
      "<strong>426</strong> unrelated people — about <strong>139 per 1,000</strong> "
      "African-ancestry participants. 15 of the 29 proteins are most common in African-ancestry "
      "participants.</p>")
    A('<div class="flag good"><span class="label">Pre-empt the obvious objection</span>'
      "<p>“TAP1 tops the list only because nobody catalogued TAP1.” That is measurable, and it was "
      "measured independently: the share of a gene’s haplotypes whose protein is <em>already</em> "
      "absent from IPD-IMGT is <strong>TAP1 8.0%, TAP2 6.9%, HLA-A/B/C 0.1–0.2%</strong>. The "
      "genes with the most novel alleles are the genes with the biggest catalogue gap — two "
      "independent estimates agreeing, not one confounding the other.</p></div>")

    # ---------------------------------------------------------------- Block B
    A('<h2>Block B · Quality check, twins and phasing<span class="block-tag">ask A2</span></h2>')
    A("<h3>How relatedness is computed — and who computes it</h3>")
    A("<p>We do <strong>not</strong> compute it. All of Us ships a genome-wide table with a "
      "kinship coefficient φ per pair, estimated from hundreds of thousands of SNPs using the "
      "<strong>KING</strong> estimator:</p>")
    A('<div class="eq">φ = [ N(Aa,Aa) − 2·N(AA,aa) ] / [ N(Aa)<sub>i</sub> + N(Aa)<sub>j</sub> ]</div>')
    A("<p>N(Aa,Aa) counts sites where both people are heterozygous; N(AA,aa) counts opposite "
      "homozygotes. Opposite homozygotes are the signal — close relatives almost never have them. "
      "The thresholds are the standard KING cut-points, each a factor of √2 apart:</p>")
    A('<div class="tablewrap"><table><thead><tr><th>φ</th><th>Relationship</th></tr></thead><tbody>'
      '<tr><td class="num">&gt; 0.354</td><td>duplicate or monozygotic twin</td></tr>'
      '<tr><td class="num">0.177 – 0.354</td><td>first degree (parent–child, full sibling)</td></tr>'
      '<tr><td class="num">0.0884 – 0.177</td><td>second degree</td></tr>'
      '<tr><td class="num">0.0442 – 0.0884</td><td>third degree</td></tr>'
      "</tbody></table></div>")
    A('<div class="flag warn"><span class="label">How sure are we they are twins?</span>'
      "<p>φ above 0.354 means essentially <em>identical genomes genome-wide</em>. What it cannot "
      "distinguish is an identical twin pair from <strong>the same person sequenced twice</strong> "
      "— which is why the report says “duplicate or MZ twin”, never just “twins”. For our purpose "
      "the distinction does not matter: either way they are a <strong>technical replicate</strong>, "
      "two independent sequencing and assembly runs of the same genome, so any HLA difference "
      "between them is <em>our</em> error, not biology.</p>"
      "<p>That count is below 20, so it stays suppressed in anything published — do not say the "
      "number out loud in a public document.</p></div>")
    A("<p>They differ at <strong>121 bases out of 306,033</strong> compared → <strong>Q34, about "
      "one error per 2,500 bases</strong>. Note that is <em>worse</em> than the Q46 the indirect "
      "method implied. Lead with that honesty.</p>")

    A("<h3>Where the 5% disagreement actually goes</h3>")
    A(img("26_qc_relatives_v2/decomposition_by_gene.png",
          "Every gene comparison between relatives, sorted into exactly one explanation"))
    A("<p><strong>87.0%</strong> byte-identical · 6.5% gene called in only one of the pair · "
      "3.2% dropout candidates · 2.3% single-base error · <strong>0.8% genuinely "
      "unexplained</strong>. The old “95% concordance” was one number hiding five different "
      "phenomena.</p>")

    A("<h3>The switch test, rebuilt so it could fail</h3>")
    A(img("26_qc_relatives_v2/switch_power.png",
          "Detection power for phase switches, measured by injecting synthetic switches"))
    A("<p>The previous version excluded disagreeing genes <em>before</em> counting switches and "
      "never reported a denominator — “0 switches in 545/545 pairs” could not distinguish “no "
      "switches” from “no test”. The rebuilt version reports <strong>3,021 testable "
      "transitions</strong>, detects injected synthetic switches <strong>100%</strong> of the "
      "time, and still finds <strong>zero</strong> real ones.</p>")

    A("<h3>Two supporting signals</h3>")
    A(img("26_qc_relatives_v2/sharing_histogram.png",
          "HLA sharing across all kinship-selected relative pairs — the pair universe was chosen "
          "on genome-wide kinship alone, never on HLA"))
    A(img("26_qc_relatives_v2/homozygosity_obs_exp.png",
          "Observed vs expected homozygosity — a dropout detector independent of relatives"))

    # ---------------------------------------------------------------- contigs
    A('<h2>Block B2 · What “same contig” actually means<span class="block-tag">your question</span></h2>')
    A("<p>Your instinct was half right. We are <strong>not</strong> looking at raw reads, and we "
      "are not looking at collapsed genotypes either. Each person’s PacBio HiFi reads were "
      "assembled <em>de novo</em> into two haplotype assemblies (hap1 / hap2), trimmed to the "
      "chr6 MHC region. Immuannot annotates genes onto those assembly <strong>contigs</strong>.</p>")
    A("<ul>"
      "<li><strong>Files:</strong> <code>~/pipeline_outputs/people/&lt;person_id&gt;/"
      "immuannot_output/*.gtf.gz</code>, parsed into <code>hla_calls_rich.tsv</code>.</li>"
      "<li><strong>The key column:</strong> GTF column 1 is the contig ID. Every earlier parser in "
      "this project threw it away — capturing it is the only reason this analysis exists.</li>"
      "<li><strong>What it proves:</strong> two genes on one continuous assembled sequence are "
      "physically on the same chromosome copy. No statistical phasing, no reference panel, no "
      "population inference.</li></ul>")
    A(img("29_hla_ld/fig_phasing_yield.png",
          "Physical phasing yield per interval — the fraction of assemblies carrying both genes "
          "that had them on a single contig"))
    A('<div class="tablewrap"><table><thead><tr><th>Pair</th><th>Assemblies with both genes</th>'
      "<th>On one contig</th><th>%</th></tr></thead><tbody>"
      '<tr><td>DPA1~DPB1</td><td class="num">22,601</td><td class="num">21,913</td><td class="num"><strong>96.96</strong></td></tr>'
      '<tr><td>DQA1~DQB1</td><td class="num">22,481</td><td class="num">21,652</td><td class="num"><strong>96.31</strong></td></tr>'
      '<tr><td>DRB1~DQB1</td><td class="num">22,055</td><td class="num">20,229</td><td class="num">91.72</td></tr>'
      '<tr><td>B~C</td><td class="num">22,416</td><td class="num">19,474</td><td class="num">86.88</td></tr>'
      '<tr><td>A~B</td><td class="num">22,187</td><td class="num">6,709</td><td class="num"><strong>30.24</strong></td></tr>'
      "</tbody></table></div>")
    A("<p>Denominator is assemblies: ~11,856 unrelated people × 2 haplotypes. <strong>The "
      "gradient is the answer to “does everyone have this?”</strong> — DPA1 and DPB1 sit ~10 kb "
      "apart so nearly every assembly spans them; HLA-A and HLA-B are 1.4 Mb apart and only 30% "
      "of assemblies reach that far. When a pair is not on one contig we <strong>drop it</strong> "
      "rather than guess.</p>")

    # ---------------------------------------------------------------- Block C LD
    A('<h2>Block C · Linkage disequilibrium<span class="block-tag">ask A3</span></h2>')
    A("<h3>What LD is, in one paragraph</h3>")
    A("<p>Two genes sit near each other on a chromosome. If their alleles were inherited "
      "independently, knowing which allele you carry at gene A would tell you nothing about which "
      "allele sits beside it at gene B <em>on the same chromosome copy</em>. LD is the departure "
      "from that independence — alleles that travel together more often than chance.</p>")
    A('<div class="eq">D   = p(AB) − p(A)·p(B)          the raw excess\n'
      "r²  = D² / [ p(A)(1−p(A))·p(B)(1−p(B)) ]   the squared correlation — what Cole asked for\n"
      "D′  = D / D(max)                          D rescaled by the most it could be</div>")
    A('<div class="flag warn"><span class="label">Why I did not just hand him r²</span>'
      "<p>r² is <strong>bounded by the allele frequencies</strong>. An allele at 2% frequency "
      "physically cannot reach a high r², even in perfect linkage. So an ancestry with rarer "
      "alleles shows lower r² for reasons that have nothing to do with linkage. On top of that, "
      "HLA genes are multi-allelic (dozens of alleles each) while r² is defined for two-allele "
      "markers. So we also compute <strong>Hedrick’s multi-allelic D′</strong>, and compare "
      "ancestries only at a <strong>common rarefied sample size</strong>.</p></div>")
    A("<h3>Why they care</h3>")
    A("<p>DQA1 and DQB1 encode the two halves of one molecule — the α and β chain of the same "
      "heterodimer. Which α pairs with which β determines what peptides that person can present. "
      "<strong>The haplotype is the functional unit, not the individual allele.</strong> Cole’s "
      "hypothesis was that these pairings differ between populations beyond just differing "
      "frequencies.</p>")
    A(img("29_hla_ld/fig_ld_multiallelic.png",
          "Multi-allelic D′ per ancestry, each gene pair rarefied to a common haplotype count"))
    A('<div class="flag good"><span class="label">Read the controls first</span>'
      "<p>HLA-B~HLA-C, 90 kb apart and known to be tightly linked → D′ <strong>0.87–0.92</strong>. "
      "HLA-A~HLA-B, 1.4 Mb apart → lowest at <strong>0.58–0.74</strong>. The method finds strong "
      "linkage where it is known to exist and weak linkage where it is not. Only then read the "
      "result.</p></div>")
    A('<div class="tablewrap"><table><thead><tr><th>Pair</th><th>African ancestry</th>'
      "<th>European ancestry</th></tr></thead><tbody>"
      '<tr><td>DQA1~DQB1</td><td class="num">0.908 [0.890–0.925]</td><td class="num">0.974 [0.964–0.983]</td></tr>'
      '<tr><td>DRB1~DQB1</td><td class="num">0.880 [0.860–0.899]</td><td class="num">0.947 [0.932–0.960]</td></tr>'
      "</tbody></table></div>")
    A("<p>Non-overlapping intervals at equal sample size: African-ancestry haplotypes carry "
      "<strong>weaker class II linkage</strong> — more distinct DQ and DR–DQ haplotype "
      "combinations, consistent with older effective population size and more accumulated "
      "recombination.</p>")
    A("<p><strong>DPA1~DPB1 does not follow the pattern.</strong> Worth saying out loud: Brandt "
      "et al. 2018 independently singled out the DP locus as the exception among HLA genes, "
      "arguing it is under directional rather than balancing selection. Our data behave the same "
      "way.</p>")
    A(img("29_hla_ld/fig_r2_DQA1_DQB1.png",
          "The literal r² Cole asked for: allele-by-allele heatmaps, one panel per ancestry"))
    A('<div class="flag warn"><span class="label">Scope it</span>'
      "<p>AFR-vs-EUR is established. “African ancestry is lowest of all six” is <strong>not</strong> "
      "— the African and East Asian intervals touch for DRB1~DQB1.</p></div>")

    # ---------------------------------------------------------------- Block D SV
    A('<h2>Block D · Deletions, duplications and KIR<span class="block-tag">asks A4, A5</span></h2>')
    A("<p>A gene missing from an assembly could be (a) really deleted, (b) the assembly broke "
      "there, or (c) detection failed. So a deletion is called <strong>only</strong> when one "
      "contig carries a gene on <em>each side</em> of the missing gene — a “bridged absence”, "
      "meaning the sequence was assembled straight through that position and the gene was not "
      "there.</p>")
    A(img("30_hla_sv/fig_deletion_rates.png",
          "Deletion rate per gene, ordered by physical position. Green = known copy-number genes "
          "(positive control); grey = genes that should never be deleted (false-positive rate)"))
    A('<div class="flag good"><span class="label">The positive control is the thing to show</span>'
      "<p>Which second DRB locus you carry is determined by your DRB1 type — textbook "
      "immunogenetics, known independently of our data. DR52 (DRB1*03/11/12/13/14) → DRB3; DR53 "
      "(*04/07/09) → DRB4; DR51 (*15/16) → DRB5; DR1/DR8/DR10 → none. "
      "<strong>We recover it in 92–98% of haplotypes across all 13 DRB1 groups</strong>, including "
      "the groups whose correct answer is “no second DRB locus at all”. The negative control gives "
      "the false-positive rate: <strong>0.04–2.0%</strong>.</p></div>")
    A('<div class="tablewrap"><table><thead><tr><th>Gene</th><th>% haplotypes lacking it</th>'
      "<th>Reading</th></tr></thead><tbody>"
      '<tr><td>DRB5 / DRB4 / DRB3</td><td class="num">83.1 / 69.6 / 49.3</td><td class="w">known DR haplotype-group structure</td></tr>'
      '<tr><td>C4B / C4A</td><td class="num">19.6 / 11.0</td><td class="w">real copy-number variation</td></tr>'
      '<tr><td>MICA / MICB</td><td class="num">4.2 / 2.9</td><td class="w">known deletion polymorphism</td></tr>'
      "</tbody></table></div>")
    A(img("30_hla_sv/fig_c4_copy_number.png",
          "C4 copy number per haplotype — 54.1% carry one C4A and one C4B, 14.8% A-only, 10.5% B-only"))
    A("<p>This is phased copy number on <em>individual haplotypes</em> — the thing short reads do "
      "badly and long reads do well.</p>")
    A('<div class="flag warn"><span class="label">A live lead, offered as a lead</span>'
      "<p>Four class I pseudogenes — H (12.15%), K (12.42%), T (12.25%), U (12.37%) — have "
      "near-identical deletion rates while every other pseudogene sits under 1.4%. Four "
      "independent genes agreeing to within 0.3 points is not what independent deletions look "
      "like. Either they sit in one segment deleted as a block on ~12% of haplotypes — a real "
      "finding — or they share a detection artifact. Deciding needs one short follow-up.</p></div>")
    A('<div class="flag say"><span class="label">KIR — the answer he did not expect</span>'
      "<p><strong>Zero KIR calls, and that is correct.</strong> The KIR cluster is on chromosome "
      "19, outside the chromosome 6 window these assemblies were trimmed to. Cole expected it to "
      "“just be called” — it cannot be. Typing KIR needs a separate extraction from the original "
      "BAMs. That is an opportunity, not an oversight.</p></div>")

    # ---------------------------------------------------------------- Block E groove
    A('<h2>Block E · Diversity in the peptide groove<span class="block-tag">asks A8, A10</span></h2>')
    A("<h3>How we decide what is inside the groove</h3>")
    A("<p>Two definitions, deliberately — one coarse and safe, one sharp and structural.</p>")
    A("<ol><li><strong>By exon (what the headline numbers use).</strong> IMGT’s own annotation "
      "says which bases belong to which exon. Class I exons 2+3 encode the α1/α2 domains that "
      "form the groove; class II exon 2 encodes β1/α1. This comes from data we ship and "
      "<em>cannot</em> be off by a signal-peptide length.</li>"
      "<li><strong>By structure (the sharper one).</strong> The published ARS residue tables could "
      "not be verified from primary sources — Parham 1988 returned 403, Bondinas 2007 is "
      "paywalled. So we downloaded six crystal structures with peptides bound (1HHK, 1A1M, 4NT6, "
      "1DLH, 1JK8, 3LQZ) and computed <strong>every MHC residue with an atom within 4.5 Å of the "
      "bound peptide</strong> — 206 residues, mapped onto our numbering by alignment, so no "
      "leader length is ever assumed.</li></ol>")
    A(img("31_aa_diversity/fig_protein_track_DRB1.png",
          "Amino-acid diversity per residue along HLA-DRB1; the shaded band is the groove-encoding "
          "exon 2"))
    A('<div class="flag good"><span class="label">The validation to lead with</span>'
      "<p>Ranking residues purely by diversity — with nothing told about disease, about KIR, or "
      "about any published paper — the top positions are <strong>DRB1 β11/13/71/74</strong> (the "
      "rheumatoid-arthritis shared epitope), <strong>DQB1 β57</strong> (type 1 diabetes, celiac), "
      "and <strong>HLA-B 77/80</strong> (the Bw4/Bw6 epitope that determines KIR binding) plus 116 "
      "(a principal peptide anchor).</p>"
      "<p>Be precise: it is <strong>8 of 10, 6 of 10 and 8 of 10</strong> of each gene’s top ten "
      "landing on canonical functional positions — a strong result, not a clean sweep.</p></div>")
    A(img("31_aa_diversity/fig_groove_enrichment.png",
          "Groove vs non-groove diversity ratio per gene; grey bars are the conserved control genes"))
    A("<p>Groove exons are <strong>2.2×–8.0×</strong> more diverse in eight genes, all "
      "p ≤ 0.0015 by position permutation. <strong>The controls show nothing</strong> — DRA 0.0×, "
      "HLA-F 0.04× — which is what makes the rest credible. HLA-C (1.3×) and DPA1 (2.0×) are not "
      "significant and are reported that way.</p>")
    A('<div class="flag warn"><span class="label">One control misbehaves — disclose it</span>'
      "<p>HLA-E, which is also a designated conserved control, shows the <em>largest</em> raw "
      "ratio in the table (20×). It is not significant (p = 0.39) because HLA-E’s mean diversity "
      "is 0.0015, so the ratio divides two near-zero quantities. But a reader is entitled to see "
      "it and to ask whether HLA-E — which presents a very restricted peptide repertoire — is the "
      "right negative control at all.</p></div>")
    A(img("31_aa_diversity/fig_allele_differentiation.png",
          "Between-ancestry differentiation per gene, on Hedrick's standardised G′st"))
    A('<div class="flag say"><span class="label">Where the data disagree with Cole</span>'
      "<p>He expected class II to be the most differentiated between ancestries. On Hedrick’s "
      "standardised G′st the leader is <strong>HLA-B (0.518)</strong>, a class I gene, then DPB1 "
      "(0.483), DRB1 (0.444), HLA-A (0.417). Class I and class II interleave.</p>"
      "<p>There is a methods point worth making here: <strong>raw F<sub>st</sub> would have given "
      "a completely different ranking</strong> — DRB5, DRB4 and DPA1 on top — purely because those "
      "genes have lower within-population heterozygosity, which mechanically permits a larger "
      "F<sub>st</sub>. At HLA-B, within-population heterozygosity is 0.954, so raw F<sub>st</sub> "
      "cannot exceed about 0.05 however different the populations are. That is exactly the "
      "artifact Brandt et al. 2018 documented — and we can demonstrate it in our own data rather "
      "than cite it.</p></div>")

    # ---------------------------------------------------------------- asks
    A("<h2>What to ask them for</h2>")
    A("<ol>"
      "<li><strong>Reframe the first section around the non-classical MHC?</strong> The evidence "
      "says yes: classical genes essentially catalogued, discovery concentrated in TAP/MIC/DM/DO.</li>"
      "<li><strong>Headline: 231 recurrent novel proteins, or 1,404 distinct?</strong> "
      "Recommendation: 231.</li>"
      "<li><strong>Restart the short-read validation.</strong> He parked it — “I wouldn’t do that "
      "right now” — but the 38-allele list now exists and every one of those people has "
      "high-coverage short reads. This is the orthogonal confirmation the whole novelty claim "
      "currently lacks.</li>"
      "<li><strong>Is HLA-E the right negative control</strong> for groove diversity?</li>"
      "<li><strong>Scope KIR.</strong> It needs a separate extraction from the original BAMs — is "
      "that worth a workstream?</li>"
      "</ol>")

    A("<h2>What is still open</h2>")
    A('<div class="tablewrap"><table><thead><tr><th>Item</th><th>Needs VM?</th><th>Size</th></tr>'
      "</thead><tbody>"
      '<tr><td class="w">Rerun scripts 06 + 10 at strict 98% admixture, commit the frequency tables</td><td>yes + small code change</td><td>S</td></tr>'
      '<tr><td class="w">Regenerate the spectrum (panel d) at the corrected novelty definition</td><td>yes</td><td>S</td></tr>'
      '<tr><td class="w">Per-count histogram of novel-allele recurrence (1, 2, 3 … 19)</td><td>yes</td><td>S</td></tr>'
      '<tr><td class="w">H/K/T/U co-deletion: same haplotypes or not?</td><td>yes</td><td>S</td></tr>'
      '<tr><td class="w">Contact-vs-non-contact diversity inside the groove (function is written and verified)</td><td>yes</td><td>S</td></tr>'
      '<tr><td class="w">Supplement figure dump</td><td>yes</td><td>M</td></tr>'
      '<tr><td class="w">Short-read validation of the 38 callouts</td><td>yes</td><td>L</td></tr>'
      '<tr><td class="w">HLA × TCR/BCR join with Aleix</td><td>yes</td><td>XL</td></tr>'
      "</tbody></table></div>")

    A("<footer>Every figure above names its file under <code>reports/hla_popgen/</code>. "
      "Method detail is in each folder’s <code>README.md</code>; the sprint boards and journals "
      "are under <code>sprints/</code>. Branch <code>fig1-drafts-and-research-map</code> is "
      "committed locally and not yet pushed.</footer>")
    A("</div>")
    return "\n".join(h)


open(OUT, "w").write(build())
print("wrote", OUT, os.path.getsize(OUT) // 1024, "KB")
