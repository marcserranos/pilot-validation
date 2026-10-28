"""117 S08 full report: are the 'novel' HLA alleles real? (HTML, inline SVG; overview first, then evidence / method / caveat drop-downs.)

  python3 scripts/hla_popgen/117_s08_full_report.py --out reports/hla_popgen/91_novel_recurrence/S08_report.html [--review notes.html]

Numbers come from sprints/S08_case_browser/LOG.md and the aggregate arrays kept in 94_s08_report.py (VM runs of scripts 90-96; counts of 1-19 are never shown).
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
R94 = _load("r94", "94_s08_report.py")
fig, ev, badge, cards, table, section = F.fig, F.ev, F.badge, F.cards, F.table, F.section
ANC = R94.ANC

# ===================================================================== data (all from LOG.md / 94_s08_report.py)
FUNNEL = [("Novel HLA calls vs IPD 3.55", 8015, "a"), ("…novel-protein calls", 7635, "a"), ("…substitution-only and rebuildable", 675, "b"), ("…distinct signatures", 441, "b"),
          ("…seen in two or more unrelated people", 49, "c"), ("…recurrent, unnamed and carried forward here", 2, "warn", "both HLA-F")]
SHORTLIST = [("Shortlist (protein-altering, substitution-only)", 337, "a"), ("…not an exact IPD 3.66 match", 279, "b"), ("…and no 3.66 partial-allele hit", 270, "b"), ("…plus tier-A short-read support, ≤15 edits, plain reference", 98, "c")]
EDIT_TYPES = [("Novel-protein calls (7,635)", [86, 8, 6]), ("Novel-synonymous (227)", [0, 0, 100])]
SR_SUPPORT = [("Indels", [4, 79, 8, 9]), ("Protein-altering substitutions", [60, 18, 0, 22]), ("Synonymous substitutions", [68, 0, 0, 32])]
# 92: indel-bearing and substitution-only novel-protein calls per person, classical genes, by platform
REVIO_INDEL, SEQ2E_INDEL, REVIO_SNV, SEQ2E_SNV = R94.REVIO_INDEL, R94.SEQ2E_INDEL, R94.REVIO_SNV, R94.SEQ2E_SNV
ALL_SNV, ALL_NOVEL, PERSONS = R94.ALL_SNV, R94.ALL_NOVEL, R94.PERSONS
HAP_NAMES = ["F T237I (unnamed)", "F P205T (unnamed)", "F Y106C (= F*01:25, known)", "G synonymous (= G*01:01:33, known)"]
HAP_CARRIER, HAP_CONTROL = R94.HAP_CARRIER, R94.HAP_CONTROL
TD_LOW, TD_MID, TD_HIGH = R94.TD_LOW, R94.TD_MID, R94.TD_HIGH
ASE = [("HLA-F T237I (unnamed)", 0.54), ("HLA-F P205T (unnamed)", 0.65), ("HLA-F Y106C = F*01:25 (known, control)", 0.58)]
SIGS = [("HLA-G synonymous (on G*01:01:03:01)", "88 (EAS 82, others <20)", "Named", "= G*01:01:33"), ("HLA-F T237I (on F*01:01:01:01)", "27 (AFR-predominant, some AMR)", "UNNAMED", "no IPD 3.66 allele has this CDS"),
        ("HLA-F P205T (on F*01:01:02:01)", "<20 (AFR-predominant)", "UNNAMED", "no IPD 3.66 allele has this CDS"), ("HLA-F Y106C", "<20 (AMR only)", "Named", "= F*01:25"),
        ("DQB1 A57D (Asp57)", "<20 (AFR + AMR)", "Named (literature)", "= DQB1*02:180:02 (T1D-relevant residue); from a literature check, not yet re-derived locally"), ("HLA-E K195R", "<20 (AFR + AMR)", "Named", "= E*01:152"),
        ("HLA-G D130H / G107S / K170E", "<20 each", "Named", "= G*01:58 / 01:64 / 01:66"), ("DRB1 R262T; many-edit B/C substitutions", "<20", "Artefact", "3 different references; truncated or suffix-flagged templates"),
        ("DPB1 T16M", "<20 (EAS)", "Named", "CpG transition; base already in older IPD")]

LEDGER = [
    ("S08-1", "Low-distance 'novel' calls are mostly consensus indel errors", "Strong", "86% of novel-protein calls indel-only; 4% of indels supported by short reads; ~2× more on Sequel II/IIe than Revio", "a2"),
    ("S08-2", "A sizeable minority of 'novelty' is database age (IPD 3.55 vs 3.66)", "Strong", "62/441 (14%) substitution signatures and 58/337 (17%) shortlist calls are already named in 3.66; 16 of the 49 recurrent signatures (33%). The larger effect is indel error (S08-1)", "a3"),
    ("S08-3", "Recurrent signatures are inherited, haplotype-linked variants of under-represented ancestries", "Supported", "each on one reference allele, ancestry-structured. Rules out random noise, not a systematic error tied to one haplotype or template", "a4"),
    ("S08-4", "HLA-F T237I and P205T are probably real alleles (inheritance not formally tested); unnamed in IPD 3.66", "Supported", "one flanking haplotype, expected platform/centre counts, relatives agree (counts withheld), variant reads in RNA. A haplotype-specific assembly or template error would give the same single-background pattern; raw reads (which S09 shows in every examined carrier) and an independent trio are what settle it; the trio is still missing", "a5"),
    ("S08-5", "…and they are transcribed", "Supported", "variant haplotype 54% / 65% of informative HLA-F RNA reads (known-allele control 58%); fewer than 20 carriers per signature, 24 analysable runs, no confidence interval", "a5"),
    ("S08-6", "…and are ready for an IPD submission", "Open", "raw reads support both in every examined carrier (S09 script 105, up to 8 carriers, relatedness not re-checked); still missing an independent trio and phased full-length exons", "a6"),
    ("S08-7", "Trio transmission test", "Not feasible yet", "only 827 of 12,233 long-read people have a relative in the cohort", "a6"),
    ("S08-8", "Reference-bias proxy separates novel alleles from artefacts", "Inconclusive", "burden rises with distance from reference: expected for both real alleles and mapping artefacts", "a5"),
]


def dt(h, r):
    return table(h, r)


def sec_link(i, n):
    return '<a href="#%s">section %s</a>' % (i, n)


# ===================================================================== sections
def overview(review_html):
    imp = """<div class="box imp"><h3>Most important or most novel findings</h3><ol>
<li><b>Most of the "novel" HLA calls are not new alleles.</b> {} <span class="why">Of 8,015 calls labelled novel, 86% of the novel-protein ones are indel-only (mostly one-base errors in repeat runs); only 4% of indel calls are supported by short reads, and indel "novelty" is 1.6 to 3.2 times as high on the older Sequel II/IIe platform as on Revio. See {}.</span></li>
<li><b>A sizeable minority of the remaining "novelty" is simply the age of the reference database.</b> {} <span class="why">Re-checking against the current IPD release (3.66, eleven releases newer than the 3.55 used for the calls) already names 62 of 441 substitution signatures (14%) and 16 of the 49 that recur in two or more people (33%), including the big East-Asian HLA-G cluster (= G*01:01:33). Indel error (item 1) is the bigger effect. See {}.</span></li>
<li><b>Two recurrent signatures survive every test we could run: HLA-F T237I and HLA-F P205T.</b> {} <span class="why">Unnamed in IPD 3.66, each sits on a single flanking haplotype (as the known controls do), shows no platform or centre skew, is shared by relatives where a relative exists, and is transcribed in blood (variant haplotype carries 54% and 65% of HLA-F RNA reads; a known-allele control gives 58%). They are <b>probably real alleles</b>, not proven: a haplotype-specific assembly error would look the same. See {}.</span></li>
<li><b>The recurrent variants follow ancestry.</b> {} <span class="why">That rules out random noise and points to inherited, haplotype-linked variants in groups that IPD under-represents, so "novel vs 3.55" partly measures reference bias. It cannot by itself rule out a systematic error tied to one haplotype. See {}.</span></li></ol></div>""".format(
        badge("strong"), sec_link("a2", 2), badge("strong"), sec_link("a3", 3), badge("lead", "Supported"), sec_link("a5", 5), badge("lead", "Supported"), sec_link("a4", 4))
    prom = """<div class="box prom"><h3>Most promising lines of research (in order)</h3><ol>
<li><b>Build the IPD submission package for HLA-F T237I and P205T.</b> The only candidate new alleles in the whole exercise. <span class="why">Not yet submittable. The S09 raw-read pile-up (script 105) already supports both variants in every examined carrier (up to 8; whether they are unrelated was not re-checked). Still missing: an independent trio and phased full-length exons 2–3. A trio is what would turn "probably real" into "real".</span></li>
<li><b>Find out what the two HLA-F variants do.</b> They sit in exon 4 (α3 domain, outside the peptide groove, so no effect on peptide presentation is expected; function is unknown). Open conformers of HLA-F bind the KIR3DS1 receptor; whether these residues matter is untested. <span class="why">Power is the limit: carrier counts are small (the S09 power check estimated under 50), so disease associations cannot be tested yet.</span></li>
<li><b>Adopt the filter that makes "novel" mean novel.</b> Match against the current IPD release, drop suffix-flagged reference alleles and calls with more than ~15 edits, and require short-read support. <span class="why">This took 337 shortlisted calls to 98 defensible ones (mostly HLA-F, HLA-E and DQA1).</span></li>
<li><b>Get frequencies and transmission from outside this cohort.</b> Population frequencies (gnomAD/ALFA) for the two HLA-F variants and trio data from the S06 yak tables; both are open items in the S08 log.</li></ol></div>"""
    under = """<div class="box under"><h3>Underwhelming or negative (still informative)</h3><ul>
<li><b>The transmission test is unusable here:</b> only 827 of 12,233 long-read people have a relative in the cohort, so carrier children with two long-read parents are under 20 per signature. {}</li>
<li><b>The reference-bias proxy is inconclusive:</b> novel-substitution burden rises with how far a person's haplotypes are from the reference, but that is what both real alleles and mapping artefacts would do. {}</li>
<li><b>No phenotype link yet:</b> the health-record path was never probed, and rare alleles have almost no power anyway.</li>
<li><b>Several recurrent "novel" calls turned out known or artefactual</b> (DQB1 A57D = DQB1*02:180, HLA-E K195R = E*01:152, DRB1 R262T an artefact). The shrinking of the novel list is itself the main result of this sprint, not a failure of it.</li></ul></div>""".format(badge("null", "Not feasible"), badge("null", "Inconclusive"))
    cheat = """<details class="bg"><summary><span class="k">Cheat-sheet</span> Terms used in this report (30 seconds)</summary><div class="in"><dl class="cheat">
<dt>Novel call</dt><dd>An allele that Immuannot could not match to any named allele of IPD-IMGT/HLA release 3.55.</dd>
<dt>Indel / substitution</dt><dd>An inserted-or-deleted base versus a changed base. One-base indels in repeat runs are the classic long-read error.</dd>
<dt>Signature</dt><dd>One specific change on one reference allele (for example "F T237I on F*01:01:01:01"); recurrence means the same signature in several people.</dd>
<dt>IPD / 3.55 / 3.66</dt><dd>The reference database of named HLA alleles and its release numbers.</dd>
<dt>Haplotype</dt><dd>A stretch of one chromosome inherited as a block; a real inherited allele sits on one or two backgrounds, a recurring artefact on many.</dd>
<dt>Carrier / control</dt><dd>A person with the signature / matched people without it.</dd>
<dt>Immuannot</dt><dd>The annotation tool that matches each assembled HLA/KIR gene to named IPD alleles and labels the rest "new".</dd>
<dt>Burden</dt><dd>Average number of calls of a given kind per person in a group.</dd>
<dt>Homopolymer</dt><dd>A run of identical bases (AAAAA); long reads often get its length wrong by one, which creates a fake indel.</dd>
<dt>Rebuildable</dt><dd>The call's alignment can be reconstructed from the public reference allele and Immuannot's edit string, so its exact changes are known.</dd>
<dt>Tier-A / suffix-flagged</dt><dd>Tier-A: the short-read check supports the call. Suffix-flagged: the reference allele carries an N/Q/S null or questionable suffix, so calls made against it are unreliable.</dd>
<dt>Effective backgrounds</dt><dd>How many genuinely different neighbouring haplotypes the carriers sit on (1 = all the same).</dd>
<dt>Bias-corrected log2</dt><dd>RNA-read ratio of the variant copy to the other copy, corrected with simulated reads for how well each copy aligns; 0 = balanced.</dd>
<dt>Counts &lt;20</dt><dd>Group sizes below 20 are never shown (written "&lt;20") to protect participants.</dd></dl></div></details>"""
    return """<h2 id="glance">At a glance</h2>
<p class="sub">What this sprint asked: 8,015 HLA calls that the annotation pipeline labelled "new" in about 12,200 long-read All of Us genomes: are they real alleles, and do they matter? Approach: try to break each claim with raw reads, short reads, platforms, ancestry, relatives, RNA and the current reference database.</p>
<div class="glance three">{}{}{}</div>{}{}""".format(imp, prom, under, cheat, review_html)


def a1():
    flow = F.flow([("Assemblies", "~12.2k people|14 HLA genes", "a"), ("Call vs 3.55", "'new' = no match|in the reference", "a"), ("Rebuild", "each alignment|from public CDS", "b"),
                   ("Classify", "indel / substitution|/ synonymous", "b"), ("Check", "short reads, platform,|ancestry, relatives", "c"), ("Re-test 3.66", "47,419 named|sequences", "c"), ("Biology", "carriers, RNA,|literature", "d")], title="How each 'novel' call was tested", h=100)
    return section("a1", "1. What 'novel' means and how it was tested", badge("known", "Method"),
                   "<b>'Novel' only means 'no match in the reference snapshot used'.</b> Every claim in this report comes from taking Immuannot's own alignment of each novel call, rebuilding it from the public reference allele, classifying the edits, and then asking independent questions: do short reads agree, does it follow the sequencing platform, does it follow ancestry, is it already named in the current database, do relatives share it, is it transcribed.",
                   cards([("8,015", "novel HLA calls vs IPD 3.55 in ~12,200 genomes"), ("14", "HLA genes annotated"), ("47,419", "named coding sequences in IPD 3.66 used for the re-check"), ("<20", "group sizes below this are never shown")]),
                   fig(flow, "Each stage uses only Immuannot's alignments plus public IPD sequence. Person-level files stay in the Workbench; only suppressed aggregates leave it.", "scripts 90–96 (VM), LOG.md of S08") +
                   ev("What is an HLA allele name, and why does a database release matter?", "<p>HLA genes are the most variable in the genome; each version is an <b>allele</b> named with up to four fields (first two change the protein). IPD-IMGT/HLA is the reference list, updated several times a year. A variant that is rare in the European donors who dominated early releases can be perfectly real and still look 'novel' until a later release names it.</p><figure>%s</figure>" % BG.fig_resolution(), "bg") +
                   ev("What raw reads versus an assembly mean", "<p>Reads are short pieces of sequence; software stitches them into an <b>assembly</b>, and alleles are called from the assembly. An error in the assembly (most often a missing or extra base in a run of identical bases) creates a 'novel allele' no read supports. Checking the raw reads is the cleanest test.</p><figure>%s</figure>" % BG.fig_reads(), "bg") +
                   ev("Scripts and data", "<p><code>90_case_browser*.py</code> (viewer and per-case alignments), <code>91_novel_recurrence.py</code> (signatures, ancestry, relatives, matching against IPD), <code>92_novel_burden.py</code> (platform × ancestry), <code>93_shortlist_ipd_release.py</code>, <code>95_signature_context.py</code> (haplotype origin, centre skew), <code>96_f_signature_ase.py</code> (RNA allele-specific expression). VM runs used one core at low priority from a throwaway clone that was removed afterwards; raw tables and carrier lists are mode 600 on the VM and not in the repository.</p>", "how"))


def a2():
    f_edit = fig(F.stack100([(a, v) for a, v in EDIT_TYPES], [("indel only", "warn"), ("indel + substitution", "b"), ("substitution only", "c")], "What kind of edit makes a call 'novel'?", label_w=250),
                 "86% of the 7,635 novel-protein calls consist only of insertions or deletions (93% carry a frameshift; at one edit, 91% are indel-only and 72% of those are in homopolymer or tandem runs). Only 6% (449 calls) are pure substitutions. The 227 novel-synonymous calls are all substitution-only.",
                 "script 90 stats pass on the VM (LOG 10-21)",
                 dt(["Class", "indel only %", "indel + substitution %", "substitution only %"], [(a, v[0], v[1], v[2]) for a, v in EDIT_TYPES]))
    f_sr = fig(F.stack100([(a, v) for a, v in SR_SUPPORT], [("supported", "c"), ("not seen in short reads", "warn"), ("contradicted", "b"), ("other / half the reads", "mut")], "Do independent short reads agree with the call?", label_w=300),
               "Short reads realigned to the called haplotype (567 cases, 1,097 sites). Only 4% of indel calls are supported (79% not seen, 8% contradicted); 60% of protein-altering substitutions and 68% of synonymous ones are supported. 'Other' is the remainder not classified in the record (half the reads or other haplotype disagrees).",
               "script 90 short-read panel (LOG 10-21)",
               dt(["Class", "supported %", "not seen %", "contradicted %", "other %"], [(a, v[0], v[1], v[2], v[3]) for a, v in SR_SUPPORT]))
    f_plat = fig(F.grouped(ANC, [("Revio", REVIO_INDEL, "a"), ("Sequel II/IIe", SEQ2E_INDEL, "b")], 1.2, "Indel-bearing novel-protein calls per person, classical genes", fmt="{:.2f}", height=300),
                 "If 'novel' means error it should follow the sequencing platform. Indel novelty is 1.6× (AFR) to 3.2× (SAS) higher on Sequel II/IIe (the older chemistry) than on Revio, at every ancestry. Sequel II/IIe has about 0.9k people in total, so some ancestry groups are small (≥20 persons shown).",
                 "script 92 (12,233 people; groups ≥20 persons)",
                 dt(["Ancestry", "Revio", "Sequel II/IIe"], [(a, "%.3f" % b, "%.3f" % c) for a, b, c in zip(ANC, REVIO_INDEL, SEQ2E_INDEL)]))
    f_snv = fig(F.grouped(ANC, [("Revio", REVIO_SNV, "a"), ("Sequel II/IIe", SEQ2E_SNV, "b")], 0.04, "Substitution-only novel-protein calls per person, classical genes", fmt="{:.3f}", height=300),
                "Substitution-only calls show no consistent platform effect (higher on Sequel II/IIe in four ancestries, lower in two, on small groups), unlike indels. (The Sequel II/IIe MID bar rests on 22 people, so it is noise.)",
                "script 92",
                dt(["Ancestry", "Revio", "Sequel II/IIe"], [(a, "%.3f" % b, "%.3f" % c) for a, b, c in zip(ANC, REVIO_SNV, SEQ2E_SNV)]))
    return section("a2", "2. Most 'novel' calls are indel errors", badge("strong"),
                   "<b>Low-distance 'novel' calls are mostly consensus errors, not new alleles.</b> 86% of novel-protein calls are indels only, almost all one base long and often in runs of identical bases; only 4% of them are supported by independent short reads; and they are about 1.6 to 3.2 times as common on the older sequencing platform (Sequel II/IIe, about 0.9k of the 12,233 people). Substitution-only calls show no consistent platform pattern.",
                   cards([("86%", "of 7,635 novel-protein calls are indel-only"), ("4%", "of indel calls supported by short reads"), ("1.6–3.2×", "more indel 'novelty' on Sequel II/IIe than Revio (AFR 1.6× … SAS 3.2×)"), ("449", "pure-substitution novel-protein calls: the real candidates")]),
                   f_edit + f_sr + f_plat + f_snv +
                   ev("Detail from the case browser", "<ul><li>All 249 ladder alignments replay Immuannot's own alignment string exactly.</li><li>Indels do not recur at fixed positions (462 distinct sites among 486 events, 23 recurrent), are spread uniformly along the coding sequence, and 41 of 140 indel-carrying cases restore the reading frame with a second indel: the pattern of sequencing error, not of mutation.</li><li>Mendelian check (30 sampled mismatches): 17 were within 6 edits of a parent haplotype, and the child's own call was 0 edits from the parent in all 17.</li></ul>", "ev") +
                   ev("What could still be wrong", "<ul><li>The 762 browser cases are stratified examples, not population rates; the whole-cohort statistics (8,015 calls) come from the separate stats pass.</li><li>Short-read verdicts can't check DRB3/DRB4 (no primary-assembly locus).</li><li>Platform effects were tested for the classical genes only.</li></ul>", "cav"))


def a3():
    f_fun = fig(F.hbars(FUNNEL, vmax=9000, title="From 8,015 novel calls to the defensible recurrent set", label_w=370, right=190),
                "'Rebuildable' means the alignment can be reconstructed from the public reference allele and Immuannot's cs string. Most signatures are seen in one person only; of the 49 seen in at least two people, 15 are exact matches to an IPD 3.66 allele and one a partial match; two (HLA-F T237I, P205T) are carried forward here, and the rest fail a filter or have too few carriers to examine. Four small clusters remain unnamed and were not pursued in this report (HLA-G A235V, HLA-G E190A on a Q-suffix reference, HLA-F T279I, one HLA-F synonymous).",
                "script 91 with --ipd-nuc hla_nuc.fasta (VM); LOG 10-23",
                dt(["Step", "count"], [(r[0], r[1]) for r in FUNNEL]))
    f_sl = fig(F.hbars(SHORTLIST, vmax=380, title="The 337-call shortlist against IPD 3.66", label_w=430, right=130),
               "58 shortlist calls are exact matches to a named allele in IPD 3.66 and 9 contain a partial entry. A further filter (tier-A short-read support, ≤15 edits, reference allele not suffix-flagged) leaves 98, mostly HLA-F (36), HLA-E (15) and DQA1 (12); 28 are in the peptide groove. These 98 are a filtered follow-up list, not a recurrence-tested set (most are single calls), which is why only 2 appear in the funnel above.",
               "script 93 (local) on the 337-case shortlist",
               dt(["Step", "count"], [(r[0], r[1]) for r in SHORTLIST]))
    return section("a3", "3. Much of the 'novelty' is the age of the reference", badge("strong"),
                   "<b>Database age explains a sizeable minority of the novelty, not most of it.</b> The pipeline used IPD release 3.55; the current release is 3.66. Re-matching names 62 of the 441 substitution signatures (14%; 54 exact, 8 partial) and 58 of the 337 shortlisted calls (17%); among signatures seen in two or more unrelated people, 16 of 49 (33%) are named. The larger share of the original 8,015 is indel error (section 2). Two named cases are worth remembering: the large East-Asian HLA-G cluster is the known allele G*01:01:33, and DQB1 A57D is reported as DQB1*02:180:02 (an Asp57 variant relevant to type 1 diabetes; from a literature check, not re-derived here).",
                   cards([("62 / 441 (14%)", "substitution signatures already named in IPD 3.66"), ("16 / 49 (33%)", "of the recurrent signatures (≥2 people) already named"), ("58 / 337 (17%)", "shortlist calls that are exact 3.66 matches"), ("11", "IPD releases between 3.55 and 3.66")]),
                   f_fun + f_sl + ev("Which signatures were named, and which were not", dt(["Signature", "Unrelated carriers", "Status vs IPD 3.66", "Note"], SIGS), "ev", True) +
                   ev("Denominators: why 8,015, 7,635, 682 and 675 differ", "<ul><li><b>8,015</b> novel calls with a coding-sequence alignment = <b>7,635</b> novel-protein + <b>227</b> novel-synonymous + <b>153</b> beyond the coding sequence (137 further calls could not be rebuilt from the edit string).</li><li><b>682</b> substitution-only calls come from the case-browser statistics pass (script 90); <b>675</b> substitution-only and rebuildable come from the recurrence script (91). They count slightly differently; the 7-call gap is consistent with the rebuildability requirement but was not traced call by call.</li><li><b>441</b> distinct signatures = distinct (gene, reference allele, set of substitutions) among the 675.</li><li><b>337</b> shortlist = protein-altering substitution-only calls kept for manual QC.</li></ul>", "how") + ev("What could still be wrong", "<ul><li>'*' in a partial IPD entry means an unsequenced position, which was not resolved.</li><li>The 3.66 re-check was run from the public nucleotide file; names are the same ones IPD would assign for exact coding matches, but the WHO committee assigns official names.</li><li>The DQB1 A57D = 02:180 naming rests on a literature check that has not been re-derived locally.</li></ul>", "cav"))


def a4():
    f_anc = fig(F.grouped(ANC, [("per person", ALL_SNV, "c")], 0.026, "Substitution-only novel-protein calls per person, by ancestry", fmt="{:.3f}", height=270, legend=False),
                "Counts per person of substitution-only novel-protein calls in the classical genes (the 'burden'). Higher in AFR, EAS and MID (about 1.5×) than in EUR, AMR and SAS: the pattern expected when the reference is built mostly from European donors.",
                "script 92 (12,233 people)",
                dt(["Ancestry", "persons", "calls per person"], [(a, p, "%.3f" % v) for a, p, v in zip(ANC, PERSONS, ALL_SNV)]))
    f_all = fig(F.grouped(ANC, [("per person", ALL_NOVEL, "d")], 4.3, "Novel calls of any kind per person (incl. synonymous / beyond the coding sequence)", fmt="{:.2f}", height=270, legend=False),
                "Ordering AFR ≈ EAS > SAS > MID > AMR > EUR. Caution: this count includes indel errors, so unlike the substitution-only panel above it is not a clean reference-bias measure.",
                "script 92",
                dt(["Ancestry", "novel calls per person"], [(a, "%.2f" % v) for a, v in zip(ANC, ALL_NOVEL)]))
    return section("a4", "4. The recurrent signatures follow ancestry", badge("lead", "Supported"),
                   "<b>Inherited, haplotype-linked variants of under-represented groups, not discoveries of new biology.</b> Recurrent signatures are concentrated in a few ancestries (for example the HLA-G cluster is almost entirely East Asian; HLA-F T237I and P205T are mostly African-ancestry), and each sits on a single reference allele. Per person, novelty is highest in AFR and EAS and lowest in EUR.",
                   cards([("49 of 441", "signatures seen in ≥2 unrelated people"), ("1 reference", "allele under every recurrent signature"), ("AFR, EAS", "highest novelty per person")]),
                   f_anc + f_all + ev("Why ancestry structure points to real variants", "<p>A random error would appear in every ancestry at about the same rate and on many different backgrounds, so ancestry structure and a single reference allele rule out <i>random</i> noise. They do not rule out a <i>systematic</i> error tied to one haplotype or to a poorly represented reference template, which would also be ancestry-structured and sit on one background. The pattern therefore supports 'inherited, haplotype-linked' but not, by itself, 'real allele'.</p>", "bg") +
                   ev("What could still be wrong", "<ul><li>Ancestry labels are AoU predictions.</li><li>Unrelated carriers are counted after removing relatives; group sizes under 20 are not shown.</li></ul>", "cav"))


def a5():
    f_hap = fig(F.grouped(HAP_NAMES, [("carrier haplotypes", HAP_CARRIER, "b"), ("matched non-carrier haplotypes (median of 500 draws)", HAP_CONTROL, "a")], 10, "Effective number of different backgrounds (HLA-A + HLA-G/F neighbours)", fmt="{:.1f}", height=310),
                "A recurrent mutation or a recurring artefact would sit on many different backgrounds; an inherited allele sits on one. Both unnamed HLA-F signatures are on essentially a single background (1.1 and 1.4), just like the two known controls (1.6, 1.0), while matched non-carriers of the same lineage have 5–8.",
                "script 95 (VM)",
                dt(["Signature", "carrier haplotypes", "non-carrier median"], [(n, "%.2f" % a, "%.2f" % b) for n, a, b in zip(HAP_NAMES, HAP_CARRIER, HAP_CONTROL)]))
    f_ase = fig(F.hbars([(a, b, "main" if "unnamed" in a else "ref") for a, b in ASE], vmax=1.0, title="Share of HLA-F RNA reads coming from the variant haplotype (blood RNA-seq)", label_w=330, fmt="{:.2f}", right=100),
                "In a heterozygote about half of the RNA should come from each copy. Both unnamed alleles are transcribed (0.54 and 0.65), like the known F*01:25 control (0.58). After correcting for reference bias with simulated reads the pattern in the controls and unnamed alleles is the same, so no allele-specific expression is claimed.",
                "script 96 on personal two-haplotype references (VM; fewer than 20 carriers per signature, 24 analysable runs). Bias-corrected median log2(variant/other): +0.26 (T237I), +0.56 (P205T), +0.47 (F*01:25 control); no confidence intervals were computed",
                dt(["Signature", "variant-haplotype read fraction", "bias-corrected median log2(variant/other)"], [(a, "%.2f" % b, c) for (a, b), c in zip(ASE, ["+0.26", "+0.56", "+0.47"])]))
    f_td = fig(F.grouped(ANC, [("low divergence", TD_LOW, "a"), ("middle", TD_MID, "c"), ("high divergence", TD_HIGH, "b")], 0.17, "Substitution-only novel signatures per person by distance from the reference", fmt="{:.2f}", height=310),
               "Within every ancestry the burden rises with how far a person's haplotypes sit from the reference. That fits real alleles living on poorly represented haplotypes, but a mapping artefact would do the same: the proxy cannot separate them (inconclusive by design).",
               "script 95",
               dt(["Ancestry", "low", "middle", "high"], [(a, "%.3f" % b, "%.3f" % c, "%.3f" % d) for a, b, c, d in zip(ANC, TD_LOW, TD_MID, TD_HIGH)]))
    return section("a5", "5. HLA-F T237I and P205T: the two alleles that survive", badge("strong"),
                   "<b>Both are probably real inherited alleles, unnamed in IPD and transcribed; this is supported, not proven.</b> Each sits on a single flanking haplotype, matches the ancestry-expected number of carriers per sequencing platform and centre, is shared by relatives wherever a relative exists, and is expressed in blood RNA (54% and 65% of HLA-F reads from the variant copy; a known allele control gives 58%). Both lie in exon 4 (α3 domain), outside the peptide-binding groove, so no effect on peptide presentation is expected and their function is unknown. A haplotype-specific assembly or template error would also give a single background, which is why raw-read confirmation in unrelated carriers and a trio are still needed.",
                   cards([("1.1 and 1.4", "effective backgrounds among carriers (non-carriers: 5–8)"), ("22 vs 22.1", "Revio carriers observed vs expected (T237I)"), ("0.54 / 0.65", "variant share of HLA-F RNA reads (control 0.58)"), ("exon 4", "α3 domain, outside the groove")]),
                   f_hap + ev("Batch effects and relatives", "<p><b>Platform and sequencing centre:</b> observed carriers per platform and per centre match the ancestry-matched expectation for all four signatures (for example Revio 22 observed vs 22.1 expected for HLA-F T237I; the Broad centre 86 vs 85.7 for the G control). One centre is mildly above expectation for T237I, not significant.</p><p><b>Relatives:</b> where a relative is in the long-read cohort the relative carries the same signature in almost every case (counts under 20 are withheld). Consistent with inheritance but too few for a formal transmission rate.</p>", "ev") + f_ase + f_td +
                   ev("Quality checks on each signature", "<ul><li>Each signature is always against a single reference allele, strands are balanced, mapping quality is ≈60 and there is no CpG context (a hypermutable site type). Short-read support is described in the log as 'mostly' clean; no per-signature count was recorded.</li><li>'Matched non-carrier' haplotypes are drawn from people of the same HLA-F lineage and ancestry (500 draws); 'expected' carriers per platform/centre come from the ancestry-matched share of that platform/centre.</li><li>HLA-F and HLA-E alleles are absent from the only local IPD FASTA (3.38), so cross-allele checks used the public 3.66 alignment files.</li></ul>", "ev") +
                   ev("What could still be wrong", "<ul><li>Still missing for a submission: an independent trio and phased full-length exons. (The S09 raw-read pile-up supports both variants in every examined carrier, up to 8; relatedness among them was not re-checked.)</li><li>Fewer than 20 carriers per signature have RNA, so the expression result rests on few people (exact counts withheld).</li></ul>", "cav"))


def a6():
    return section("a6", "6. What it would take to claim a new allele", badge("open"),
                   "<b>IPD-IMGT/HLA requires phased, full-length exons (2–3 for class I; 2 for class II), the new allele sequenced apart from its partner, and an official name assigned by the WHO Nomenclature Committee.</b> The transmission evidence that would normally help is not available in this cohort.",
                   cards([("827 / 12,233", "long-read people with a relative in the cohort"), ("~41%", "of assembled HLA/KIR sequences novel vs IPD 3.52 in the Immuannot paper (mostly non-classical)"), ("trio", "an independent parent-child trio is the main missing piece")]),
                   ev("The think-tank consensus (what we would do next, in order)", "<ol><li><b>Transmission oracle across platforms</b> (child long-read vs parental short-read, because same-pipeline errors transmit).</li><li><b>Simulated reference-bias null.</b></li><li><b>RNA-seq allele-specific expression</b> on personal references (done, section 5).</li><li><b>Primate outgroup</b> (can only confirm antiquity).</li><li><b>Per-run error barcode</b> for hidden batch effects.</li><li><b>Mark-recapture per-platform error</b> estimate.</li></ol><p>Check trio × carrier × RNA overlap before committing to any of them: power may be near zero for HLA-F T237I/P205T.</p>", "ev") +
                   ev("Honest limits", "<ul><li>Novelty was defined against IPD 3.55; 3.66 changes the answer and 3.67 will again.</li><li>Signatures are keyed on one reference allele; ancestry labels are AoU predictions; relatedness data is sparse.</li><li>No phenotype comparison: the health-record path was not probed and rare alleles have little power.</li></ul>", "cav"))


def ledger_section():
    rows = [(a, b, "<span class='badge %s'>%s</span>" % ({"Strong": "strong", "Supported": "lead", "Open": "open", "Not feasible yet": "null", "Inconclusive": "null"}[c], c), d, '<a href="#%s">%s</a>' % (e, e[1:])) for a, b, c, d, e in LEDGER]
    return '<h2 id="ledger">Ledger of every hypothesis</h2>' + dt(["ID", "Hypothesis", "Status", "Evidence", "Section"], rows)


def build(review_html=""):
    toc = [("glance", "At a glance"), ("a1", "1 Method"), ("a2", "2 Indel errors"), ("a3", "3 Reference age"), ("a4", "4 Ancestry"), ("a5", "5 HLA-F alleles"), ("a6", "6 Next steps"), ("ledger", "Ledger")]
    lede = '<p class="sub">Omni-HLA · sprint S08 · rebuilt 28 Oct 2026 · aggregate results only (counts of 1–19 written &lt;20) · companion to the S09 report</p>'
    body = overview(review_html) + a1() + a2() + a3() + a4() + a5() + a6() + ledger_section()
    return F.page("S08 novel HLA alleles", "Omni-HLA · S08", "Are the 'novel' HLA alleles real, and do they matter?", lede, toc, body)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--review", default="")
    a = ap.parse_args(argv)
    rv = open(a.review).read() if a.review and os.path.exists(a.review) else ""
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w").write(build(rv))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
