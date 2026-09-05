#!/usr/bin/env python3
"""Synthetic fixtures reproducing the REAL Immuannot output layout, faithfully enough that every
script in scripts/hla_popgen/ can be developed and tested without VM access.

Why this exists: the production data lives on a Verily Workbench VM behind a VPC-SC perimeter that
blocks SSH entirely (ENVIRONMENT.md quirk #28). Every round-trip to real data costs a human
paste-cycle and real money. This project's own history is emphatic that fixtures catch real bugs
before handoff (5 bugs caught this way in scripts/production_analysis/, plus the merge-bug
reproduction in rebuild_immuannot_calls.py) -- so nothing gets handed to Marc until it runs clean
against these.

Format fidelity is sourced from a full read of Immuannot's upstream source, archived at
reference/IMMUANNOT_GTF_SPEC.md. The details that matter and are easy to get wrong:
  - template_distance is UNQUOTED in the GTF; every other attribute is quoted.
  - "new" is spliced INTO the consensus allele string at a variable field depth (2/3/4), which
    encodes novelty severity -- it is not a separate flag attribute.
  - cds_distance / cds_mut / template_warning are CONDITIONAL and legitimately absent.
  - column 1 is the real contig name, and ONE hap file can hold MULTIPLE contigs -- the cis-phasing
    trap this whole sub-project exists to fix.
  - the #-prefixed header carries copy-number and per-gene contig lists.

Usage:
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures
"""

import argparse
import gzip
import hashlib
import json
import os
import random

# Deterministic: a fixture that shifts between runs is worthless for regression testing.
SEED = 20260902

GENE_CLASSES = {
    "HLA-A": "classical_I", "HLA-B": "classical_I", "HLA-C": "classical_I",
    "HLA-DPA1": "classical_II", "HLA-DPB1": "classical_II", "HLA-DQA1": "classical_II",
    "HLA-DQB1": "classical_II", "HLA-DRB1": "classical_II",
    "HLA-DRA": "class_II_accessory", "HLA-DMA": "class_II_accessory",
    "HLA-DRB3": "class_II_paralog", "HLA-DRB4": "class_II_paralog",
    "HLA-DRB5": "class_II_paralog", "HLA-DQA2": "class_II_paralog",
    "HLA-E": "nonclassical_I", "HLA-F": "nonclassical_I", "HLA-G": "nonclassical_I",
    "HLA-H": "pseudogene_I", "HLA-J": "pseudogene_I",
    "MICA": "mic_tap", "MICB": "mic_tap", "TAP1": "mic_tap", "TAP2": "mic_tap",
    "C4A": "complement", "C4B": "complement",
}

# Realistic-ish allele pools. Not exhaustive -- enough that frequency spectra and ancestry
# stratification have something to bite on.
ALLELE_POOL = {
    "HLA-A": ["01:01:01:01", "02:01:01:01", "03:01:01:01", "23:01:01:01", "24:02:01:01",
              "30:01:01:01", "33:03:01:01", "68:02:01:01"],
    "HLA-B": ["07:02:01:01", "08:01:01:01", "15:03:01:01", "35:01:01:01", "42:01:01:01",
              "53:01:01:01", "57:01:01:01", "58:01:01:01"],
    "HLA-C": ["04:01:01:01", "06:02:01:01", "07:01:01:01", "07:02:01:01", "17:01:01:01"],
    "HLA-DPA1": ["01:03:01:01", "02:01:01:01", "02:02:02:01"],
    "HLA-DPB1": ["01:01:01:01", "02:01:02:01", "04:01:01:01", "04:02:01:01", "17:01:01:01"],
    "HLA-DQA1": ["01:02:01:01", "02:01:01:01", "03:01:01:01", "04:01:01:01", "05:01:01:01"],
    "HLA-DQB1": ["02:01:01:01", "02:02:01:01", "03:01:01:01", "05:01:01:01", "06:02:01:01"],
    "HLA-DRB1": ["03:01:01:01", "04:01:01:01", "07:01:01:01", "11:01:01:01", "13:02:01:01",
                 "15:01:01:01", "15:03:01:01"],
    "HLA-DRA": ["01:01:01:01", "01:02:02:01"],
    "HLA-DMA": ["01:01:01:01", "01:02:01:01"],
    "HLA-DRB3": ["01:01:02:01", "02:02:01:01"],
    "HLA-DRB4": ["01:03:01:01"],
    "HLA-DRB5": ["01:01:01:01", "02:02:01:01"],
    "HLA-DQA2": ["01:01:01:01"],
    "HLA-E": ["01:01:01:01", "01:03:02:01"],
    "HLA-F": ["01:01:01:01"],
    "HLA-G": ["01:01:01:01", "01:04:01:01"],
    "HLA-H": ["01:01:01:01"],
    "HLA-J": ["01:01:01:01"],
    "MICA": ["002:01:01", "008:01:01", "010:01:01"],
    "MICB": ["002:01:01", "005:02:01"],
    "TAP1": ["01:01:01", "02:01:01"],
    "TAP2": ["01:01:01", "01:02:01"],
}

# Ancestry-skewed sampling weights, so novel-allele-by-ancestry and frequency-by-ancestry figures
# have real signal to detect rather than uniform noise. Deliberately crude.
ANCESTRY_BIAS = {
    "afr": {"HLA-B": ["15:03:01:01", "53:01:01:01", "58:01:01:01"],
            "HLA-A": ["23:01:01:01", "30:01:01:01", "33:03:01:01"]},
    "eur": {"HLA-B": ["07:02:01:01", "08:01:01:01", "57:01:01:01"],
            "HLA-A": ["01:01:01:01", "02:01:01:01", "03:01:01:01"]},
    "eas": {"HLA-B": ["35:01:01:01", "15:03:01:01"], "HLA-A": ["24:02:01:01", "02:01:01:01"]},
    "amr": {"HLA-B": ["35:01:01:01", "07:02:01:01"], "HLA-A": ["02:01:01:01", "24:02:01:01"]},
    "sas": {"HLA-B": ["07:02:01:01", "42:01:01:01"], "HLA-A": ["01:01:01:01", "33:03:01:01"]},
    "mid": {"HLA-B": ["35:01:01:01", "08:01:01:01"], "HLA-A": ["02:01:01:01", "01:01:01:01"]},
}
ANCESTRIES = ["afr", "amr", "eas", "eur", "mid", "sas"]
# Deliberately non-uniform, roughly echoing All of Us's own long-read composition.
ANCESTRY_WEIGHTS = [0.30, 0.20, 0.08, 0.34, 0.04, 0.04]

BASES = "ACGT"


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def _build_novel_variants():
    """Pre-generate a FINITE pool of novel alleles per gene, each with a fixed observed CDS
    sequence and an ancestry affinity.

    This is essential, not cosmetic. An earlier version of this generator drew a fresh random CDS
    for every novel call, which meant no novel allele was ever seen twice -- so cross-person
    recurrence never occurred, `passes_qc` (which requires >=2 unrelated carriers, the primary
    defence against mistaking an assembly artifact for a real allele) could never fire, and Chao2
    hit its degenerate Q2=0 regime on every gene. A finite pool with reuse is also simply what
    reality looks like: a real novel allele segregates in a population and recurs.

    Ancestry affinity makes the pool non-uniform so that novel alleles CLUSTER by ancestry -- the
    exact claim the paper wants to make, and therefore the exact thing the analysis must be able
    to detect on data where we know the ground truth.
    """
    rng = random.Random(SEED + 1)
    variants = {}
    for gene, pool in ALLELE_POOL.items():
        n_variants = max(4, len(pool))
        gene_variants = []
        for _ in range(n_variants):
            base_allele = rng.choice(pool)
            fields = base_allele.split(":")
            depth = rng.choices([2, 3, 4], weights=[0.3, 0.35, 0.35])[0]
            novel_allele = ":".join(fields[: depth - 1] + ["new"])
            td = rng.randint(1, 12)
            # Affinity: most novel alleles are near-private to one ancestry, a few are shared.
            affinity = rng.choice(ANCESTRIES) if rng.random() < 0.75 else None
            # Depth 1 ("HLA-A*new" -- even the first field undetermined) is rare but REAL: the
            # production recon found 3 such calls in a 50-person sample. It arises when
            # consensusCall()'s commonprefix truncation collapses tied candidates that disagree at
            # the very first field. Downstream code must classify it, not crash or mislabel it.
            if rng.random() < 0.03:
                gene_variants.append({
                    "consensus": f"{gene}*new", "td": rng.randint(1, 40), "cds_dist": None,
                    "cds_mut": None, "warning": "partial_CDS", "affinity": affinity,
                    "cds_seq": rand_seq(rng, rng.randrange(810, 1101, 3)),
                })
                continue
            if depth == 4:
                # Difference lies outside the CDS (intron/UTR) -- no CDS-level diff to report.
                cds_dist, cds_mut, warning = 0, None, None
            else:
                cds_dist = rng.randint(1, 4)
                if depth == 3:
                    # Synonymous: codon changed, amino acid did not.
                    muts = ":".join(f"K(AA{rng.choice('AG')})<K(AA{rng.choice('AG')})"
                                    for _ in range(cds_dist))
                else:
                    muts = ":".join(f"{rng.choice('KRDEG')}({rand_seq(rng, 3)})<"
                                    f"{rng.choice('LVIST')}({rand_seq(rng, 3)})"
                                    for _ in range(cds_dist))
                cs = ":" + str(rng.randint(100, 900)) + "*" + rand_seq(rng, 2).lower()
                cds_mut = f"{gene}*{base_allele}|{cs}|{muts}"
                warning = rng.choice([None, None, None, "partial_CDS", "no-stop_codon"])
            gene_variants.append({
                "consensus": f"{gene}*{novel_allele}", "td": td, "cds_dist": cds_dist,
                "cds_mut": cds_mut, "warning": warning, "affinity": affinity,
                # The observed CDS is fixed per variant, so the same novel allele in two different
                # people yields byte-identical sequence and therefore clusters correctly.
                "cds_seq": rand_seq(rng, rng.randrange(810, 1101, 3)),
            })
        variants[gene] = gene_variants
    return variants


NOVEL_VARIANTS = _build_novel_variants()


def make_consensus(rng, gene, ancestry):
    """Return (consensus, template_distance, cds_distance, cds_mut, warning, cds_seq).

    Reproduces Immuannot's real novelty encoding: "new" replaces a field at depth 2 (protein
    altering), 3 (synonymous) or 4 (beyond-CDS), and cds_distance/cds_mut are absent entirely when
    the gene-level template was already a perfect match. `cds_seq` is None for non-novel calls.
    """
    pool = ALLELE_POOL[gene]
    biased = ANCESTRY_BIAS.get(ancestry, {}).get(gene)
    allele = rng.choice(biased) if (biased and rng.random() < 0.6) else rng.choice(pool)

    # Novel rate is deliberately ancestry-skewed: the IPD reference is European-biased, so
    # non-European haplotypes should throw more novel calls. This is the signal the real analysis
    # is hunting for -- the fixture must contain it so we can prove the code detects it.
    novel_rate = {"afr": 0.22, "amr": 0.13, "eas": 0.12, "sas": 0.12,
                  "mid": 0.11, "eur": 0.05}[ancestry]

    if rng.random() >= novel_rate:
        # Exact match to a documented allele: template_distance 0, and Immuannot never runs a
        # CDS-level search at all, so cds_distance/cds_mut are legitimately ABSENT (not zero).
        return f"{gene}*{allele}", 0, None, None, None, None

    # Draw from the finite pool, strongly preferring variants matching this person's ancestry so
    # novel alleles cluster by ancestry the way real ones do.
    cands = NOVEL_VARIANTS[gene]
    matched = [v for v in cands if v["affinity"] == ancestry]
    v = rng.choice(matched) if (matched and rng.random() < 0.8) else rng.choice(cands)
    return (v["consensus"], v["td"], v["cds_dist"], v["cds_mut"], v["warning"], v["cds_seq"])


def write_gtf(path, contigs_genes, rng):
    """contigs_genes: {contig_name: [(gene, copy_index, ancestry), ...]}"""
    body = []
    per_gene_contigs = {}
    copynum = {}
    novel_cds = {}

    for contig, entries in contigs_genes.items():
        pos = 1000
        for gene, copy_index, ancestry in entries:
            consensus, td, cds_dist, cds_mut, warning, cds_seq = make_consensus(
                rng, gene, ancestry)
            # THE "NA" TRAP. Immuannot writes template_warning "NA" to mean "no warning" far more
            # often than it omits the attribute -- measured on the production cohort (2026-09-03,
            # 00b_warning_census.py, 200 people): 57.4% literal "NA", 4.6% attribute absent, and
            # only ~38% a real warning token. Any consumer doing bool(template_warning) or
            # splitting on "," without excluding "NA" will score ~62% of all calls as warned when
            # they are clean. The existing scripts/production_orchestrator/rebuild_immuannot_calls.py
            # gets this right (`warn_val not in {"", "NA"}`); a first pass of 01/03 here did not.
            #
            # So the fixture emits BOTH clean forms in roughly the observed ratio, and real tokens
            # at the real ~38% rate. A consumer that mishandles "NA" now fails loudly here.
            #
            # Token mix also mirrors reality: partial_CDS dominates, no-start/no-stop are
            # concentrated in pseudogenes (correct biology -- HLA-P/T/W genuinely lack valid
            # codons), and inframe_stop was observed ZERO times in 200 people.
            # Warning rate is GENE-CLASS DEPENDENT in reality, not uniform. Measured per-gene
            # clean rates (200 people): the 8 classical genes are 96-99% clean (97.6% overall,
            # partial_CDS the only meaningful token), while pseudogenes are almost never clean --
            # HLA-N and HLA-S are 100% partial_CDS, HLA-P/T/W are dominated by paired
            # no-start_codon/no-stop_codon. Modelling this as one flat 38% would make a
            # partial_CDS-disqualifying QC gate look far more destructive on fixtures than it is
            # on real data, and would mask the fact that warnings are correct biology for the
            # loci that carry them.
            if warning is None:
                gclass = GENE_CLASSES.get(gene, "other")
                if gclass in ("pseudogene_I", "class_II_paralog"):
                    warn_rate, pseudo = 0.95, True
                elif gclass in ("classical_I", "classical_II"):
                    warn_rate, pseudo = 0.024, False
                else:
                    warn_rate, pseudo = 0.15, False
                roll = rng.random()
                if roll < warn_rate:
                    warning = rng.choices(
                        ["partial_CDS", "no-start_codon,no-stop_codon"],
                        weights=[0.6, 0.4])[0] if pseudo else "partial_CDS"
                elif roll < warn_rate + (1 - warn_rate) * 0.92:
                    warning = "NA"      # literal string meaning "no warning"
                # else: attribute omitted entirely (the other clean form)
            template = f"{gene}*{rng.choice(ALLELE_POOL[gene])}"
            glen = rng.randint(3000, 15000)
            start, end = pos, pos + glen
            pos = end + rng.randint(500, 3000)
            strand = rng.choice("+-")
            sfx = "" if copy_index == 1 else f".{copy_index}"
            gid = f"IAG{abs(hash((contig, gene, copy_index))) % 900000 + 100000}{sfx}"
            tid = f"IAT{abs(hash((contig, gene, copy_index, 't'))) % 900000 + 100000}{sfx}"

            per_gene_contigs.setdefault(gene, set()).add(contig)
            copynum[gene] = max(copynum.get(gene, 0), copy_index)

            # gene row -- note template_distance is emitted UNQUOTED, matching upstream exactly.
            body.append((contig, start, "gene", start, end, strand,
                         f'gene_id "{gid}"; template_allele "{template}"; '
                         f'template_distance {td}; gene_name "{gene}";'))

            # transcript row -- conditional attributes appear only when applicable.
            attrs = [f'gene_id "{gid}"', f'transcript_id "{tid}"', f'gene_name "{gene}"',
                     f'consensus "{consensus}"']
            tied = [template] if rng.random() < 0.75 else [
                template, f"{gene}*{rng.choice(ALLELE_POOL[gene])}"]
            attrs.append(f'alleles "{",".join(tied)}"')
            if warning:
                attrs.append(f'template_warning "{warning}"')
            if cds_dist is not None:
                attrs.append(f"cds_distance {cds_dist}")
            if cds_mut:
                attrs.append(f'cds_mut "{cds_mut}"')
            body.append((contig, start, "transcript", start, end, strand,
                         "; ".join(attrs) + ";"))

            # A few exon/CDS rows so coordinate-based extraction has something to chew on.
            epos = start + 200
            for en in range(1, rng.randint(3, 6)):
                elen = rng.randint(150, 400)
                base = (f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gene}"')
                body.append((contig, epos, "exon", epos, epos + elen, strand,
                             base + f'; exon_number "{en}";'))
                body.append((contig, epos, "CDS", epos, epos + elen, strand, base + ";"))
                epos += elen + rng.randint(200, 900)

            if cds_seq is not None:
                # Fixed per novel variant, so the same novel allele carried by two different
                # people produces byte-identical sequence and clusters into one Table 3 row.
                novel_cds[f"{contig}_{gene}_{copy_index}"] = cds_seq

    # combine.py sorts the final GTF by (contig, start) -- reproduce that, because any consumer
    # that assumes detection order instead of coordinate order is wrong against real data.
    body.sort(key=lambda r: (r[0], r[3]))

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write("## format: gtf\n## date: 2026-08-08\n")
        zero = sorted(g for g in GENE_CLASSES if g not in copynum)
        one = sorted(g for g, c in copynum.items() if c == 1)
        many = sorted(g for g, c in copynum.items() if c > 1)
        f.write(f"## gene (copy num = 0): {','.join(zero)}\n")
        f.write(f"## gene (copy num = 1): {','.join(one)}\n")
        if many:
            f.write(f"## gene (copy num > 1): {','.join(many)}\n")
        for gene in sorted(per_gene_contigs):
            f.write(f"## contigs for {gene}: {','.join(sorted(per_gene_contigs[gene]))}\n")
        for contig, _, feat, start, end, strand, attrs in body:
            f.write(f"{contig}\tIPD-IMGT/HLA-V3.55.0\t{feat}\t{start}\t{end}\t.\t"
                    f"{strand}\t.\t{attrs}\n")
    return novel_cds


def build(outroot, n_people):
    rng = random.Random(SEED)
    outroot = os.path.expanduser(outroot)
    os.makedirs(outroot, exist_ok=True)

    core = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
    extra = ["HLA-DRA", "HLA-DMA", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J",
             "MICA", "MICB", "TAP1", "TAP2", "HLA-DQA2"]
    drb_paralogs = ["HLA-DRB3", "HLA-DRB4", "HLA-DRB5"]

    cohort_rows, ancestry_rows = [], []
    people = [f"{1000000 + i}" for i in range(n_people)]

    for i, pid in enumerate(people):
        ancestry = rng.choices(ANCESTRIES, weights=ANCESTRY_WEIGHTS)[0]
        # Person-id directories live under people/, not directly at <outroot> -- matches the real
        # production layout after the 2026-09-04 move (RUNBOOK.md "Step 1c"): ~12,000 top-level
        # person_id directories made the Workbench Jupyter file browser unusably slow, so person
        # dirs were moved one level deeper into people/, leaving only the aggregate .tsv files (and
        # this one named subfolder) at the top level. Aggregate files below (immuannot_cohort_full.tsv,
        # ancestry_preds.tsv, hla_genotypes.tsv, FIXTURE_MANIFEST.json) intentionally stay directly
        # under <outroot> -- only the person-directory placement moves.
        pdir = os.path.join(outroot, "people", pid, "immuannot_output")

        # ~8% of people get NO output at all -- the real run had ~92% "any output".
        if rng.random() < 0.08:
            cohort_rows.append((pid, rng.choice(["revio", "sequel2e"]),
                                "paf_region", ancestry, 0))
            ancestry_rows.append((pid, ancestry))
            continue

        for hap in ("hap1", "hap2"):
            genes = list(core)
            genes += [g for g in extra if rng.random() < 0.7]
            # DRB paralog presence is haplotype-dependent in reality -- model that, since it is a
            # documented source of copy-number confusion in the MHC.
            genes += [g for g in drb_paralogs if rng.random() < 0.35]

            # THE CIS TRAP: genes in the same hap file are NOT necessarily in cis. Any script that
            # pairs DQA1~DQB1 without checking the contig will silently produce wrong haplotypes.
            # Rate calibrated to the REAL production cohort: 00_recon_vm.py measured 80% of
            # haplotypes spanning >1 contig over a 50-person sample (2026-09-03) -- the MHC is
            # genuinely fragmented across assembly contigs far more often than not. An earlier
            # guess of 18% here was badly optimistic and under-exercised the pairing logic.
            if rng.random() < 0.80:
                split = rng.randint(2, max(3, len(genes) - 2))
                contigs = {f"{hap}_ctg_{pid}_a": [], f"{hap}_ctg_{pid}_b": []}
                names = list(contigs)
                for j, g in enumerate(genes):
                    contigs[names[0 if j < split else 1]].append((g, 1, ancestry))
            else:
                contigs = {f"{hap}_ctg_{pid}_a": [(g, 1, ancestry) for g in genes]}

            # Occasional genuine copy number >1 (segmental duplication).
            if rng.random() < 0.06:
                c = rng.choice(list(contigs))
                contigs[c].append((rng.choice(drb_paralogs + ["HLA-DQA2"]), 2, ancestry))

            novel_cds = write_gtf(os.path.join(pdir, f"{hap}.gtf.gz"), contigs, rng)

            # cds.fa.gz -- Immuannot writes observed CDS for EVERY detected gene copy, in coding
            # orientation, and its own cleanup does NOT delete this file. This is what makes
            # novel-allele sequence work possible with no re-run.
            os.makedirs(os.path.join(pdir, hap), exist_ok=True)
            with gzip.open(os.path.join(pdir, hap, "cds.fa.gz"), "wt") as f:
                for key, seq in novel_cds.items():
                    f.write(f">{key} intron_bnd=GT..AG\n")
                    for k in range(0, len(seq), 60):
                        f.write(seq[k:k + 60] + "\n")

        cohort_rows.append((pid, rng.choice(["revio", "sequel2e"]), "paf_region", ancestry, 16))
        ancestry_rows.append((pid, ancestry))

    with open(os.path.join(outroot, "immuannot_cohort_full.tsv"), "w") as f:
        f.write("person_id\tplatform\ttrim_tier\tancestry_pred\tn_rows\n")
        for r in cohort_rows:
            r = list(r)
            # A small fraction of real rows have a blank/malformed ancestry_pred -- confirmed on
            # the production cohort (ENVIRONMENT.md quirk #30: ~24 non-lowercase rows out of
            # ~13,252, mostly blank, one a literal stray header string), attributed to a minor
            # checkpoint-concatenation artifact in build_immuannot_cohort.py. This is exactly what
            # crashed 06_figures_structure.py and 07_figures_crosscohort.py at real scale
            # (`pd.NA in list` raises TypeError instead of returning False) -- both now fixed, but
            # a 300-person fixture at the real ~0.18% rate would only average ~0.5 such rows, too
            # rare to reliably re-exercise the fix. Using ~2% here (deliberately higher than
            # reality, for reliable test coverage, not realism) guarantees several every run.
            if rng.random() < 0.02:
                r[3] = ""
            f.write("\t".join(map(str, r)) + "\n")

    # ancestry_preds.tsv mirrors the real AoU schema: a hard label plus a 6-way probability array
    # in AFR/AMR/EAS/EUR/MID/SAS order. Continuous admixture lives here.
    with open(os.path.join(outroot, "ancestry_preds.tsv"), "w") as f:
        f.write("research_id\tancestry_pred\tprobabilities\n")
        for pid, anc in ancestry_rows:
            probs = [rng.random() * 0.12 for _ in ANCESTRIES]
            probs[ANCESTRIES.index(anc)] += rng.uniform(1.2, 3.0)
            tot = sum(probs)
            arr = "[" + ", ".join(f"{p / tot:.4f}" for p in probs) + "]"
            f.write(f"{pid}\t{anc}\t{arr}\n")

    # AoU-native short-read calls: bigger N (superset of the LR cohort), 2-field, unphased.
    #
    # This is NOT random noise, and it must not be -- the headline cross-cohort figure (same-allele
    # frequency in SR vs LR, evidencing short-read reference bias) is only testable if the fixture
    # encodes a bias for it to find. So SR is generated as a DEGRADED COPY of the long-read truth:
    #   - 2-field truncation (real AoU calls are 2-3 field)
    #   - a miscall rate that is ANCESTRY-DEPENDENT, higher for non-European haplotypes, because
    #     short-read HLA callers are trained/benchmarked on European-biased references
    #   - miscalls are biased TOWARD common European alleles, which is the specific direction real
    #     reference bias pushes (the reference allele's frequency gets overestimated)
    # Analyses should therefore recover: SR and LR agree well in EUR, diverge in AFR, and the
    # divergence is directional, not symmetric.
    sr_path = os.path.join(outroot, "hla_genotypes.tsv")
    sr_genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    sr_miscall_rate = {"afr": 0.28, "amr": 0.16, "eas": 0.15, "sas": 0.15,
                       "mid": 0.14, "eur": 0.04}

    def two_field(a):
        return ":".join(a.split(":")[:2])

    def sr_call(gene, ancestry, truth_allele):
        """Truth -> observed short-read call, with ancestry-dependent, EUR-directional error."""
        if truth_allele is not None and rng.random() >= sr_miscall_rate[ancestry]:
            return two_field(truth_allele)
        eur_pool = ANCESTRY_BIAS["eur"].get(gene)
        pool = eur_pool if (eur_pool and rng.random() < 0.7) else ALLELE_POOL[gene]
        return two_field(rng.choice(pool))

    with open(sr_path, "w") as f:
        cols = [f"{g}_{i}" for g in sr_genes for i in (1, 2)]
        f.write("research_id\t" + "\t".join(cols) + "\n")
        for pid, anc in ancestry_rows:
            vals = []
            for g in sr_genes:
                gene = "HLA-" + g
                # Draw the person's true diploid genotype the same ancestry-biased way the
                # long-read side did, then degrade it.
                truth = []
                for _ in (1, 2):
                    biased = ANCESTRY_BIAS.get(anc, {}).get(gene)
                    truth.append(rng.choice(biased) if (biased and rng.random() < 0.6)
                                 else rng.choice(ALLELE_POOL[gene]))
                for t in truth:
                    vals.append(f"{g}*{sr_call(gene, anc, t)}")
            f.write(pid + "\t" + "\t".join(vals) + "\n")
        # SR-only people (no long-read data) -- the real SR cohort is ~40x larger than the LR one.
        for i in range(n_people * 3):
            pid = f"{9000000 + i}"
            anc = rng.choices(ANCESTRIES, weights=ANCESTRY_WEIGHTS)[0]
            vals = []
            for g in sr_genes:
                gene = "HLA-" + g
                for _ in (1, 2):
                    biased = ANCESTRY_BIAS.get(anc, {}).get(gene)
                    t = (rng.choice(biased) if (biased and rng.random() < 0.6)
                         else rng.choice(ALLELE_POOL[gene]))
                    vals.append(f"{g}*{sr_call(gene, anc, t)}")
            f.write(pid + "\t" + "\t".join(vals) + "\n")

    manifest = {"seed": SEED, "n_people": n_people, "outroot": outroot,
                "notes": "Synthetic. Contains: multi-contig haps (~18%), copy_index>1 (~6%), "
                         "no-output people (~8%), ancestry-skewed novel rate (afr .22 / eur .05), "
                         "conditional cds_distance/cds_mut, unquoted template_distance."}
    with open(os.path.join(outroot, "FIXTURE_MANIFEST.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2))
    return outroot


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default="/tmp/hla_fixtures")
    ap.add_argument("-n", "--n-people", type=int, default=200)
    a = ap.parse_args()
    build(a.outroot, a.n_people)
