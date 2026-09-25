#!/usr/bin/env python3
"""S04 WS-D -- KIR-HLA ligand co-occurrence by ancestry (seed-list item: ORCHESTRATOR_HANDOFF.md
sect2 WS-D, "KIR-HLA ligand co-occurrence by ancestry"). Aggregate-only outputs.

## Biology (why this test)

Three KIR inhibitory/activating receptors recognise specific HLA class-I epitope groups, not
individual alleles -- the functional "ligand" is the group, not the two-field type:
  - **KIR2DL1** (inhibitory) recognises the **HLA-C2** group.
  - **KIR2DL2 / KIR2DL3** (inhibitory) recognise the **HLA-C1** group. KIR2DL2 binds C1 (and
    weakly C2) more promiscuously than KIR2DL3 (Moesta et al. 2008 J Immunol; Parham 2005 Nat
    Rev Immunol review).
  - **KIR3DL1** (inhibitory) recognises **HLA-Bw4** (all HLA-B alleles carrying the Bw4 public
    epitope, plus the HLA-A Bw4 alleles A*23/A*24/A*32), with allotype-dependent avidity split
    by the Bw4 **position-80 dimorphism** (80I = higher avidity, 80T = lower) -- Cella et al.
    1994; Gumperz et al. 1995 J Exp Med.
  - **KIR3DS1** (activating, no confirmed HLA ligand by direct binding assay, but epidemiologic
    association specifically with **Bw4-80I**) -- Martin et al. 2002 Nat Genet (HIV progression).

Because a person's genotype fixes *both* halves of a pair in trans (unlike a linked haplotype),
co-carriage of a KIR gene and its ligand across a cohort should, under independence (no epistasis,
no shared population-structure confound), equal the product of the two marginal frequencies. This
script tests observed vs. expected co-carriage per ancestry with a 2x2 contingency table
(Fisher's exact + chi2), an odds ratio, and a person-label permutation baseline. **A real
enrichment or depletion here is NOT expected** under simple population genetics -- HLA (chr6) and
KIR (chr19) segregate independently, so any signal is either (a) genuine linkage disequilibrium
with a third factor (ancestry is the obvious confound -- both loci have ancestry-stratified allele
frequencies, so within-ancestry stratification is the whole point of this script), (b) selection
acting on specific KIR-HLA combinations (a real, published phenomenon -- e.g. HIV/pre-eclampsia
literature), or (c) an assay/miscall artifact. This script cannot distinguish those on its own; it
flags what is worth a closer look.

## Epitope group assignment -- sequence-derived, lookup table as fallback + cross-check

**HLA-C1/C2** is determined by the residue at *mature-protein* (leader-peptide-stripped) position
80, together with position 77 (Colonna et al. 1993 PNAS; Winter & Long 1997 J Immunol): **C1 =
Asn80 with Ser77**, **C2 = Lys80**. Verified against known reference alleles: C*01:02, C*03:04,
C*07:01 are C1 (mature residue 80 = Asn/N); C*02:02, C*04:01, C*05:01, C*06:02 are C2 (mature
residue 80 = Lys/K) -- see `TestClassifyC1C2Seq` in the test file for synthetic-sequence
assertions pinned to exactly this rule.

**HLA-Bw4/Bw6** is determined from the alpha-1-helix residues 77-83 of the mature protein (Gumperz
et al. 1995 J Exp Med; Cella et al. 1994; Parham reviews, e.g. Parham 2005 Nat Rev Immunol): the
**Bw6** reference pattern is Ser77-Asn80-Leu81-Arg82-Gly83; **Bw4** is defined by **Arg83 together
with Ile80 or Thr80** (the Bw4-80I/80T avidity dimorphism that KIR3DL1/3DS1 binding strength
tracks -- 80I = higher avidity). Both sub-types are reported (`n_bw4_80I` / `n_bw4_80T` in
`ligand_lookup_qc.tsv`).

**How the sequence assignment works (`--cds-dir`, default `~/tools/Immuannot_refdata/CDSseq`):**
each `<gene>.fa.gz` (IPD-IMGT/HLA CDSseq, headers like
`>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp`) holds one nucleotide CDS per reference allele. For
each Table 1 allele call this script (1) maps the call to a reference CDS record -- an exact match
at the call's own field resolution, else the alphabetically-first reference allele sharing the
call's two-field prefix (`mapping_level` = `exact` / `two_field_fallback`; recorded per call); (2)
translates the CDS respecting the header's `frame=` offset, stopping at the first in-frame stop
codon; (3) strips the 24-aa class-I leader peptide (HLA-A/B/C signal peptide length) to get the
mature protein; (4) reads residues 77-83 directly. **Novel calls** (Table 1's spliced `new` field
token) cannot be resolved this way -- their protein at 77-83 is not guaranteed to match any
reference allele's -- and are reported `unresolved` (counted, never guessed).

**The two-field lookup tables below (`C1C2_TABLE`, `BW4_B_GROUPS`, `BW4_A_GROUPS`) are now a
fallback (used only when sequence resolution fails: missing CDS file, novel call, or no reference
allele shares the two-field prefix) and a cross-check** (`ligand_seq_vs_lookup_crosscheck.tsv`
compares the sequence-derived label against the lookup-table label for every two-field group where
sequence *did* resolve, at the two-field-group level -- catalogue facts, not per-person data, so
this table is never suppressed, but per instruction is never printed next to a person-level carrier
count). `ligand_lookup_qc.tsv` reports `n_two_field_groups_seq_vs_lookup_{compared,agree,disagree}`
and lists disagreeing two-field group names (no counts) in the crosscheck TSV.

**Known Bw4 exception not applied by default:** most B*15 alleles are Bw6, but a documented subset
(B*15:13, B*15:16, B*15:17, B*15:24 and a few others) are Bw4. The *lookup fallback* treats all
B*15 as Bw6 unless `--b15-as-bw4-list` supplies the exception 2-field alleles; the *sequence path*
is unaffected by this (it reads the actual mature-protein residues, so a true Bw4 B*15 allele is
called correctly whenever its CDS record is available and resolvable).

## KIR receptor definitions

Presence/absence of KIR2DL1, KIR2DL2, KIR2DL3, KIR3DL1, KIR3DS1 per person (either haplotype
carries the gene at all -- copy-number and allele-level avidity differences are out of scope for
this pass), from the same `parse_hap_gtf`/`KIR_GENES` logic as 41_kir_pilot.py /
43_kir_full_aggregate.py.

## Method

Per ancestry (ANCESTRY_ORDER, unrelated subset, same kinship threshold as 41/43): for each of the
5 functional pairs (2DL1xC2, 2DL2xC1, 2DL3xC1, 3DL1xBw4, 3DS1xBw4), build the person-level 2x2
table (KIR present/absent x ligand-carrier present/absent), report counts (disclosure-masked),
odds ratio with a 2x2 Haldane-Anscombe-corrected OR when a cell is 0, Fisher's exact p-value,
chi2 p-value, and a >=20-shuffle permutation null of the OR (KIR presence label shuffled within
ancestry, ligand fixed). Disclosure: any 2x2 cell in [1,19] masks the OR/CI/both p-values for that
row (the raw masked cell counts are still reported as `<20`, never blank, never 0 unless truly 0).

## Running (VM, real mode)

    cd ~/s04 && PYTHONPATH=~/s04:~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 48_kir_hla_ligand_cooccurrence.py \\
        --table1 ~/pipeline_outputs/hla_calls_rich.tsv \\
        --cohort-membership ~/pipeline_outputs/cohort_membership.tsv \\
        --relatedness-table ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv \\
        --kir-outroot ~/pipeline_outputs_kir \\
        --kir-pilot-script ~/s03/41_kir_pilot.py \\
        --cds-dir ~/tools/Immuannot_refdata/CDSseq \\
        --n-perms 1000 \\
        --out-dir ~/s04/results/48

## Local synthetic dry run

    python3 scripts/hla_popgen/48_kir_hla_ligand_cooccurrence.py --synthetic --n-people 400 \\
        --out-dir /tmp/ligand_synthetic
"""
import argparse
import gzip
import importlib.util
import os
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]

# ---------------------------------------------------------------------------
# Ligand-group lookup tables -- FALLBACK + cross-check only now (see module docstring). Used
# when sequence resolution fails: missing CDS file, novel call, or no reference allele shares the
# called allele's two-field prefix.
# ---------------------------------------------------------------------------
# HLA-C two-field groups -> C1/C2, per the position-80 dimorphism (documented, common
# low-resolution assignment reproduced in KIR-ligand review tables).
C1C2_TABLE = {
    # C2 group (Lys80)
    "02:02": "C2", "04:01": "C2", "05:01": "C2", "06:02": "C2", "07:04": "C2",
    "08:02": "C2", "12:03": "C2", "15:02": "C2", "16:02": "C2", "17:01": "C2", "18:01": "C2",
    # C1 group (Asn80 + Ser77)
    "01:02": "C1", "03:02": "C1", "03:03": "C1", "03:04": "C1", "07:01": "C1", "07:02": "C1",
    "08:01": "C1", "12:02": "C1", "14:02": "C1", "16:01": "C1",
}
BW4_B_GROUPS = {"13", "27", "37", "38", "44", "47", "49", "51", "52", "53", "57", "58", "59",
                "63", "77"}
BW4_A_GROUPS = {"23", "24", "32"}

# ---------------------------------------------------------------------------
# Sequence-derived assignment: translation + class-I leader stripping + residue rules.
# ---------------------------------------------------------------------------
CLASS_I_LEADER_LEN = 24  # HLA-A/B/C signal peptide length.

_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_CODONS, _AMINO_ACIDS))

GENE_CDS_FILENAME = {"HLA-A": "A.fa.gz", "HLA-B": "B.fa.gz", "HLA-C": "C.fa.gz"}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def suppressed(n):
    n = int(n)
    if n == 0:
        return "0"
    return "<%d" % SUPPRESS_BELOW if n < SUPPRESS_BELOW else str(n)


def hla_two_field(consensus):
    if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
        return None
    s = str(consensus).strip()
    if not s or s.upper() in {"UNDETERMINED", "NA", ""}:
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = []
    for f in s.split(":"):
        if f == "":
            continue
        if f.strip().lower() == "new":
            break
        fields.append(f)
    if len(fields) < 2:
        return None
    return ":".join(fields[:2])


def classify_c1c2(two_field):
    """Lookup-table fallback (see module docstring) -- NOT the primary classification anymore."""
    return C1C2_TABLE.get(two_field, "unclassified")


def classify_bw4(gene, two_field, b15_as_bw4=frozenset()):
    """Lookup-table fallback (see module docstring) -- NOT the primary classification anymore."""
    if two_field is None:
        return "unclassified"
    group = two_field.split(":", 1)[0]
    if gene == "HLA-B":
        if group == "15":
            return "Bw4" if two_field in b15_as_bw4 else "Bw6"
        return "Bw4" if group in BW4_B_GROUPS else "Bw6"
    if gene == "HLA-A":
        return "Bw4" if group in BW4_A_GROUPS else "not_bw4_locus_A"
    return "unclassified"


def load_kir_module(script_path):
    spec = importlib.util.spec_from_file_location("kir_pilot_48", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# CDS translation + mature-protein residue rules.
# ---------------------------------------------------------------------------
def translate_cds(nt_seq, frame=1):
    """Translate a nucleotide CDS starting at the 1-based reading frame given by the IPD-IMGT/HLA
    CDSseq FASTA header's `frame=` field, stopping at the first in-frame stop codon (or end of
    sequence, whichever comes first). Non-ACGT/ambiguous codons translate to 'X' (never silently
    dropped or guessed)."""
    seq = nt_seq.upper().replace("U", "T")
    start = (frame - 1) if frame in (1, 2, 3) else 0
    protein = []
    for i in range(start, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def mature_protein_from_cds(nt_seq, frame, leader_len=CLASS_I_LEADER_LEN):
    """Translate + strip the class-I leader peptide. None if the translated protein doesn't even
    reach past the leader (malformed/truncated record -- never silently returns a wrong protein)."""
    protein = translate_cds(nt_seq, frame)
    if len(protein) <= leader_len:
        return None
    return protein[leader_len:]


def classify_c1c2_seq(mature_protein):
    """HLA-C ligand group from the mature (leader-stripped) protein, using the Ser77/Asn80 (C1)
    vs Lys80 (C2) dimorphism (Colonna et al. 1993 PNAS; Winter & Long 1997 J Immunol). Verified
    against known reference alleles: C*01:02, C*03:04, C*07:01 are C1 (mature residue 80 = Asn/N,
    with Ser/S at 77); C*02:02, C*04:01, C*05:01, C*06:02 are C2 (mature residue 80 = Lys/K).
    1-based mature-protein numbering (positions 77/80 -> 0-based indices 76/79).
    Returns 'C1' / 'C2' / 'other' (neither canonical pattern -- reported, never silently folded
    into a group) / None if the sequence doesn't reach position 80."""
    if mature_protein is None or len(mature_protein) < 80:
        return None
    r77, r80 = mature_protein[76], mature_protein[79]
    if r80 == "K":
        return "C2"
    if r80 == "N" and r77 == "S":
        return "C1"
    return "other"


def classify_bw4_seq(mature_protein):
    """Bw4/Bw6 public epitope from the mature protein's alpha-1-helix residues 77-83 (Gumperz
    et al. 1995 J Exp Med; Cella et al. 1994; Parham reviews, e.g. Parham 2005 Nat Rev Immunol).
    Bw6 reference pattern: Ser77-Asn80-Leu81-Arg82-Gly83. Bw4 is defined by Arg83 together with
    Ile80 or Thr80 (the Bw4-80I/80T avidity dimorphism KIR3DL1/3DS1 binding strength tracks; 80I
    = higher avidity, 80T = lower -- Cella et al. 1994; Gumperz et al. 1995).
    Returns (label, subtype) with label in {'Bw4','Bw6','other',None} and subtype in
    {'80I','80T',None}; None label if the sequence doesn't reach position 83."""
    if mature_protein is None or len(mature_protein) < 83:
        return None, None
    r77, r80, r81, r82, r83 = (mature_protein[76], mature_protein[79], mature_protein[80],
                               mature_protein[81], mature_protein[82])
    if r83 == "R" and r80 in ("I", "T"):
        return "Bw4", ("80I" if r80 == "I" else "80T")
    if r77 == "S" and r80 == "N" and r81 == "L" and r82 == "R" and r83 == "G":
        return "Bw6", None
    return "other", None


# ---------------------------------------------------------------------------
# CDSseq reference loading + allele-to-reference mapping.
# ---------------------------------------------------------------------------
def parse_cds_header(line):
    """'>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp' -> ('HLA-A*01:01:01:01', 1)."""
    parts = line[1:].strip().split()
    allele = parts[0] if parts else ""
    frame = 1
    for p in parts[1:]:
        if p.startswith("frame="):
            try:
                frame = int(p.split("=", 1)[1])
            except ValueError:
                frame = 1
    return allele, frame


def load_cds_fasta_gz(path):
    """Parse one IPD-IMGT/HLA CDSseq reference FASTA (gzip or plain text): one nucleotide CDS per
    allele. Returns {allele_full_name: (nt_seq, frame)}. Missing file -> {} (callers fall back to
    the lookup table for that gene, never crash)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    opener = gzip.open if path.endswith(".gz") else open
    header, frame, chunks = None, 1, []

    def flush():
        if header:
            out[header] = ("".join(chunks).upper(), frame)

    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header, frame = parse_cds_header(line)
                chunks = []
            else:
                chunks.append(line.strip())
    flush()
    return out


def load_all_cds(cds_dir):
    """{gene: {allele_full_name: (nt_seq, frame)}} for HLA-A/B/C; empty per-gene dict if the file
    is missing or cds_dir is falsy (falls back to the lookup table for that gene)."""
    out = {}
    cds_dir = os.path.expanduser(cds_dir) if cds_dir else ""
    for gene, fname in GENE_CDS_FILENAME.items():
        out[gene] = load_cds_fasta_gz(os.path.join(cds_dir, fname)) if cds_dir else {}
    return out


def two_field_from_full(full_allele_name):
    """'HLA-A*01:01:01:01' -> '01:01' (same field-cleaning rule as hla_two_field())."""
    return hla_two_field(full_allele_name)


def build_cds_two_field_index(cds_index):
    """Group a gene's CDS allele names by two-field prefix: {twofield: [sorted allele names]}.
    The mapping fallback for a called allele whose own field resolution isn't present in the
    reference verbatim uses the alphabetically-first member of the matching group as the
    representative reference sequence -- documented here rather than picked silently."""
    groups = defaultdict(list)
    for allele in cds_index:
        tf = two_field_from_full(allele)
        if tf:
            groups[tf].append(allele)
    return {tf: sorted(v) for tf, v in groups.items()}


def normalize_full_allele(consensus):
    """Table1 'consensus' -> the resolved-fields-only allele string ('HLA-A*01:01:01:01'),
    truncated before any spliced 'new' token (novel-call encoding, SCHEMA.md). None if
    unresolvable (missing/undetermined/no '*')."""
    if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
        return None
    s = str(consensus).strip()
    if not s or s.upper() in {"UNDETERMINED", "NA", ""} or "*" not in s:
        return None
    prefix, rest = s.split("*", 1)
    fields = []
    for f in rest.split(":"):
        if f == "":
            continue
        if f.strip().lower() == "new":
            break
        fields.append(f)
    if not fields:
        return None
    return f"{prefix}*{':'.join(fields)}"


def is_novel_consensus(consensus):
    """True iff consensus has a spliced 'new' token at any field depth (SCHEMA.md novelty
    encoding) -- such a call's protein at 77-83 is not guaranteed to match any reference allele,
    so it cannot be resolved from the reference alone."""
    if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
        return False
    s = str(consensus).strip()
    if "*" not in s:
        return False
    rest = s.split("*", 1)[1]
    return any(f.strip().lower() == "new" for f in rest.split(":"))


def map_called_allele_to_cds(consensus, cds_index, cds_twofield_index):
    """Map a Table1 allele call to a reference CDS record. Returns (mature_protein_or_None,
    mapping_level) with mapping_level in:
      - 'exact'             : the call's resolved fields matched a CDSseq header verbatim.
      - 'two_field_fallback': no exact match -- fell back to the alphabetically-first reference
                              allele sharing the call's two-field prefix.
      - 'unresolved'        : novel call, no CDS file for this gene, or no reference allele shares
                              the two-field prefix at all.
    """
    if not cds_index:
        return None, "unresolved"
    if is_novel_consensus(consensus):
        return None, "unresolved"
    full = normalize_full_allele(consensus)
    if full is None:
        return None, "unresolved"
    if full in cds_index:
        nt_seq, frame = cds_index[full]
        return mature_protein_from_cds(nt_seq, frame), "exact"
    tf = hla_two_field(consensus)
    reps = cds_twofield_index.get(tf) if tf else None
    if reps:
        nt_seq, frame = cds_index[reps[0]]
        return mature_protein_from_cds(nt_seq, frame), "two_field_fallback"
    return None, "unresolved"


def classify_c1c2_with_crosscheck(consensus, cds_index, cds_twofield_index):
    """Sequence-derived C1/C2 first; lookup table as fallback + cross-check. Returns dict:
    label ('C1'/'C2'/'other'/'unclassified'), source ('sequence'/'lookup_fallback'),
    mapping_level, lookup_label (always computed, for the cross-check table even when sequence
    resolves), two_field."""
    tf = hla_two_field(consensus)
    lookup_label = classify_c1c2(tf) if tf else "unclassified"
    mature, level = map_called_allele_to_cds(consensus, cds_index, cds_twofield_index)
    seq_label = classify_c1c2_seq(mature) if mature is not None else None
    if seq_label in ("C1", "C2", "other"):
        return {"label": seq_label, "source": "sequence", "mapping_level": level,
                "lookup_label": lookup_label, "two_field": tf}
    return {"label": lookup_label, "source": "lookup_fallback", "mapping_level": level,
            "lookup_label": lookup_label, "two_field": tf}


def classify_bw4_with_crosscheck(gene, consensus, cds_index, cds_twofield_index,
                                 b15_as_bw4=frozenset()):
    """Sequence-derived Bw4/Bw6 (+80I/80T subtype) first; lookup table as fallback + cross-check.
    Returns dict: label ('Bw4'/'Bw6'/'other'/'unclassified'/'not_bw4_locus_A'), subtype
    ('80I'/'80T'/None), source, mapping_level, lookup_label, two_field."""
    tf = hla_two_field(consensus)
    lookup_label = classify_bw4(gene, tf, b15_as_bw4) if tf else "unclassified"
    mature, level = map_called_allele_to_cds(consensus, cds_index, cds_twofield_index)
    seq_label, seq_subtype = (None, None)
    if mature is not None:
        seq_label, seq_subtype = classify_bw4_seq(mature)
    if seq_label in ("Bw4", "Bw6", "other"):
        return {"label": seq_label, "subtype": seq_subtype, "source": "sequence",
                "mapping_level": level, "lookup_label": lookup_label, "two_field": tf}
    return {"label": lookup_label, "subtype": None, "source": "lookup_fallback",
            "mapping_level": level, "lookup_label": lookup_label, "two_field": tf}


# ---------------------------------------------------------------------------
# Per-person epitope + receptor carriage.
# ---------------------------------------------------------------------------
def build_person_epitopes(table1_df, b15_as_bw4=frozenset(), cds_by_gene=None):
    """Returns (epitopes_df, qc, crosscheck_rows).
    epitopes_df: indexed by person_id, columns has_C1, has_C2, has_Bw4 (sequence-derived when a
    reference CDS record resolves, lookup-table fallback otherwise).
    qc: dict of QC counters (n_unclassified_{C,Bw}_allele_calls, n_other_{C,Bw}_seq_pattern,
    n_{C,Bw}_source_{sequence,lookup_fallback}, n_bw4_80I, n_bw4_80T).
    crosscheck_rows: one row per distinct (gene, two_field) allele group where sequence actually
    resolved, comparing the sequence-derived label against the lookup-table label -- aggregate/
    catalogue-level (allele-group names, not person data), deduped so no per-person counts are
    attached."""
    cds_by_gene = cds_by_gene or {}
    cds_idx = {g: cds_by_gene.get(g, {}) for g in ("HLA-A", "HLA-B", "HLA-C")}
    cds_2f = {g: build_cds_two_field_index(cds_idx[g]) for g in cds_idx}

    df = table1_df[table1_df["gene"].isin(["HLA-B", "HLA-C", "HLA-A"])].copy()

    people = sorted(table1_df["person_id"].unique())
    has_c1 = pd.Series(False, index=people)
    has_c2 = pd.Series(False, index=people)
    has_bw4 = pd.Series(False, index=people)

    qc = defaultdict(int)
    # Always-present metrics (0 is a real, reportable value here -- never omit the key just
    # because this run/synthetic-cohort happened not to hit that path).
    for _k in ("n_unclassified_C_allele_calls", "n_other_C_seq_pattern",
              "n_C_source_sequence", "n_C_source_lookup_fallback",
              "n_unclassified_Bw_allele_calls", "n_other_Bw_seq_pattern",
              "n_Bw_source_sequence", "n_Bw_source_lookup_fallback",
              "n_bw4_80I", "n_bw4_80T"):
        qc[_k] = 0
    seen_groups = {}  # (gene, two_field) -> last classification dict seen (catalogue-level dedupe)

    c_rows = df[df["gene"] == "HLA-C"]
    for pid, grp in c_rows.groupby("person_id"):
        for consensus in grp["consensus"].dropna():
            r = classify_c1c2_with_crosscheck(consensus, cds_idx["HLA-C"], cds_2f["HLA-C"])
            if r["label"] == "C1":
                has_c1.loc[pid] = True
            elif r["label"] == "C2":
                has_c2.loc[pid] = True
            elif r["label"] == "other":
                qc["n_other_C_seq_pattern"] += 1
            else:
                qc["n_unclassified_C_allele_calls"] += 1
            qc["n_C_source_%s" % r["source"]] += 1
            if r["two_field"]:
                seen_groups[("HLA-C", r["two_field"])] = r

    b_rows = df[df["gene"].isin(["HLA-B", "HLA-A"])]
    for pid, grp in b_rows.groupby("person_id"):
        for gene, consensus in zip(grp["gene"], grp["consensus"]):
            if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
                continue
            r = classify_bw4_with_crosscheck(gene, consensus, cds_idx[gene], cds_2f[gene],
                                             b15_as_bw4)
            if r["label"] == "Bw4":
                has_bw4.loc[pid] = True
                if r["subtype"]:
                    qc["n_bw4_%s" % r["subtype"]] += 1
            elif r["label"] == "unclassified":
                qc["n_unclassified_Bw_allele_calls"] += 1
            elif r["label"] == "other":
                qc["n_other_Bw_seq_pattern"] += 1
            qc["n_Bw_source_%s" % r["source"]] += 1
            if r["two_field"]:
                seen_groups[(gene, r["two_field"])] = r

    crosscheck_rows = []
    for (gene, tf), r in sorted(seen_groups.items()):
        if r["source"] != "sequence":
            continue  # cross-check is only meaningful where sequence actually resolved
        seq_label, lookup_label = r["label"], r["lookup_label"]
        comparable = seq_label != "other" and lookup_label not in ("unclassified",
                                                                    "not_bw4_locus_A")
        crosscheck_rows.append({
            "gene": gene, "two_field": tf, "seq_label": seq_label, "lookup_label": lookup_label,
            "mapping_level": r["mapping_level"],
            "agree": (seq_label == lookup_label) if comparable else "",
        })

    out = pd.DataFrame({"has_C1": has_c1, "has_C2": has_c2, "has_Bw4": has_bw4})
    out.index.name = "person_id"
    return out, dict(qc), crosscheck_rows


def build_kir_receptor_carriage(pids, kir_outroot, kir_mod, receptor_genes):
    presence = pd.DataFrame(False, index=pids, columns=receptor_genes)
    for pid in pids:
        genes_this_person = set()
        for hap in ("hap1", "hap2"):
            gtf_path = os.path.join(kir_outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
            if not os.path.exists(gtf_path):
                continue
            try:
                rows = kir_mod.parse_hap_gtf(gtf_path)
            except (OSError, EOFError):
                continue
            genes_this_person.update(r["gene"] for r in rows if r["gene"] in receptor_genes)
        for g in genes_this_person:
            presence.at[pid, g] = True
    presence.index.name = "person_id"
    return presence


# ---------------------------------------------------------------------------
# 2x2 contingency stats.
# ---------------------------------------------------------------------------
def two_by_two_stats(kir_present, ligand_present, n_perms=100, seed=0):
    """kir_present/ligand_present: boolean arrays, same length (one ancestry's people).
    Returns a dict with masked cell counts + OR/CI/p-values (blanked if any cell in [1,19])."""
    a = int(np.sum(kir_present & ligand_present))            # KIR+ ligand+
    b = int(np.sum(kir_present & ~ligand_present))           # KIR+ ligand-
    c = int(np.sum(~kir_present & ligand_present))           # KIR- ligand+
    d = int(np.sum(~kir_present & ~ligand_present))          # KIR- ligand-
    n = a + b + c + d

    cells = [a, b, c, d]
    disclosive = any(0 < x < SUPPRESS_BELOW for x in cells)

    result = {
        "n_a_kir_and_ligand": suppressed(a), "n_b_kir_only": suppressed(b),
        "n_c_ligand_only": suppressed(c), "n_d_neither": suppressed(d), "n_total": n,
    }

    if n == 0 or disclosive:
        result.update({"odds_ratio": "", "or_ci_lo": "", "or_ci_hi": "",
                       "fisher_p": "", "chi2_p": "", "perm_p": "", "n_perms": 0})
        return result

    # Haldane-Anscombe correction only when a zero cell would make the OR undefined.
    aa, bb, cc, dd = a, b, c, d
    if 0 in cells:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds_ratio = (aa * dd) / (bb * cc)

    # Fisher exact (exact, small-sample safe) + chi2 (large-sample cross-check).
    fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]])
    try:
        chi2, chi2_p, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=True)
    except ValueError:
        chi2_p = float("nan")

    # log-OR CI (Woolf's method).
    se_log_or = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    log_or = np.log(odds_ratio)
    ci_lo = float(np.exp(log_or - 1.96 * se_log_or))
    ci_hi = float(np.exp(log_or + 1.96 * se_log_or))

    # Permutation null: shuffle the KIR-presence label within this ancestry, ligand fixed.
    rng = np.random.default_rng(seed)
    perm_ors = []
    for _ in range(n_perms):
        kp = rng.permutation(kir_present)
        pa = int(np.sum(kp & ligand_present))
        pb = int(np.sum(kp & ~ligand_present))
        pc = int(np.sum(~kp & ligand_present))
        pd_ = int(np.sum(~kp & ~ligand_present))
        if 0 in (pa, pb, pc, pd_):
            pa, pb, pc, pd_ = pa + 0.5, pb + 0.5, pc + 0.5, pd_ + 0.5
        perm_ors.append((pa * pd_) / (pb * pc))
    perm_arr = np.array(perm_ors)
    perm_p = float((np.sum(np.abs(np.log(perm_arr)) >= abs(log_or)) + 1) / (len(perm_arr) + 1))

    result.update({
        "odds_ratio": round(float(odds_ratio), 3),
        "or_ci_lo": round(ci_lo, 3), "or_ci_hi": round(ci_hi, 3),
        "fisher_p": float(fisher_p), "chi2_p": float(chi2_p) if chi2_p == chi2_p else "",
        "perm_p": perm_p, "n_perms": n_perms,
    })
    return result


PAIRS = [
    ("KIR2DL1", "has_C2", "C2"),
    ("KIR2DL2", "has_C1", "C1"),
    ("KIR2DL3", "has_C1", "C1"),
    ("KIR3DL1", "has_Bw4", "Bw4"),
    ("KIR3DS1", "has_Bw4", "Bw4"),
]
RECEPTOR_GENES = sorted({p[0] for p in PAIRS})


def run_pipeline(table1_df, epitopes_df, kir_presence_df, ancestry_of, out_dir, n_perms,
                 b15_as_bw4, status_path, qc=None, crosscheck_rows=None):
    os.makedirs(out_dir, exist_ok=True)
    qc = qc or {}
    crosscheck_rows = crosscheck_rows or []
    people = sorted(set(epitopes_df.index) & set(kir_presence_df.index))
    epitopes_df = epitopes_df.loc[people]
    kir_presence_df = kir_presence_df.loc[people]
    anc = np.array([ancestry_of.get(p, "UNASSIGNED") for p in people])

    rows = []
    for ancestry in ANCESTRY_ORDER:
        mask = anc == ancestry
        n_people_anc = int(mask.sum())
        if n_people_anc < SUPPRESS_BELOW:
            log(f"[48] ancestry={ancestry}: skipped (n={n_people_anc} < 20)")
            continue
        for kir_gene, ligand_col, ligand_name in PAIRS:
            kir_present = kir_presence_df.loc[mask, kir_gene].values
            ligand_present = epitopes_df.loc[mask, ligand_col].values
            r = two_by_two_stats(kir_present, ligand_present, n_perms=n_perms)
            r["ancestry"] = ancestry
            r["kir_gene"] = kir_gene
            r["ligand"] = ligand_name
            rows.append(r)

    cols = ["ancestry", "kir_gene", "ligand", "n_total", "n_a_kir_and_ligand", "n_b_kir_only",
            "n_c_ligand_only", "n_d_neither", "odds_ratio", "or_ci_lo", "or_ci_hi",
            "fisher_p", "chi2_p", "perm_p", "n_perms"]
    pd.DataFrame(rows)[cols].to_csv(
        os.path.join(out_dir, "kir_hla_ligand_cooccurrence.tsv"), sep="\t", index=False)

    # Epitope frequencies per ancestry (Cheap, feeds the ligand analysis directly).
    freq_rows = []
    for ancestry in ANCESTRY_ORDER:
        mask = anc == ancestry
        n_people_anc = int(mask.sum())
        if n_people_anc < SUPPRESS_BELOW:
            continue
        for col, label in (("has_C1", "C1"), ("has_C2", "C2"), ("has_Bw4", "Bw4")):
            n_carriers = int(epitopes_df.loc[mask, col].sum())
            n_str, d_str = suppressed(n_carriers), n_people_anc
            pct = (round(100.0 * n_carriers / n_people_anc, 1)
                   if not (0 < n_carriers < SUPPRESS_BELOW) else "")
            freq_rows.append({"ancestry": ancestry, "epitope": label,
                              "n_carriers": n_str, "n_total": d_str, "pct": pct})
    pd.DataFrame(freq_rows).to_csv(
        os.path.join(out_dir, "epitope_freq_by_ancestry.tsv"), sep="\t", index=False)

    # QC: unclassified-allele coverage + sequence-vs-lookup source mix.
    n_agree = sum(1 for r in crosscheck_rows if r["agree"] is True)
    n_disagree = sum(1 for r in crosscheck_rows if r["agree"] is False)
    qc_rows = [{"metric": k, "value": suppressed(v)} for k, v in sorted(qc.items())]
    qc_rows += [
        {"metric": "n_people_in_analysis", "value": len(people)},
        {"metric": "c1c2_table_size", "value": len(C1C2_TABLE)},
        {"metric": "bw4_b_groups_size", "value": len(BW4_B_GROUPS)},
        {"metric": "b15_as_bw4_exceptions_supplied", "value": len(b15_as_bw4)},
        {"metric": "n_two_field_groups_seq_vs_lookup_compared", "value": n_agree + n_disagree},
        {"metric": "n_two_field_groups_seq_vs_lookup_agree", "value": n_agree},
        {"metric": "n_two_field_groups_seq_vs_lookup_disagree", "value": n_disagree},
    ]
    pd.DataFrame(qc_rows).to_csv(os.path.join(out_dir, "ligand_lookup_qc.tsv"),
                                 sep="\t", index=False)

    # Cross-check table: catalogue-level (allele two-field group names), never paired with a
    # person-level carrier count -- see module docstring.
    pd.DataFrame(crosscheck_rows,
                columns=["gene", "two_field", "seq_label", "lookup_label", "mapping_level",
                         "agree"]).to_csv(
        os.path.join(out_dir, "ligand_seq_vs_lookup_crosscheck.tsv"), sep="\t", index=False)

    with open(status_path, "w") as f:
        f.write(f"done n_people={len(people)} n_pairs_tested={len(rows)} "
                f"n_unclassified_C={qc.get('n_unclassified_C_allele_calls', 0)} "
                f"n_unclassified_Bw={qc.get('n_unclassified_Bw_allele_calls', 0)}\n")
    log(f"[48] wrote {len(rows)} pair x ancestry rows -> {out_dir}")


# ---------------------------------------------------------------------------
# Synthetic cohort generator.
# ---------------------------------------------------------------------------
def make_synthetic(n_people=400, seed=20260926, planted_enrichment=True, cds_by_gene=None):
    rng = np.random.default_rng(seed)
    pids = [f"synthP{i:05d}" for i in range(n_people)]
    ancestries = rng.choice(ANCESTRY_ORDER, size=n_people, p=[0.2, 0.15, 0.15, 0.3, 0.1, 0.1])
    ancestry_of = dict(zip(pids, ancestries))

    c_alleles = list(C1C2_TABLE.keys())
    b_bw4 = ["57:01", "51:01"]  # Bw4
    b_bw6 = ["07:02", "08:01"]  # Bw6
    table1_rows = []
    for pid in pids:
        for hap in ("hap1", "hap2"):
            c_allele = rng.choice(c_alleles)
            table1_rows.append({"person_id": pid, "hap": hap, "gene": "HLA-C",
                                "consensus": f"HLA-C*{c_allele}"})
            b_allele = rng.choice(b_bw4 if rng.random() < 0.35 else b_bw6)
            table1_rows.append({"person_id": pid, "hap": hap, "gene": "HLA-B",
                                "consensus": f"HLA-B*{b_allele}"})
    table1_df = pd.DataFrame(table1_rows)

    epitopes_df, qc, crosscheck_rows = build_person_epitopes(table1_df, cds_by_gene=cds_by_gene)

    kir_presence = pd.DataFrame(False, index=pids, columns=RECEPTOR_GENES)
    for pid in pids:
        has_c2 = epitopes_df.at[pid, "has_C2"]
        has_bw4 = epitopes_df.at[pid, "has_Bw4"]
        for gene in RECEPTOR_GENES:
            base_p = 0.5
            if planted_enrichment:
                if gene == "KIR2DL1" and has_c2:
                    base_p = 0.8  # planted enrichment: 2DL1 co-occurs with its ligand C2
                if gene == "KIR3DL1" and has_bw4:
                    base_p = 0.8
            kir_presence.at[pid, gene] = rng.random() < base_p
    kir_presence.index.name = "person_id"

    return table1_df, epitopes_df, kir_presence, ancestry_of, qc, crosscheck_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--cohort-membership",
                    default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table",
                    default=os.path.expanduser(
                        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/"
                        "aux/relatedness/samples_relatedness.tsv"))
    ap.add_argument("--kir-outroot", default=os.path.expanduser("~/pipeline_outputs_kir"))
    ap.add_argument("--kir-pilot-script", default=os.path.expanduser("~/s03/41_kir_pilot.py"))
    ap.add_argument("--cds-dir", default=os.path.expanduser("~/tools/Immuannot_refdata/CDSseq"),
                    help="Directory with IPD-IMGT/HLA CDSseq <gene>.fa.gz reference files "
                         "(A.fa.gz, B.fa.gz, C.fa.gz). Used for sequence-derived C1/C2/Bw4/Bw6 "
                         "assignment; missing files fall back to the two-field lookup table.")
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s04/results/48"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--n-perms", type=int, default=100)
    ap.add_argument("--b15-as-bw4-list", default="",
                    help="Comma-separated 2-field B*15 alleles to treat as Bw4 in the LOOKUP "
                         "FALLBACK only (documented exceptions to the default B*15=Bw6 rule; "
                         "the sequence path is unaffected -- it reads the real residues). Empty "
                         "until a VM operator supplies the authoritative IPD-IMGT/HLA exception "
                         "list.")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--n-people", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20260926)
    args = ap.parse_args()

    t0 = time.perf_counter()
    os.makedirs(args.out_dir, exist_ok=True)
    status_path = os.path.join(args.out_dir, "STATUS.txt")
    b15_as_bw4 = frozenset(x.strip() for x in args.b15_as_bw4_list.split(",") if x.strip())

    if args.synthetic:
        (table1_df, epitopes_df, kir_presence, ancestry_of,
         qc, crosscheck_rows) = make_synthetic(args.n_people, args.seed)
        run_pipeline(table1_df, epitopes_df, kir_presence, ancestry_of, args.out_dir,
                    args.n_perms, b15_as_bw4, status_path, qc, crosscheck_rows)
        log(f"[48] synthetic run done in {time.perf_counter()-t0:.0f}s")
        return

    from _viz_common import load_table1, load_cohort_membership
    table1_df = load_table1(args.table1)
    if args.limit:
        keep = sorted(table1_df["person_id"].unique())[:args.limit]
        table1_df = table1_df[table1_df["person_id"].isin(keep)]
    cohort_df = load_cohort_membership(args.cohort_membership)
    ancestry_of = dict(zip(cohort_df["person_id"].astype(str),
                          cohort_df["ancestry_pred"].astype(str).str.upper()))

    kir_mod = load_kir_module(args.kir_pilot_script)
    pairs = kir_mod.load_relatedness_pairs(args.relatedness_table)
    all_pids = sorted(set(table1_df["person_id"]) | set(cohort_df["person_id"].astype(str)))
    kept, removed = kir_mod.greedy_unrelated(all_pids, pairs, kin_min=kir_mod.KIN_MIN)
    log(f"[48] unrelated subset: {len(kept)}/{len(all_pids)} ({len(removed)} relatives dropped)")
    table1_df = table1_df[table1_df["person_id"].isin(kept)]

    cds_by_gene = load_all_cds(args.cds_dir)
    for gene, idx in cds_by_gene.items():
        log(f"[48] CDS reference loaded for {gene}: {len(idx)} alleles"
            if idx else f"[48] CDS reference NOT FOUND for {gene} (--cds-dir={args.cds_dir}) "
                        f"-- falling back to the two-field lookup table for this gene")

    epitopes_df, qc, crosscheck_rows = build_person_epitopes(table1_df, b15_as_bw4, cds_by_gene)
    log(f"[48] unclassified HLA-C calls: {qc.get('n_unclassified_C_allele_calls', 0)}; "
        f"unclassified Bw calls: {qc.get('n_unclassified_Bw_allele_calls', 0)}")
    unrelated_pids = sorted(kept & set(table1_df["person_id"]))
    kir_presence = build_kir_receptor_carriage(unrelated_pids, args.kir_outroot, kir_mod,
                                               RECEPTOR_GENES)

    run_pipeline(table1_df, epitopes_df, kir_presence, ancestry_of, args.out_dir, args.n_perms,
                b15_as_bw4, status_path, qc, crosscheck_rows)
    log(f"[48] done in {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
