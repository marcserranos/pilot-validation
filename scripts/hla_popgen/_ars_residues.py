#!/usr/bin/env python3
"""Peptide-contact (ARS) residues, derived from crystal structures rather than copied from a
published table.

Why this exists
---------------
The classical HLA "antigen recognition site" residue lists (Hughes & Nei 1988 Table 1; Parham
1988; Bondinas 2007) are cited constantly in the HLA population-genetics literature, but this
project could not verify any of the three from a primary source: Hughes & Nei's Table 1 numbering
predates a stable IMGT/HLA numbering convention, Parham 1988 is a review whose residue list is not
reproduced in its abstract/methods, and Bondinas 2007's supplementary table was not accessible from
here. Rather than propagate an unverifiable table, this module recomputes peptide contacts directly
from solved crystal structures: any MHC heavy atom within a distance cutoff of any peptide heavy
atom is a peptide-contact residue, full stop. That is auditable and reproducible from PDB files
alone.

What this module does NOT try to do
------------------------------------
It does not use Biopython (repo convention: avoid the dependency for a well-scoped parsing job) and
it does not implement a general PDB parser -- only the ATOM/SEQRES records needed here. It does not
attempt anything beyond a simple global (Needleman-Wunsch) alignment for transferring structure
contacts onto our own translated sequences -- see `align_and_transfer` below for why that is the
right amount of sophistication for >80%-identical sequences.

Layout
------
  parse_pdb_text        -- ATOM/SEQRES parser (no third-party deps)
  chain_sequence        -- one-letter sequence + author resnums for one chain
  find_contacts         -- distance-based contact residues between an MHC chain set and a peptide
  needleman_wunsch       -- simple global alignment (match +1 / mismatch -1 / gap -2)
  align_and_transfer    -- maps structure contact positions onto 0-based indices of our own sequence
  compute_contact_rows  -- glues the above into TSV-ready rows for reference/ars_peptide_contacts.tsv

CLI usage (one structure/chain-set/cutoff per call; append rows to a TSV):

    python3 scripts/hla_popgen/_ars_residues.py \\
        --pdb 1hhk.pdb --gene HLA-A --pdb-id 1HHK --mhc-chains A --peptide-chain C \\
        --cutoff 4.5 --source "RCSB 1HHK, this project's contact recomputation" \\
        --out reference/ars_peptide_contacts.tsv
"""
import argparse
import os
import sys

# ---------------------------------------------------------------------------
# amino acid tables
# ---------------------------------------------------------------------------
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",  # selenomethionine, common in crystal structures
}


def three_to_one(resname):
    return THREE_TO_ONE.get(resname.strip().upper(), "X")


# ---------------------------------------------------------------------------
# PDB parsing (ATOM + SEQRES only; fixed-column per the PDB format spec)
# ---------------------------------------------------------------------------
def parse_pdb_text(text):
    """Parse plain-text PDB content into {'atoms': [...], 'seqres': {chain: [resname,...]}}.

    Only ATOM records are kept (HETATM -- waters, ions, modified residues -- are deliberately
    excluded; a bound peptide in these reference structures is always given as ATOM records).
    Alternate conformers other than blank/'A' are skipped so a residue isn't double-counted.
    """
    atoms = []
    seqres = {}
    for line in text.splitlines():
        rec = line[0:6].strip() if len(line) >= 6 else line.strip()
        if rec == "ATOM":
            if len(line) < 54:
                continue
            try:
                name = line[12:16].strip()
                altloc = line[16]
                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = int(line[22:26])
                icode = line[26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            if altloc not in (" ", "A", ""):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                stripped = name.lstrip("0123456789")
                element = stripped[0] if stripped else name[0]
            atoms.append({
                "chain": chain, "resseq": resseq, "icode": icode, "resname": resname,
                "name": name, "element": element.upper(), "x": x, "y": y, "z": z,
            })
        elif rec == "SEQRES":
            parts = line.split()
            if len(parts) >= 5:
                chain = parts[2]
                seqres.setdefault(chain, []).extend(parts[4:])
    return {"atoms": atoms, "seqres": seqres}


def _atom_residues_in_order(atoms, chain):
    """Unique (resseq, icode) residues for one chain, in file order."""
    out = []
    seen = set()
    for a in atoms:
        if a["chain"] != chain:
            continue
        key = (a["resseq"], a["icode"])
        if key not in seen:
            seen.add(key)
            out.append({"resseq": a["resseq"], "icode": a["icode"], "resname": a["resname"]})
    return out


def chain_sequence(parsed, chain):
    """One-letter sequence for `chain`, plus a parallel list of (resseq, icode) author numbers.

    Prefers SEQRES (the full construct, including residues never resolved in the density) and
    walks it in lockstep with the ATOM-derived residue order, assigning an author resnum wherever
    the SEQRES entry's resname matches the next unclaimed ATOM residue and `None` where it does
    not (residue present in the construct but not resolved -- it cannot be a reported contact
    since it has no coordinates anyway). Falls back to the ATOM-derived sequence entirely when the
    chain has no SEQRES record.
    """
    atom_residues = _atom_residues_in_order(parsed["atoms"], chain)
    seqres = parsed.get("seqres", {}).get(chain)
    if not seqres:
        seq = "".join(three_to_one(r["resname"]) for r in atom_residues)
        resnums = [(r["resseq"], r["icode"]) for r in atom_residues]
        return seq, resnums

    seq_chars = []
    resnums = []
    ai = 0
    for resname in seqres:
        seq_chars.append(three_to_one(resname))
        if ai < len(atom_residues) and atom_residues[ai]["resname"] == resname:
            resnums.append((atom_residues[ai]["resseq"], atom_residues[ai]["icode"]))
            ai += 1
        else:
            resnums.append(None)
    return "".join(seq_chars), resnums


# ---------------------------------------------------------------------------
# contacts
# ---------------------------------------------------------------------------
def _group_residues(atoms, chains):
    groups = {}
    order = []
    for a in atoms:
        if a["chain"] not in chains or a["element"] == "H":
            continue
        key = (a["chain"], a["resseq"], a["icode"])
        if key not in groups:
            groups[key] = {"resname": a["resname"], "coords": []}
            order.append(key)
        groups[key]["coords"].append((a["x"], a["y"], a["z"]))
    return groups, order


def find_contacts(atoms, mhc_chains, peptide_chain, cutoff=4.5):
    """Every MHC residue (chain in `mhc_chains`) with a heavy atom within `cutoff` angstrom of any
    peptide (chain `peptide_chain`) heavy atom. Returns a list of dicts sorted by (chain, resseq,
    icode), each with resname/aa/author resnum/the observed minimum distance.
    """
    mhc_chains = set(mhc_chains) if not isinstance(mhc_chains, str) else {mhc_chains}
    mhc_groups, mhc_order = _group_residues(atoms, mhc_chains)
    pep_groups, _ = _group_residues(atoms, {peptide_chain})
    pep_coords = [xyz for g in pep_groups.values() for xyz in g["coords"]]

    cutoff2 = cutoff * cutoff
    contacts = []
    for key in mhc_order:
        g = mhc_groups[key]
        mind2 = None
        for (x1, y1, z1) in g["coords"]:
            for (x2, y2, z2) in pep_coords:
                d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                if mind2 is None or d2 < mind2:
                    mind2 = d2
        if mind2 is not None and mind2 <= cutoff2:
            chain, resseq, icode = key
            contacts.append({
                "chain": chain, "resseq": resseq, "icode": icode,
                "resname": g["resname"], "aa": three_to_one(g["resname"]),
                "mindist": mind2 ** 0.5,
            })
    contacts.sort(key=lambda c: (c["chain"], c["resseq"], c["icode"]))
    return contacts


# ---------------------------------------------------------------------------
# alignment: transfer structure contact positions onto our own sequence
# ---------------------------------------------------------------------------
def needleman_wunsch(seq_a, seq_b, match=1, mismatch=-1, gap=-2):
    """Simple global alignment. Deliberately not BLOSUM-scored: for the >80%-identical sequences
    this is used on (our translated protein vs. the structure's own chain, both real HLA alleles),
    a flat match/mismatch/gap scheme is enough to get the gap placement right, which is all that
    matters for transferring residue positions.
    """
    n, m = len(seq_a), len(seq_b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
    for j in range(1, m + 1):
        score[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            score[i][j] = max(score[i - 1][j - 1] + s, score[i - 1][j] + gap, score[i][j - 1] + gap)

    aligned_a, aligned_b = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            if score[i][j] == score[i - 1][j - 1] + s:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append(seq_b[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and score[i][j] == score[i - 1][j] + gap:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
            continue
        aligned_a.append("-")
        aligned_b.append(seq_b[j - 1])
        j -= 1
    aligned_a.reverse()
    aligned_b.reverse()
    return "".join(aligned_a), "".join(aligned_b)


def align_and_transfer(query_seq, struct_seq, struct_contact_indices):
    """Align `query_seq` (our own translated mature-protein sequence) to `struct_seq` (the
    structure chain's sequence) and map `struct_contact_indices` (0-based indices into
    `struct_seq`) onto 0-based indices into `query_seq`.

    This is the step that removes any need to guess a signal-peptide length: whatever offset or
    internal indel separates the two sequences, the alignment places it correctly. A structure
    contact that falls in a gap relative to the query (no corresponding query residue) is dropped
    rather than guessed at.
    """
    aligned_q, aligned_s = needleman_wunsch(query_seq, struct_seq)
    qi = -1
    si = -1
    struct_to_query = {}
    for qc, sc in zip(aligned_q, aligned_s):
        if qc != "-":
            qi += 1
        if sc != "-":
            si += 1
        if qc != "-" and sc != "-":
            struct_to_query[si] = qi
    return {struct_to_query[idx] for idx in struct_contact_indices if idx in struct_to_query}


# ---------------------------------------------------------------------------
# glue: PDB text -> TSV rows
# ---------------------------------------------------------------------------
TSV_COLUMNS = ["gene", "pdb_id", "mhc_chain", "peptide_chain", "cutoff_angstrom",
               "author_resnum", "aa", "source"]


def compute_contact_rows(pdb_text, gene, pdb_id, mhc_chains, peptide_chain, cutoff, source):
    parsed = parse_pdb_text(pdb_text)
    contacts = find_contacts(parsed["atoms"], mhc_chains, peptide_chain, cutoff)
    rows = []
    for c in contacts:
        rows.append({
            "gene": gene,
            "pdb_id": pdb_id,
            "mhc_chain": c["chain"],
            "peptide_chain": peptide_chain,
            "cutoff_angstrom": cutoff,
            "author_resnum": c["resseq"],
            "aa": c["aa"],
            "source": source,
        })
    return rows


def _write_tsv_rows(rows, out_path):
    file_exists = os.path.exists(out_path) and os.path.getsize(out_path) > 0
    with open(out_path, "a") as fh:
        if not file_exists:
            fh.write("\t".join(TSV_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in TSV_COLUMNS) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="path to a PDB-format coordinate file")
    ap.add_argument("--gene", required=True, help="gene label to record, e.g. HLA-A or DRB1")
    ap.add_argument("--pdb-id", required=True)
    ap.add_argument("--mhc-chains", required=True, help="comma-separated chain id(s), e.g. A")
    ap.add_argument("--peptide-chain", required=True)
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--source", required=True, help="free-text provenance note for the TSV")
    ap.add_argument("--out", required=True, help="TSV path to append rows to")
    args = ap.parse_args()

    with open(args.pdb) as fh:
        pdb_text = fh.read()
    mhc_chains = [c.strip() for c in args.mhc_chains.split(",") if c.strip()]
    rows = compute_contact_rows(pdb_text, args.gene, args.pdb_id, mhc_chains,
                                 args.peptide_chain, args.cutoff, args.source)
    _write_tsv_rows(rows, args.out)
    print(f"{args.pdb_id} chain(s) {mhc_chains} vs peptide {args.peptide_chain} "
          f"@ {args.cutoff}A: {len(rows)} contact residues -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
