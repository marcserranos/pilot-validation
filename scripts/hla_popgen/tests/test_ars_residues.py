#!/usr/bin/env python3
"""Unit tests for _ars_residues.py.

The whole point of that module is "peptide-contact residues come from measuring a structure, not
from copying a table nobody could verify." So the tests check the two places that claim could go
wrong silently:

  1. The distance geometry: a residue genuinely within the cutoff of the peptide must be called a
     contact, and one clearly outside it must not be -- on tiny synthetic PDB text we build by
     hand, so the expected answer is known exactly.
  2. The alignment transfer: given a structure contact position, mapping it onto OUR sequence must
     survive (a) a leading offset (the signal-peptide case this whole alignment step exists to
     avoid guessing at) and (b) an internal gap, without silently shifting or fabricating a match.

Run: python3 scripts/hla_popgen/tests/test_ars_residues.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import importlib.util


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load_module("_ars_residues.py", "ars_residues")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def _atom_line(serial, name, resname, chain, resseq, x, y, z, icode=" "):
    # Fixed-column PDB ATOM record, per the format spec columns used by parse_pdb_text.
    return (
        f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}{icode}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00           {name[0]:>2}"
    )


# ---------------------------------------------------------------------------
# PDB parsing + distance contacts
# ---------------------------------------------------------------------------
def _tiny_pdb():
    lines = [
        _atom_line(1, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "ALA", "A", 2, 20.0, 20.0, 20.0),
        _atom_line(3, "CA", "SER", "P", 1, 2.0, 2.0, 2.0),
    ]
    return "\n".join(lines) + "\n"


def test_parse_pdb_extracts_atoms():
    parsed = m.parse_pdb_text(_tiny_pdb())
    check("three ATOM lines parsed", len(parsed["atoms"]) == 3, str(parsed["atoms"]))
    chains = {a["chain"] for a in parsed["atoms"]}
    check("two chains present", chains == {"A", "P"}, str(chains))


def test_contact_within_cutoff_is_found():
    parsed = m.parse_pdb_text(_tiny_pdb())
    contacts = m.find_contacts(parsed["atoms"], mhc_chains=["A"], peptide_chain="P", cutoff=4.5)
    resnums = {c["resseq"] for c in contacts}
    check("residue 1 (distance ~3.46A) is called a contact", 1 in resnums, str(resnums))
    check("only one contact found at this cutoff", resnums == {1}, str(resnums))
    c1 = [c for c in contacts if c["resseq"] == 1][0]
    check("the amino acid is reported correctly", c1["aa"] == "G", str(c1))


def test_contact_outside_cutoff_is_excluded():
    parsed = m.parse_pdb_text(_tiny_pdb())
    contacts = m.find_contacts(parsed["atoms"], mhc_chains=["A"], peptide_chain="P", cutoff=4.5)
    resnums = {c["resseq"] for c in contacts}
    check("residue 2 (distance ~31A) is NOT called a contact", 2 not in resnums, str(resnums))


def test_cutoff_is_sensitive_as_expected():
    parsed = m.parse_pdb_text(_tiny_pdb())
    tight = m.find_contacts(parsed["atoms"], mhc_chains=["A"], peptide_chain="P", cutoff=3.0)
    check("tightening the cutoff below the true distance drops the contact",
          len(tight) == 0, str(tight))
    loose = m.find_contacts(parsed["atoms"], mhc_chains=["A"], peptide_chain="P", cutoff=40.0)
    check("loosening the cutoff picks up both residues",
          {c["resseq"] for c in loose} == {1, 2}, str(loose))


def test_chain_sequence_from_atoms_fallback():
    parsed = m.parse_pdb_text(_tiny_pdb())
    seq, resnums = m.chain_sequence(parsed, "A")
    check("sequence built from ATOM records when SEQRES is absent", seq == "GA", seq)
    check("resnums line up with the sequence", resnums == [(1, ""), (2, "")], str(resnums))


def test_chain_sequence_prefers_seqres():
    pdb_text = _tiny_pdb() + "SEQRES   1 A    2  GLY ALA\n"
    parsed = m.parse_pdb_text(pdb_text)
    check("SEQRES parsed for chain A", parsed["seqres"].get("A") == ["GLY", "ALA"],
          str(parsed["seqres"]))
    seq, resnums = m.chain_sequence(parsed, "A")
    check("sequence matches SEQRES order", seq == "GA", seq)
    check("resnums still map back to author numbers", resnums == [(1, ""), (2, "")], str(resnums))


# ---------------------------------------------------------------------------
# alignment transfer
# ---------------------------------------------------------------------------
def test_needleman_wunsch_identical_sequences():
    a, b = m.needleman_wunsch("GATTACA", "GATTACA")
    check("identical sequences align with no gaps", a == b == "GATTACA", f"{a} / {b}")


def test_transfer_with_leading_offset():
    """Simulates a signal peptide: the structure's own chain sequence starts right at the mature
    protein, but OUR sequence (query) has an extra leading stretch. A struct contact at index 0
    must land at query index 2, not query index 0 -- getting this backwards is exactly the
    off-by-signal-peptide-length bug this alignment step exists to prevent."""
    query_seq = "XXGA"      # 2-residue leader + mature "GA"
    struct_seq = "GA"       # structure resolved starting at the mature sequence
    transferred = m.align_and_transfer(query_seq, struct_seq, {0})
    check("a leading offset shifts the transferred index by the offset length",
          transferred == {2}, str(transferred))
    transferred2 = m.align_and_transfer(query_seq, struct_seq, {1})
    check("the second structure residue also shifts correctly",
          transferred2 == {3}, str(transferred2))


def test_transfer_with_internal_gap():
    """struct_seq has an extra residue ('A') relative to query_seq. A contact on that extra
    residue has no corresponding query position and must be dropped, not misassigned to a
    neighbour; a contact on a residue that IS shared must still transfer correctly."""
    query_seq = "GT"
    struct_seq = "GAT"
    aligned_q, aligned_s = m.needleman_wunsch(query_seq, struct_seq)
    check("alignment places a gap in the query opposite the extra residue",
          "-" in aligned_q, f"{aligned_q} / {aligned_s}")
    transferred = m.align_and_transfer(query_seq, struct_seq, {1, 2})
    check("the gapped residue contributes no extra transferred index "
          "(only the shared trailing residue's does)",
          len(transferred) == 1, str(transferred))
    check("the contact on the shared trailing residue transfers to query index 1",
          transferred == {1}, str(transferred))


def test_align_and_transfer_ignores_out_of_range_indices():
    transferred = m.align_and_transfer("GA", "GA", {5})
    check("an index with no alignment counterpart is simply dropped, not an error",
          transferred == set(), str(transferred))


# ---------------------------------------------------------------------------
# end-to-end glue
# ---------------------------------------------------------------------------
def test_compute_contact_rows_shape():
    rows = m.compute_contact_rows(_tiny_pdb(), gene="TEST", pdb_id="0XYZ", mhc_chains=["A"],
                                   peptide_chain="P", cutoff=4.5, source="unit test")
    check("one contact row produced", len(rows) == 1, str(rows))
    row = rows[0]
    check("all TSV columns present", set(m.TSV_COLUMNS) <= set(row.keys()), str(row))
    check("author_resnum matches the geometric contact", row["author_resnum"] == 1, str(row))
    check("aa is the one-letter code", row["aa"] == "G", str(row))


def main():
    print("test_ars_residues.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
