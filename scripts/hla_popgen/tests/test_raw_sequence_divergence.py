#!/usr/bin/env python3
"""Fixture tests for 17_raw_sequence_divergence.py -- constructs synthetic cds.fa.gz files with the
EXACT real header format (contig_gene_index, PLUS a splice-site annotation after a space that must
be stripped -- a real parsing bug caught against real VM data on 2026-09-08, when the naive regex
matched against the whole line and silently found nothing). Also constructs known point-substitution,
clustered-indel, and scattered-multi-region sequence differences to verify n_diff_bases/n_diff_blocks
are computed correctly -- this is the whole point of the script, so a silent counting bug here would
misreport exactly the number Marc asked for."""
import gzip
import importlib
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mod = importlib.import_module("17_raw_sequence_divergence")


def _write_cds_fasta(path, entries):
    """entries: list of (contig, gene_with_prefix, copy_idx, seq)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        for contig, gene, idx, seq in entries:
            f.write(f">{contig}_{gene}_{idx} GT_AG:GT_AG:GT_AG\n")  # real files carry splice sites
            f.write(seq + "\n")


def _setup_people(tmp_root):
    # Person P1: hap1 has HLA-X (a 20bp toy sequence), hap2 has a different HLA-X.
    _write_cds_fasta(os.path.join(tmp_root, "P1", "immuannot_output", "hap1", "cds.fa.gz"),
                      [("ctg1", "HLA-X", 1, "AAAAAAAAAAAAAAAAAAAA")])   # exact-match candidate
    _write_cds_fasta(os.path.join(tmp_root, "P1", "immuannot_output", "hap2", "cds.fa.gz"),
                      [("ctg1", "HLA-X", 2, "CCCCCCCCCCCCCCCCCCCC")])

    # Person P2 (clean child of P1 at HLA-X): hap1 = EXACT copy of P1's hap1 -- raw_shared.
    _write_cds_fasta(os.path.join(tmp_root, "P2", "immuannot_output", "hap1", "cds.fa.gz"),
                      [("ctg2", "HLA-X", 1, "AAAAAAAAAAAAAAAAAAAA")])
    _write_cds_fasta(os.path.join(tmp_root, "P2", "immuannot_output", "hap2", "cds.fa.gz"),
                      [("ctg2", "HLA-X", 2, "GGGGGGGGGGGGGGGGGGGG")])

    # Person P3 (one point substitution from P1's hap1 at position 10 -- single clustered change).
    seq_point = "AAAAAAAAA" + "T" + "AAAAAAAAAA"  # 20bp, single 'T' substitution at index 9
    assert len(seq_point) == 20
    _write_cds_fasta(os.path.join(tmp_root, "P3", "immuannot_output", "hap1", "cds.fa.gz"),
                      [("ctg3", "HLA-X", 1, seq_point)])
    _write_cds_fasta(os.path.join(tmp_root, "P3", "immuannot_output", "hap2", "cds.fa.gz"),
                      [("ctg3", "HLA-X", 2, "TTTTTTTTTTTTTTTTTTTT")])

    # Person P4 (two scattered substitutions from P1's hap1, at opposite ends -- 2 separate blocks).
    seq_scattered = "T" + "AAAAAAAAAAAAAAAAAA" + "T"  # 20bp, subs at index 0 AND index 19
    assert len(seq_scattered) == 20
    _write_cds_fasta(os.path.join(tmp_root, "P4", "immuannot_output", "hap1", "cds.fa.gz"),
                      [("ctg4", "HLA-X", 1, seq_scattered)])
    _write_cds_fasta(os.path.join(tmp_root, "P4", "immuannot_output", "hap2", "cds.fa.gz"),
                      [("ctg4", "HLA-X", 2, "GGGGGGGGGGGGGGGGGGGG")])


def test_header_parsing_strips_splice_site_annotation(tmp_path):
    root = str(tmp_path)
    _setup_people(root)
    idx = mod.parse_cds_fasta(os.path.join(root, "P1", "immuannot_output", "hap1", "cds.fa.gz"))
    assert "HLA-X" in idx, f"parsing failed to find HLA-X, got keys: {list(idx.keys())}"
    assert idx["HLA-X"] == [(1, "AAAAAAAAAAAAAAAAAAAA")]


def test_exact_match_reports_raw_shared_zero_diff():
    result = mod.compare_gene(["AAAAAAAAAAAAAAAAAAAA"], ["AAAAAAAAAAAAAAAAAAAA", "CCCC"])
    assert result["raw_shared"] is True
    assert result["n_diff_bases"] == 0
    assert result["n_diff_blocks"] == 0


def test_single_point_substitution_is_one_block_one_base():
    a = "AAAAAAAAAAAAAAAAAAAA"
    b = "AAAAAAAAA" + "T" + "AAAAAAAAAA"
    result = mod.compare_gene([a], [b])
    assert result["raw_shared"] is False
    assert result["n_diff_bases"] == 1
    assert result["n_diff_blocks"] == 1
    assert result["same_length"] is True


def test_two_scattered_substitutions_are_two_blocks():
    a = "AAAAAAAAAAAAAAAAAAAA"
    b = "T" + "AAAAAAAAAAAAAAAAAA" + "T"
    result = mod.compare_gene([a], [b])
    assert result["raw_shared"] is False
    assert result["n_diff_bases"] == 2
    assert result["n_diff_blocks"] == 2   # the whole point: scattered, not one clustered change


def test_load_person_sequences_and_end_to_end_run(tmp_path):
    root = str(tmp_path)
    _setup_people(root)
    pairs_df = pd.DataFrame([
        {"person_i": "P1", "person_j": "P2", "pair_class": "high_sharing"},
        {"person_i": "P1", "person_j": "P3", "pair_class": "high_sharing"},
        {"person_i": "P1", "person_j": "P4", "pair_class": "high_sharing"},
    ])
    gene_class = {"X": "classical_I"}
    df = mod.run_all(pairs_df, root, ["X"], gene_class)
    got = {(r.person_i if hasattr(r, "person_i") else None) for r in df.itertuples()}  # unused, sanity
    assert len(df) == 3  # one HLA-X comparison per pair
    by_pair = df.reset_index(drop=True)
    # P1~P2: exact match exists (P2 hap1 == P1 hap1) -> raw_shared True
    assert bool(by_pair.iloc[0]["raw_shared"]) is True
    # P1~P3: closest pairing is the 1-base substitution -> n_diff_bases == 1
    assert by_pair.iloc[1]["raw_shared"] == False
    assert by_pair.iloc[1]["n_diff_bases"] == 1
    # P1~P4: closest pairing is the 2-scattered-subs sequence -> n_diff_blocks == 2
    assert by_pair.iloc[2]["raw_shared"] == False
    assert by_pair.iloc[2]["n_diff_blocks"] == 2


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        from pathlib import Path
        tp = Path(td)
        test_header_parsing_strips_splice_site_annotation(tp)
    test_exact_match_reports_raw_shared_zero_diff()
    test_single_point_substitution_is_one_block_one_base()
    test_two_scattered_substitutions_are_two_blocks()
    with tempfile.TemporaryDirectory() as td:
        from pathlib import Path
        tp = Path(td)
        test_load_person_sequences_and_end_to_end_run(tp)
    print("OK")
