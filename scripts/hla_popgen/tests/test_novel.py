#!/usr/bin/env python3
"""Real assertions against 03_novel_alleles.py and 04_allele_saturation.py, per SCHEMA.md hard
rule #6 ("every script here must run against the synthetic fixtures before it is handed to Marc").

Three groups of tests:
  1. Unit tests of the Chao2/ACE/jackknife formulas against hand-computed toy examples (arithmetic
     worked out in comments below -- verifiable by any reader without external tooling).
  2. Unit tests of 03's homopolymer-indel QC heuristic and clustering/passes_qc logic against
     hand-built matched-row fixtures (the real fixtures never produce genuine cross-person
     recurrence -- see make_fixtures.py's per-call independent-random-sequence design -- so the
     recurrence gate has to be exercised synthetically here, not via the full pipeline).
  3. An integration test against tests/make_fixtures.py's output, run through 01/02/03, asserting
     the ancestry-skewed novel-rate gradient (afr 0.22 vs eur 0.05) is recovered end to end -- the
     key assertion the task spec calls out explicitly.

Run:
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures_test -n 300
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures_test/people \\
        --out-dir /tmp/hla_fixtures_test --sample
    python3 scripts/hla_popgen/tests/test_novel.py --fixtures /tmp/hla_fixtures_test
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


extract = _load_module("01_extract_rich.py", "hla_popgen_extract_rich")
cohorts = _load_module("02_build_cohorts.py", "hla_popgen_build_cohorts")
novel = _load_module("03_novel_alleles.py", "hla_popgen_novel_alleles")
sat = _load_module("04_allele_saturation.py", "hla_popgen_saturation")

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        msg = f"  FAIL  {name}" + (f" -- {detail}" if detail else "")
        print(msg)
        FAILURES.append(name)


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol


# ---------------------------------------------------------------------------
# 1. Richness estimator formulas -- hand-computed toy examples
# ---------------------------------------------------------------------------
def test_chao2_hand_example():
    # 5 sampling units. Species incidence:
    #   sp1: units {1}          -> occurs in 1 unit  (unique)
    #   sp2: units {1,2}        -> occurs in 2 units  (duplicate)
    #   sp3: units {2,3}        -> occurs in 2 units  (duplicate)
    #   sp4: units {3}          -> occurs in 1 unit   (unique)
    #   sp5: units {4}          -> occurs in 1 unit   (unique)
    # S_obs = 5, m = 5, Q1 = 3 (sp1,sp4,sp5), Q2 = 2 (sp2,sp3).
    # Chao2 = S_obs + (m-1)/m * Q1^2/(2*Q2) = 5 + (4/5)*(9/4) = 5 + 1.8 = 6.8
    unit_sets = [
        {"sp1", "sp2"},
        {"sp2", "sp3"},
        {"sp3"},
        {"sp4"},
        set(),
    ]
    # sp4 above was meant as a distinct singleton in unit 4 -- fix the toy set to match the
    # narrative exactly (unit indices 0..4 correspond to units 1..5 above):
    unit_sets = [
        {"sp1", "sp2"},   # unit 1
        {"sp2", "sp3"},   # unit 2
        {"sp3"},          # unit 3
        {"sp4"},          # unit 4 (renamed from sp4 above to avoid confusion with sp4/sp5 clash)
        {"sp5"},          # unit 5
    ]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    check("hand example: S_obs == 5", S_obs == 5, detail=str(S_obs))
    check("hand example: m == 5", m == 5)
    check("hand example: Q1 == 3", Q.get(1, 0) == 3, detail=str(dict(Q)))
    check("hand example: Q2 == 2", Q.get(2, 0) == 2, detail=str(dict(Q)))
    expected_chao2 = 5 + (4 / 5) * (3 ** 2) / (2 * 2)
    got = sat.chao2(S_obs, Q, m)
    check("Chao2 matches hand calculation (6.8)", close(got, expected_chao2, 1e-9),
          detail=f"got {got}, expected {expected_chao2}")


def test_chao2_bias_corrected_when_q2_zero():
    # 3 units, 3 all-singleton species -> Q1=3, Q2=0. Bias-corrected form:
    # Chao2 = S_obs + (m-1)/m * Q1*(Q1-1)/2 = 3 + (2/3)*(3*2/2) = 3 + 2 = 5
    unit_sets = [{"a"}, {"b"}, {"c"}]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    expected = 3 + (2 / 3) * (3 * 2 / 2)
    got = sat.chao2(S_obs, Q, m)
    check("Chao2 bias-corrected form (Q2=0) matches hand calc (5.0)", close(got, expected, 1e-9),
          detail=f"got {got}, expected {expected}")


def test_jackknife1_hand_example():
    # Same 5-unit example as above: Jack1 = S_obs + Q1*(m-1)/m = 5 + 3*(4/5) = 5 + 2.4 = 7.4
    unit_sets = [{"sp1", "sp2"}, {"sp2", "sp3"}, {"sp3"}, {"sp4"}, {"sp5"}]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    expected = 5 + 3 * (4 / 5)
    got = sat.jackknife1_incidence(S_obs, Q, m)
    check("Jackknife1 matches hand calculation (7.4)", close(got, expected, 1e-9),
          detail=f"got {got}, expected {expected}")


def test_jackknife2_hand_example():
    # Jack2 = S_obs + Q1*(2m-3)/m - Q2*(m-2)^2/(m*(m-1))
    #       = 5 + 3*(10-3)/5 - 2*(3)^2/(5*4) = 5 + 4.2 - 0.9 = 8.3
    unit_sets = [{"sp1", "sp2"}, {"sp2", "sp3"}, {"sp3"}, {"sp4"}, {"sp5"}]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    expected = 5 + 3 * (2 * 5 - 3) / 5 - 2 * ((5 - 2) ** 2) / (5 * 4)
    got = sat.jackknife2_incidence(S_obs, Q, m)
    check("Jackknife2 matches hand calculation (8.3)", close(got, expected, 1e-9),
          detail=f"got {got}, expected {expected}")


def test_ace_reduces_sensibly_no_rare_species():
    # All species occur > kappa(=10) times -> no rare species -> ACE should just return S_obs.
    unit_sets = [{"a", "b"} for _ in range(15)]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    got = sat.ace_incidence(S_obs, Q, m, kappa=10)
    check("ACE with no rare species returns S_obs (2)", close(got, 2.0, 1e-9), detail=str(got))


def test_ewens_watterson_expected_k_monotonic_and_solvable():
    # E[K] must increase monotonically in theta, and solving for theta from a given K must recover
    # a theta whose forward evaluation reproduces K almost exactly (self-consistency check).
    n = 50
    thetas = [0.5, 2.0, 10.0, 50.0]
    ks = [sat.ewens_watterson_expected_k(t, n) for t in thetas]
    check("Ewens E[K] increasing in theta", all(ks[i] < ks[i + 1] for i in range(len(ks) - 1)),
          detail=str(ks))
    k_obs = 12.3
    theta_hat = sat.solve_theta_for_k(k_obs, n)
    k_check = sat.ewens_watterson_expected_k(theta_hat, n)
    check("solve_theta_for_k self-consistent (round trip within 1e-3)",
          close(k_check, k_obs, 1e-3), detail=f"theta={theta_hat}, k_check={k_check}")


def test_chao2_extrapolate_no_signal_returns_sobs():
    # No unseen-mass evidence (Q1=0, every species seen >=2x) -> extrapolation should not invent
    # growth.
    unit_sets = [{"a", "b"}, {"a", "b"}, {"a", "b"}]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    got = sat.chao2_extrapolate(S_obs, Q, m, t=30)
    check("Extrapolation with Q1=0 returns S_obs unchanged", close(got, S_obs, 1e-9),
          detail=str(got))


def test_chao2_extrapolate_grows_with_target():
    unit_sets = [{"sp1", "sp2"}, {"sp2", "sp3"}, {"sp3"}, {"sp4"}, {"sp5"}]
    S_obs, m, Q = sat.incidence_freqs(unit_sets)
    est_2x = sat.chao2_extrapolate(S_obs, Q, m, t=2 * m)
    est_10x = sat.chao2_extrapolate(S_obs, Q, m, t=10 * m)
    c2 = sat.chao2(S_obs, Q, m)
    check("Extrapolation at 2x >= S_obs", est_2x >= S_obs, detail=str(est_2x))
    check("Extrapolation at 10x >= extrapolation at 2x", est_10x >= est_2x,
          detail=f"{est_10x} vs {est_2x}")
    check("Extrapolation never exceeds the asymptotic Chao2 estimate",
          est_10x <= c2 + 1e-6, detail=f"est_10x={est_10x}, chao2={c2}")


# ---------------------------------------------------------------------------
# 2. 03_novel_alleles.py: homopolymer QC heuristic + clustering/passes_qc logic
# ---------------------------------------------------------------------------
def test_homopolymer_indel_detection():
    pure_homopolymer = "HLA-A*01:01:01:01|:100+aaaaaaa|K(AAG)<K(AAG)"
    # a pure indel of a homopolymer run ('aaaaaaa') with no substitution token -- should flag True.
    check("pure homopolymer insertion flagged True",
          novel.is_homopolymer_indel_only(pure_homopolymer) is True)

    mixed_run = "HLA-A*01:01:01:01|:100+acgtacgt|K(AAG)<K(AAG)"
    check("non-homopolymer indel run flagged False", novel.is_homopolymer_indel_only(mixed_run) is False)

    substitution_only = "HLA-A*01:01:01:01|:301*ag|K(AAG)<E(GAG)"
    check("pure substitution not flagged as homopolymer-indel",
          novel.is_homopolymer_indel_only(substitution_only) is False)

    substitution_and_homopolymer_indel = "HLA-A*01:01:01:01|:100+ttttt*ag|K(AAG)<E(GAG)"
    check("mix of substitution + homopolymer indel is NOT homopolymer-indel-ONLY",
          novel.is_homopolymer_indel_only(substitution_and_homopolymer_indel) is False)

    check("None/empty cds_mut never flagged", novel.is_homopolymer_indel_only(None) is False)
    check("no indel/substitution tokens at all -> False",
          novel.is_homopolymer_indel_only("HLA-A*01:01:01:01|:500|") is False)


def test_warning_tokens_na_excluded():
    """Regression (2026-09-03 real-data census, SCHEMA.md 'template_warning policy'): Immuannot
    writes the literal string template_warning "NA" to mean NO WARNING on 57.4% of transcript
    rows. warning_tokens() previously excluded 'nan'/'none' (pandas/Python missing-value
    spellings) but NOT 'NA', so it returned frozenset({'NA'}) for a clean call -- scoring 57.4%
    of real data as warned. A literal 'NA' (any case, with surrounding whitespace, standalone or
    comma-joined with real tokens) must never appear in the returned token set."""
    for na_spelling in ("NA", "na", "Na", " NA ", "", None, float("nan")):
        toks = novel.warning_tokens(na_spelling)
        check(f"warning_tokens({na_spelling!r}) == empty set (NA/absent means clean)",
              toks == frozenset(), detail=str(toks))

    check("warning_tokens('partial_CDS') == {'partial_CDS'}",
          novel.warning_tokens("partial_CDS") == frozenset(["partial_CDS"]))
    check("warning_tokens('NA,partial_CDS') drops the NA token, keeps the real one",
          novel.warning_tokens("NA,partial_CDS") == frozenset(["partial_CDS"]),
          detail=str(novel.warning_tokens("NA,partial_CDS")))
    check("warning_tokens('no-start_codon,no-stop_codon') keeps both real tokens",
          novel.warning_tokens("no-start_codon,no-stop_codon") ==
          frozenset(["no-start_codon", "no-stop_codon"]))


def test_build_table3_recurrence_and_passes_qc():
    """Hand-built matched rows exercising the recurrence gate the real fixtures can't reach (every
    novel CDS in make_fixtures.py is drawn independently at random, so genuine cross-person identity
    never happens by chance -- see module docstring). Three scenarios in one gene:
      - HASH_A: same sequence, 2 different persons, 1 ancestry each (AFR, AFR) -> passes_qc True.
      - HASH_B: same sequence, 1 person, 2 haplotypes (both haps of person P2) -> only 1 person ->
        passes_qc False (recurrence gate requires DIFFERENT persons, not just >1 haplotype).
      - HASH_C: singleton, 1 person, carries the default-disqualifying 'inframe_stop' warning token
        -> passes_qc False on two counts (recurrence AND token-aware warning gate).
    """
    matched_rows = [
        {"person_id": "P1", "hap": "hap1", "gene": "HLA-A", "gene_class": "classical_I",
         "cds_seq_sha1": "HASH_A" + "0" * 34, "cds_len": 1000, "nearest_allele": "HLA-A*01:01",
         "cds_distance": 2, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P2", "hap": "hap1", "gene": "HLA-A", "gene_class": "classical_I",
         "cds_seq_sha1": "HASH_A" + "0" * 34, "cds_len": 1000, "nearest_allele": "HLA-A*01:01",
         "cds_distance": 2, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P3", "hap": "hap1", "gene": "HLA-A", "gene_class": "classical_I",
         "cds_seq_sha1": "HASH_B" + "0" * 34, "cds_len": 900, "nearest_allele": "HLA-A*02:01",
         "cds_distance": 1, "n_aa_changes": 0, "novelty_class": "synonymous",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P3", "hap": "hap2", "gene": "HLA-A", "gene_class": "classical_I",
         "cds_seq_sha1": "HASH_B" + "0" * 34, "cds_len": 900, "nearest_allele": "HLA-A*02:01",
         "cds_distance": 1, "n_aa_changes": 0, "novelty_class": "synonymous",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P4", "hap": "hap1", "gene": "HLA-A", "gene_class": "classical_I",
         "cds_seq_sha1": "HASH_C" + "0" * 34, "cds_len": 800, "nearest_allele": "HLA-A*03:01",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["inframe_stop"]), "is_homopolymer_indel_only": False},
    ]
    ancestry_by_person = {"P1": "AFR", "P2": "AFR", "P3": "EUR", "P4": "AMR"}
    table3 = novel.build_table3(matched_rows, ancestry_by_person)
    check("3 distinct clusters produced", len(table3) == 3, detail=str(len(table3)))

    row_a = table3[table3["cds_seq_sha1"].str.startswith("HASH_A")].iloc[0]
    check("HASH_A: n_persons == 2", row_a["n_persons"] == 2, detail=str(row_a["n_persons"]))
    check("HASH_A: n_haplotypes == 2", row_a["n_haplotypes"] == 2)
    check("HASH_A: passes_qc True (2 unrelated persons, no warning, not homopolymer)",
          bool(row_a["passes_qc"]) is True)

    row_b = table3[table3["cds_seq_sha1"].str.startswith("HASH_B")].iloc[0]
    check("HASH_B: n_persons == 1 despite 2 haplotypes (same person)", row_b["n_persons"] == 1,
          detail=str(row_b["n_persons"]))
    check("HASH_B: n_haplotypes == 2", row_b["n_haplotypes"] == 2)
    check("HASH_B: passes_qc False (only 1 unrelated person -- the primary recurrence gate)",
          bool(row_b["passes_qc"]) is False)

    row_c = table3[table3["cds_seq_sha1"].str.startswith("HASH_C")].iloc[0]
    check("HASH_C: passes_qc False (singleton AND disqualifying 'inframe_stop' warning token)",
          bool(row_c["passes_qc"]) is False)

    check("novel_id is deterministic (sha1-based, recomputable)",
          row_a["novel_id"] == f"HLA-A_nov_{('HASH_A' + '0'*34)[:8]}")


def test_build_table3_token_aware_warning_gate():
    """Regression for Fix 2 (SCHEMA.md 'template_warning policy'): a blanket 'any warning'
    passes_qc gate would reject ~95% of real-data candidates -- and that ~95% figure was itself
    a bug (it counted the attribute's mere presence; Immuannot's literal "NA" token means "no
    warning" on 57.4% of rows, so the TRUE warning rate is ~38%). The gate must be TOKEN-AWARE:
    'no-start_codon'/'no-stop_codon' benign by default (correct pseudogene biology),
    'partial_CDS' and 'inframe_stop' disqualify by default (updated default per the 2026-09-03
    census: partial_CDS is only ~2.4% of classical-gene calls, so excluding it is cheap and a
    truncated CDS can't support a novel-allele claim), and the disqualifying set must be
    caller-configurable (the --disqualifying-warnings CLI flag)."""
    matched_rows = [
        {"person_id": "P1", "hap": "hap1", "gene": "HLA-B", "gene_class": "classical_I",
         "cds_seq_sha1": "BENIGN0" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-B*07:02",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["no-start_codon"]), "is_homopolymer_indel_only": False},
        {"person_id": "P2", "hap": "hap1", "gene": "HLA-B", "gene_class": "classical_I",
         "cds_seq_sha1": "BENIGN0" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-B*07:02",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["no-start_codon", "no-stop_codon"]),
         "is_homopolymer_indel_only": False},
        {"person_id": "P3", "hap": "hap1", "gene": "HLA-B", "gene_class": "classical_I",
         "cds_seq_sha1": "BADONE0" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-B*08:01",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["inframe_stop"]), "is_homopolymer_indel_only": False},
        {"person_id": "P4", "hap": "hap1", "gene": "HLA-B", "gene_class": "classical_I",
         "cds_seq_sha1": "BADONE0" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-B*08:01",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["inframe_stop"]), "is_homopolymer_indel_only": False},
        {"person_id": "P5", "hap": "hap1", "gene": "HLA-B", "gene_class": "classical_I",
         "cds_seq_sha1": "PARTIAL0" + "0" * 32, "cds_len": 1000, "nearest_allele": "HLA-B*09:01",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["partial_CDS"]), "is_homopolymer_indel_only": False},
        {"person_id": "P6", "hap": "hap1", "gene": "HLA-B", "gene_class": "classical_I",
         "cds_seq_sha1": "PARTIAL0" + "0" * 32, "cds_len": 1000, "nearest_allele": "HLA-B*09:01",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(["partial_CDS"]), "is_homopolymer_indel_only": False},
    ]
    ancestry_by_person = {"P1": "AFR", "P2": "EUR", "P3": "AFR", "P4": "EUR",
                           "P5": "AFR", "P6": "EUR"}

    table3_default = novel.build_table3(matched_rows, ancestry_by_person)
    benign_row = table3_default[table3_default["cds_seq_sha1"].str.startswith("BENIGN0")].iloc[0]
    bad_row = table3_default[table3_default["cds_seq_sha1"].str.startswith("BADONE0")].iloc[0]
    partial_row = table3_default[table3_default["cds_seq_sha1"].str.startswith("PARTIAL0")].iloc[0]
    check("default gate: no-start_codon/no-stop_codon do NOT disqualify (2 unrelated persons, "
          "benign tokens only) -> passes_qc True", bool(benign_row["passes_qc"]) is True)
    check("default gate: inframe_stop DOES disqualify even with 2 unrelated persons -> passes_qc "
          "False", bool(bad_row["passes_qc"]) is False)
    check("default gate: partial_CDS DOES disqualify by default (updated default, 2026-09-03) "
          "even with 2 unrelated persons -> passes_qc False",
          bool(partial_row["passes_qc"]) is False)

    table3_no_gate = novel.build_table3(matched_rows, ancestry_by_person,
                                          disqualifying_warnings=frozenset())
    bad_row_no_gate = table3_no_gate[
        table3_no_gate["cds_seq_sha1"].str.startswith("BADONE0")].iloc[0]
    partial_row_no_gate = table3_no_gate[
        table3_no_gate["cds_seq_sha1"].str.startswith("PARTIAL0")].iloc[0]
    check("--disqualifying-warnings '' (empty set) makes inframe_stop no longer disqualifying "
          "-> passes_qc True", bool(bad_row_no_gate["passes_qc"]) is True)
    check("--disqualifying-warnings '' (empty set) makes partial_CDS no longer disqualifying "
          "-> passes_qc True", bool(partial_row_no_gate["passes_qc"]) is True)

    table3_strict = novel.build_table3(
        matched_rows, ancestry_by_person,
        disqualifying_warnings=frozenset(["no-start_codon", "partial_CDS", "inframe_stop"]))
    benign_row_strict = table3_strict[
        table3_strict["cds_seq_sha1"].str.startswith("BENIGN0")].iloc[0]
    check("configuring no-start_codon as disqualifying (via the CLI-equivalent argument) actually "
          "disqualifies it -> passes_qc False", bool(benign_row_strict["passes_qc"]) is False)


def test_novelty_depth1_undetermined_classification():
    """Regression for Fix 1 (SCHEMA.md Table 1): novelty_depth==1 (even the first/gene-resolution
    field unresolved, from consensusCall()'s commonprefix truncation colliding on tied candidates)
    must classify as novelty_class == 'undetermined', not crash, not silently drop, and not get
    mis-bucketed into protein_altering/synonymous/beyond_cds."""
    n_fields, is_novel, depth, novelty_class, warn = extract.parse_consensus("HLA-A*new")
    check("depth-1 novel allele is_novel True", is_novel is True)
    check("depth-1 novelty_depth == 1", depth == 1, detail=str(depth))
    check("depth-1 novelty_class == 'undetermined'", novelty_class == "undetermined",
          detail=str(novelty_class))
    check("depth-1 does not raise an 'unexpected novelty depth' warning", warn is None,
          detail=str(warn))

    check("NOVELTY_CLASS_BY_DEPTH has an entry for every depth 1-4 (no KeyError/None risk on a "
          "dict .get(depth) elsewhere)",
          all(d in extract.NOVELTY_CLASS_BY_DEPTH for d in (1, 2, 3, 4)),
          detail=str(extract.NOVELTY_CLASS_BY_DEPTH))

    # End-to-end through build_table1_row: must not crash and must carry the classification through
    # to the Table 1 row dict, with no fabricated warning.
    row = {
        "contig": "ctg1", "gene_id": "IAG100001", "gene_name_raw": "HLA-A",
        "template_allele": "HLA-A*01:01:01:01", "template_distance": "3",
        "gene_start": 1000, "gene_end": 4000, "strand": "+",
        "consensus": "HLA-A*new", "alleles": "HLA-A*01:01:01:01,HLA-A*01:02:01:01",
        "template_warning": None, "cds_distance": None, "cds_mut": None,
    }
    t1_row, row_warn = extract.build_table1_row("P1", "hap1", row)
    check("build_table1_row does not crash on a depth-1 novel consensus", t1_row is not None)
    check("build_table1_row: novelty_depth == 1 survives to the Table 1 row",
          t1_row["novelty_depth"] == 1, detail=str(t1_row["novelty_depth"]))
    check("build_table1_row: novelty_class == 'undetermined' survives to the Table 1 row",
          t1_row["novelty_class"] == "undetermined", detail=str(t1_row["novelty_class"]))
    check("build_table1_row: is_novel True for a depth-1 call", t1_row["is_novel"] is True)
    check("build_table1_row: no spurious warning raised for a legitimate depth-1 call",
          row_warn is None, detail=str(row_warn))


def test_novelty_depth1_excluded_from_headline_but_kept_in_table3():
    """Regression for Fix 1's steer in 03_novel_alleles.py: an 'undetermined' (depth-1) novel
    cluster must still appear in Table 3 (the TSV's grain is every distinct novel sequence cluster)
    but must be excluded from the report's headline novel-allele counts and surfaced as its own
    explicitly-labelled, explicitly-counted category."""
    matched_rows = [
        {"person_id": "P1", "hap": "hap1", "gene": "HLA-C", "gene_class": "classical_I",
         "cds_seq_sha1": "UNDET00" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-C*01:02",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "undetermined",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P2", "hap": "hap1", "gene": "HLA-C", "gene_class": "classical_I",
         "cds_seq_sha1": "UNDET00" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-C*01:02",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "undetermined",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P3", "hap": "hap1", "gene": "HLA-C", "gene_class": "classical_I",
         "cds_seq_sha1": "RESOLVD" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-C*02:02",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
        {"person_id": "P4", "hap": "hap1", "gene": "HLA-C", "gene_class": "classical_I",
         "cds_seq_sha1": "RESOLVD" + "0" * 33, "cds_len": 1000, "nearest_allele": "HLA-C*02:02",
         "cds_distance": 1, "n_aa_changes": 1, "novelty_class": "protein_altering",
         "warning_tokens": frozenset(), "is_homopolymer_indel_only": False},
    ]
    ancestry_by_person = {"P1": "AFR", "P2": "EUR", "P3": "AFR", "P4": "EUR"}
    table3 = novel.build_table3(matched_rows, ancestry_by_person)
    check("Table 3 keeps the undetermined cluster (2 clusters total, not silently dropped)",
          len(table3) == 2, detail=str(len(table3)))
    undet_row = table3[table3["cds_seq_sha1"].str.startswith("UNDET00")].iloc[0]
    check("undetermined cluster still gets passes_qc computed like any other (recurrence-gated)",
          bool(undet_row["passes_qc"]) is True)
    check("undetermined cluster's novelty_class is preserved verbatim in Table 3",
          undet_row["novelty_class"] == "undetermined")

    for c in novel.TABLE3_COLUMNS:
        if c not in table3.columns:
            table3[c] = None
    md_path = "/tmp/hla_test_novel_report.md"
    md_text = novel.write_report(
        md_path, table3,
        {"n_novel_table1_rows": 4, "n_matched": 4, "n_ambiguous_copy": 0, "n_missing_cds": 0},
        __import__("pandas").DataFrame(
            {"novel_rate": [0.1], "n_calls": [10], "n_novel": [1]}, index=["AFR"]),
        {"table3": "novel_alleles.tsv"}, matched_rows)
    check("report explicitly states the undetermined cluster count",
          "undetermined` (gene-level-unresolved) clusters: 1" in md_text, detail=md_text[:2000])
    check("report's headline 'resolved identity' cluster count is 1, not 2 (excludes the "
          "undetermined cluster)",
          "resolved identity (`novel_id`) | 1 |" in md_text)


def test_ambiguous_copy_index_excluded():
    """A gene with copy_index>1 on the same contig must be counted, not silently matched -- per
    reference/IMMUANNOT_GTF_SPEC.md part D's documented ambiguity."""
    import pandas as pd
    table1 = pd.DataFrame([
        {"person_id": "P1", "hap": "hap1", "contig": "ctgA", "gene": "HLA-DRB3", "copy_index": 1,
         "gene_class": "class_II_paralog", "is_novel": True, "novelty_class": "synonymous",
         "cds_distance": 1, "n_aa_changes": 0, "template_warning": None, "has_warning": False,
         "alleles": "HLA-DRB3*01:01", "template_allele": "HLA-DRB3*01:01", "cds_mut": None},
        {"person_id": "P1", "hap": "hap1", "contig": "ctgA", "gene": "HLA-DRB3", "copy_index": 2,
         "gene_class": "class_II_paralog", "is_novel": False, "novelty_class": None,
         "cds_distance": None, "n_aa_changes": None, "template_warning": None, "has_warning": False,
         "alleles": "HLA-DRB3*02:02", "template_allele": "HLA-DRB3*02:02", "cds_mut": None},
    ])
    matched_rows, seqs, stats = novel.match_novel_rows(table1, "/nonexistent/outroot")
    check("copy_index>1 novel row counted as ambiguous, not matched",
          stats["n_ambiguous_copy"] == 1 and stats["n_matched"] == 0, detail=str(stats))


# ---------------------------------------------------------------------------
# 3. Integration: ancestry-skewed novel rate must survive the full pipeline (01 -> 02 -> 03)
# ---------------------------------------------------------------------------
def test_ancestry_gradient_integration(fixtures_dir):
    import subprocess
    outroot = fixtures_dir
    py = sys.executable

    def run(args):
        r = subprocess.run([py] + args, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout, file=sys.stderr)
            print(r.stderr, file=sys.stderr)
        return r

    t1 = os.path.join(outroot, "hla_calls_rich.sample.tsv")
    t4 = os.path.join(outroot, "cohort_membership.sample.tsv")
    reports_dir = os.path.join(outroot, "reports")
    # Person-id directories live under <outroot>/people/ (RUNBOOK.md "Step 1c" / DEFAULT_OUTROOT in
    # 01/03_*.py) -- the aggregate .tsv files (hla_calls_rich.tsv, cohort_membership.tsv, ...) stay
    # directly at <outroot>, so 01 needs an explicit --out-dir pointing back there.
    people_root = os.path.join(outroot, "people")

    r1 = run([os.path.join(HLA_POPGEN_DIR, "01_extract_rich.py"), "--outroot", people_root,
              "--out-dir", outroot, "--sample"])
    check("01_extract_rich.py exits 0 against fixtures", r1.returncode == 0, detail=r1.stderr[-500:])

    r2 = run([os.path.join(HLA_POPGEN_DIR, "02_build_cohorts.py"), "--outroot", outroot,
              "--table1", t1,
              "--cohort-full", os.path.join(outroot, "immuannot_cohort_full.tsv"),
              "--ancestry-preds", os.path.join(outroot, "ancestry_preds.tsv"),
              "--hla-genotypes", os.path.join(outroot, "hla_genotypes.tsv"),
              "--skip-mount-check", "--sample"])
    check("02_build_cohorts.py exits 0 against fixtures", r2.returncode == 0, detail=r2.stderr[-500:])

    seqs_path = os.path.join(reports_dir, "novel_alleles_seqs.sample.fa")
    r3 = run([os.path.join(HLA_POPGEN_DIR, "03_novel_alleles.py"), "--outroot", people_root,
              "--table1", t1, "--cohort-membership", t4, "--out-dir", reports_dir,
              "--seqs-path", seqs_path, "--sample"])
    check("03_novel_alleles.py exits 0 against fixtures", r3.returncode == 0, detail=r3.stderr[-500:])

    import pandas as pd
    t1_df = pd.read_csv(t1, sep="\t", dtype={"person_id": str})
    t4_df = pd.read_csv(t4, sep="\t", dtype={"person_id": str})
    t4_df["ancestry_pred"] = t4_df["ancestry_pred"].astype(str).str.upper()
    ancestry_by_person = dict(zip(t4_df["person_id"], t4_df["ancestry_pred"]))

    rate_df = novel.novel_rate_by_ancestry(t1_df, ancestry_by_person)
    check("AFR novel_rate present", "AFR" in rate_df.index, detail=str(rate_df.index.tolist()))
    check("EUR novel_rate present", "EUR" in rate_df.index, detail=str(rate_df.index.tolist()))
    if "AFR" in rate_df.index and "EUR" in rate_df.index:
        afr_rate = rate_df.loc["AFR", "novel_rate"]
        eur_rate = rate_df.loc["EUR", "novel_rate"]
        check(f"AFR novel_rate ({afr_rate:.3f}) > EUR novel_rate ({eur_rate:.3f}) -- the key "
              f"ancestry-gradient assertion from the fixture design (afr .22 vs eur .05)",
              afr_rate > eur_rate, detail=f"afr={afr_rate}, eur={eur_rate}")
        check("AFR novel_rate roughly in the fixture's designed neighborhood (0.15-0.30)",
              0.15 <= afr_rate <= 0.30, detail=str(afr_rate))
        check("EUR novel_rate roughly in the fixture's designed neighborhood (0.02-0.10)",
              0.02 <= eur_rate <= 0.10, detail=str(eur_rate))

    table3_path = os.path.join(reports_dir, "novel_alleles.sample.tsv")
    check("Table 3 (novel_alleles.tsv) written", os.path.exists(table3_path))
    if os.path.exists(table3_path):
        t3 = pd.read_csv(table3_path, sep="\t")
        for col in novel.TABLE3_COLUMNS:
            check(f"Table 3 has column '{col}'", col in t3.columns)
        check("Table 3 novel_id values are unique (grain: one row per cluster)",
              t3["novel_id"].is_unique, detail=str(t3["novel_id"].duplicated().sum()))

    # novel_alleles_seqs.fa is written wherever --seqs-path points (default: a DEFAULT_DATA_ROOT
    # sibling of the other aggregate .tsv files on a real run) -- pointed at the hermetic
    # reports_dir above so this test run never touches a real ~/pipeline_outputs on whatever
    # machine happens to run the suite.
    check("VM-only sequence FASTA written (real nucleotide sequences, never in a report path "
          "that would otherwise be shared/committed)", os.path.exists(seqs_path), detail=seqs_path)
    if os.path.exists(table3_path):
        with open(table3_path) as f:
            content = f.read()
        check("Table 3 file contains no raw nucleotide-looking sequence column header leak",
              "cds_seq\t" not in content and "sequence\t" not in content)


# ---------------------------------------------------------------------------
# 4. Regression: 04_allele_saturation.py's pct_discovered / CI coherence
#
# Real bug found once make_fixtures.py was rewritten to give novel alleles genuine cross-person
# recurrence (finite per-gene pool instead of an independent-random sequence per call): a raw
# bootstrap percentile CI on Chao2 is not guaranteed to bracket the point estimate computed on the
# full sample -- this bit hardest exactly when Q2_duplicates==0 (Chao2's bias-corrected form blows
# up quadratically on the full sample, but most bootstrap replicates don't reproduce Q2==0 that
# cleanly). Concretely: HLA-B/novel/AFR had chao2=11.85 with a raw CI of [3.00, 10.85] -- the point
# estimate sat OUTSIDE its own interval -- which in turn made the naive pct_discovered inversion
# produce ci_lo (55.3%) > the point estimate (50.6%), and other rows had ci_hi over 300%. Fixed by
# (a) widening the chao2 CI to guarantee it brackets the point estimate whenever it doesn't, and (b)
# clamping pct_discovered and both of its CI bounds into [0, 100] (a monotonic operation, so it
# cannot re-invert an already-correct ordering). This test asserts both invariants hold for every
# row this project's actual `analyze()` produces against real (fixture) data, not just against a
# hand-picked example.
# ---------------------------------------------------------------------------
def test_saturation_ci_bounds_are_coherent(fixtures_dir):
    import pandas as pd
    t1_path = os.path.join(fixtures_dir, "hla_calls_rich.sample.tsv")
    t4_path = os.path.join(fixtures_dir, "cohort_membership.sample.tsv")
    if not (os.path.exists(t1_path) and os.path.exists(t4_path)):
        check("saturation CI coherence test has its prerequisite Table 1 / Table 4 files", False,
              detail=f"missing {t1_path!r} or {t4_path!r} -- run the integration test first")
        return

    table1 = sat.load_table1(t1_path)
    ancestry_by_person = sat.load_ancestry(t4_path)
    matched_rows, _seqs, _stats = novel.match_novel_rows(table1, os.path.join(fixtures_dir, "people"))
    table1_ident = sat.build_allele_identity_table(table1, matched_rows)
    genes = sorted(table1["gene"].unique())

    richness_df, _extrap_df, _curves, _neutral_df = sat.analyze(
        table1_ident, ancestry_by_person, genes, n_bootstrap=150, n_permutations=20)
    check("richness_df is non-empty (test has real data to check)", len(richness_df) > 0)
    if richness_df.empty:
        return

    def _first_few(df, cols):
        return df[cols].head(5).to_dict("records")

    over_100 = richness_df[
        (richness_df["pct_discovered"] > 100 + 1e-9) |
        (richness_df["pct_discovered_ci_lo"].fillna(0) > 100 + 1e-9) |
        (richness_df["pct_discovered_ci_hi"].fillna(0) > 100 + 1e-9)]
    check("no pct_discovered value (point or either CI bound) exceeds 100% in any row",
          len(over_100) == 0,
          detail=_first_few(over_100, ["gene", "ancestry", "category", "pct_discovered",
                                        "pct_discovered_ci_lo", "pct_discovered_ci_hi"]))

    under_0 = richness_df[
        (richness_df["pct_discovered"] < -1e-9) |
        (richness_df["pct_discovered_ci_lo"].fillna(0) < -1e-9) |
        (richness_df["pct_discovered_ci_hi"].fillna(0) < -1e-9)]
    check("no pct_discovered value is negative in any row", len(under_0) == 0,
          detail=_first_few(under_0, ["gene", "ancestry", "category", "pct_discovered",
                                       "pct_discovered_ci_lo", "pct_discovered_ci_hi"]))

    has_pct_ci = richness_df["pct_discovered_ci_lo"].notna() & richness_df["pct_discovered_ci_hi"].notna()
    pct_bracket_bad = richness_df[has_pct_ci & ~(
        (richness_df["pct_discovered_ci_lo"] <= richness_df["pct_discovered"] + 1e-6) &
        (richness_df["pct_discovered"] <= richness_df["pct_discovered_ci_hi"] + 1e-6))]
    check("pct_discovered's own CI brackets its point estimate (ci_lo <= point <= ci_hi) in every "
          "row that has a CI", len(pct_bracket_bad) == 0,
          detail=_first_few(pct_bracket_bad, ["gene", "ancestry", "category", "pct_discovered",
                                               "pct_discovered_ci_lo", "pct_discovered_ci_hi"]))

    has_c2_ci = richness_df["chao2_ci_lo"].notna() & richness_df["chao2_ci_hi"].notna()
    c2_bracket_bad = richness_df[has_c2_ci & ~(
        (richness_df["chao2_ci_lo"] <= richness_df["chao2"] + 1e-6) &
        (richness_df["chao2"] <= richness_df["chao2_ci_hi"] + 1e-6))]
    check("chao2's own CI brackets its point estimate in every row that has a CI (the underlying "
          "fix pct_discovered's coherence depends on)", len(c2_bracket_bad) == 0,
          detail=_first_few(c2_bracket_bad, ["gene", "ancestry", "category", "chao2",
                                              "chao2_ci_lo", "chao2_ci_hi"]))

    check("at least one row actually exercised the Q2=0 low-confidence regime on this fixture "
          "(sanity check that this test isn't vacuously passing)",
          bool(richness_df["low_confidence_chao2"].any()),
          detail="no low_confidence_chao2 rows found -- fixtures may have changed shape again")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", default="/tmp/hla_fixtures_test",
                    help="Path to fixtures built by tests/make_fixtures.py. Built fresh if missing.")
    args = ap.parse_args()

    if not os.path.isdir(args.fixtures):
        print(f"Building fixtures at {args.fixtures!r} ...", file=sys.stderr)
        spec = importlib.util.spec_from_file_location(
            "hla_popgen_make_fixtures", os.path.join(HERE, "make_fixtures.py"))
        make_fixtures = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(make_fixtures)
        make_fixtures.build(args.fixtures, 300)

    print("\n-- Unit tests: richness estimator formulas (hand-computed examples) --")
    test_chao2_hand_example()
    test_chao2_bias_corrected_when_q2_zero()
    test_jackknife1_hand_example()
    test_jackknife2_hand_example()
    test_ace_reduces_sensibly_no_rare_species()
    test_ewens_watterson_expected_k_monotonic_and_solvable()
    test_chao2_extrapolate_no_signal_returns_sobs()
    test_chao2_extrapolate_grows_with_target()

    print("\n-- Unit tests: 03_novel_alleles.py QC / clustering logic --")
    test_homopolymer_indel_detection()
    test_warning_tokens_na_excluded()
    test_build_table3_recurrence_and_passes_qc()
    test_build_table3_token_aware_warning_gate()
    test_novelty_depth1_undetermined_classification()
    test_novelty_depth1_excluded_from_headline_but_kept_in_table3()
    test_ambiguous_copy_index_excluded()

    print("\n-- Integration test: ancestry gradient survives 01 -> 02 -> 03 --")
    test_ancestry_gradient_integration(args.fixtures)

    print("\n-- Regression: 04_allele_saturation.py pct_discovered / CI coherence --")
    test_saturation_ci_bounds_are_coherent(args.fixtures)

    print(f"\n{'ALL PASS' if not FAILURES else f'{len(FAILURES)} FAILURE(S)'}")
    if FAILURES:
        for name in FAILURES:
            print(f"  - {name}")
        sys.exit(1)


if __name__ == "__main__":
    main()
