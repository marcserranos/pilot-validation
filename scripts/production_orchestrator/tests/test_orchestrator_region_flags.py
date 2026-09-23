#!/usr/bin/env python3
"""Unit tests for run_production_orchestrator.py's --region/--pad/--out-suffix flags (added for
S03 WS5's KIR full run, 2026-09-23). The one hard requirement: **backward compatibility** -- an
invocation with none of these flags set must behave byte-for-byte like before this change
(DEFAULT_REGION/DEFAULT_PAD unchanged, canonical output filenames unchanged). Also covers the new
safety guard (a non-default --region without --out-suffix must refuse to start, not silently
merge into the real HLA canonical files -- see the inline comment at that die() call for the
ENVIRONMENT.md quirk #29 precedent this guards against).

Run: python3 scripts/production_orchestrator/tests/test_orchestrator_region_flags.py
"""
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ORCH_DIR = os.path.dirname(HERE)
ORCH_PATH = os.path.join(ORCH_DIR, "run_production_orchestrator.py")


def _load_module():
    spec = importlib.util.spec_from_file_location("orch_test_target", ORCH_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


orch = _load_module()


class TestDefaultsUnchanged(unittest.TestCase):
    """The literal HLA defaults this whole project has depended on since 2026-08-04 -- a typo
    here would silently change the standing HLA window for every future default-args run."""

    def test_default_region_is_hla_window(self):
        self.assertEqual(orch.DEFAULT_REGION, "chr6:29500000-33500000")

    def test_default_pad_unchanged(self):
        self.assertEqual(orch.DEFAULT_PAD, 100_000)


class TestMergeFragmentsBackwardCompatible(unittest.TestCase):
    """out_suffix="" (the default) must reproduce the exact prior filenames/glob pattern."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.outroot = self.tmpdir.name
        pd.DataFrame([{"person_id": "1", "gene": "HLA-A", "immuannot_1": "A*01:01",
                       "immuannot_2": "A*02:01"}]).to_csv(
            os.path.join(self.outroot, "immuannot_calls.1.tsv"), sep="\t", index=False)
        pd.DataFrame([{"person_id": "1", "hap": "hap1", "immuannot_seconds": "10"}]).to_csv(
            os.path.join(self.outroot, "immuannot_timing.1.tsv"), sep="\t", index=False)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_default_suffix_writes_unsuffixed_canonical_files(self):
        orch.merge_fragments(self.outroot)  # out_suffix defaults to ""
        self.assertTrue(os.path.exists(os.path.join(self.outroot, "immuannot_calls.tsv")))
        self.assertTrue(os.path.exists(os.path.join(self.outroot, "immuannot_timing.tsv")))
        # Fragments are consumed (deleted) once merged -- same as before this change.
        self.assertFalse(os.path.exists(os.path.join(self.outroot, "immuannot_calls.1.tsv")))

    def test_out_suffix_keeps_a_distinct_canonical_file(self):
        # A KIR-style fragment, isolated by a DIFFERENT out_suffix, must not appear in or be
        # consumed by a default (HLA) merge, and vice versa.
        pd.DataFrame([{"person_id": "2", "gene": "KIR3DL2", "immuannot_1": "KIR3DL2*001",
                       "immuannot_2": "KIR3DL2*002"}]).to_csv(
            os.path.join(self.outroot, "immuannot_calls.kir.2.tsv"), sep="\t", index=False)
        orch.merge_fragments(self.outroot, out_suffix="")
        orch.merge_fragments(self.outroot, out_suffix=".kir")
        self.assertTrue(os.path.exists(os.path.join(self.outroot, "immuannot_calls.tsv")))
        self.assertTrue(os.path.exists(os.path.join(self.outroot, "immuannot_calls.kir.tsv")))
        hla_df = pd.read_csv(os.path.join(self.outroot, "immuannot_calls.tsv"), sep="\t", dtype=str)
        kir_df = pd.read_csv(os.path.join(self.outroot, "immuannot_calls.kir.tsv"), sep="\t",
                              dtype=str)
        self.assertEqual(set(hla_df["gene"]), {"HLA-A"})
        self.assertEqual(set(kir_df["gene"]), {"KIR3DL2"})


class TestRegionOutSuffixGuard(unittest.TestCase):
    """Invoking the real CLI: a non-default --region without --out-suffix must die loudly before
    touching the mount/cohort (fast, no fixture files needed) -- not silently proceed toward a
    canonical-file collision with any concurrent/prior HLA run sharing the same --outroot."""

    def test_non_default_region_without_out_suffix_dies(self):
        proc = subprocess.run(
            [sys.executable, ORCH_PATH, "--region", "chr19:54600000-54920000"],
            capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("out-suffix", proc.stderr)
        self.assertIn("FATAL", proc.stderr)

    def test_non_default_region_with_out_suffix_passes_the_guard(self):
        # Should get past the region/out-suffix guard and fail later for an unrelated reason
        # (no real mount/cohort in this test environment) -- confirms the guard itself doesn't
        # block a properly-flagged KIR-style invocation.
        proc = subprocess.run(
            [sys.executable, ORCH_PATH, "--region", "chr19:54600000-54920000",
             "--out-suffix", ".kir", "--concurrency", "1", "--threads-per-person", "1"],
            capture_output=True, text=True)
        self.assertNotIn("out-suffix", proc.stderr.split("mount")[0] if "mount" in proc.stderr
                          else proc.stderr)


if __name__ == "__main__":
    unittest.main()
