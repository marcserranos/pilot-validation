#!/usr/bin/env python3
"""Unit tests for the layout linter (`check_layout`/`mark_marginal`/`save_fig` strict mode) added
to `_viz_common.py` in S04 WS-C phase 2 (reference/FIGURE_STYLE.md, CRITIC_WSC.md). Plants each
fault the critic flagged in `fig_dq_g1g2_committed_MAIN_POOLED.png` -- overlapping text, text
clipped off the figure, a marginal axes not aligned with its main axes -- and checks
`check_layout()` catches each one, plus a clean figure that passes with zero errors.

Run: python3 -m pytest scripts/hla_popgen/tests/test_layout_linter.py -q
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

spec = importlib.util.spec_from_file_location(
    "_viz_common", os.path.join(HLA_POPGEN_DIR, "_viz_common.py"))
vc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vc)


def _err_types(violations):
    return {v["type"] for v in violations if v["severity"] == "error"}


def test_clean_figure_has_no_error_violations():
    """A plain, well-spaced single-axes figure must pass with zero error-severity violations."""
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1, 2], [0, 1, 0])
    ax.set_xlabel("x label")
    ax.set_ylabel("y label")
    ax.set_title("a clean title")
    violations = vc.check_layout(fig)
    assert _err_types(violations) == set(), violations
    plt.close(fig)


def test_overlapping_text_detected():
    """Two text artists placed on top of each other must trigger a text_overlap error."""
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.text(0.5, 0.5, "first label", transform=ax.transAxes, fontsize=14)
    ax.text(0.5, 0.5, "second label", transform=ax.transAxes, fontsize=14)
    violations = vc.check_layout(fig)
    assert "text_overlap" in _err_types(violations), violations
    plt.close(fig)


def test_clipped_text_detected():
    """Text placed far outside the figure canvas (in figure-fraction coords, clip_on=False) must
    trigger a text_clipped error."""
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.text(5.0, 5.0, "way off the canvas", fontsize=10)
    violations = vc.check_layout(fig)
    assert "text_clipped" in _err_types(violations), violations
    plt.close(fig)


def test_misaligned_marginal_detected():
    """A marginal axes registered via mark_marginal() whose xlim doesn't match the main axes'
    (e.g. built with an independent gridspec/width rather than sharex) must trigger
    marginal_misaligned -- this is exactly CRITIC_WSC.md's top marginal / heatmap-column fault in
    fig_dq_g1g2_committed_MAIN_POOLED.png."""
    fig, (ax_top, ax_main) = plt.subplots(2, 1, figsize=(4, 4))
    ax_main.imshow([[0, 1, 2], [1, 2, 0]], aspect="auto")
    ax_main.set_xlim(-0.5, 2.5)
    # Marginal built with its OWN unrelated xlim -- not sharex'd, not width-matched -- the exact
    # bug pattern that produced the misaligned bars.
    ax_top.bar([0, 1, 2, 3], [1, 2, 1, 3])
    ax_top.set_xlim(-1, 4)
    vc.mark_marginal(ax_top, ax_main, axis="x")
    violations = vc.check_layout(fig)
    assert "marginal_misaligned" in _err_types(violations), violations
    plt.close(fig)


def test_aligned_marginal_via_sharex_passes():
    """The fix pattern (sharex, matching gridspec) must NOT trigger marginal_misaligned."""
    fig = plt.figure(figsize=(4, 4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 3])
    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_main.imshow([[0, 1, 2], [1, 2, 0]], aspect="auto")
    ax_main.set_xlim(-0.5, 2.5)
    ax_top.bar([0, 1, 2], [1, 2, 1])
    vc.mark_marginal(ax_top, ax_main, axis="x")
    violations = vc.check_layout(fig)
    assert "marginal_misaligned" not in _err_types(violations), violations
    plt.close(fig)


def test_save_fig_strict_raises_on_violation(tmp_path):
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.text(0.5, 0.5, "a", transform=ax.transAxes, fontsize=14)
    ax.text(0.5, 0.5, "b", transform=ax.transAxes, fontsize=14)
    stem = str(tmp_path / "broken")
    with pytest.raises(RuntimeError):
        vc.save_fig(fig, stem)
    assert not os.path.exists(stem + ".png")


def test_save_fig_strict_false_requires_reason(tmp_path):
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.text(0.5, 0.5, "a", transform=ax.transAxes, fontsize=14)
    ax.text(0.5, 0.5, "b", transform=ax.transAxes, fontsize=14)
    stem = str(tmp_path / "broken2")
    with pytest.raises(ValueError):
        vc.save_fig(fig, stem, strict=False)


def test_save_fig_strict_false_with_reason_saves():
    import tempfile
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.text(0.5, 0.5, "a", transform=ax.transAxes, fontsize=14)
    ax.text(0.5, 0.5, "b", transform=ax.transAxes, fontsize=14)
    with tempfile.TemporaryDirectory() as d:
        stem = os.path.join(d, "override")
        pdf_path, png_path = vc.save_fig(fig, stem, strict=False, reason="deliberate test override")
        assert os.path.exists(pdf_path)
        assert os.path.exists(png_path)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
