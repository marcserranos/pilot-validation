#!/usr/bin/env python3
"""Unit tests for the figure-drawing side of the S04 disclosure layer (_viz_common.py's
add_internal_watermark / mark_internal_lt20 / mark_review_flag) -- INTERNAL figures must carry the
"do not export" watermark, PUBLIC figures must not, and check_layout(strict=True) must still pass
for a figure that uses the review-flag/lt20 markers (task item 3's "Figures must still pass
check_layout(strict=True)").

Run: pytest scripts/hla_popgen/tests/test_disclosure_figures.py -q
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
import matplotlib.text as mtext
import pytest

import _viz_common as vc


def _make_simple_fig():
    fig, ax = plt.subplots(figsize=(vc.mm(89), vc.mm(60)))
    ax.plot([0, 1, 2], [0, 1, 0.5])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return fig, ax


def _has_watermark_text(fig, text=vc.INTERNAL_WATERMARK_TEXT):
    return any(t.get_text() == text for t in fig.findobj(mtext.Text) if t.get_visible())


def test_watermark_present_after_add_internal_watermark():
    fig, _ax = _make_simple_fig()
    vc.add_internal_watermark(fig)
    assert _has_watermark_text(fig)
    plt.close(fig)


def test_watermark_absent_without_add_internal_watermark():
    fig, _ax = _make_simple_fig()
    assert not _has_watermark_text(fig)
    plt.close(fig)


def test_watermark_does_not_itself_trigger_a_layout_error():
    """The watermark is expected to visually cross real data/text by design -- it must never be
    the thing that makes check_layout(strict=True) refuse to save a figure."""
    fig, _ax = _make_simple_fig()
    vc.add_internal_watermark(fig)
    violations = vc.check_layout(fig)
    errors = [v for v in violations if v["severity"] == "error"]
    assert errors == [], f"watermark introduced layout error(s): {errors}"
    plt.close(fig)


def test_mark_internal_lt20_adds_patch_and_tag_without_layout_error():
    fig, ax = _make_simple_fig()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    rect, txt = vc.mark_internal_lt20(ax, 1, 1, 2, 2)
    assert rect in ax.patches
    assert txt.get_text() == "n<20"
    violations = vc.check_layout(fig)
    errors = [v for v in violations if v["severity"] == "error"]
    assert errors == [], f"mark_internal_lt20 introduced layout error(s): {errors}"
    plt.close(fig)


def test_mark_review_flag_adds_unfilled_patch_without_layout_error():
    fig, ax = _make_simple_fig()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    rect = vc.mark_review_flag(ax, 4, 4, 2, 2)
    assert rect in ax.patches
    assert rect.get_facecolor()[3] == 0.0 or rect.get_fill() is False or \
        rect.get_facecolor() == (0.0, 0.0, 0.0, 0.0)
    violations = vc.check_layout(fig)
    errors = [v for v in violations if v["severity"] == "error"]
    assert errors == [], f"mark_review_flag introduced layout error(s): {errors}"
    plt.close(fig)


def test_mark_internal_lt20_and_review_flag_are_visually_distinct_from_hatch_suppressed():
    """hatch_suppressed() (PUBLIC 'value hidden') must remain visually distinct from
    mark_review_flag() (PUBLIC 'value shown, please review') and mark_internal_lt20() (INTERNAL
    'value shown, was <20') -- different fill/edge encodings, not just different call sites."""
    fig, ax = _make_simple_fig()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    hatched = vc.hatch_suppressed(ax, 0, 0, 1, 1)
    review = vc.mark_review_flag(ax, 2, 2, 1, 1)
    internal_rect, _tag = vc.mark_internal_lt20(ax, 4, 4, 1, 1)

    assert hatched.get_hatch() == "////"
    assert review.get_hatch() is None
    assert internal_rect.get_hatch() is None
    assert hatched.get_facecolor()[:3] != review.get_facecolor()[:3] or \
        hatched.get_facecolor()[3] != review.get_facecolor()[3]
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
