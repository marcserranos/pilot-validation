#!/usr/bin/env python3
"""Unit tests for the cnsplots-derived style port in _viz_common.py (S04 WS-C, reference/
FIGURE_STYLE.md). Covers: nature_style() actually sets the rcParams the style guideline requires
(no grid, top/right spines off, font sizes, legend frameon=False) -- so a future edit that
accidentally drops one of these regresses silently instead of loudly, and JOURNAL_PALETTES /
ACCENT_COLOR exist with the expected shape.

Run: python3 scripts/hla_popgen/tests/test_viz_common_style.py
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import matplotlib

spec = importlib.util.spec_from_file_location(
    "_viz_common", os.path.join(HLA_POPGEN_DIR, "_viz_common.py"))
vc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vc)


def test_nature_style_context_sets_rcparams():
    """Inside the `with nature_style():` block, every cnsplots-derived rcParam must be active;
    outside it, rcParams must be restored (rc_context's whole point) -- prevents a future figure
    accidentally leaking style into the next one."""
    baseline_grid = matplotlib.rcParams["axes.grid"]
    with vc.nature_style():
        assert matplotlib.rcParams["axes.grid"] is False, "cnsplots: axes_grid=False"
        assert matplotlib.rcParams["axes.spines.top"] is False
        assert matplotlib.rcParams["axes.spines.right"] is False
        assert matplotlib.rcParams["legend.frameon"] is False, "cnsplots: legend_frameon=False"
        assert matplotlib.rcParams["font.size"] == 7
        assert matplotlib.rcParams["axes.titlesize"] == 7
        assert matplotlib.rcParams["axes.labelsize"] == 7
        assert matplotlib.rcParams["xtick.labelsize"] == 5
        assert matplotlib.rcParams["ytick.labelsize"] == 5
        assert matplotlib.rcParams["legend.fontsize"] == 6
        assert matplotlib.rcParams["xtick.major.size"] == 2, "cnsplots: xtick_major_size=2"
        assert matplotlib.rcParams["ytick.major.size"] == 2
        assert matplotlib.rcParams["mathtext.fontset"] == "custom"
        assert matplotlib.rcParams["pdf.fonttype"] == 42
        assert matplotlib.rcParams["ps.fonttype"] == 42
        assert matplotlib.rcParams["savefig.dpi"] == 600
        sans = matplotlib.rcParams["font.sans-serif"]
        for fallback in ("Helvetica Neue", "Nimbus Sans", "Liberation Sans"):
            assert fallback in sans, f"expected {fallback!r} in font.sans-serif fallback list"
    # restored on exit
    assert matplotlib.rcParams["axes.grid"] == baseline_grid


def test_nature_style_apply_true_mutates_global_rcparams():
    """apply=True pushes onto global rcParams with no auto-restore (documented escape hatch)."""
    try:
        vc.nature_style(apply=True)
        assert matplotlib.rcParams["legend.frameon"] is False
        assert matplotlib.rcParams["axes.grid"] is False
    finally:
        matplotlib.rcdefaults()


def test_legend_kw_matches_cnsplots_small_glyph_spec():
    assert vc.LEGEND_KW["frameon"] is False
    assert vc.LEGEND_KW["markerscale"] == 0.5
    assert vc.LEGEND_KW["handlelength"] == 0.7
    assert vc.LEGEND_KW["handletextpad"] == 0.3


def test_journal_palettes_present_and_well_formed():
    assert set(vc.JOURNAL_PALETTES) >= {"nature", "science", "lancet", "nejm", "cell"}
    for name, hexes in vc.JOURNAL_PALETTES.items():
        assert len(hexes) >= 5, f"{name} palette too short"
        for h in hexes:
            assert h.startswith("#") and len(h) == 7, f"{name} has a malformed hex {h!r}"
    # Accent picked FROM a journal palette, not from the ancestry palette (must stay distinct).
    assert vc.ACCENT_COLOR in vc.JOURNAL_PALETTES["nature"]
    assert vc.ACCENT_COLOR not in vc.ANCESTRY_COLORS.values()


def test_backward_compatible_names_still_exist():
    """Old scripts import these by name -- the port must not rename or remove any of them."""
    for name in ("nature_style", "mm", "save_fig", "panel_letter", "SUPPRESSED_COLOR",
                 "hatch_suppressed", "diverging_cmap", "diverging_norm", "NATURE_RC",
                 "NATURE_SINGLE_COL_MM", "NATURE_DOUBLE_COL_MM", "ANCESTRY_COLORS"):
        assert hasattr(vc, name), f"backward-compat break: _viz_common.{name} missing"


def _run_all():
    tests = [v for k, v in globals().items() if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {t.__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    _run_all()
