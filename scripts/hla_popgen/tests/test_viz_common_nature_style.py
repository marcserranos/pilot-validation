#!/usr/bin/env python3
"""Unit tests for the Nature-style publication additions to _viz_common.py (S03 Task A):
`nature_style()`, `mm()`, `save_fig()`, `panel_letter()`, the diverging colormap/norm helpers,
and `hatch_suppressed()`.

The one thing worth actually checking here is that `save_fig()` writes BOTH a .pdf and a .png
for a given path stem, since that's the contract every future figure script will depend on.

Run: python3 scripts/hla_popgen/tests/test_viz_common_nature_style.py
"""
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


vc = _load_module("_viz_common.py", "hla_popgen_viz_common_nature_style")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def test_module_imports_and_exposes_the_new_names():
    for name in ("nature_style", "mm", "save_fig", "panel_letter", "SUPPRESSED_COLOR",
                 "diverging_cmap", "diverging_norm", "hatch_suppressed", "NATURE_RC",
                 "NATURE_SINGLE_COL_MM", "NATURE_DOUBLE_COL_MM", "PANEL_LETTER_FONTSIZE"):
        check(f"_viz_common exposes {name}", hasattr(vc, name))


def test_mm_conversion():
    check("mm(25.4) == 1.0 inch", abs(vc.mm(25.4) - 1.0) < 1e-9)
    check("mm(NATURE_SINGLE_COL_MM) is a sane figure width in inches",
          3.0 < vc.mm(vc.NATURE_SINGLE_COL_MM) < 4.0)


def test_nature_style_context_manager_restores_rcparams():
    import matplotlib
    before = matplotlib.rcParams["font.size"]
    with vc.nature_style():
        inside = matplotlib.rcParams["font.size"]
        check("nature_style() applies NATURE_RC font.size inside the block",
              inside == vc.NATURE_RC["font.size"])
    after = matplotlib.rcParams["font.size"]
    check("nature_style() restores rcParams on exit", after == before, f"{after} != {before}")


def test_save_fig_writes_both_pdf_and_png():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with tempfile.TemporaryDirectory() as tmp:
        stem = os.path.join(tmp, "nested", "dir", "test_panel")
        with vc.nature_style():
            fig, ax = plt.subplots(figsize=(vc.mm(vc.NATURE_SINGLE_COL_MM), vc.mm(50)))
            ax.plot([0, 1, 2], [0, 1, 0])
            vc.panel_letter(ax, "A")  # should lowercase to "a"
        pdf_path, png_path = vc.save_fig(fig, stem)

        check("save_fig returns the expected .pdf path", pdf_path == stem + ".pdf")
        check("save_fig returns the expected .png path", png_path == stem + ".png")
        check("save_fig wrote a non-empty .pdf",
              os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 0)
        check("save_fig wrote a non-empty .png",
              os.path.exists(png_path) and os.path.getsize(png_path) > 0)


def test_diverging_cmap_and_norm_centered_at_zero():
    cmap = vc.diverging_cmap()
    norm = vc.diverging_norm(vmin=-2.0, vmax=3.0)
    check("diverging_norm centers at 0", norm.vcenter == 0.0)
    lo = cmap(norm(-2.0))
    mid = cmap(norm(0.0))
    hi = cmap(norm(3.0))
    check("diverging cmap: low end is not equal to the midpoint color", lo != mid)
    check("diverging cmap: high end is not equal to the midpoint color", hi != mid)
    check("diverging cmap: low and high ends differ (blue vs red)", lo != hi)


def test_suppressed_color_is_reserved_and_distinct():
    check("SUPPRESSED_COLOR is a hex string", isinstance(vc.SUPPRESSED_COLOR, str)
          and vc.SUPPRESSED_COLOR.startswith("#"))
    check("SUPPRESSED_COLOR is not any ANCESTRY_COLORS value",
          vc.SUPPRESSED_COLOR not in vc.ANCESTRY_COLORS.values())


def test_hatch_suppressed_adds_a_patch():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    n_before = len(ax.patches)
    patch = vc.hatch_suppressed(ax, 0, 0, 1, 1)
    check("hatch_suppressed adds exactly one patch", len(ax.patches) == n_before + 1)
    check("hatch_suppressed patch uses SUPPRESSED_COLOR",
          patch.get_facecolor()[:3] == matplotlib.colors.to_rgb(vc.SUPPRESSED_COLOR))
    plt.close(fig)


def main():
    print("test_viz_common_nature_style.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
