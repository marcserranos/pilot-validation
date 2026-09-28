#!/usr/bin/env python3
"""Single disclosure-policy module for the S04 KIR/HLA pipeline (scripts 43/43b/44/45/46).

Source of truth: `reference/AOU_SMALL_CELL_POLICY.md` (the research) and the newest
`context/DECISIONS.md` entry, "Resolved (2026-09-28, Marc): two disclosure versions, split by
location" (authoritative). Summary of that decision:

  - INTERNAL (`--disclosure internal`): exact values everywhere, n<20 cells visibly marked
    (never hidden), written ONLY under `~/s04/internal/<NN>/` on the VM -- never pulled,
    never committed, never written into this repo's working tree.
  - PUBLIC (`--disclosure public`, the default): the existing hard `<20` masking rule for
    participant/carrier counts and small-back-calculable rates, but LOOSENED for counts that
    are properties of alleles/QC outcomes rather than of people (distinct-allele counts,
    richness, QC tallies, recurrence-class counts) -- these are reported exact, with a
    `review_flag` on recurrence-class cells whose underlying stratum has few total carriers,
    so a human (not this module) decides publishability case by case.

This module is the ONLY place that decides "mask or not" for a given (value, count_type, mode)
triple -- callers (43/43b/44/45/46) must route every disclosure decision through `mask()` /
`mask_rate()` rather than reimplementing the `<20` string literal locally, so the two disclosure
versions can never silently drift apart from this policy.
"""
import dataclasses
import math
import os

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))

# ---------------------------------------------------------------------------
# Count-type / mode vocabulary
# ---------------------------------------------------------------------------
PARTICIPANT = "participant"
CARRIER_NAMED_ALLELE = "carrier_named_allele"
RATE_NUM_DENOM = "rate_numerator/denominator"
ALLELE_DISTINCT = "allele_distinct"
RICHNESS = "richness"
QC_TALLY = "qc_tally"
RECURRENCE_CLASS = "recurrence_class"

COUNT_TYPES = {PARTICIPANT, CARRIER_NAMED_ALLELE, RATE_NUM_DENOM, ALLELE_DISTINCT, RICHNESS,
               QC_TALLY, RECURRENCE_CLASS}

# Count types that stay masked (1-19 -> "<20") in PUBLIC mode -- the non-negotiable "any count of
# PARTICIPANTS" clause of the AoU policy (AOU_SMALL_CELL_POLICY.md "Proposed rule for this
# project"). Everything else in COUNT_TYPES is reported exact in PUBLIC mode.
_PUBLIC_MASKED_TYPES = {PARTICIPANT, CARRIER_NAMED_ALLELE, RATE_NUM_DENOM}
_PUBLIC_EXACT_TYPES = {ALLELE_DISTINCT, RICHNESS, QC_TALLY, RECURRENCE_CLASS}

PUBLIC = "public"
INTERNAL = "internal"
MODES = {PUBLIC, INTERNAL}

SUPPRESS_BELOW = 20
# Default review-flag factor for recurrence_class cells (task instruction: "some documented
# factor; choose a defensible default, e.g. total unrelated carriers of the gene in that
# ancestry < 100"). 100 = 5x SUPPRESS_BELOW -- chosen so a recurrence-class cell is flagged
# whenever the WHOLE stratum it comes from is itself thin enough that even a handful of
# eq1/eq2 alleles could plausibly account for most of its carriers (AOU_SMALL_CELL_POLICY.md
# answer (2)'s "apply clause 2 caution per-locus" caveat), not just when the individual cell
# value is small.
REVIEW_TOTAL_CARRIERS_FLOOR = 100

INTERNAL_HEADER_COMMENT = "# INTERNAL — contains n<20 cells, do not export"


@dataclasses.dataclass(frozen=True)
class MaskResult:
    """`display` is what a caller should write into the output cell (a string, possibly empty).
    `lt20` is True iff the true underlying value is in the disclosive 1-19 band (populated in
    BOTH modes -- INTERNAL uses it as the visible `lt20` flag column; PUBLIC callers can use it
    to decide hatching in a figure). `review_flag` is only ever True for `recurrence_class` cells
    in PUBLIC mode whose stratum total fell below the review floor; it is always False otherwise."""
    display: str
    lt20: bool
    review_flag: bool = False


def _is_missing(value):
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def mask(value, count_type, mode, *, total_carriers=None,
         review_floor=REVIEW_TOTAL_CARRIERS_FLOOR):
    """Core policy function. `value`: a single non-negative integer count (or None/NaN for "not
    applicable", which always renders as an empty string in both modes -- never guessed at).
    `count_type`: one of the COUNT_TYPES constants above. `mode`: 'public' or 'internal'.

    `total_carriers`: for `count_type='recurrence_class'` only, the total number of unrelated
    carriers of the gene in that locus/ancestry stratum (i.e. the denominator the review-flag
    threshold is judged against) -- required to compute `review_flag` in PUBLIC mode; if omitted,
    `review_flag` conservatively defaults to True (flag for human review rather than silently
    assume the stratum is well-powered) for a `recurrence_class` cell.

    Returns a MaskResult. Never raises for a normal (missing or in-range) value; raises
    ValueError only for a genuinely malformed count_type/mode/negative value, since those are
    caller bugs, not disclosure edge cases.
    """
    if count_type not in COUNT_TYPES:
        raise ValueError(f"_disclosure.mask: unknown count_type {count_type!r}; expected one of "
                          f"{sorted(COUNT_TYPES)}")
    if mode not in MODES:
        raise ValueError(f"_disclosure.mask: unknown mode {mode!r}; expected one of {sorted(MODES)}")

    if _is_missing(value):
        return MaskResult(display="", lt20=False, review_flag=False)

    n = int(value)
    if n < 0:
        raise ValueError(f"_disclosure.mask: count must be >= 0, got {n}")
    lt20 = 0 < n < SUPPRESS_BELOW

    if mode == INTERNAL:
        # Exact everywhere in INTERNAL mode; lt20 flag column tells the reader which cells
        # PUBLIC would have masked. review_flag is a PUBLIC-only concept (it steers what a
        # supervisor should double check before something becomes public) -- always False here.
        return MaskResult(display=str(n), lt20=lt20, review_flag=False)

    # mode == PUBLIC
    if count_type in _PUBLIC_EXACT_TYPES:
        review_flag = False
        if count_type == RECURRENCE_CLASS:
            if total_carriers is None:
                review_flag = True  # conservative default -- see docstring
            else:
                review_flag = int(total_carriers) < review_floor
        return MaskResult(display=str(n), lt20=lt20, review_flag=review_flag)

    # _PUBLIC_MASKED_TYPES: participant / carrier_named_allele / rate_numerator/denominator
    if n == 0:
        # A true zero is never disclosive (feedback_suppressed_counts_are_not_zero.md) and must
        # never collapse into "<20" alongside a real 1-19 count.
        return MaskResult(display="0", lt20=False, review_flag=False)
    if lt20:
        return MaskResult(display=f"<{SUPPRESS_BELOW}", lt20=True, review_flag=False)
    return MaskResult(display=str(n), lt20=False, review_flag=False)


def suppressed(n, mode=PUBLIC):
    """Backward-compatible thin wrapper matching 43_kir_full_aggregate.suppressed()'s old
    signature/semantics for a bare participant/carrier count (count_type=PARTICIPANT). Prefer
    calling mask() directly with the correct count_type in new code."""
    return mask(n, PARTICIPANT, mode).display


def wilson_ci(count, nobs, z=1.96):
    """Wilson score interval, scalar-friendly -- identical formula to _viz_common.wilson_ci /
    43_kir_full_aggregate.wilson_ci, reproduced here (not imported) to keep this module
    dependency-light (numpy only)."""
    if nobs <= 0:
        return float("nan"), float("nan"), float("nan")
    p = count / nobs
    denom = 1 + z ** 2 / nobs
    center = (p + z ** 2 / (2 * nobs)) / denom
    half = (z * np.sqrt(p * (1 - p) / nobs + z ** 2 / (4 * nobs ** 2))) / denom
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return p, lo, hi


def mask_rate(n, d, mode, pct_decimals=1):
    """Rate/frequency disclosure, count_type=RATE_NUM_DENOM implicitly (numerator AND
    denominator are both back-calculable small-cell risks per AOU_SMALL_CELL_POLICY.md answer
    (4)). Returns a dict: n, d (display strings), pct, lo, hi (formatted percentage strings, or
    "" when blanked/undefined), lt20_n, lt20_d (bool flags, populated in both modes).

    PUBLIC: point estimate + CI are blanked (not just the raw n/d) whenever d==0 or either n or d
    is in the disclosive 1-19 band -- matches 43_kir_full_aggregate.rate_row()'s existing
    behaviour exactly, so this is a drop-in replacement.
    INTERNAL: n/d are always exact; pct/lo/hi are computed whenever d>0 (never blanked -- that's
    the point of INTERNAL), and lt20_n/lt20_d mark which side(s) PUBLIC would have masked.
    """
    n, d = int(n), int(d)
    n_res = mask(n, RATE_NUM_DENOM, mode)
    d_res = mask(d, RATE_NUM_DENOM, mode)
    fmt = f"%.{pct_decimals}f"

    if mode == PUBLIC and (d == 0 or (0 < n < SUPPRESS_BELOW) or (0 < d < SUPPRESS_BELOW)):
        return {"n": n_res.display, "d": d_res.display, "pct": "", "lo": "", "hi": "",
                "lt20_n": n_res.lt20, "lt20_d": d_res.lt20}
    if d == 0:
        return {"n": n_res.display, "d": d_res.display, "pct": "", "lo": "", "hi": "",
                "lt20_n": n_res.lt20, "lt20_d": d_res.lt20}
    p, lo, hi = wilson_ci(n, d)
    return {"n": n_res.display, "d": d_res.display,
            "pct": fmt % (100 * p), "lo": fmt % (100 * lo), "hi": fmt % (100 * hi),
            "lt20_n": n_res.lt20, "lt20_d": d_res.lt20}


# ---------------------------------------------------------------------------
# Output-path / TSV-writing routing (task item 2).
# ---------------------------------------------------------------------------
def default_out_dir(mode, script_num, public_default):
    """`script_num` like '44' or '43b'. PUBLIC keeps whatever the caller already defaults to
    (`public_default`, unchanged); INTERNAL defaults to ~/s04/internal/<script_num>/ regardless
    of what public_default was, per the two-disclosure-versions decision."""
    if mode not in MODES:
        raise ValueError(f"_disclosure.default_out_dir: unknown mode {mode!r}")
    if mode == PUBLIC:
        return public_default
    return os.path.expanduser(os.path.join("~/s04/internal", str(script_num)))


def assert_internal_path_allowed(path):
    """Raises ValueError if `path` is inside this repo's working tree, or looks like a
    `reports/` path anywhere (even outside the repo, since a report is by definition meant to be
    read/shared) -- INTERNAL output must never land in either place (task item 2 / DECISIONS.md:
    INTERNAL "is generated and kept ONLY on the VM under ~/s04/internal/ ... never pulled or
    committed"). No-op for a path that passes.
    """
    abspath = os.path.abspath(os.path.expanduser(str(path)))
    repo_root_abs = os.path.abspath(_REPO_ROOT)
    if abspath == repo_root_abs or abspath.startswith(repo_root_abs + os.sep):
        raise ValueError(
            f"_disclosure: refusing to write INTERNAL disclosure output to {path!r} -- it is "
            f"inside the repo working tree ({repo_root_abs}). INTERNAL (uncensored, n<20-visible) "
            f"output must live only under ~/s04/internal/ on the VM (context/DECISIONS.md, "
            f"'two disclosure versions'), never in the repo.")
    parts = abspath.replace("\\", "/").split("/")
    if "reports" in parts:
        raise ValueError(
            f"_disclosure: refusing to write INTERNAL disclosure output to {path!r} -- it looks "
            f"like a reports/ path. INTERNAL output must never be written into any reports/ "
            f"directory, in or out of the repo.")


def write_tsv(df, path, mode):
    """Writes `df` to `path` as a TSV. PUBLIC: identical to plain `df.to_csv(path, sep='\\t',
    index=False)` (no behaviour change for existing PUBLIC callers). INTERNAL: first calls
    `assert_internal_path_allowed(path)` (raises rather than writes on a disallowed path), then
    prepends `INTERNAL_HEADER_COMMENT` as the file's first line before the TSV body."""
    if mode not in MODES:
        raise ValueError(f"_disclosure.write_tsv: unknown mode {mode!r}")
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    if mode == INTERNAL:
        assert_internal_path_allowed(path)
        with open(path, "w") as f:
            f.write(INTERNAL_HEADER_COMMENT + "\n")
            df.to_csv(f, sep="\t", index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def add_disclosure_arg(ap, default=PUBLIC):
    """Shared argparse wiring: `--disclosure {public,internal}` (default 'public'), for 43/43b/
    44/45/46's main() to call on their existing ArgumentParser."""
    ap.add_argument("--disclosure", choices=sorted(MODES), default=default,
                     help="'public' (default): existing <20 masking for participant/carrier "
                          "counts and small-back-calculable rates, exact allele-distinct/"
                          "richness/QC-tally/recurrence-class counts (with a review_flag on thin "
                          "recurrence-class strata). 'internal': everything exact, with an lt20 "
                          "flag column on every cell PUBLIC would mask; writes ONLY under "
                          "~/s04/internal/<script>/, never into the repo or reports/.")
    return ap
