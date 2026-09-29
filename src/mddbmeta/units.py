"""Unit conversions and comparison tolerances.

Internal canonical units: time in ps, temperature in K, lengths in Angstrom.
MDDB inputs use fs for the integration timestep and ns for the frame step.
"""

from __future__ import annotations

import math

PS_PER_NS = 1000.0
FS_PER_PS = 1000.0

# Continuity tolerance floor (ps): a gap smaller than this is always healthy.
CONTINUITY_FLOOR_PS = 0.1

TRUNCATED_OCTAHEDRON_ANGLE = 109.4712206  # degrees, arccos(-1/3)


def ps_to_fs(ps: float) -> float:
    return clean_float(ps * FS_PER_PS)


def ps_to_ns(ps: float) -> float:
    return clean_float(ps / PS_PER_NS)


def clean_float(x: float, digits: int = 10) -> float:
    """Round away binary noise (0.1 * 3 -> 0.3) without losing real precision."""
    if x == 0 or not math.isfinite(x):
        return x
    return float(f"{x:.{digits}g}")


def rel_close(a: float, b: float, rel: float = 1e-6, abs_tol: float = 1e-12) -> bool:
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_tol)


def continuity_tolerance(frame_interval_ps: float | None) -> float:
    """Tolerance for a hand-off gap: a small floor, or half a frame interval."""
    tol = CONTINUITY_FLOOR_PS
    if frame_interval_ps is not None and frame_interval_ps > 0:
        tol = max(tol, 0.5 * frame_interval_ps)
    return tol
