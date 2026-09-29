"""Periodic box type from the topology's box flag and the box angles/lengths."""

from __future__ import annotations

from mddbmeta.model import SystemInfo
from mddbmeta.provenance import Confidence, Mined, Source
from mddbmeta.units import TRUNCATED_OCTAHEDRON_ANGLE

# Labels written to MDDB. Keep them in one place so they can be aligned with the
# vocabulary already used in the database ("check existing values" in the template).
CUBIC = "Cubic"
RECTANGULAR = "Rectangular"
TRUNC_OCT = "Truncated Octahedron"
TRICLINIC = "Triclinic"
DODECAHEDRON = "Dodecahedron"

ANGLE_TOL = 0.05
LENGTH_REL_TOL = 1e-3


def _close(a, b) -> bool:
    return all(abs(x - y) <= ANGLE_TOL * 20 for x, y in zip(a, b, strict=True))


def classify_box(
    lengths: tuple[float, float, float] | None,
    angles: tuple[float, float, float] | None,
    flag: int | None,
) -> str | None:
    if flag == 0:
        return None
    if angles is not None and all(abs(a - TRUNCATED_OCTAHEDRON_ANGLE) <= ANGLE_TOL for a in angles):
        return TRUNC_OCT
    if angles is not None:
        srt = sorted(angles)
        # GROMACS triclinic forms of the same cells
        if _close(srt, (60.0, 60.0, 90.0)):  # rhombic dodecahedron (xy-square)
            return DODECAHEDRON
        if _close(srt, (60.0, 60.0, 60.0)):  # rhombic dodecahedron (xy-hexagon)
            return DODECAHEDRON
        if _close(
            srt,
            (
                180.0 - TRUNCATED_OCTAHEDRON_ANGLE,
                180.0 - TRUNCATED_OCTAHEDRON_ANGLE,
                TRUNCATED_OCTAHEDRON_ANGLE,
            ),
        ):
            return TRUNC_OCT
    if flag == 2:
        return TRUNC_OCT
    if angles is not None and all(abs(a - 90.0) <= ANGLE_TOL for a in angles):
        if lengths and max(lengths) - min(lengths) <= LENGTH_REL_TOL * max(lengths):
            return CUBIC
        return RECTANGULAR
    if angles is not None:
        return TRICLINIC
    if flag == 1:
        return RECTANGULAR
    return None


def mine_boxtype(
    system: SystemInfo | None,
    coords_box: tuple[tuple[float, float, float], tuple[float, float, float], str] | None = None,
) -> Mined[str]:
    """``coords_box`` is (lengths, angles, source path) from a restart or trajectory frame."""
    if system is None:
        return Mined.missing("no topology to inspect")
    flag = system.periodic_box_type
    if flag == 0:
        return Mined(
            value=None,
            method="non-periodic",
            confidence=Confidence.EXACT,
            sources=(Source(system.source, "POINTERS IFBOX"),),
            note="non-periodic system: no box",
        )
    lengths, angles = system.box_lengths, system.box_angles
    sources = [Source(system.source, "box flag and dimensions")]
    if coords_box is not None:
        lengths, angles, path = coords_box
        sources.append(Source(path, "box line"))
    label = classify_box(lengths, angles, flag)
    if label is None:
        return Mined.missing("box shape could not be determined")
    return Mined(
        value=label,
        sources=tuple(sources),
        method="derived: box angles/lengths",
        confidence=Confidence.DERIVED,
    )
