"""AMBER ASCII coordinate files: restart/inpcrd (rst7) and trajectory (mdcrd).

rst7:  line 1 title; line 2 ``NATOM [TIME]``; coordinates as 6F12.7 (two
       atoms per line); optional velocities (same layout); optional final
       box line (a, b, c, alpha, beta, gamma).
mdcrd: line 1 title; then frames of 10F8.3 values, each frame optionally
       followed by one box line of 3F8.3.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from mddbmeta.io.safe import open_text

_FIELD12 = 12


class Rst7FormatError(ValueError):
    pass


@dataclass
class Rst7:
    path: str
    title: str
    natom: int
    time_ps: float | None
    has_velocities: bool
    box: tuple[float, float, float, float, float, float] | None


def _fixed_floats(line: str, width: int) -> list[float]:
    line = line.rstrip("\n")
    vals = []
    for i in range(0, len(line), width):
        chunk = line[i : i + width].strip()
        if chunk:
            vals.append(float(chunk))
    return vals


def read_rst7(path: str | Path) -> Rst7:
    with open_text(path) as fh:
        title = fh.readline().rstrip("\n")
        header = fh.readline()
        parts = header.split()
        if not parts:
            raise Rst7FormatError("missing NATOM line")
        try:
            natom = int(parts[0])
        except ValueError as exc:
            raise Rst7FormatError(f"line 2 is not 'NATOM [TIME]': {header.strip()!r}") from exc
        time_ps = None
        if len(parts) > 1:
            try:
                time_ps = float(parts[1].replace("D", "E"))
            except ValueError:
                time_ps = None
        n_data = 0
        last = ""
        for line in fh:
            if line.strip():
                n_data += 1
                last = line
    coord_lines = math.ceil(3 * natom / 6)
    has_vel = n_data >= 2 * coord_lines
    expected = coord_lines * (2 if has_vel else 1)
    box = None
    if n_data == expected + 1:
        vals = _fixed_floats(last, _FIELD12)
        if len(vals) != 6:
            vals = [float(x) for x in last.split()]
        if len(vals) >= 3:
            vals = (vals + [90.0, 90.0, 90.0])[:6]
            box = tuple(vals)  # type: ignore[assignment]
    elif n_data < coord_lines:
        raise Rst7FormatError(
            f"file ends early: {n_data} data lines, {coord_lines} needed for {natom} atoms"
        )
    return Rst7(
        path=str(path),
        title=title.strip(),
        natom=natom,
        time_ps=time_ps,
        has_velocities=has_vel,
        box=box,  # type: ignore[arg-type]
    )


@dataclass
class AsciiTraj:
    path: str
    title: str
    n_data_lines: int
    fields_histogram: dict[int, int]  # number of 8-wide fields on a line -> line count


def read_ascii_traj(path: str | Path) -> AsciiTraj:
    hist: dict[int, int] = {}
    n = 0
    with open_text(path) as fh:
        title = fh.readline().rstrip("\n")
        for line in fh:
            body = line.rstrip()
            if not body:
                continue
            n += 1
            k = math.ceil(len(body) / 8)
            hist[k] = hist.get(k, 0) + 1
    return AsciiTraj(path=str(path), title=title.strip(), n_data_lines=n, fields_histogram=hist)


def _frame_pattern(natom: int, box: bool) -> dict[int, int]:
    full, rem = divmod(3 * natom, 10)
    pattern: dict[int, int] = {}
    if full:
        pattern[10] = full
    if rem:
        pattern[rem] = pattern.get(rem, 0) + 1
    if box:
        pattern[3] = pattern.get(3, 0) + 1
    return pattern


def ascii_traj_frames(traj: AsciiTraj, natom: int) -> tuple[int | None, bool | None]:
    """Frame count for a given atom count; also whether each frame carries a box line.

    Both hypotheses (with/without a box line per frame) are checked against the
    exact per-line field-count histogram, so they cannot be confused.
    """
    if natom <= 0:
        return None, None
    for box in (True, False):
        pattern = _frame_pattern(natom, box)
        per_frame = sum(pattern.values())
        if traj.n_data_lines % per_frame:
            continue
        frames = traj.n_data_lines // per_frame
        if {k: v * frames for k, v in pattern.items()} == traj.fields_histogram:
            return frames, box
    return None, None


_INT_ONLY = re.compile(r"^\s*\d+\s*$")


def sniff_ascii_coords(head: bytes) -> str | None:
    """'rst7', 'mdcrd' or None, from the second line's shape."""
    text = head.decode("latin-1", errors="replace")
    lines = text.splitlines()
    if len(lines) < 2:
        return None
    second = lines[1]
    parts = second.split()
    if not parts:
        return None
    if _INT_ONLY.match(parts[0]) and len(parts) <= 2:
        if len(parts) == 2:
            try:
                float(parts[1].replace("D", "E"))
            except ValueError:
                return None
        if len(lines) > 2 and lines[2].strip():
            try:
                vals = _fixed_floats(lines[2], _FIELD12)
            except ValueError:
                return None
            if 1 <= len(vals) <= 6:
                return "rst7"
            return None
        return "rst7"
    try:
        vals = _fixed_floats(second, 8)
    except ValueError:
        return None
    if len(vals) >= 3 and all(
        "." in second[i : i + 8] for i in range(0, min(len(second.rstrip()), 24), 8)
    ):
        return "mdcrd"
    return None
