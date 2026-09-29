"""GROMACS text formats: .gro coordinates and .top topologies.

gro: title line; atom count; one fixed-width line per atom
     (resid 5, resname 5, atom name 5, atom id 5, x y z in nm, optional velocities);
     box line of 3 (rectangular) or 9 values (v1x v2y v3z v1y v1z v2x v2z v3x v3y), nm.
top: C-preprocessed text; ``#include`` lines name force-field / water files,
     ``[ molecules ]`` lists the molecule counts.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

NM_TO_A = 10.0


class GroFormatError(ValueError):
    pass


@dataclass
class Gro:
    path: str
    title: str
    natom: int
    time_ps: float | None
    residue_counts: dict[str, int] = field(default_factory=dict)
    residue_atoms: dict[str, int] = field(default_factory=dict)
    residue_atom_names: dict[str, list[str]] = field(default_factory=dict)
    nres: int = 0
    box_vectors_A: tuple[tuple[float, float, float], ...] | None = None

    def box_lengths_angles(
        self,
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
        if self.box_vectors_A is None:
            return None
        return lengths_angles(self.box_vectors_A)


def lengths_angles(vectors: tuple[tuple[float, float, float], ...]):
    a, b, c = vectors

    def norm(v):
        return math.sqrt(sum(x * x for x in v))

    def angle(u, v):
        nu, nv = norm(u), norm(v)
        if nu == 0 or nv == 0:
            return 90.0
        cos = max(-1.0, min(1.0, sum(x * y for x, y in zip(u, v, strict=True)) / (nu * nv)))
        return round(math.degrees(math.acos(cos)), 4)

    return (norm(a), norm(b), norm(c)), (angle(b, c), angle(a, c), angle(a, b))


_TIME_IN_TITLE = re.compile(r"\bt=\s*(-?[\d.]+)")


def read_gro(path: str | Path) -> Gro:
    with open_text(path) as fh:
        title = fh.readline().rstrip("\n")
        count_line = fh.readline()
        try:
            natom = int(count_line.strip())
        except ValueError as exc:
            raise GroFormatError(f"line 2 is not an atom count: {count_line.strip()!r}") from exc
        m = _TIME_IN_TITLE.search(title)
        gro = Gro(
            path=str(path),
            title=title.strip(),
            natom=natom,
            time_ps=float(m.group(1)) if m else None,
        )
        prev_key = None
        cur_names: list[str] = []
        cur_name = None

        def close_residue() -> None:
            if cur_name is None:
                return
            gro.residue_counts[cur_name] = gro.residue_counts.get(cur_name, 0) + 1
            if cur_name not in gro.residue_atoms:
                gro.residue_atoms[cur_name] = len(cur_names)
                gro.residue_atom_names[cur_name] = list(cur_names)
            gro.nres += 1

        for i in range(natom):
            line = fh.readline()
            if not line:
                raise GroFormatError(f"file ends after {i} of {natom} atoms")
            key = (line[0:5], line[5:10])
            if key != prev_key:
                close_residue()
                cur_name = line[5:10].strip()
                cur_names = []
                prev_key = key
            cur_names.append(line[10:15].strip())
        close_residue()
        box_line = fh.readline().split()
    if len(box_line) >= 3:
        v = [float(x) * NM_TO_A for x in box_line[:9]] + [0.0] * (9 - min(len(box_line), 9))
        gro.box_vectors_A = ((v[0], v[3], v[4]), (v[5], v[1], v[6]), (v[7], v[8], v[2]))
    return gro


_GRO_ATOM = re.compile(r"^[ \d]{5}.{5}.{5}[ \d]{5}\s*-?\d+\.\d{3}")


def looks_like_gro(head: bytes) -> bool:
    lines = head.decode("latin-1", errors="replace").splitlines()
    if len(lines) < 3 or not lines[1].strip().isdigit():
        return False
    return bool(_GRO_ATOM.match(lines[2]))


# ---------------------------------------------------------------- .top

_INCLUDE = re.compile(r'^\s*#include\s+"([^"]+)"', re.MULTILINE)
_SECTION = re.compile(r"^\s*\[\s*(\w+)\s*\]")


@dataclass
class Top:
    path: str
    includes: list[str] = field(default_factory=list)
    molecules: list[tuple[str, int]] = field(default_factory=list)
    system: str | None = None
    moleculetypes: list[str] = field(default_factory=list)


def read_top(path: str | Path) -> Top:
    top = Top(path=str(path))
    section = None
    with open_text(path) as fh:
        text = fh.read()
    top.includes = _INCLUDE.findall(text)
    expect_moltype = False
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        m = _SECTION.match(line)
        if m:
            section = m.group(1).lower()
            expect_moltype = section == "moleculetype"
            continue
        if section == "molecules":
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit():
                top.molecules.append((parts[0], int(parts[1])))
        elif section == "system" and top.system is None:
            top.system = line
        elif expect_moltype:
            top.moleculetypes.append(line.split()[0])
            expect_moltype = False
    return top


def looks_like_top(head: bytes, name: str) -> bool:
    text = head.decode("latin-1", errors="replace")
    if text.startswith("%VERSION") or text.startswith("%FLAG"):
        return False
    low = name.lower()
    if not low.endswith(".top"):
        return False
    return (
        "#include" in text or "[ moleculetype" in text or "[ system" in text or "[ defaults" in text
    )


# force-field directory stem -> label
FORCE_FIELD_DIRS = {
    "amber99sb-ildn": "Amber ff99SB-ILDN", "amber99sb": "Amber ff99SB", "amber99": "Amber ff99",
    "amber03": "Amber ff03", "amber94": "Amber ff94", "amber96": "Amber ff96", "ambergs": "Amber GS",
    "amber14sb": "Amber ff14SB", "amber14sb_parmbsc1": "Amber ff14SB + parmbsc1", "amber19sb": "Amber ff19SB",
    "charmm27": "CHARMM27", "charmm36": "CHARMM36", "charmm36m": "CHARMM36m",
    "oplsaa": "OPLS-AA/L", "gromos43a1": "GROMOS 43A1", "gromos43a2": "GROMOS 43A2",
    "gromos45a3": "GROMOS 45A3", "gromos53a5": "GROMOS 53A5", "gromos53a6": "GROMOS 53A6",
    "gromos54a7": "GROMOS 54A7", "martini_v3.0.0": "Martini 3", "martini3": "Martini 3",
    "martini_v2.2": "Martini 2.2",
}  # fmt: skip
WATER_ITPS = {
    "tip3p": "TIP3P", "tips3p": "TIPS3P", "tip4p": "TIP4P", "tip4pew": "TIP4P-Ew", "tip4p2005": "TIP4P/2005",
    "tip5p": "TIP5P", "spc": "SPC", "spce": "SPC/E", "opc": "OPC", "opc3": "OPC3", "tip3p-fb": "TIP3P-FB",
    "tip4p-fb": "TIP4P-FB", "tip4pd": "TIP4P-D",
}  # fmt: skip


def force_field_label(ff_dir: str) -> str:
    stem = ff_dir.lower().removesuffix(".ff")
    if stem in FORCE_FIELD_DIRS:
        return FORCE_FIELD_DIRS[stem]
    for key in sorted(FORCE_FIELD_DIRS, key=len, reverse=True):
        if stem.startswith(key):
            variant = ff_dir[len(key) :].strip("-_. ").removesuffix(".ff")
            return f"{FORCE_FIELD_DIRS[key]} ({variant})" if variant else FORCE_FIELD_DIRS[key]
    return ff_dir.removesuffix(".ff")


def evidence_from_includes(includes: list[str]) -> tuple[list[str], list[str]]:
    """(force-field labels, water models) named by #include lines."""
    ffs: list[str] = []
    waters: list[str] = []
    for inc in includes:
        parts = inc.replace("\\", "/").split("/")
        base = parts[-1].lower()
        stem = base.rsplit(".", 1)[0]
        if len(parts) >= 2 and parts[-2].lower().endswith(".ff") and base == "forcefield.itp":
            label = force_field_label(parts[-2])
            if label not in ffs:
                ffs.append(label)
        if stem in WATER_ITPS and WATER_ITPS[stem] not in waters:
            waters.append(WATER_ITPS[stem])
    return ffs, waters
