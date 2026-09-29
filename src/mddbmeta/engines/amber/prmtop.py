"""AMBER parameter/topology (prmtop, "new" %FLAG format) reader.

Format: a ``%VERSION`` line, then sections introduced by ``%FLAG <NAME>``
followed by ``%FORMAT(<count><type><width>[.<decimals>])`` and fixed-width
data lines. Only the sections mddbmeta needs are materialized.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

_FORMAT_RE = re.compile(r"%FORMAT\s*\(\s*(\d*)\s*([aAiIeEfF])\s*(\d+)(?:\.(\d+))?\s*\)")
_COMPOUND_FORMAT_RE = re.compile(r"%FORMAT\s*\(.*,.*\)")

# Index of each POINTERS entry (0-based), per the AMBER file-format spec.
POINTER_NAMES = [
    "NATOM", "NTYPES", "NBONH", "MBONA", "NTHETH", "MTHETA", "NPHIH", "MPHIA",
    "NHPARM", "NPARM", "NNB", "NRES", "NBONA", "NTHETA", "NPHIA", "NUMBND",
    "NUMANG", "NPTRA", "NATYP", "NPHB", "IFPERT", "NBPER", "NGPER", "NDPER",
    "MBPER", "MGPER", "MDPER", "IFBOX", "NMXRS", "IFCAP", "NUMEXTRA", "NCOPY",
]  # fmt: skip

DEFAULT_SECTIONS = frozenset(
    {
        "TITLE",
        "CTITLE",
        "POINTERS",
        "ATOM_NAME",
        "RESIDUE_LABEL",
        "RESIDUE_POINTER",
        "BOX_DIMENSIONS",
        "SOLVENT_POINTERS",
        "FORCE_FIELD_TYPE",
    }
)


class PrmtopFormatError(ValueError):
    pass


@dataclass
class Prmtop:
    path: str
    version_stamp: str | None = None
    title: str | None = None
    pointers: dict[str, int] = field(default_factory=dict)
    atom_names: list[str] = field(default_factory=list)
    residue_labels: list[str] = field(default_factory=list)
    residue_pointers: list[int] = field(default_factory=list)  # 1-based first atom of each residue
    box: tuple[float, float, float, float] | None = None  # (beta, a, b, c)
    solvent_pointers: tuple[int, int, int] | None = None  # (IPTRES, NSPM, NSPSOL)
    force_field_type: str | None = None  # present in CHAMBER (CHARMM-converted) topologies
    is_chamber: bool = False
    flags_seen: list[str] = field(default_factory=list)

    @property
    def natom(self) -> int | None:
        return self.pointers.get("NATOM")

    @property
    def nres(self) -> int | None:
        return self.pointers.get("NRES")

    @property
    def ifbox(self) -> int | None:
        return self.pointers.get("IFBOX")

    @property
    def numextra(self) -> int | None:
        return self.pointers.get("NUMEXTRA")

    def residue_sizes(self) -> list[int]:
        """Atom count of every residue, from RESIDUE_POINTER."""
        if not self.residue_pointers or self.natom is None:
            return []
        starts = self.residue_pointers
        ends = starts[1:] + [self.natom + 1]
        return [e - s for s, e in zip(starts, ends, strict=True)]

    def residue_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for name in self.residue_labels:
            counts[name] = counts.get(name, 0) + 1
        return counts

    def residue_templates(self) -> tuple[dict[str, int], dict[str, list[str]]]:
        """Atoms per residue and atom names for the first instance of each residue name."""
        sizes = self.residue_sizes()
        atoms: dict[str, int] = {}
        names: dict[str, list[str]] = {}
        for i, (label, size) in enumerate(zip(self.residue_labels, sizes, strict=False)):
            if label in atoms:
                continue
            atoms[label] = size
            start = self.residue_pointers[i] - 1
            names[label] = self.atom_names[start : start + size] if self.atom_names else []
        return atoms, names


def _split_fixed(line: str, width: int, count: int) -> list[str]:
    line = line.rstrip("\n").rstrip("\r")
    out = []
    for i in range(count):
        chunk = line[i * width : (i + 1) * width]
        if not chunk:
            break
        out.append(chunk)
    return out


def _convert(chunks: list[str], typ: str) -> list:
    typ = typ.lower()
    if typ == "a":
        return [c.strip() for c in chunks]
    values = []
    for c in chunks:
        c = c.strip()
        if not c:
            continue
        if typ == "i":
            values.append(int(c))
        else:
            values.append(float(c.replace("D", "E").replace("d", "e")))
    return values


def read_prmtop(path: str | Path, sections: frozenset[str] | None = DEFAULT_SECTIONS) -> Prmtop:
    """Read a prmtop. ``sections=None`` materializes every section."""
    top = Prmtop(path=str(path))
    data: dict[str, list] = {}
    current: str | None = None
    fmt: tuple[int, str, int] | None = None
    first = True
    with open_text(path) as fh:
        for line in fh:
            if first:
                first = False
                if line.startswith("%VERSION"):
                    top.version_stamp = line[len("%VERSION") :].strip()
                    continue
                if not line.startswith("%"):
                    raise PrmtopFormatError(
                        "no %VERSION/%FLAG header: old-format (pre-AMBER 7) topologies are unsupported"
                    )
            if line.startswith("%FLAG"):
                current = line[5:].strip().split()[0] if line[5:].strip() else ""
                top.flags_seen.append(current)
                fmt = None
                continue
            if line.startswith("%FORMAT"):
                m = _FORMAT_RE.match(line)
                if m:
                    fmt = (int(m.group(1) or 1), m.group(2), int(m.group(3)))
                elif _COMPOUND_FORMAT_RE.match(line):
                    fmt = (1, "a", 10**6)  # e.g. (i2,a78): keep each line as one string
                else:
                    raise PrmtopFormatError(f"unreadable format line for {current}: {line.strip()}")
                continue
            if line.startswith("%COMMENT") or line.startswith("%VERSION"):
                continue
            if current is None or (sections is not None and current not in sections):
                continue
            if fmt is None:
                raise PrmtopFormatError(f"section {current} has data before its %FORMAT line")
            count, typ, width = fmt
            try:
                data.setdefault(current, []).extend(_convert(_split_fixed(line, width, count), typ))
            except ValueError as exc:
                raise PrmtopFormatError(f"bad value in section {current}: {exc}") from exc

    if "POINTERS" not in data and "POINTERS" not in top.flags_seen:
        raise PrmtopFormatError("no %FLAG POINTERS section; not an AMBER topology")
    pointers = data.get("POINTERS", [])
    top.pointers = {
        name: int(v) for name, v in zip(POINTER_NAMES, pointers, strict=False)
    }  # old files carry fewer
    title = data.get("TITLE") or data.get("CTITLE")
    if title:
        top.title = "".join(title).strip() or None
    top.is_chamber = "CTITLE" in top.flags_seen or "FORCE_FIELD_TYPE" in top.flags_seen
    if data.get("FORCE_FIELD_TYPE"):
        top.force_field_type = " ".join(str(x) for x in data["FORCE_FIELD_TYPE"]).strip() or None
    top.atom_names = [str(x) for x in data.get("ATOM_NAME", [])]
    top.residue_labels = [str(x) for x in data.get("RESIDUE_LABEL", [])]
    top.residue_pointers = [int(x) for x in data.get("RESIDUE_POINTER", [])]
    box = data.get("BOX_DIMENSIONS")
    if box and len(box) >= 4:
        top.box = (float(box[0]), float(box[1]), float(box[2]), float(box[3]))
    solv = data.get("SOLVENT_POINTERS")
    if solv and len(solv) >= 3:
        top.solvent_pointers = (int(solv[0]), int(solv[1]), int(solv[2]))
    return top


def looks_like_prmtop(head: bytes) -> bool:
    text = head[:400].decode("latin-1", errors="replace")
    return text.startswith("%VERSION") or (text.startswith("%FLAG") and "%FORMAT" in text)
