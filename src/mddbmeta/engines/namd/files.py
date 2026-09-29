"""NAMD / CHARMM data files, read with the stdlib only.

dcd:  Fortran unformatted records. Header: 'CORD', NSET, ISTART (first frame
      step), NSAVC (steps between frames), NSTEP, ..., DELTA (timestep in AKMA
      units, 1 AKMA = 48.88821 fs), ..., flag[10] (unit cell present), version.
      Then a title record, the atom count, and per frame an optional 6-double
      unit-cell record plus X, Y and Z float records.
coor: NAMD binary coordinates/velocities: int32 atom count + 3N doubles.
xsc:  text; '#$LABELS step a_x a_y a_z b_x ... o_x o_y o_z ...' then a data line.
psf:  '!NATOM' section: id segid resid resname name type charge mass.
pdb:  ATOM/HETATM records and an optional CRYST1 box.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

AKMA_FS = 48.88821


class NamdFileError(ValueError):
    pass


def _endian(first4: bytes, expect: int) -> str | None:
    if len(first4) < 4:
        return None
    if struct.unpack("<i", first4)[0] == expect:
        return "<"
    if struct.unpack(">i", first4)[0] == expect:
        return ">"
    return None


def is_dcd(head: bytes) -> bool:
    return _endian(head[:4], 84) is not None and head[4:8] == b"CORD"


@dataclass
class Dcd:
    path: str
    natom: int
    n_frames: int
    first_step: int
    step_interval: int
    timestep_fs: float | None
    has_cell: bool
    truncated: bool
    cell_A: tuple[float, float, float] | None = None

    def time_ps(self, frame: int) -> float | None:
        if self.timestep_fs is None:
            return None
        return round((self.first_step + frame * self.step_interval) * self.timestep_fs / 1000.0, 6)


def read_dcd(path: str | Path) -> Dcd:
    size = Path(path).stat().st_size
    with open(path, "rb") as fh:
        head = fh.read(92)
        e = _endian(head[:4], 84)
        if e is None or head[4:8] != b"CORD":
            raise NamdFileError("not a DCD file")
        icntrl = struct.unpack(f"{e}20i", head[8:88])
        delta = struct.unpack(f"{e}f", head[44:48])[0]
        nset, istart, nsavc = icntrl[0], icntrl[1], icntrl[2]
        has_cell = icntrl[10] == 1
        # title record
        tlen = struct.unpack(f"{e}i", fh.read(4))[0]
        fh.seek(tlen + 4, 1)
        natom_rec = fh.read(12)
        natom = struct.unpack(f"{e}iii", natom_rec)[1]
        header_size = fh.tell()
        cell = None
        if has_cell and size > header_size + 56:
            rec = fh.read(56)
            v = struct.unpack(f"{e}6d", rec[4:52])
            cell = (v[0], v[2], v[5])  # A, B, C (angles interleaved)
    frame_size = (56 if has_cell else 0) + 3 * (8 + 4 * natom)
    body = size - header_size
    n_frames, rem = divmod(body, frame_size) if frame_size else (0, 0)
    timestep_fs = round(delta * AKMA_FS, 6) if delta > 0 else None
    return Dcd(
        path=str(path),
        natom=natom,
        n_frames=n_frames,
        first_step=istart,
        step_interval=nsavc,
        timestep_fs=timestep_fs,
        has_cell=has_cell,
        truncated=rem != 0 or (nset not in (0, n_frames)),
        cell_A=cell,
    )


def binary_natom(head: bytes, size: int) -> int | None:
    """Atom count of a NAMD binary .coor/.vel file, if the size matches."""
    for e in ("<", ">"):
        if len(head) < 4:
            return None
        n = struct.unpack(f"{e}i", head[:4])[0]
        if n > 0 and size == 4 + 24 * n:
            return n
    return None


@dataclass
class NamdBinary:
    path: str
    natom: int
    velocities: bool


@dataclass
class Xsc:
    path: str
    step: int
    vectors_A: tuple[tuple[float, float, float], ...] | None
    origin: tuple[float, float, float] | None


def read_xsc(path: str | Path) -> Xsc:
    labels: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#$LABELS"):
                labels = s.split()[1:]
            elif s and not s.startswith("#"):
                vals = s.split()
                data = dict(zip(labels, vals, strict=False)) if labels else {}
                step = int(float(data.get("step", vals[0])))
                try:
                    vec = tuple(tuple(float(data[f"{v}_{c}"]) for c in "xyz") for v in "abc")
                    origin = tuple(float(data[f"o_{c}"]) for c in "xyz")
                except KeyError:
                    vec, origin = None, None
                return Xsc(path=str(path), step=step, vectors_A=vec, origin=origin)  # type: ignore[arg-type]
    raise NamdFileError("no data line in the extended-system file")


def looks_like_xsc(head: bytes) -> bool:
    return head.startswith(b"# NAMD extended system")


@dataclass
class Psf:
    path: str
    natom: int
    residue_counts: dict[str, int] = field(default_factory=dict)
    residue_atoms: dict[str, int] = field(default_factory=dict)
    residue_atom_names: dict[str, list[str]] = field(default_factory=dict)
    nres: int = 0
    segments: list[str] = field(default_factory=list)


def read_psf(path: str | Path) -> Psf:
    with open_text(path) as fh:
        first = fh.readline()
        if not first.startswith("PSF"):
            raise NamdFileError("not a PSF file")
        natom = None
        for line in fh:
            if "!NATOM" in line:
                natom = int(line.split()[0])
                break
        if natom is None:
            raise NamdFileError("no !NATOM section")
        psf = Psf(path=str(path), natom=natom)
        prev = None
        names: list[str] = []
        cur = None

        def close() -> None:
            if cur is None:
                return
            psf.residue_counts[cur] = psf.residue_counts.get(cur, 0) + 1
            if cur not in psf.residue_atoms:
                psf.residue_atoms[cur] = len(names)
                psf.residue_atom_names[cur] = list(names)
            psf.nres += 1

        for _ in range(natom):
            parts = fh.readline().split()
            if len(parts) < 5:
                raise NamdFileError("atom section ends early")
            seg, resid, resname, name = parts[1], parts[2], parts[3], parts[4]
            if seg not in psf.segments:
                psf.segments.append(seg)
            key = (seg, resid, resname)
            if key != prev:
                close()
                cur, names, prev = resname, [], key
            names.append(name)
        close()
    return psf


@dataclass
class Pdb:
    path: str
    natom: int
    cryst1_A: tuple[float, float, float, float, float, float] | None
    residue_counts: dict[str, int] = field(default_factory=dict)
    residue_atoms: dict[str, int] = field(default_factory=dict)
    residue_atom_names: dict[str, list[str]] = field(default_factory=dict)
    nres: int = 0


def read_pdb(path: str | Path) -> Pdb:
    n = 0
    cryst = None
    counts: dict[str, int] = {}
    atoms: dict[str, int] = {}
    names: dict[str, list[str]] = {}
    prev = None
    cur: list[str] = []
    cur_name: str | None = None
    nres = 0

    def close() -> None:
        nonlocal nres
        if cur_name is None:
            return
        counts[cur_name] = counts.get(cur_name, 0) + 1
        if cur_name not in atoms:
            atoms[cur_name] = len(cur)
            names[cur_name] = list(cur)
        nres += 1

    with open_text(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                n += 1
                key = (line[21:22], line[22:27], line[17:21].strip())
                if key != prev:
                    close()
                    prev, cur_name, cur = key, line[17:21].strip(), []
                cur.append(line[12:16].strip())
            elif line.startswith("CRYST1") and cryst is None:
                try:
                    cryst = tuple(
                        float(line[i : i + w])
                        for i, w in ((6, 9), (15, 9), (24, 9), (33, 7), (40, 7), (47, 7))
                    )
                except ValueError:
                    cryst = None
            elif line.startswith("ENDMDL"):
                break
    close()
    return Pdb(path=str(path), natom=n, cryst1_A=cryst, residue_counts=counts,  # type: ignore[arg-type]
               residue_atoms=atoms, residue_atom_names=names, nres=nres)  # fmt: skip


def looks_like_pdb(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    lines = text.splitlines()[:40]
    return any(line.startswith(("ATOM  ", "HETATM", "CRYST1")) for line in lines)
