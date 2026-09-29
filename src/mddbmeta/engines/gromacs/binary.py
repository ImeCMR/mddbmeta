"""GROMACS binary files, read with the stdlib only (big-endian XDR).

xtc: frames of  magic(1995) natoms step time box[9]  natoms, then either
     natoms*3 floats (natoms <= 9) or precision, minint[3], maxint[3],
     smallidx, nbytes and nbytes of compressed data padded to 4 bytes.
     We walk the frame headers and skip the coordinates.
trr: frames of  magic(1993), a version string, 13 block sizes, then t and
     lambda (float or double) and the data blocks.
tpr: a version string ("VERSION 2024.2"), precision, file version,
     generation, a file tag ("release") and the atom count.
cpt: checkpoint magic 171817 (recognized; state is read from the log instead).
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

XTC_MAGIC = 1995
TRR_MAGIC = 1993
CPT_MAGIC = 171817
EDR_MAGIC = -55555

NM_TO_A = 10.0


class GmxBinaryError(ValueError):
    pass


def first_int(head: bytes) -> int | None:
    if len(head) < 4:
        return None
    return struct.unpack(">i", head[:4])[0]


def _pad4(n: int) -> int:
    return (n + 3) // 4 * 4


@dataclass
class Frames:
    path: str
    format: str  # "xtc" | "trr"
    natom: int | None = None
    n_frames: int = 0
    first_step: int | None = None
    last_step: int | None = None
    times: list[float] = field(default_factory=list)
    box_vectors_A: tuple[tuple[float, float, float], ...] | None = None
    truncated: bool = False

    @property
    def first_time_ps(self) -> float | None:
        return self.times[0] if self.times else None

    @property
    def last_time_ps(self) -> float | None:
        return self.times[-1] if self.times else None

    @property
    def median_interval_ps(self) -> float | None:
        deltas = sorted(b - a for a, b in zip(self.times, self.times[1:], strict=False))
        if not deltas:
            return None
        mid = len(deltas) // 2
        return deltas[mid] if len(deltas) % 2 else (deltas[mid - 1] + deltas[mid]) / 2


def _box(values: tuple[float, ...]) -> tuple[tuple[float, float, float], ...]:
    v = [x * NM_TO_A for x in values]
    return ((v[0], v[1], v[2]), (v[3], v[4], v[5]), (v[6], v[7], v[8]))


def read_xtc(path: str | Path) -> Frames:
    out = Frames(path=str(path), format="xtc")
    size = Path(path).stat().st_size
    with open(path, "rb") as fh:
        while True:
            head = fh.read(52)
            if not head:
                break
            if len(head) < 52:
                out.truncated = True
                break
            magic, natoms, step, time = struct.unpack(">iiif", head[:16])
            if magic != XTC_MAGIC:
                if out.n_frames == 0:
                    raise GmxBinaryError(f"not an xtc file (magic {magic})")
                out.truncated = True
                break
            if out.box_vectors_A is None:
                out.box_vectors_A = _box(struct.unpack(">9f", head[16:52]))
            n2 = fh.read(4)
            if len(n2) < 4:
                out.truncated = True
                break
            if natoms <= 9:
                skip = natoms * 12
            else:
                meta = fh.read(36)  # precision, minint[3], maxint[3], smallidx, nbytes
                if len(meta) < 36:
                    out.truncated = True
                    break
                nbytes = struct.unpack(">i", meta[32:36])[0]
                skip = _pad4(nbytes)
            if fh.tell() + skip > size:  # last frame cut short
                out.truncated = True
                break
            fh.seek(skip, 1)
            out.natom = natoms
            out.n_frames += 1
            out.times.append(round(float(time), 6))
            if out.first_step is None:
                out.first_step = step
            out.last_step = step
    return out


def read_trr(path: str | Path) -> Frames:
    out = Frames(path=str(path), format="trr")
    size = Path(path).stat().st_size
    with open(path, "rb") as fh:
        while True:
            head = fh.read(12)
            if not head:
                break
            if len(head) < 12:
                out.truncated = True
                break
            magic, _slen, strlen = struct.unpack(">iii", head)
            if magic != TRR_MAGIC:
                if out.n_frames == 0:
                    raise GmxBinaryError(f"not a trr file (magic {magic})")
                out.truncated = True
                break
            fh.seek(_pad4(strlen), 1)
            sizes_raw = fh.read(13 * 4)
            if len(sizes_raw) < 52:
                out.truncated = True
                break
            (ir, e, box, vir, pres, top, sym, x, v, f, natoms, step, _nre) = struct.unpack(
                ">13i", sizes_raw
            )
            real = 8 if (box and box // 9 == 8) or (x and natoms and x // (natoms * 3) == 8) else 4
            tl = fh.read(2 * real)
            if len(tl) < 2 * real:
                out.truncated = True
                break
            time = struct.unpack(">d" if real == 8 else ">f", tl[:real])[0]
            if box and out.box_vectors_A is None:
                raw = fh.read(box)
                out.box_vectors_A = _box(
                    struct.unpack(f">9{'d' if real == 8 else 'f'}", raw[: 9 * real])
                )
                rest = ir + e + vir + pres + top + sym + x + v + f
            else:
                rest = ir + e + box + vir + pres + top + sym + x + v + f
            if fh.tell() + rest > size:
                out.truncated = True
                break
            fh.seek(rest, 1)
            out.natom = natoms
            out.n_frames += 1
            out.times.append(round(float(time), 6))
            if out.first_step is None:
                out.first_step = step
            out.last_step = step
    return out


@dataclass
class TprHeader:
    path: str
    version: str  # e.g. "2024.2" (from "VERSION 2024.2")
    precision: int
    file_version: int
    natom: int | None


def _xdr_string(buf: bytes, off: int) -> tuple[str, int]:
    """GROMACS writes strings as int(size incl. NUL) + XDR string (int len + bytes)."""
    _size, length = struct.unpack(">ii", buf[off : off + 8])
    if not 0 <= length < 4096:
        raise GmxBinaryError("bad string length")
    start = off + 8
    text = buf[start : start + length].decode("latin-1").rstrip("\x00")
    return text, start + _pad4(length)


def read_tpr_header(path: str | Path) -> TprHeader:
    with open(path, "rb") as fh:
        buf = fh.read(512)
    text, off = _xdr_string(buf, 0)
    if not text.startswith("VERSION"):
        raise GmxBinaryError("not a tpr file (no VERSION string)")
    version = text[len("VERSION") :].strip()
    precision, file_version = struct.unpack(">ii", buf[off : off + 8])
    off += 8
    natom = None
    try:
        # file version >= 77: generation, then the tag string, then natoms
        _gen = struct.unpack(">i", buf[off : off + 4])[0]
        tag, off2 = _xdr_string(buf, off + 4)
        if tag.isprintable() and tag:
            natom = struct.unpack(">i", buf[off2 : off2 + 4])[0]
        else:
            raise GmxBinaryError("no tag")
    except (GmxBinaryError, struct.error):
        try:
            tag, off2 = _xdr_string(buf, off)
            natom = struct.unpack(">i", buf[off2 + 4 : off2 + 8])[0]
        except (GmxBinaryError, struct.error):
            natom = None
    if natom is not None and not 0 < natom < 10**9:
        natom = None
    return TprHeader(
        path=str(path), version=version, precision=precision, file_version=file_version, natom=natom
    )


def looks_like_tpr(head: bytes) -> bool:
    return len(head) > 16 and head[8:15] == b"VERSION"
