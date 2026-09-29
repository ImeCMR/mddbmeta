"""NAMD log (stdout) reader.

NAMD prints its resolved configuration as ``Info:`` lines (``TIMESTEP 2``,
``FIRST TIMESTEP 1500``, ``BINARY COORDINATES md_1.coor``, ``LANGEVIN
TEMPERATURE 300``...), Tcl actions as ``TCL:`` lines (``Running for 1500
steps``, ``Setting parameter langevinTemp to 100``), and energies as
``ETITLE:`` / ``ENERGY:`` rows. A failed run still ends with ``WallClock``
and ``End of program``, so ``FATAL ERROR`` lines must be checked.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

_VERSION = re.compile(r"^Info: NAMD (\S+) for (\S+)")
_INFO_KV = re.compile(r"^Info: ([A-Z][A-Z0-9 ()/_\-.]*?[A-Z0-9)])\s{2,}(\S.*?)\s*$")
_INFO_FLAG = re.compile(r"^Info: ([A-Z][A-Z0-9 ()/_\-.]*?(?:ACTIVE|ON|USED|DATA))\s*$")
_ATOMS = re.compile(r"^Info: (\d+) ATOMS\s*$")
_TCL_RUN = re.compile(r"^TCL: (Running|Minimizing) for (\d+) steps")
_TCL_SET = re.compile(r"^TCL: Setting parameter (\S+) to (.+?)\s*$")
_WALL = re.compile(r"^WallClock:\s*([\d.]+)")


@dataclass
class EnergyRow:
    step: int
    values: dict[str, str]


@dataclass
class NamdLog:
    path: str
    version: str | None = None
    platform: str | None = None
    info: dict[str, list[str]] = field(default_factory=dict)
    flags: set[str] = field(default_factory=set)
    natom: int | None = None
    runs: list[tuple[str, int]] = field(default_factory=list)
    parameter_changes: list[tuple[str, str]] = field(default_factory=list)
    first: EnergyRow | None = None
    last: EnergyRow | None = None
    n_energy: int = 0
    fatal: list[str] = field(default_factory=list)
    wallclock_s: float | None = None
    wrote_output: bool = False
    gpu: bool = False

    def one(self, key: str) -> str | None:
        values = self.info.get(key)
        return values[0] if values else None

    def all(self, key: str) -> list[str]:
        return self.info.get(key, [])

    @property
    def completed(self) -> bool:
        return not self.fatal and self.wallclock_s is not None

    @property
    def dynamics_steps(self) -> int:
        return sum(n for kind, n in self.runs if kind == "Running")

    @property
    def minimization_steps(self) -> int:
        return sum(n for kind, n in self.runs if kind == "Minimizing")


def read_log(path: str | Path) -> NamdLog:
    log = NamdLog(path=str(path))
    titles: list[str] | None = None
    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("Info:"):
                if log.version is None:
                    m = _VERSION.match(line)
                    if m:
                        log.version, log.platform = m.group(1), m.group(2)
                        continue
                m = _ATOMS.match(line)
                if m and log.natom is None:
                    log.natom = int(m.group(1))
                    continue
                m = _INFO_KV.match(line)
                if m:
                    log.info.setdefault(m.group(1), []).append(m.group(2))
                    continue
                m = _INFO_FLAG.match(line)
                if m:
                    log.flags.add(m.group(1))
                if "CUDA" in line or "GPU" in line:
                    log.gpu = True
                continue
            if line.startswith("TCL:"):
                m = _TCL_RUN.match(line)
                if m:
                    log.runs.append((m.group(1), int(m.group(2))))
                    continue
                m = _TCL_SET.match(line)
                if m:
                    log.parameter_changes.append((m.group(1), m.group(2)))
                continue
            if line.startswith("ETITLE:"):
                titles = line.split()[1:]
                continue
            if line.startswith("ENERGY:"):
                # NAMD prints the first rows before any ETITLE header; label them later
                parts = line.split()[1:]
                try:
                    step = int(parts[0])
                except (ValueError, IndexError):
                    continue
                names = titles[1:] if titles else [f"#{i}" for i in range(1, len(parts))]
                row = EnergyRow(step, dict(zip(names, parts[1:], strict=False)))
                if log.first is None:
                    log.first = row
                log.last = row
                log.n_energy += 1
                continue
            if "FATAL ERROR" in line:
                log.fatal.append(line.strip())
            elif line.startswith("WallClock:"):
                m = _WALL.match(line)
                if m:
                    log.wallclock_s = float(m.group(1))
            elif line.startswith("WRITING COORDINATES TO OUTPUT FILE"):
                log.wrote_output = True
    if titles:
        for row in (log.first, log.last):
            if row is not None and any(k.startswith("#") for k in row.values):
                row.values = {titles[1 + int(k[1:]) - 1] if k.startswith("#") and int(k[1:]) < len(titles) else k: v
                              for k, v in row.values.items()}  # fmt: skip
    return log


def looks_like_log(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    return "Info: NAMD" in text or ("Charm++" in text and "NAMD" in text)
