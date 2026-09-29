"""GROMACS mdrun log (md.log) reader.

A log may hold several *sessions*: ``mdrun -cpi`` in append mode writes a
new banner into the same file (after truncating it to the checkpoint). Each
session records the command line, the GROMACS version, the resolved run
parameters ("Input Parameters", including grpopts such as ``ref-t``), the
checkpoint it continued from, and the energy records ("Step Time" blocks
followed by an "Energies" table).
"""

from __future__ import annotations

import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta.engines.gromacs.mdp import norm_key, parse_value
from mddbmeta.io.safe import open_text

_BANNER = re.compile(r":-\)\s+GROMACS - gmx (\w+), (\S+)")
_VERSION = re.compile(r"^GROMACS version:\s+(\S+)")
_PARAM = re.compile(r"^ {3}([A-Za-z][\w\-]*)\s*=\s*(.*?)\s*$")  # top-level keys only
_GRPOPT = re.compile(r"^\s*([A-Za-z][\w\-]*):\s+(.*?)\s*$")
_STEP_TIME = re.compile(r"^\s+(-?\d+)\s+(-?[\d.]+(?:e[-+]?\d+)?)\s*$", re.IGNORECASE)
_CPT_STEP = re.compile(r"^\s*step:\s+(-?\d+)")
_CPT_TIME = re.compile(r"^\s*time:\s+(-?[\d.eE+-]+)")
_CPT_PART = re.compile(r"simulation part #:\s+(\d+)")
_NUM = re.compile(r"[-+]?\d+\.\d+e[-+]\d+|[-+]?\d+\.\d*|[-+]?\d+", re.IGNORECASE)

GRPOPTS = {"nrdf", "ref-t", "tau-t", "annealing", "annealing-npoints"}
MIN_INTEGRATORS = {"steep", "cg", "l-bfgs", "nm", "tpi", "tpic"}
# energy terms that describe the physical state; used to recognise continuations
STATE_TERMS = ("Bond", "Angle", "Proper Dih.", "Ryckaert-Bell.", "Improper Dih.", "Per. Imp. Dih.",
               "LJ-14", "Coulomb-14", "LJ (SR)", "Coulomb (SR)", "Coul. recip.", "Potential",
               "Kinetic En.", "Total Energy", "Temperature")  # fmt: skip


@dataclass
class EnergyRecord:
    step: int
    time_ps: float
    energies: dict[str, str] = field(default_factory=dict)  # term -> printed value


@dataclass
class Session:
    tool: str | None = None
    version: str | None = None
    executable: str | None = None
    command: list[str] = field(default_factory=list)
    params: dict[str, Any] = field(default_factory=dict)
    checkpoint_file: str | None = None
    checkpoint_step: int | None = None
    checkpoint_time_ps: float | None = None
    simulation_part: int | None = None
    appending: bool = False
    first: EnergyRecord | None = None
    last: EnergyRecord | None = None
    n_records: int = 0
    finished: bool = False
    converged: bool | None = None  # minimizers
    performance_ns_day: float | None = None

    def option(self, flag: str) -> str | None:
        """Value of a command-line option (e.g. '-deffnm'), or None."""
        for i, tok in enumerate(self.command):
            if (
                tok == flag
                and i + 1 < len(self.command)
                and not self.command[i + 1].startswith("-")
            ):
                return self.command[i + 1]
        return None

    def has_flag(self, flag: str) -> bool:
        return flag in self.command


@dataclass
class GmxLog:
    path: str
    sessions: list[Session] = field(default_factory=list)

    @property
    def first(self) -> Session:
        return self.sessions[0]

    @property
    def last(self) -> Session:
        return self.sessions[-1]

    @property
    def params(self) -> dict[str, Any]:
        """Resolved run parameters; later sessions win (convert-tpr -extend changes nsteps)."""
        merged: dict[str, Any] = {}
        for s in self.sessions:
            merged.update(s.params)
        return merged

    @property
    def minimization(self) -> bool:
        return str(self.params.get("integrator", "")).lower() in MIN_INTEGRATORS

    def records(self) -> tuple[EnergyRecord | None, EnergyRecord | None]:
        firsts = [s.first for s in self.sessions if s.first is not None]
        lasts = [s.last for s in self.sessions if s.last is not None]
        return (firsts[0] if firsts else None, lasts[-1] if lasts else None)


def _parse_energy_block(lines: list[str]) -> dict[str, str]:
    """Header lines hold 15-char-wide term names; value lines hold numbers."""
    out: dict[str, str] = {}
    for header, values in zip(lines[0::2], lines[1::2], strict=False):
        names = [header[i : i + 15].strip() for i in range(0, len(header.rstrip()), 15)]
        nums = _NUM.findall(values)
        for name, num in zip(names, nums, strict=False):
            if name:
                out[name] = num
    return out


def read_log(path: str | Path) -> GmxLog:
    log = GmxLog(path=str(path))
    cur: Session | None = None
    mode = None  # "command", "params", "grpopts", "checkpoint"
    pending_step: tuple[int, float] | None = None
    expect_step = False
    energy_lines: list[str] = []
    in_energy = False
    in_averages = False
    next_appending = False
    carry: tuple | None = None

    def flush_energy() -> None:
        nonlocal energy_lines, in_energy, pending_step
        if cur is not None and pending_step is not None and not in_averages:
            rec = EnergyRecord(pending_step[0], pending_step[1], _parse_energy_block(energy_lines))
            if cur.first is None:
                cur.first = rec
            cur.last = rec
            cur.n_records += 1
        energy_lines, in_energy, pending_step = [], False, None

    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            m = _BANNER.search(line)
            if m:
                if in_energy:
                    flush_energy()
                cur = Session(tool=m.group(1), version=m.group(2), appending=next_appending)
                if next_appending and carry is not None:
                    (cur.checkpoint_file, cur.checkpoint_step, cur.checkpoint_time_ps,
                     cur.simulation_part) = carry  # fmt: skip
                next_appending, carry = False, None
                log.sessions.append(cur)
                mode, in_averages = None, False
                continue
            if cur is None:
                continue
            stripped = line.strip()

            if in_energy:
                if stripped == "" and energy_lines:
                    flush_energy()
                elif stripped and not stripped.startswith("Energies"):
                    energy_lines.append(line)
                continue

            if mode == "command":
                if stripped:
                    try:
                        cur.command = shlex.split(stripped)
                    except ValueError:
                        cur.command = stripped.split()
                mode = None
                continue
            if mode == "params":
                # the block ends at the first blank line; matrices and sub-sections are skipped
                if stripped == "":
                    mode = None
                elif stripped == "grpopts:":
                    mode = "grpopts"
                else:
                    pm = _PARAM.match(line)
                    if pm:
                        cur.params.setdefault(norm_key(pm.group(1)), parse_value(pm.group(2)))
                continue
            if mode == "grpopts":
                gm = _GRPOPT.match(line)
                if stripped == "" or gm is None:
                    mode = None
                elif norm_key(gm.group(1)) in GRPOPTS:
                    cur.params[norm_key(gm.group(1))] = parse_value(gm.group(2))
                if mode is not None:
                    continue
            if mode == "checkpoint":
                sm = _CPT_STEP.match(line)
                if sm:
                    cur.checkpoint_step = int(sm.group(1))
                tm = _CPT_TIME.match(line)
                if tm:
                    cur.checkpoint_time_ps = float(tm.group(1))
                    mode = None
                pm2 = _CPT_PART.search(line)
                if pm2:
                    cur.simulation_part = int(pm2.group(1))
                continue

            if line.startswith("Executable:"):
                cur.executable = line.split(":", 1)[1].strip()
            elif line.startswith("Command line:"):
                mode = "command"
            elif _VERSION.match(line):
                cur.version = _VERSION.match(line).group(1)  # type: ignore[union-attr]
            elif stripped == "Input Parameters:":
                mode = "params"
            elif stripped.startswith("Reading checkpoint file"):
                cur.checkpoint_file = stripped.split("Reading checkpoint file", 1)[1].split()[0]
                mode = "checkpoint"
            elif "appending to previous log file" in line:
                next_appending = True
                carry = (
                    cur.checkpoint_file,
                    cur.checkpoint_step,
                    cur.checkpoint_time_ps,
                    cur.simulation_part,
                )
                cur.checkpoint_file = cur.checkpoint_step = cur.checkpoint_time_ps = None
                cur.simulation_part = None
            elif "A V E R A G E S" in line:
                in_averages = True
            elif (
                stripped.startswith("Step")
                and "Time" in stripped
                and stripped.split() == ["Step", "Time"]
            ):
                expect_step = True
            elif expect_step:
                sm2 = _STEP_TIME.match(line)
                expect_step = False
                if sm2 and not in_averages:
                    pending_step = (int(sm2.group(1)), float(sm2.group(2)))
            elif stripped.startswith("Energies (") and pending_step is not None:
                in_energy = True
                energy_lines = []
            elif stripped.startswith("Finished mdrun"):
                cur.finished = True
            elif " converged to " in line:
                cur.converged = True
            elif "did not converge" in line:
                cur.converged = False
            elif stripped.startswith("Performance:"):
                nums = _NUM.findall(stripped)
                if nums:
                    cur.performance_ns_day = float(nums[0])
        if in_energy:
            flush_energy()
    return log


def looks_like_log(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    return bool(_BANNER.search(text)) and "gmx mdrun" in text
