"""AMBER output log (mdout) scanner, single pass, constant memory.

Extracts: engine identity (SANDER / PMEMD, release, executable), the
"File Assignments" block (which files the run actually read and wrote), the
echoed input, the "CONTROL DATA FOR THE RUN" values, atom count, box type,
the begin time read from the input coordinates, the first/last energy
records (outside the averages blocks), and whether the run finished.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

_ASSIGN_RE = re.compile(r"^\|?\s*([A-Z]+)\s*:\s*(\S+)\s*$")
_ENGINE_BANNER = re.compile(r"Amber\s+(\d+[\w.]*)\s+(SANDER|PMEMD)\b", re.IGNORECASE)
_RELEASE = re.compile(r"(PMEMD|SANDER)\b.*Release\s+(\d+[\w.]*)", re.IGNORECASE)
_EXEC_PATH = re.compile(r"Executable path:\s*(\S+)")
_NSTEP_TIME = re.compile(
    r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-\d.Ee+]+)(?:\s+TEMP\(K\)\s*=\s*([-\d.Ee+*]+))?"
)
_KV = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)")
_NATOM = re.compile(r"NATOM\s*=\s*(\d+)")
_BEGIN_TIME = re.compile(
    r"begin time read from input coords\s*=\s*([-\d.Ee+]+)\s*ps", re.IGNORECASE
)
_WALL = re.compile(r"(?:Master\s+)?Total wall time:\s*([\d.]+)\s*seconds", re.IGNORECASE)
_BOX_TYPE = re.compile(r"BOX TYPE:\s*(.+?)\s*$")
_SECTION = re.compile(r"^\s*(\d)\.\s+([A-Z][A-Z .&]+?):?\s*$")

_THERMOSTAT_HEADERS = (
    ("langevin", 3),
    ("berendsen", 1),
    ("anderson", 2),
    ("andersen", 2),
    ("bussi", 11),
    ("isokinetic", 9),
)

AVERAGE_MARKERS = ("A V E R A G E S", "R M S  F L U C T U A T I O N S", "F L U C T U A T I O N S")


@dataclass
class EnergyRecord:
    nstep: int
    time_ps: float | None
    temp_K: float | None


@dataclass
class Mdout:
    path: str
    program: str | None = None  # "PMEMD" | "SANDER"
    version: str | None = None  # e.g. "22"
    executable: str | None = None  # e.g. "pmemd.cuda_SPFP"
    gpu: bool = False
    mpi: bool = False
    assignments: dict[str, str] = field(
        default_factory=dict
    )  # MDIN, INPCRD, PARM, RESTRT, MDCRD...
    input_echo: str | None = None
    control: dict[str, float | int] = field(default_factory=dict)
    natom: int | None = None
    box_type: str | None = None
    begin_time_ps: float | None = None
    first_record: EnergyRecord | None = None
    last_record: EnergyRecord | None = None
    n_records: int = 0
    minimization: bool = False
    last_min_step: int | None = None
    completed: bool = False
    wall_time_s: float | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def variant(self) -> str | None:
        if self.executable:
            return Path(self.executable).name
        if self.program == "PMEMD" and self.gpu:
            return "pmemd.cuda"
        return None


def _num(text: str) -> float | int:
    t = text.replace("D", "E").replace("d", "e")
    try:
        return int(t)
    except ValueError:
        return float(t)


def read_mdout(path: str | Path) -> Mdout:
    out = Mdout(path=str(path))
    in_assign = False
    in_echo = False
    echo_lines: list[str] = []
    section: int | None = None
    skip_next_nstep = False
    expect_min_values = False
    with open_text(path) as fh:
        for line in fh:
            stripped = line.strip()

            if out.program is None:
                m = _ENGINE_BANNER.search(line)
                if m:
                    out.version = m.group(1)
                    out.program = m.group(2).upper()
            m = _RELEASE.search(line)
            if m and out.version is None:
                out.version = m.group(2)
                out.program = out.program or m.group(1).upper()
            if "Executable path:" in line:
                m = _EXEC_PATH.search(line)
                if m:
                    out.executable = m.group(1)
            if "CUDA" in line or "GPU DEVICE INFO" in line:
                out.gpu = True
            if "AMBER/MPI" in line or "Running AMBER/MPI" in line:
                out.mpi = True

            if stripped.startswith("File Assignments"):
                in_assign = True
                continue
            if in_assign:
                m = _ASSIGN_RE.match(stripped)
                if m:
                    out.assignments[m.group(1)] = m.group(2)
                    continue
                if stripped == "" and not out.assignments:
                    continue
                in_assign = False

            if stripped.startswith("Here is the input file:"):
                in_echo = True
                continue
            if in_echo:
                if _SECTION.match(line) or stripped.startswith("-----"):
                    in_echo = False
                    out.input_echo = "\n".join(echo_lines).strip() or None
                else:
                    echo_lines.append(line.rstrip("\n"))
                    continue

            sm = _SECTION.match(line)
            if sm:
                section = int(sm.group(1))
                continue

            if section == 1 or section is None or section == 2:
                if out.natom is None and "NATOM" in line:
                    m = _NATOM.search(line)
                    if m:
                        out.natom = int(m.group(1))
                if "BOX TYPE:" in line:
                    m = _BOX_TYPE.search(line)
                    if m:
                        out.box_type = m.group(1)
            if section == 2:
                for key, value in _KV.findall(line):
                    out.control.setdefault(key.lower(), _num(value))
                # pmemd reports the thermostat and MC barostat as section headers, not keys
                low = line.lower()
                if "temperature" in low and ("regulation" in low or "randomization" in low):
                    for word, ntt in _THERMOSTAT_HEADERS:
                        if word in low:
                            out.control.setdefault("ntt", ntt)
                            break
                if "monte-carlo barostat" in low or "monte carlo barostat" in low:
                    out.control.setdefault("barostat", 2)
            if "begin time read from input coords" in line:
                m = _BEGIN_TIME.search(line)
                if m:
                    out.begin_time_ps = float(m.group(1))

            if any(marker in line for marker in AVERAGE_MARKERS):
                skip_next_nstep = True
                continue
            if "NSTEP" in line:
                m = _NSTEP_TIME.search(line)
                if m:
                    if skip_next_nstep:
                        skip_next_nstep = False
                        continue
                    temp = m.group(3)
                    rec = EnergyRecord(
                        nstep=int(m.group(1)),
                        time_ps=float(m.group(2)),
                        temp_K=float(temp) if temp and "*" not in temp else None,
                    )
                    out.n_records += 1
                    if out.first_record is None:
                        out.first_record = rec
                    out.last_record = rec
                    continue
                if "ENERGY" in line and "RMS" in line:
                    out.minimization = True
                    expect_min_values = True
                    continue
            if expect_min_values and stripped:
                expect_min_values = False
                first = stripped.split()[0]
                if first.isdigit():
                    out.last_min_step = int(first)
                    out.n_records += 1
                continue

            if "wall time" in line.lower():
                m = _WALL.search(line)
                if m:
                    out.completed = True
                    out.wall_time_s = float(m.group(1))
            if "Final Performance Info" in line:
                out.completed = True
    if in_echo and echo_lines:
        out.input_echo = "\n".join(echo_lines).strip() or None
    return out


def looks_like_mdout(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    return bool(_ENGINE_BANNER.search(text) or _RELEASE.search(text)) or (
        "File Assignments" in text and "MDIN" in text
    )
