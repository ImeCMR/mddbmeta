"""OpenMM outputs: StateDataReporter logs and XmlSerializer files.

StateDataReporter writes a header of quoted column names (``#"Step","Time
(ps)","Potential Energy (kJ/mole)",...``) and one row per report. Rows start
one report interval after the run began. The log may be a file or stdout
captured together with other output.

XmlSerializer writes ``<State>`` (openmmVersion, stepCount, time in ps, box
vectors in nm, parameters), ``<System>`` (box, particles, forces with their
parameters) and ``<Integrator>`` (type, stepSize in ps, temperature, friction).
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

NM_TO_A = 10.0
_HEADER = re.compile(r'^#"Step"|^#"Time|^#"Progress')


class OpenmmFileError(ValueError):
    pass


@dataclass
class StateLog:
    path: str
    columns: list[str]
    n_rows: int = 0
    first: dict[str, float] = field(default_factory=dict)
    last: dict[str, float] = field(default_factory=dict)
    temp_sum: float = 0.0
    volume_min: float | None = None
    volume_max: float | None = None
    report_interval: int | None = None  # steps between rows (from the first two rows)

    @property
    def mean_temperature(self) -> float | None:
        return (
            self.temp_sum / self.n_rows
            if self.n_rows and "Temperature (K)" in self.columns
            else None
        )


def looks_like_state_log(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    return any(_HEADER.match(line) for line in text.splitlines()[:200])


def read_state_log(path: str | Path, sep: str | None = None) -> StateLog:
    log: StateLog | None = None
    second: dict[str, float] | None = None
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if log is None:
                if _HEADER.match(s):
                    header = s[1:]
                    sep = sep or ("," if "," in header else "\t")
                    log = StateLog(
                        path=str(path), columns=[c.strip().strip('"') for c in header.split(sep)]
                    )
                continue
            if not s or s.startswith("#"):
                continue
            parts = s.split(sep)
            if len(parts) != len(log.columns):
                continue  # other program output mixed into stdout
            try:
                row = {
                    c: float(v)
                    for c, v in zip(log.columns, parts, strict=True)
                    if v not in ("--", "")
                }
            except ValueError:
                continue
            if not log.first:
                log.first = row
            elif second is None:
                second = row
            log.last = row
            log.n_rows += 1
            if "Temperature (K)" in row:
                log.temp_sum += row["Temperature (K)"]
            vol = row.get("Box Volume (nm^3)")
            if vol is not None:
                log.volume_min = vol if log.volume_min is None else min(log.volume_min, vol)
                log.volume_max = vol if log.volume_max is None else max(log.volume_max, vol)
    if log is None:
        raise OpenmmFileError("no StateDataReporter header")
    if second is not None and "Step" in log.first and "Step" in second:
        log.report_interval = int(second["Step"] - log.first["Step"])
    return log


@dataclass
class XmlState:
    path: str
    kind: str  # "State" | "System" | "Integrator"
    version: str | None
    step: int | None = None
    time_ps: float | None = None
    box_vectors_A: tuple[tuple[float, float, float], ...] | None = None
    natom: int | None = None
    forces: list[dict[str, str]] = field(default_factory=list)  # System: attributes of each <Force>
    integrator: dict[str, str] = field(default_factory=dict)  # Integrator: its attributes
    parameters: dict[str, str] = field(default_factory=dict)  # State: <Parameters>


def looks_like_xml_state(head: bytes) -> bool:
    text = head[:400].decode("latin-1", errors="replace")
    return bool(re.search(r"<(State|System|Integrator) [^>]*\b(openmmVersion|type)=", text))


def _vectors(el: ET.Element | None):
    if el is None:
        return None
    try:
        return tuple(
            tuple(float(el.find(n).get(c)) * NM_TO_A for c in "xyz")  # type: ignore[union-attr]
            for n in "ABC"
        )
    except (AttributeError, TypeError, ValueError):
        return None


def read_xml_state(path: str | Path) -> XmlState:
    """Parse the header parts only; positions/velocities are counted, not kept."""
    kind = None
    out: XmlState | None = None
    natom = 0
    for event, el in ET.iterparse(str(path), events=("start", "end")):
        if event == "start" and kind is None:
            kind = el.tag
            if kind not in ("State", "System", "Integrator"):
                raise OpenmmFileError(f"not an OpenMM XML file (<{kind}>)")
            out = XmlState(path=str(path), kind=kind, version=el.get("openmmVersion"))
            if kind == "State":
                out.step = int(el.get("stepCount")) if el.get("stepCount") else None
                out.time_ps = float(el.get("time")) if el.get("time") else None
            if kind == "Integrator":
                out.integrator = dict(el.attrib)
            continue
        if event != "end" or out is None:
            continue
        tag = el.tag
        if tag == "PeriodicBoxVectors" and out.box_vectors_A is None:
            out.box_vectors_A = _vectors(el)
        elif tag == "Parameters" and kind == "State":
            out.parameters = dict(el.attrib)
        elif tag == "Force" and kind == "System":
            out.forces.append({k: v for k, v in el.attrib.items()})
            el.clear()
        elif tag in ("Position", "Particle"):
            if tag == ("Position" if kind == "State" else "Particle"):
                natom += 1
            el.clear()
        elif tag in ("Velocity", "Bond", "Angle", "Torsion", "Exception", "Constraint"):
            el.clear()
    if out is None:
        raise OpenmmFileError("empty XML")
    out.natom = natom or None
    return out
