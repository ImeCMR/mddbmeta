"""NAMD configuration files (Tcl-flavoured ``keyword value`` lines).

We read what a config states without running Tcl: ``set`` variables are
substituted (``$var`` / ``${var}``), ``source`` includes are followed, and
top-level ``run`` / ``minimize`` commands are counted. Anything inside Tcl
blocks (``for``, ``if``, ``proc``...) or computed with ``[...]`` is recorded
as unresolved -- the log records the values NAMD actually used.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

KEYWORDS = {
    "structure", "coordinates", "bincoordinates", "binvelocities", "extendedsystem", "parameters",
    "paratypecharmm", "outputname", "restartname", "dcdfile", "timestep", "firsttimestep",
    "numsteps", "run", "minimize", "langevin", "langevintemp", "langevinpiston", "temperature",
    "cellbasisvector1", "pme", "cutoff", "rigidbonds", "dcdfreq", "restartfreq", "amber", "parmfile",
    "ambercoor", "outputenergies", "switching", "exclude", "wrapall",
}  # fmt: skip
TCL_BLOCKS = {"for", "foreach", "while", "if", "proc", "else", "elseif"}
_VAR = re.compile(r"\$\{(\w+)\}|\$(\w+)")


class NamdConfigError(ValueError):
    pass


@dataclass
class NamdConfig:
    path: str
    values: dict[str, str] = field(default_factory=dict)  # lower-case keyword -> last stated value
    parameters: list[str] = field(default_factory=list)
    runs: list[tuple[str, int]] = field(default_factory=list)  # ("run" | "minimize", steps)
    variables: dict[str, str] = field(default_factory=dict)
    unresolved: set[str] = field(default_factory=set)  # keywords whose value needs Tcl
    has_tcl_blocks: bool = False
    sources: list[str] = field(default_factory=list)

    def get(self, key: str) -> str | None:
        return self.values.get(key.lower())


def _strip(value: str) -> str:
    v = value.strip()
    if len(v) >= 2 and v[0] == v[-1] == '"' or (v.startswith("{") and v.endswith("}")):
        return v[1:-1].strip()
    return v


def _substitute(text: str, variables: dict[str, str]) -> tuple[str, bool]:
    resolved = True

    def repl(m: re.Match) -> str:
        nonlocal resolved
        name = m.group(1) or m.group(2)
        if name in variables:
            return variables[name]
        resolved = False
        return m.group(0)

    out = _VAR.sub(repl, text)
    return out, resolved and "[" not in out


def _logical_lines(text: str) -> list[str]:
    lines, buf = [], ""
    for raw in text.splitlines():
        line = raw.rstrip()
        if line.endswith("\\"):
            buf += line[:-1] + " "
            continue
        lines.append(buf + line)
        buf = ""
    if buf:
        lines.append(buf)
    return lines


def _uncomment(line: str) -> str:
    s = line.strip()
    if s.startswith("#"):
        return ""
    # ';#' starts an inline comment in Tcl; a bare '#' after whitespace is commonly used too
    for marker in (";#", " #", "\t#"):
        i = s.find(marker)
        if i >= 0:
            s = s[:i]
    return s.strip()


def read_config(path: str | Path, _depth: int = 0, _cfg: NamdConfig | None = None) -> NamdConfig:
    path = Path(path)
    cfg = _cfg or NamdConfig(path=str(path))
    with open_text(path) as fh:
        text = fh.read()
    depth = 0
    for line in _logical_lines(text):
        s = _uncomment(line)
        if not s:
            continue
        opening, closing = s.count("{"), s.count("}")
        if depth > 0:
            depth += opening - closing
            continue
        parts = s.split(None, 1)
        key = parts[0].lower()
        rest = parts[1] if len(parts) > 1 else ""
        if key in TCL_BLOCKS:
            cfg.has_tcl_blocks = True
            depth += opening - closing
            continue
        if key == "set":
            name, _, value = rest.partition(" ")
            value, ok = _substitute(_strip(value), cfg.variables)
            if ok:
                cfg.variables[name.strip()] = value
            continue
        if key == "source":
            target, ok = _substitute(_strip(rest), cfg.variables)
            inc = (path.parent / target) if ok else None
            if inc is not None and inc.is_file() and _depth < 5:
                cfg.sources.append(str(inc))
                read_config(inc, _depth + 1, cfg)
            continue
        value, ok = _substitute(_strip(rest), cfg.variables)
        if key in ("run", "minimize"):
            try:
                cfg.runs.append((key, int(float(value))))
            except ValueError:
                cfg.unresolved.add(key)
            continue
        if not ok:
            cfg.unresolved.add(key)
            continue
        if key == "parameters":
            cfg.parameters.append(value)
        cfg.values[key] = value
        depth += max(0, opening - closing)
    if _depth == 0 and not (cfg.values or cfg.runs):
        raise NamdConfigError("no NAMD keywords found")
    return cfg


_KEYLINE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9]*)\s+\S", re.MULTILINE)


def looks_like_config(head: bytes, name: str) -> bool:
    text = head.decode("latin-1", errors="replace")
    if "&cntrl" in text.lower() or text.startswith("%VERSION") or "Info: NAMD" in text:
        return False
    keys = {k.lower() for k in _KEYLINE.findall(text)}
    hits = len(keys & KEYWORDS)
    return hits >= 4 or (name.lower().endswith((".namd", ".conf")) and hits >= 2)
