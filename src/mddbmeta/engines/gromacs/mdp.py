"""GROMACS run parameters (.mdp): ``key = value ; comment`` lines.

Keys are normalized the way GROMACS treats them: case-insensitive and with
``_`` equivalent to ``-`` (``ref_t`` == ``ref-t``). Only stated keys are
returned; defaults live in :mod:`mddbmeta.engines.gromacs.semantics`.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from mddbmeta.io.safe import open_text

# pre-2016 names -> current names
ALIASES = {
    "nstxtcout": "nstxout-compressed",
    "xtc-precision": "compressed-x-precision",
    "xtc-grps": "compressed-x-grps",
}
KNOWN_KEYS = {
    "integrator", "dt", "nsteps", "tinit", "init-step", "tcoupl", "pcoupl", "ref-t", "tc-grps",
    "nstxout", "nstxout-compressed", "nstenergy", "nstlog", "constraints", "coulombtype",
    "rcoulomb", "rvdw", "cutoff-scheme", "gen-vel", "continuation", "pbc", "emtol", "define",
}  # fmt: skip


class MdpFormatError(ValueError):
    pass


def norm_key(key: str) -> str:
    k = key.strip().lower().replace("_", "-")
    return ALIASES.get(k, k)


def parse_value(text: str) -> Any:
    tokens = text.split()
    if not tokens:
        return ""
    values = [_scalar(t) for t in tokens]
    return values[0] if len(values) == 1 else values


def _scalar(token: str) -> Any:
    try:
        return int(token)
    except ValueError:
        pass
    try:
        return float(token)
    except ValueError:
        return token


def parse_mdp_text(text: str) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for line in text.splitlines():
        body = line.split(";", 1)[0].strip()
        if not body or "=" not in body:
            continue
        key, _, value = body.partition("=")
        key = norm_key(key)
        if not key:
            continue
        # define keeps its raw text (-DPOSRES -DFLEXIBLE)
        values[key] = value.strip() if key in ("define", "include") else parse_value(value)
    if not values:
        raise MdpFormatError("no 'key = value' lines")
    return values


def read_mdp(path: str | Path) -> dict[str, Any]:
    with open_text(path) as fh:
        return parse_mdp_text(fh.read())


_LINE = re.compile(r"^\s*([A-Za-z][\w\-]*)\s*=", re.MULTILINE)


def looks_like_mdp(head: bytes, name: str) -> bool:
    text = head.decode("latin-1", errors="replace")
    if "&cntrl" in text.lower() or text.startswith("%VERSION"):
        return False
    keys = {norm_key(k) for k in _LINE.findall(text)}
    hits = len(keys & KNOWN_KEYS)
    return hits >= 3 or (name.lower().endswith(".mdp") and hits >= 1)
