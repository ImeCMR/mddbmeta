"""AMBER control input (mdin) reader: title plus Fortran namelists.

Only values actually written in the file are returned; defaults are applied
later, explicitly, in :mod:`mddbmeta.engines.amber.semantics`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta.io.safe import open_text

_INT_RE = re.compile(r"^[+-]?\d+$")
_FLOAT_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")
_NAMELIST_START = re.compile(r"^\s*[&$]\s*([A-Za-z_]\w*)", re.MULTILINE)


class MdinFormatError(ValueError):
    pass


@dataclass
class Mdin:
    path: str | None
    title: str = ""
    namelists: list[tuple[str, dict[str, Any]]] = field(default_factory=list)
    trailer: str = ""  # free text after the namelists (DISANG=, LISTIN=, group restraints...)

    def first(self, name: str) -> dict[str, Any] | None:
        for nl_name, values in self.namelists:
            if nl_name == name:
                return values
        return None

    def all(self, name: str) -> list[dict[str, Any]]:
        return [values for nl_name, values in self.namelists if nl_name == name]

    @property
    def cntrl(self) -> dict[str, Any]:
        return self.first("cntrl") or {}

    @property
    def n_cntrl(self) -> int:
        return len(self.all("cntrl"))


def parse_value(token: str) -> Any:
    t = token.strip()
    if len(t) >= 2 and t[0] == t[-1] and t[0] in "'\"":
        return t[1:-1]
    low = t.lower()
    if low in (".true.", ".t.", "t", "true"):
        return True
    if low in (".false.", ".f.", "f", "false"):
        return False
    if _INT_RE.match(t):
        return int(t)
    if _FLOAT_RE.match(t):
        return float(t.replace("d", "e").replace("D", "E"))
    return t


def _strip_comments(text: str) -> str:
    """Remove '!' comments that are not inside quotes."""
    out_lines = []
    for line in text.splitlines():
        quote = None
        cut = len(line)
        for i, ch in enumerate(line):
            if quote:
                if ch == quote:
                    quote = None
            elif ch in "'\"":
                quote = ch
            elif ch == "!":
                cut = i
                break
        out_lines.append(line[:cut])
    return "\n".join(out_lines)


def _tokenize(body: str) -> tuple[list[str], int]:
    """Tokens of a namelist body up to its terminator; returns (tokens, consumed chars)."""
    tokens: list[str] = []
    i = 0
    n = len(body)
    while i < n:
        ch = body[i]
        if ch in " \t\r\n,":
            i += 1
            continue
        if ch in "'\"":
            j = body.find(ch, i + 1)
            if j < 0:
                raise MdinFormatError("unterminated string in namelist")
            tokens.append(body[i : j + 1])
            i = j + 1
            continue
        if ch == "/":
            return tokens, i + 1
        if ch in "&$":
            m = re.match(r"[&$]\s*end\b", body[i:], re.IGNORECASE)
            if m:
                return tokens, i + m.end()
            # a new namelist started without closing the previous one
            return tokens, i
        if ch == "=":
            tokens.append("=")
            i += 1
            continue
        j = i
        while j < n and body[j] not in " \t\r\n,=/'\"":
            j += 1
        tokens.append(body[i:j])
        i = j
    return tokens, n


def _assign(tokens: list[str]) -> dict[str, Any]:
    values: dict[str, Any] = {}
    i = 0
    while i < len(tokens):
        key = tokens[i]
        if i + 1 < len(tokens) and tokens[i + 1] == "=":
            items = []
            j = i + 2
            while j < len(tokens):
                if j + 1 < len(tokens) and tokens[j + 1] == "=":
                    break
                if tokens[j] != "=":
                    items.append(parse_value(tokens[j]))
                j += 1
            name = key.lower()
            if not items:
                values[name] = None
            else:
                values[name] = items[0] if len(items) == 1 else items
            i = j
        else:
            i += 1  # stray token; tolerated
    return values


def parse_mdin_text(text: str, path: str | None = None) -> Mdin:
    clean = _strip_comments(text)
    mdin = Mdin(path=path)
    m = _NAMELIST_START.search(clean)
    if m is None:
        raise MdinFormatError("no namelist (&cntrl ...) found")
    mdin.title = clean[: m.start()].strip()
    pos = last_end = m.start()
    while (m := _NAMELIST_START.search(clean, pos)) is not None:
        name = m.group(1).lower()
        body_start = m.end()
        if name == "end":
            pos = last_end = body_start
            continue
        tokens, consumed = _tokenize(clean[body_start:])
        mdin.namelists.append((name, _assign(tokens)))
        pos = last_end = body_start + max(consumed, 1)
    mdin.trailer = clean[last_end:].strip()
    return mdin


def read_mdin(path: str | Path) -> Mdin:
    with open_text(path) as fh:
        return parse_mdin_text(fh.read(), path=str(path))


def looks_like_mdin(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    return re.search(r"^\s*[&$]\s*cntrl\b", text, re.IGNORECASE | re.MULTILINE) is not None and (
        "Amber" not in text[:600] and "PMEMD" not in text[:600]
    )
