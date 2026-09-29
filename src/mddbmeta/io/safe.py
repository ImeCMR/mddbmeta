"""File access that never raises: problems come back as FileLoadError."""

from __future__ import annotations

import os
from pathlib import Path

from mddbmeta.provenance import FileLoadError, classify_exception

HEAD_BYTES = 4096


def read_head(path: str | Path, n: int = HEAD_BYTES) -> bytes | FileLoadError:
    try:
        with open(path, "rb") as fh:
            return fh.read(n)
    except OSError as exc:
        return FileLoadError("unknown", str(path), classify_exception(exc), str(exc))


def open_text(path: str | Path):
    """Open a text file tolerant of stray non-UTF-8 bytes (engines write latin-1 at times)."""
    return open(path, encoding="utf-8", errors="replace", newline=None)


def read_text(path: str | Path, kind: str = "unknown") -> str | FileLoadError:
    try:
        with open(path, "rb") as fh:
            raw = fh.read()
    except OSError as exc:
        return FileLoadError(kind, str(path), classify_exception(exc), str(exc))
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def is_readable_file(path: str | Path) -> bool:
    return os.path.isfile(path) and os.access(path, os.R_OK)
