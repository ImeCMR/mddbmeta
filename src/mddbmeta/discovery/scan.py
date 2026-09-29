"""Walk a project directory, sniff every file, parse what an adapter claims."""

from __future__ import annotations

import fnmatch
import os
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.engines.base import EngineAdapter, ParsedFile
from mddbmeta.engines.registry import sniff_all
from mddbmeta.findings import finding
from mddbmeta.io.safe import read_head
from mddbmeta.provenance import FileLoadError, Finding

SKIP_DIRS = {"__pycache__", "node_modules", ".git", ".svn", ".hg"}
# MDDB-workflow writes these into a project; they are outputs, not raw inputs.
MWF_CACHE_NAMES = (".mwf_cache.json", ".cache.json")
MWF_OUTPUT_PATTERNS = ("topology.*", "structure.pdb", "trajectory.xtc", "mdf.*", "mda.*")


@dataclass
class ScanResult:
    root: Path
    files: list[ParsedFile] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)
    load_errors: list[FileLoadError] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    engines: dict[str, list[str]] = field(default_factory=dict)  # engine -> claimed files


def _is_mwf_dir(path: Path) -> bool:
    return any((path / name).exists() for name in MWF_CACHE_NAMES)


def scan(
    root: Path,
    adapters: list[EngineAdapter],
    max_depth: int = 4,
    exclude: list[str] | None = None,
) -> tuple[ScanResult, dict[str, EngineAdapter]]:
    """Return parsed files plus which adapter claimed each (by relative path)."""
    root = root.expanduser().resolve()
    result = ScanResult(root=root)
    owner: dict[str, EngineAdapter] = {}
    exclude = list(exclude or [])
    seen_real: set[str] = set()

    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        current = Path(dirpath)
        real = os.path.realpath(dirpath)
        if real in seen_real:  # symlink loop guard
            dirnames[:] = []
            continue
        seen_real.add(real)
        rel_dir = current.relative_to(root).as_posix()
        depth = 0 if rel_dir == "." else len(Path(rel_dir).parts)
        dirnames[:] = sorted(
            d
            for d in dirnames
            if not d.startswith(".")
            and d not in SKIP_DIRS
            and not _excluded(_join(rel_dir, d), exclude)
            and depth < max_depth
        )
        mwf_here = _is_mwf_dir(current)
        for name in sorted(filenames):
            if name.startswith("."):
                continue
            rel = _join(rel_dir, name)
            if _excluded(rel, exclude):
                result.skipped.append(rel)
                continue
            if mwf_here and any(fnmatch.fnmatch(name, p) for p in MWF_OUTPUT_PATTERNS):
                result.skipped.append(rel)
                continue
            path = current / name
            if not path.is_file():
                continue
            head = read_head(path)
            if isinstance(head, FileLoadError):
                err = FileLoadError("unknown", rel, head.error_type, head.detail)
                result.load_errors.append(err)
                continue
            claim = sniff_all(path, head, adapters)
            if claim is None:
                continue
            adapter, sniffed = claim
            parsed = adapter.parse(path, rel, sniffed)
            result.files.append(parsed)
            owner[rel] = adapter
            result.engines.setdefault(adapter.name, []).append(rel)
            result.findings.extend(parsed.findings)
            for err in parsed.errors:
                result.load_errors.append(err)
                result.findings.append(
                    finding(
                        "LOAD001",
                        f"{rel}: {err.error_type} ({err.detail})",
                        subject=rel,
                        paths=[rel],
                    )
                )
    return result, owner


def _join(rel_dir: str, name: str) -> str:
    return name if rel_dir in ("", ".") else f"{rel_dir}/{name}"


def _excluded(rel: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(Path(rel).name, p) for p in patterns)
