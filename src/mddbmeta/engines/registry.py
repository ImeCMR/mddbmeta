"""Adapter registry: built-ins plus the ``mddbmeta.engines`` entry-point group."""

from __future__ import annotations

from importlib.metadata import entry_points
from pathlib import Path

from mddbmeta.engines.base import EngineAdapter, SniffResult


def _builtin() -> dict[str, type[EngineAdapter]]:
    from mddbmeta.engines.amber import AmberAdapter
    from mddbmeta.engines.gromacs import GromacsAdapter
    from mddbmeta.engines.meld import MeldAdapter
    from mddbmeta.engines.namd import NamdAdapter
    from mddbmeta.engines.openmm import OpenmmAdapter

    return {
        a.name: a for a in (AmberAdapter, GromacsAdapter, NamdAdapter, OpenmmAdapter, MeldAdapter)
    }


def available_adapters() -> dict[str, type[EngineAdapter]]:
    adapters = _builtin()
    try:
        eps = entry_points(group="mddbmeta.engines")
    except TypeError:  # pragma: no cover - very old importlib.metadata
        eps = entry_points().get("mddbmeta.engines", [])  # type: ignore[attr-defined]
    for ep in eps:
        if ep.name in adapters:
            continue
        try:
            cls = ep.load()
        except Exception:  # a broken third-party plugin must not break discovery
            continue
        adapters[ep.name] = cls
    return adapters


def get_adapter(name: str) -> EngineAdapter:
    adapters = available_adapters()
    if name not in adapters:
        raise KeyError(f"unknown engine '{name}' (available: {', '.join(sorted(adapters))})")
    return adapters[name]()


def sniff_all(
    path: Path, head: bytes, adapters: list[EngineAdapter]
) -> tuple[EngineAdapter, SniffResult] | None:
    """Ask every adapter; the highest-scoring claim wins."""
    best: tuple[EngineAdapter, SniffResult] | None = None
    for adapter in adapters:
        result = adapter.sniff(path, head)
        if result is None:
            continue
        if best is None or result.score > best[1].score:
            best = (adapter, result)
    return best
