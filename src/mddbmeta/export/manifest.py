"""project.yaml: a snapshot of discovery plus the user's overrides.

The snapshot is regenerated from the run files on every load (so it can
never go stale); only the ``discovery`` options and the ``overrides`` block
are read back. Overrides:

* an MDDB field:              ``wat: OPC``, ``ff: [Amber ff14SB]``
* a per-MD field:             ``mds.rep1.temp: 310``
* a step role:                ``steps.rep1/prod_0003.role: equilibration``
* a replica's merge order:    ``replicas.rep1.production_chain: [rep1/prod_1, rep1/prod_2]``
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from mddbmeta import __version__
from mddbmeta.errors import MddbmetaError
from mddbmeta.io import yamlio
from mddbmeta.model import Project

SCHEMA_VERSION = 1

HEADER = (
    "mddbmeta project manifest.\n"
    "Only `discovery` and `overrides` are read back; `snapshot` is regenerated from the run\n"
    "files on every load. Examples of overrides:\n"
    "  wat: OPC\n"
    "  mds.rep1.temp: 310\n"
    "  steps.rep1/prod_0003.role: equilibration\n"
    "  replicas.rep1.production_chain: [rep1/prod_0001, rep1/prod_0002]"
)


def write_manifest(
    project: Project,
    path: str | Path,
    discovery: dict[str, Any] | None = None,
    overrides: dict[str, Any] | None = None,
) -> None:
    path = Path(path).resolve()
    root = Path(project.root)
    try:
        root_ref = os.path.relpath(root, path.parent)
    except ValueError:
        root_ref = str(root)
    payload = {
        "mddbmeta": {"schema_version": SCHEMA_VERSION, "version": __version__},
        "discovery": {"root": root_ref, "engine": project.engine, **(discovery or {})},
        "overrides": dict(overrides or project.overrides or {}),
        "snapshot": project.to_dict(),
    }
    yamlio.dump(payload, path, header=HEADER)


def read_manifest(path: str | Path) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    """(absolute root, discovery options, overrides)."""
    path = Path(path).resolve()
    data = yamlio.load(path)
    if not isinstance(data, dict) or "mddbmeta" not in data or "discovery" not in data:
        raise MddbmetaError(
            f"{path} is not an mddbmeta manifest (write one with `mddbmeta discover DIR --write`)"
        )
    disc = dict(data.get("discovery") or {})
    root_ref = disc.pop("root", ".")
    root = (path.parent / root_ref).resolve() if not os.path.isabs(root_ref) else Path(root_ref)
    overrides = data.get("overrides") or {}
    if not isinstance(overrides, dict):
        raise MddbmetaError(f"{path}: overrides must be a mapping")
    return root, disc, overrides


def split_overrides(
    overrides: dict[str, Any],
) -> tuple[dict[str, str], dict[str, list[str]], dict[str, Any]]:
    """(step roles, replica chains, field overrides)."""
    roles: dict[str, str] = {}
    chains: dict[str, list[str]] = {}
    fields: dict[str, Any] = {}
    for key, value in overrides.items():
        if key.startswith("steps.") and key.endswith(".role"):
            roles[key[len("steps.") : -len(".role")]] = str(value)
        elif key.startswith("replicas.") and key.endswith(".production_chain"):
            chains[key[len("replicas.") : -len(".production_chain")]] = [
                str(v) for v in value or []
            ]
        else:
            fields[key] = value
    return roles, chains, fields
