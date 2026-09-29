"""YAML/JSON reading and writing with PyYAML as an optional extra.

Without PyYAML, writing falls back to JSON syntax (which is valid YAML, so
mwf still reads it); reading a YAML file requires PyYAML.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mddbmeta.errors import MddbmetaError, MissingDependencyError

try:  # pragma: no cover - exercised implicitly
    import yaml as _yaml
except ImportError:  # pragma: no cover
    _yaml = None


def have_yaml() -> bool:
    return _yaml is not None


def load(path: str | Path) -> Any:
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        return json.loads(text)
    if _yaml is None:
        try:
            return json.loads(text)
        except json.JSONDecodeError as exc:
            raise MissingDependencyError(f"Reading {path.name}", "yaml") from exc
    try:
        return _yaml.safe_load(text)
    except _yaml.YAMLError as exc:
        raise MddbmetaError(f"{path}: invalid YAML: {exc}") from exc


def dumps_yaml(data: Any, header: str | None = None) -> str:
    """Serialize to YAML (block style, insertion order kept)."""
    if _yaml is None:
        body = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    else:
        body = _yaml.safe_dump(
            data, sort_keys=False, default_flow_style=False, allow_unicode=True, width=100
        )
    if header:
        lines = "".join(f"# {line}\n" if line else "#\n" for line in header.splitlines())
        return lines + body
    return body


def dump(data: Any, path: str | Path, header: str | None = None) -> None:
    path = Path(path)
    if path.suffix.lower() == ".json":
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    else:
        path.write_text(dumps_yaml(data, header=header), encoding="utf-8")


def dump_scalar(value: Any) -> str:
    """Render one value as an inline YAML/JSON-compatible literal (for Jinja templates)."""
    if value is None:
        return "null"
    return json.dumps(value, ensure_ascii=False)
