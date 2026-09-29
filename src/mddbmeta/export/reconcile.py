"""Fill an existing inputs file without overwriting user values, and compare the two.

``fill`` only writes fields that are missing, null or empty; it never
reorders a user's trajectory list. ``reconcile`` is read-only: it expands the
user's ``input_trajectory_filepaths`` the way MDDB-workflow does (globs, then
the MD directory, the project directory) and checks the resulting merge
order against the continuity chain.
"""

from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta.export.mddb import MINED_KEYS, TEMPLATE_ORDER, Export
from mddbmeta.findings import finding
from mddbmeta.model import Project
from mddbmeta.provenance import Finding, Severity

GLOB_CHARS = re.compile(r"[*?\[]")


def _empty(value: Any) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _match_md(existing: list[dict[str, Any]], entry: dict[str, Any]) -> dict[str, Any] | None:
    for key in ("name", "mdir"):
        for md in existing:
            if (
                isinstance(md, dict)
                and entry.get(key) is not None
                and md.get(key) == entry.get(key)
            ):
                return md
    return None


@dataclass
class FillResult:
    inputs: dict[str, Any]
    filled: list[str] = field(default_factory=list)
    kept: list[str] = field(default_factory=list)
    added_mds: list[str] = field(default_factory=list)
    skipped_mds: list[str] = field(default_factory=list)


def fill(existing: dict[str, Any], export: Export, add_mds: bool = False) -> FillResult:
    merged = dict(existing or {})
    res = FillResult(inputs=merged)
    for key, value in export.inputs.items():
        if key in ("mds", "mdref"):
            continue
        if _empty(merged.get(key)):
            merged[key] = value
            res.filled.append(key)
        elif merged.get(key) != value:
            res.kept.append(key)
    user_mds = merged.get("mds")
    had_mds = isinstance(user_mds, list) and len(user_mds) > 0
    mds: list[dict[str, Any]] = [dict(m) if isinstance(m, dict) else m for m in (user_mds or [])]
    for entry in export.inputs.get("mds", []):
        target = _match_md(mds, entry) if had_mds else None
        if target is None:
            if had_mds and not add_mds:
                res.skipped_mds.append(entry.get("name", "?"))
                continue
            mds.append(dict(entry))
            res.added_mds.append(entry.get("name", "?"))
            continue
        for key, value in entry.items():
            if _empty(target.get(key)):
                target[key] = value
                res.filled.append(f"mds.{entry.get('name')}.{key}")
            elif target.get(key) != value:
                res.kept.append(f"mds.{entry.get('name')}.{key}")
    merged["mds"] = mds
    if "mdref" not in merged and mds:
        merged["mdref"] = 0
    rank = {k: i for i, k in enumerate(TEMPLATE_ORDER)}
    original = list((existing or {}).keys())
    new_keys = sorted(
        (k for k in merged if k not in original), key=lambda k: rank.get(k, len(rank))
    )
    res.inputs = {k: merged[k] for k in original + new_keys}
    return res


# ---------------------------------------------------------------- reconcile


@dataclass
class Item:
    field: str
    status: str  # agree | fill | conflict | user_only | unmined
    user: Any = None
    mined: Any = None
    detail: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if v is not None}


def _norm(value: Any) -> Any:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        return re.sub(r"[\s_\-/]+", "", value).lower()
    if isinstance(value, (list, tuple)):
        return sorted(str(_norm(v)) for v in value)
    return value


def same(user: Any, mined: Any) -> bool:
    if isinstance(user, str) and isinstance(mined, list) and len(mined) == 1:
        mined = mined[0]
    if isinstance(mined, str) and isinstance(user, list) and len(user) == 1:
        user = user[0]
    u, m = _norm(user), _norm(mined)
    if isinstance(u, float) and isinstance(m, float):
        return abs(u - m) <= 1e-6 * max(abs(u), abs(m), 1e-12)
    if isinstance(u, float) and isinstance(m, str) or isinstance(u, str) and isinstance(m, float):
        return str(user).strip() == str(mined).strip()
    return u == m


def expand_like_mwf(paths: Any, md_dir: Path, project_dir: Path) -> tuple[list[str], bool]:
    """(absolute paths in merge order, whether any glob was involved)."""
    entries = [paths] if isinstance(paths, str) else list(paths or [])
    used_glob = any(GLOB_CHARS.search(str(e)) for e in entries)
    for context in (md_dir, project_dir, Path.cwd()):
        out: list[str] = []
        for e in entries:
            pattern = str(e) if os.path.isabs(str(e)) else str(context / str(e))
            if GLOB_CHARS.search(pattern):
                out.extend(glob.glob(pattern))  # unsorted, exactly like mwf
            elif os.path.exists(pattern):
                out.append(pattern)
        if out:
            return [os.path.realpath(p) for p in out], used_glob
    return [], used_glob


def reconcile(
    existing: dict[str, Any], project: Project, export: Export
) -> tuple[list[Item], list[Finding]]:
    items: list[Item] = []
    findings: list[Finding] = []
    root = Path(project.root)
    for key in MINED_KEYS:
        m = project.mined.get(key)
        mined_value = export.inputs.get(key, m.value if m is not None else None)
        user = existing.get(key)
        if _empty(user) and mined_value is None:
            items.append(Item(key, "unmined"))
        elif _empty(user):
            items.append(Item(key, "fill", mined=mined_value))
        elif mined_value is None:
            items.append(Item(key, "user_only", user=user))
        elif same(user, mined_value):
            items.append(Item(key, "agree", user=user, mined=mined_value))
        else:
            detail = None
            if m is not None and m.confidence.value == "heuristic":
                detail = "mined value is only a heuristic guess"
            items.append(Item(key, "conflict", user=user, mined=mined_value, detail=detail))
            findings.append(
                finding(
                    "REC001",
                    f"{key}: inputs file says {user!r}, the run files say {mined_value!r}",
                    data={"field": key},
                    severity=Severity.SUGGESTION if detail else None,
                )
            )

    user_mds = existing.get("mds") or []
    for entry in export.inputs.get("mds", []):
        name = entry.get("name")
        md = _match_md(user_mds, entry) if isinstance(user_mds, list) else None
        label = f"mds.{name}.input_trajectory_filepaths"
        chain = [os.path.realpath(root / p) for p in entry.get("input_trajectory_filepaths", [])]
        if md is None:
            items.append(
                Item(f"mds.{name}", "fill", mined=entry, detail="MD not in the inputs file")
            )
            continue
        for key in ("temp", "framestep", "timestep", "ensemble"):
            if key in entry:
                user = md.get(key, existing.get(key))
                if _empty(user):
                    items.append(Item(f"mds.{name}.{key}", "fill", mined=entry[key]))
                elif not same(user, entry[key]):
                    items.append(Item(f"mds.{name}.{key}", "conflict", user=user, mined=entry[key]))
                    findings.append(
                        finding("REC001", f"mds.{name}.{key}: {user!r} vs mined {entry[key]!r}")
                    )
        user_paths = md.get(
            "input_trajectory_filepaths", existing.get("input_trajectory_filepaths")
        )
        if _empty(user_paths):
            items.append(Item(label, "fill", mined=entry.get("input_trajectory_filepaths")))
            continue
        mdir = md.get("mdir") or ""
        expanded, used_glob = expand_like_mwf(user_paths, root / str(mdir), root)
        rel = [os.path.relpath(p, root) for p in expanded]
        mined_rel = entry.get("input_trajectory_filepaths")
        if expanded == chain:
            items.append(Item(label, "agree", user=user_paths, mined=mined_rel))
            if used_glob and len(chain) > 1:
                findings.append(
                    finding(
                        "ORDER002",
                        f"{label}: the glob happens to expand in the right order here, but glob order "
                        "is filesystem-dependent; list the parts explicitly",
                        subject=name,
                        severity=Severity.WARNING,
                    )
                )
        elif sorted(expanded) == sorted(chain):
            items.append(
                Item(
                    label,
                    "conflict",
                    user=rel,
                    mined=mined_rel,
                    detail="same parts, different order",
                )
            )
            findings.append(
                finding(
                    "ORDER002",
                    f"{label}: MDDB-workflow would merge the parts as {rel}, but the continuity chain is "
                    f"{mined_rel}",
                    subject=name,
                    user=rel,
                    chain=mined_rel,
                )
            )
        else:
            missing = [os.path.relpath(p, root) for p in chain if p not in expanded]
            extra = [os.path.relpath(p, root) for p in expanded if p not in chain]
            items.append(Item(label, "conflict", user=rel, mined=mined_rel,
                              detail=f"missing {missing}, extra {extra}"))  # fmt: skip
            findings.append(
                finding(
                    "REC001",
                    f"{label}: parts differ from the production chain (missing {missing}, extra {extra})",
                    subject=name,
                )
            )
    return items, findings


def has_conflicts(items: list[Item], findings: list[Finding]) -> bool:
    return any(f.code == "ORDER002" and f.severity is Severity.ERROR for f in findings) or any(
        i.status == "conflict" and i.detail != "mined value is only a heuristic guess"
        for i in items
    )
