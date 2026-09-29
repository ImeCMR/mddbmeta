"""Plain-text rendering for the CLI (JSON output bypasses this)."""

from __future__ import annotations

from typing import Any

from mddbmeta.export.reconcile import Item
from mddbmeta.model import FileRole, Project
from mddbmeta.provenance import Finding, Mined, Severity

FIELD_ORDER = ["program", "version", "type", "timestep", "framestep", "temp", "ensemble", "boxtype",
               "wat", "ff", "input_topology_filepath"]  # fmt: skip


def _fmt_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, list):
        return ", ".join(map(str, value)) if value else "[]"
    return str(value)


def _src(m: Mined[Any]) -> str:
    if not m.sources:
        return ""
    first = m.sources[0]
    more = f" (+{len(m.sources) - 1})" if len(m.sources) > 1 else ""
    return f"{first.path}{':' + first.locator if first.locator else ''}{more}"


def findings_block(findings: list[Finding], verbose: bool = False) -> list[str]:
    shown = [f for f in findings if verbose or f.severity is not Severity.INFO]
    counts: dict[str, int] = {}
    for f in findings:
        counts[f.severity.value] = counts.get(f.severity.value, 0) + 1
    summary = ", ".join(
        f"{counts[s.value]} {s.value}" for s in reversed(Severity) if s.value in counts
    )
    lines = [
        f"Findings: {summary or 'none'}"
        + ("" if verbose or "info" not in counts else "  (info hidden; -v shows it)")
    ]
    order = {s: i for i, s in enumerate(reversed(list(Severity)))}
    for f in sorted(shown, key=lambda f: order[f.severity]):
        lines.append(f"  [{f.severity.value}] {f.code} {f.message}")
    return lines


def render_project(project: Project, verbose: bool = False) -> str:
    lines = [f"Project: {project.root}", f"Engine:  {project.engine or '-'}"]
    if project.topologies:
        tops = ", ".join(t.path for t in project.topologies)
        lines.append(f"Topologies: {tops}")
    if project.starting_structure:
        lines.append(f"Starting structure: {project.starting_structure.path}")
    shared = [s for s in project.steps.values() if s.replica is None]
    if shared:
        lines.append(f"Shared steps (no replica): {', '.join(s.id for s in shared)}")
    for rep in project.replicas:
        mdir = project.md_mined.get(rep.id, {}).get("mdir")
        lines.append("")
        lines.append(
            f"Replica {rep.id}  [dir: {rep.dir or '.'}]  -> mdir {mdir.value if mdir else '-'}"
        )
        for phase in project.phases(rep):
            for sid in phase.step_ids:
                s = project.steps[sid]
                src = s.input_coords
                came = src.ref or src.path or "?"
                traj = s.files.get(FileRole.TRAJECTORY)
                frames = s.trajectory.n_frames.value
                extra = (
                    f"  traj {traj.path}" + (f" ({frames} frames)" if frames is not None else "")
                    if traj
                    else ""
                )
                status = "" if s.status == "done" else f" [{s.status}]"
                done = s.runtime.completed.value
                mark = "" if done in (True, None) else " [incomplete]"
                lines.append(f"  {phase.role.value:<14} {s.id}{status}{mark}  <- {came}{extra}")
        trajs = [project.steps[s].files[FileRole.TRAJECTORY].path for s in rep.production_chain
                 if FileRole.TRAJECTORY in project.steps[s].files]  # fmt: skip
        lines.append(
            f"  merge order: {', '.join(trajs) if trajs else '(no production trajectories)'}"
        )
    if project.mined:
        lines.append("")
        lines.append("MDDB fields:")
        for key in FIELD_ORDER:
            m = project.mined.get(key)
            if m is None:
                continue
            value = project.overrides.get(key, m.value)
            conf = (
                "user" if key in project.overrides else (m.confidence.value if m.is_known else "")
            )
            lines.append(
                f"  {key:<24} {_fmt_value(value):<22} {conf:<10} {m.method}  {_src(m)}".rstrip()
            )
        for rid, fields in project.md_mined.items():
            per = {
                k: v.value
                for k, v in fields.items()
                if k in ("temp", "framestep", "timestep", "ensemble", "input_topology_filepath")
            }
            if per:
                lines.append(
                    f"  mds.{rid}: " + ", ".join(f"{k}={_fmt_value(v)}" for k, v in per.items())
                )
    lines.append("")
    lines += findings_block(project.findings, verbose)
    return "\n".join(lines)


def render_reconcile(items: list[Item], findings: list[Finding]) -> str:
    lines = [f"{'field':<44} {'status':<10} {'inputs file':<24} mined"]
    for it in items:
        lines.append(
            f"{it.field:<44} {it.status:<10} {_fmt_value(it.user)[:24]:<24} {_fmt_value(it.mined)[:40]}"
            + (f"  ({it.detail})" if it.detail else "")
        )
    lines.append("")
    lines += findings_block(findings, verbose=True)
    return "\n".join(lines)
