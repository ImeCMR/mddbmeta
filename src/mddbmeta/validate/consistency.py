"""Do the production steps agree with each other and with the topology?"""

from __future__ import annotations

from typing import Any

from mddbmeta.findings import finding
from mddbmeta.model import FileRole, Project, Step
from mddbmeta.provenance import Finding

PARAMS = {
    "dt_ps": "timestep",
    "traj_interval_steps": "trajectory write interval",
    "temp_K": "target temperature",
    "ensemble": "ensemble",
}


def distinct(values: list[Any]) -> list[Any]:
    out: list[Any] = []
    for v in values:
        if v is not None and v not in out:
            out.append(v)
    return out


def check_consistency(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    for rep in project.replicas:
        prod = [s for s in project.production_steps(rep) if s.status == "done"]
        for attr, label in PARAMS.items():
            values = {s.id: getattr(s.settings, attr).value for s in prod}
            if len(distinct(list(values.values()))) > 1:
                findings.append(
                    finding(
                        "CONS001",
                        f"{rep.id}: production runs use different {label}: "
                        + ", ".join(f"{k}={v}" for k, v in values.items()),
                        subject=rep.id,
                        parameter=attr,
                        values=values,
                    )
                )
    all_prod = [
        s for rep in project.replicas for s in project.production_steps(rep) if s.status == "done"
    ]
    for attr in ("program", "version"):
        values = {s.id: getattr(s.runtime, attr).value for s in all_prod}
        if len(distinct(list(values.values()))) > 1:
            findings.append(
                finding(
                    "CONS003",
                    f"production runs report different {attr}s: "
                    + ", ".join(sorted(set(map(str, values.values())))),
                    parameter=attr,
                    values=values,
                )
            )
    findings.extend(_atom_counts(project))
    return findings


def _atom_counts(project: Project) -> list[Finding]:
    findings = []
    adapter = project.adapter
    top_natom: dict[str, int | None] = {}
    for step in project.steps.values():
        counts: dict[str, int] = {}
        top = step.files.get(FileRole.TOPOLOGY)
        if top is not None and adapter is not None and top.path in project.parsed:
            if top.path not in top_natom:
                info = adapter.system_info(project.parsed[top.path])
                top_natom[top.path] = info.natom if info else None
            if top_natom[top.path]:
                counts[top.path] = int(top_natom[top.path])  # type: ignore[arg-type]
        if step.runtime.natom.is_known:
            counts[step.files[FileRole.LOG].path if FileRole.LOG in step.files else "log"] = int(
                step.runtime.natom.value
            )
        if step.trajectory.natom.is_known and step.trajectory.natom.method == "stated":
            counts[step.files[FileRole.TRAJECTORY].path] = int(step.trajectory.natom.value)
        coords = step.files.get(FileRole.INPUT_COORDS)
        if coords is not None and adapter is not None and coords.path in project.parsed:
            n, _t = adapter.coords_natom_time(project.parsed[coords.path])
            if n:
                counts[coords.path] = int(n)
        if len(set(counts.values())) > 1:
            findings.append(
                finding(
                    "CONS004",
                    f"{step.id}: atom counts disagree: "
                    + ", ".join(f"{p}={n}" for p, n in counts.items()),
                    subject=step.id,
                    paths=list(counts),
                    counts=counts,
                )
            )
    return findings


def production_values(steps: list[Step], getter) -> list[Any]:
    return distinct([getter(s) for s in steps])
