"""Hand-off continuity between chained steps, and holes in numbered run families."""

from __future__ import annotations

from itertools import pairwise

from mddbmeta.findings import finding
from mddbmeta.model import CoordSourceKind, FileRole, PhaseRole, Project, Step
from mddbmeta.provenance import Finding, Severity
from mddbmeta.units import continuity_tolerance


def restart_time(project: Project, step: Step) -> float | None:
    """Simulation time stored in the restart a step wrote (exact), else its log end time."""
    ref = step.files.get(FileRole.OUTPUT_RESTART)
    if ref is not None and project.adapter is not None and ref.path in project.parsed:
        _natom, t = project.adapter.coords_natom_time(project.parsed[ref.path])
        if t is not None:
            return float(t)
    if step.runtime.end_time_ps.is_known:
        return float(step.runtime.end_time_ps.value)
    return None


def input_time(project: Project, step: Step) -> float | None:
    if step.runtime.start_time_ps.is_known:
        return float(step.runtime.start_time_ps.value)
    ref = step.files.get(FileRole.INPUT_COORDS)
    if ref is not None and project.adapter is not None and ref.path in project.parsed:
        _natom, t = project.adapter.coords_natom_time(project.parsed[ref.path])
        return float(t) if t is not None else None
    return None


def frame_interval_ps(step: Step) -> float | None:
    s = step.settings
    if s.traj_interval_steps.is_known and s.dt_ps.is_known and s.traj_interval_steps.value:
        return float(s.traj_interval_steps.value) * float(s.dt_ps.value)
    if step.trajectory.frame_interval_ps.is_known:
        return float(step.trajectory.frame_interval_ps.value)
    return None


def _is_md(step: Step) -> bool:
    return step.settings.minimization.value is False or (
        step.settings.minimization.value is None and step.status == "orphan"
    )


def check_continuity(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    steps = project.steps
    for step in steps.values():
        if step.status == "queued":
            continue
        src = step.input_coords
        if src.kind is CoordSourceKind.UNKNOWN and src.path:
            findings.append(
                finding(
                    "CONT002",
                    f"{step.id} reads {src.path}, which is not in the project",
                    subject=step.id,
                    paths=[src.path],
                )
            )
            continue
        if src.kind is not CoordSourceKind.STEP or src.ref not in steps:
            if step.settings.restart.value is True and src.kind is CoordSourceKind.PATH:
                findings.append(
                    finding(
                        "CONT002",
                        f"{step.id} continues from {src.path}, which no step in the project wrote",
                        subject=step.id,
                        paths=[src.path] if src.path else [],
                        severity=None,
                    )
                )
            continue
        producer = steps[src.ref]
        if not _is_md(step) or not _is_md(producer):
            continue  # minimization has no clock
        both_prod = (
            step.phase_role is PhaseRole.PRODUCTION and producer.phase_role is PhaseRole.PRODUCTION
        )
        if step.settings.restart.value is False:
            if both_prod:
                findings.append(
                    finding(
                        "CONT002",
                        f"{step.id} restarts the clock and velocities (irest=0) mid production chain "
                        f"after {producer.id}",
                        subject=step.id,
                    )
                )
            continue
        prev_end, cur_start = restart_time(project, producer), input_time(project, step)
        if prev_end is None or cur_start is None:
            findings.append(
                finding(
                    "CONT004",
                    f"continuity {producer.id} -> {step.id} not measured (times unavailable)",
                    subject=step.id,
                )
            )
            continue
        gap = round(cur_start - prev_end, 6) + 0.0
        tol = continuity_tolerance(frame_interval_ps(producer))
        data = {"gap_ps": gap, "tolerance_ps": tol, "from": producer.id, "to": step.id}
        prev_start = input_time(project, producer)
        if abs(gap) > tol and prev_start is not None and cur_start <= prev_start + tol:
            # harmless between preparation stages; within production it breaks the merged time axis
            findings.append(
                finding(
                    "CONT005",
                    f"{step.id} continues {producer.id}'s state but its clock restarts at {cur_start:g} ps "
                    f"({producer.id} ended at {prev_end:g} ps)",
                    subject=step.id,
                    data=data,
                    severity=None if both_prod else Severity.INFO,
                )
            )
            continue
        if abs(gap) <= tol:
            findings.append(
                finding(
                    "CONT000",
                    f"{producer.id} -> {step.id}: gap {gap:g} ps",
                    subject=step.id,
                    data=data,
                )
            )
        elif gap > 0:
            findings.append(
                finding(
                    "CONT001",
                    f"{producer.id} -> {step.id}: {gap:g} ps missing between the runs (tolerance {tol:g} ps)",
                    subject=step.id,
                    data=data,
                )
            )
        else:
            findings.append(
                finding(
                    "CONT003",
                    f"{producer.id} -> {step.id}: runs overlap by {-gap:g} ps (tolerance {tol:g} ps)",
                    subject=step.id,
                    data=data,
                )
            )
    findings.extend(_trajectory_joins(project))
    return findings


def _trajectory_joins(project: Project) -> list[Finding]:
    """Frames at each join of the production chain: what MDDB-workflow will concatenate."""
    findings: list[Finding] = []
    for rep in project.replicas:
        chain = project.production_steps(rep)
        for prev, nxt in pairwise(chain):
            pt, nt = prev.trajectory, nxt.trajectory
            if not (pt.last_time_ps.is_known and nt.first_time_ps.is_known):
                continue
            interval = frame_interval_ps(nxt) or frame_interval_ps(prev)
            if interval is None:
                continue
            delta = float(nt.first_time_ps.value) - float(pt.last_time_ps.value)
            # GROMACS re-writes the boundary frame (expected delta 0); AMBER starts one interval in
            expected = 0.0 if nt.first_frame_at_start else interval
            off = round(delta - expected, 6) + 0.0
            tol = continuity_tolerance(interval)
            data = {
                "from": prev.id,
                "to": nxt.id,
                "frame_delta_ps": round(delta, 6),
                "frame_interval_ps": interval,
            }
            linked = (
                nxt.input_coords.kind is CoordSourceKind.STEP and nxt.input_coords.ref == prev.id
            )
            reset = (
                pt.first_time_ps.is_known
                and float(nt.first_time_ps.value) <= float(pt.first_time_ps.value) + tol
            )
            if abs(off) > tol and reset:
                if linked:  # already reported on the hand-off itself
                    continue
                findings.append(
                    finding(
                        "CONT005",
                        f"trajectory join {prev.id} -> {nxt.id}: time goes back from "
                        f"{pt.last_time_ps.value:g} to {nt.first_time_ps.value:g} ps",
                        subject=nxt.id,
                        data=data,
                    )
                )
            elif abs(off) <= tol:
                if not linked:
                    findings.append(
                        finding(
                            "CONT000",
                            f"trajectory join {prev.id} -> {nxt.id} is seamless",
                            subject=nxt.id,
                            data=data,
                        )
                    )
            elif off > 0:
                findings.append(
                    finding(
                        "CONT001",
                        f"trajectory join {prev.id} -> {nxt.id}: {off:g} ps of frames missing",
                        subject=nxt.id,
                        data=data,
                    )
                )
            else:
                findings.append(
                    finding(
                        "CONT003",
                        f"trajectory join {prev.id} -> {nxt.id}: {-off:g} ps of duplicated/overlapping frames",
                        subject=nxt.id,
                        data=data,
                    )
                )
    return findings


def check_sequences(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    groups: dict[tuple[str | None, str, str], list[Step]] = {}
    for step in project.steps.values():
        if step.family is None or step.seq_index is None:
            continue
        groups.setdefault((step.replica, step.dir, step.family), []).append(step)
    for (replica, run_dir, base), members in sorted(
        groups.items(), key=lambda kv: tuple(str(x) for x in kv[0])
    ):
        if len(members) < 2:
            continue
        indices = sorted(s.seq_index for s in members)  # type: ignore[misc]
        label = f"{run_dir + '/' if run_dir else ''}{base}"
        dupes = sorted({i for i in indices if indices.count(i) > 1})
        if dupes:
            findings.append(
                finding(
                    "SEQ002",
                    f"{label} has duplicate member(s) {', '.join(map(str, dupes))}",
                    subject=replica,
                    family=base,
                    duplicates=dupes,
                )
            )
        missing = sorted(set(range(indices[0], indices[-1] + 1)) - set(indices))
        if missing:
            findings.append(
                finding(
                    "SEQ001",
                    f"{label} sequence is missing member(s) {', '.join(map(str, missing))}",
                    subject=replica,
                    family=base,
                    missing=missing,
                    dir=run_dir,
                )
            )
    return findings
