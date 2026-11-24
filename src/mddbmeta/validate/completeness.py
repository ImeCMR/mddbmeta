"""Did every run finish, write its trajectory, and does the trajectory match the log?"""

from __future__ import annotations

from mddbmeta.findings import finding
from mddbmeta.model import CoordSourceKind, FileRole, Project
from mddbmeta.provenance import Finding, Severity
from mddbmeta.units import continuity_tolerance
from mddbmeta.validate.continuity import frame_interval_ps


def check_completeness(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    continued: set[str] = {
        s.input_coords.ref
        for s in project.steps.values()
        if s.input_coords.kind is CoordSourceKind.STEP
        and s.input_coords.ref
        and s.status != "queued"
    }
    unfinished_by_log: dict[str, list] = {}
    for step in project.steps.values():
        if step.status != "done":
            continue
        s, rt = step.settings, step.runtime
        md = s.minimization.value is False
        if (
            md
            and s.traj_interval_steps.is_known
            and s.traj_interval_steps.value
            and FileRole.TRAJECTORY not in step.files
        ):
            named = step.raw.get("unresolved_assignments", {}).get("MDCRD")
            findings.append(
                finding(
                    "CMP001",
                    f"{step.id} writes a frame every {s.traj_interval_steps.value} steps but its trajectory"
                    + (f" ({named})" if named else "")
                    + " is not in the project",
                    subject=step.id,
                )
            )
        if rt.completed.value is False:
            mid = step.id in continued
            log = step.files.get(FileRole.LOG)
            if not mid and log is not None:
                # several runs logged in one file (e.g. every replica of one REMD job): one finding per log
                unfinished_by_log.setdefault(log.path, []).append(step)
            else:
                findings.append(
                    finding(
                        "CMP002",
                        f"{step.id} did not finish"
                        + (
                            " but a later run continues from it"
                            if mid
                            else " (last run of its chain)"
                        ),
                        subject=step.id,
                        severity=Severity.ERROR if mid else None,
                    )
                )
        traj = step.trajectory
        if (
            md
            and traj.n_frames.is_known
            and rt.nsteps_completed.is_known
            and s.traj_interval_steps.is_known
            and s.traj_interval_steps.value
        ):
            expected = int(rt.nsteps_completed.value) // int(s.traj_interval_steps.value)
            if traj.first_frame_at_start:
                expected += 1
            if int(traj.n_frames.value) != expected:
                findings.append(
                    finding(
                        "CMP003",
                        f"{step.id}: trajectory has {traj.n_frames.value} frames, the log implies {expected}",
                        subject=step.id,
                        frames=traj.n_frames.value,
                        expected=expected,
                    )
                )
        interval = frame_interval_ps(step)
        if md and traj.last_time_ps.is_known and rt.end_time_ps.is_known and interval:
            tol = continuity_tolerance(interval)
            span_end = float(traj.last_time_ps.value)
            if abs(span_end - float(rt.end_time_ps.value)) > tol:
                findings.append(
                    finding(
                        "CMP004",
                        f"{step.id}: last frame at {span_end:g} ps but the log ends at {rt.end_time_ps.value:g} ps",
                        subject=step.id,
                    )
                )
    for log, steps in unfinished_by_log.items():
        note = steps[0].runtime.completed.note
        if len(steps) == 1:
            msg = f"{steps[0].id} did not finish (last run of its chain)"
        else:
            msg = f"{len(steps)} runs logged in {log} did not finish (last runs of their chains)"
        findings.append(finding("CMP002", msg + (f": {note}" if note else ""),
                                subject=steps[0].id if len(steps) == 1 else None,
                                paths=[log], data={"steps": [st.id for st in steps]}))  # fmt: skip
    return findings
