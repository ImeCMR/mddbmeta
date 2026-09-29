"""Input-coordinate resolution, lineage chains and production-chain ordering.

A step's ``input_coords`` names the coordinate file the engine read. When
that file is another step's output restart, the step *continues from* that
step -- an explicit, evidence-backed edge. Production chains are ordered by
these edges; only where no edge exists do we fall back to times, then to
sequence numbers / natural name order, and that fallback is announced
(ORDER001). Filesystem or glob order is never used.
"""

from __future__ import annotations

from collections import Counter

from mddbmeta.discovery.replicas import natural_key
from mddbmeta.findings import finding
from mddbmeta.model import CoordSource, CoordSourceKind, FileRole, PhaseRole, Project, Replica, Step
from mddbmeta.provenance import Finding


def resolve_input_coords(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    producers: dict[str, str] = {}
    for step in project.steps.values():
        ref = step.files.get(FileRole.OUTPUT_RESTART)
        if ref is None:
            continue
        if ref.path in producers:
            findings.append(
                finding(
                    "LIN001",
                    f"{ref.path} is written by both {producers[ref.path]} and {step.id}",
                    subject=step.id,
                    paths=[ref.path],
                )
            )
        producers[ref.path] = step.id

    # Starting structure: coordinates read by some step but written by none.
    consumed_unproduced = Counter(
        s.files[FileRole.INPUT_COORDS].path
        for s in project.steps.values()
        if FileRole.INPUT_COORDS in s.files
        and s.files[FileRole.INPUT_COORDS].path not in producers
        and s.input_coords.kind is not CoordSourceKind.STEP  # already explained by the adapter
    )
    starting = None
    if len(consumed_unproduced) == 1:
        starting = next(iter(consumed_unproduced))
    for step in project.steps.values():
        if step.input_coords.kind is CoordSourceKind.STEP:
            continue  # the adapter already resolved it from engine evidence
        ref = step.files.get(FileRole.INPUT_COORDS)
        if ref is not None:
            if ref.path in producers and producers[ref.path] != step.id:
                step.input_coords = CoordSource(
                    CoordSourceKind.STEP,
                    ref=producers[ref.path],
                    path=ref.path,
                    method="engine log",
                )
            elif ref.path == starting:
                step.input_coords = CoordSource(
                    CoordSourceKind.STARTING_STRUCTURE, path=ref.path, method="engine log"
                )
            else:
                step.input_coords = CoordSource(
                    CoordSourceKind.PATH, path=ref.path, method="engine log"
                )
            continue
        missing = step.raw.get("unresolved_assignments", {}).get("INPCRD")
        if missing:
            step.input_coords = CoordSource(
                CoordSourceKind.UNKNOWN, path=missing, method="named in the log but not found"
            )
    if starting is not None:
        for f in project.coords_files:
            if f.path == starting:
                project.starting_structure = f
                break

    findings.extend(_cycles(project))
    return findings


def check_forks(project: Project) -> list[Finding]:
    """Run after replicas are known: a fork inside one replica is suspicious."""
    findings: list[Finding] = []
    # Forks inside one replica are suspicious; across replicas they are how replicas start.
    consumers: dict[str, list[str]] = {}
    for step in project.steps.values():
        if step.input_coords.kind is CoordSourceKind.STEP and step.input_coords.ref:
            consumers.setdefault(step.input_coords.ref, []).append(step.id)
    for producer, cons in consumers.items():
        reps = Counter(project.steps[c].replica for c in cons)
        for rep, n in reps.items():
            if n > 1:
                findings.append(
                    finding(
                        "LIN001",
                        f"{producer} feeds {n} steps in the same replica: "
                        + ", ".join(c for c in cons if project.steps[c].replica == rep),
                        subject=producer,
                    )
                )
    return findings


def edges(project: Project) -> dict[str, str]:
    """consumer step id -> producer step id."""
    return {
        s.id: s.input_coords.ref
        for s in project.steps.values()
        if s.input_coords.kind is CoordSourceKind.STEP and s.input_coords.ref in project.steps
    }


def _cycles(project: Project) -> list[Finding]:
    findings = []
    for start in project.steps.values():
        seen = {start.id}
        cur = start
        while (
            cur.input_coords.kind is CoordSourceKind.STEP and cur.input_coords.ref in project.steps
        ):
            nxt = project.steps[cur.input_coords.ref]
            if nxt.id in seen:
                if nxt.id == start.id:
                    findings.append(
                        finding("LIN002", f"lineage cycle through {start.id}", subject=start.id)
                    )
                break
            seen.add(nxt.id)
            cur = nxt
    return findings[:1]


def start_time(step: Step) -> float | None:
    if step.runtime.start_time_ps.is_known:
        return float(step.runtime.start_time_ps.value)
    if step.trajectory.first_time_ps.is_known:
        # first frame is written one frame interval after the run starts
        first = float(step.trajectory.first_time_ps.value)
        if step.trajectory.first_frame_at_start:
            return first
        if step.trajectory.frame_interval_ps.is_known:
            return first - float(step.trajectory.frame_interval_ps.value)
        return first
    return None


def _order_key(step: Step) -> tuple:
    t = start_time(step)
    return (
        0 if t is not None else 1,
        t if t is not None else 0.0,
        step.seq_index if step.seq_index is not None else 10**12,
        natural_key(step.id),
    )


def build_chains(project: Project, replica: Replica) -> list[Finding]:
    """Fill replica.lineages and replica.production_chain."""
    findings: list[Finding] = []
    members = [project.steps[s] for s in replica.step_ids]
    ids = {s.id for s in members}

    # lineage chains: follow single consumers from each head
    consumers: dict[str, list[Step]] = {}
    for s in members:
        if s.input_coords.kind is CoordSourceKind.STEP and s.input_coords.ref in ids:
            consumers.setdefault(s.input_coords.ref, []).append(s)
    heads = [
        s
        for s in members
        if not (s.input_coords.kind is CoordSourceKind.STEP and s.input_coords.ref in ids)
    ]
    chains: list[list[str]] = []
    stack = sorted(heads, key=_order_key, reverse=True)
    while stack:
        cur = stack.pop()
        chain = [cur.id]
        while True:
            nxt = sorted(consumers.get(cur.id, []), key=_order_key)
            if len(nxt) != 1:
                stack.extend(reversed(nxt))
                break
            cur = nxt[0]
            chain.append(cur.id)
        chains.append(chain)
    replica.lineages = chains

    prod = [
        s
        for s in members
        if s.phase_role is PhaseRole.PRODUCTION and s.status in ("done", "orphan")
    ]
    replica.production_chain, how = order_production(prod)
    if how == "names":
        findings.append(
            finding(
                "ORDER001",
                f"{replica.id}: production order taken from run numbering/names "
                "(no restart links or times to prove it)",
                subject=replica.id,
                data={"order": list(replica.production_chain)},
            )
        )
    return findings


def order_production(prod: list[Step]) -> tuple[list[str], str]:
    """Topological order over restart edges; ties by start time, then numbering.

    Returns (ordered ids, how) with how in {"lineage", "time", "names", "single"}.
    """
    if len(prod) <= 1:
        return [s.id for s in prod], "single"
    ids = {s.id for s in prod}
    parent = {
        s.id: s.input_coords.ref
        for s in prod
        if s.input_coords.kind is CoordSourceKind.STEP and s.input_coords.ref in ids
    }
    children: dict[str, list[Step]] = {}
    for s in prod:
        if s.id in parent:
            children.setdefault(parent[s.id], []).append(s)
    ready = [s for s in prod if s.id not in parent]
    order: list[str] = []
    how = "lineage"
    while ready:
        ready.sort(key=_order_key)
        if len(ready) > 1:
            times = [start_time(s) for s in ready]
            if all(t is not None for t in times) and len(set(times)) == len(times):
                how = "time" if how == "lineage" else how
            else:
                how = "names"
        cur = ready.pop(0)
        order.append(cur.id)
        ready.extend(children.get(cur.id, []))
    # anything left is on a cycle; append deterministically
    for s in sorted(prod, key=_order_key):
        if s.id not in order:
            order.append(s.id)
            how = "names"
    return order, how
