"""discover(): directory -> fully resolved, validated and mined Project."""

from __future__ import annotations

import posixpath
from collections import Counter
from pathlib import Path
from typing import Any

from mddbmeta.discovery.lineage import build_chains, check_forks, edges, resolve_input_coords
from mddbmeta.discovery.replicas import ReplicaResult, infer_replicas, natural_key
from mddbmeta.discovery.roles import classify, family_index
from mddbmeta.discovery.scan import scan
from mddbmeta.engines.base import EngineAdapter, ParsedFile
from mddbmeta.engines.registry import available_adapters, get_adapter
from mddbmeta.errors import MddbmetaError
from mddbmeta.export.manifest import read_manifest, split_overrides
from mddbmeta.findings import finding
from mddbmeta.mining.fields import mine
from mddbmeta.model import FileKind, FileRole, PhaseRole, Project, Replica, Step
from mddbmeta.provenance import Confidence, Mined
from mddbmeta.validate import run_all

GENERIC_FORMATS = {"pdb", "pdbx"}


def discover(
    root: str | Path,
    engine: str | None = None,
    replica_glob: str | None = None,
    max_depth: int = 4,
    exclude: list[str] | None = None,
    overrides: dict[str, Any] | None = None,
) -> Project:
    root = Path(root).expanduser().resolve()
    roles, chains, field_overrides = split_overrides(overrides or {})
    if not root.is_dir():
        raise MddbmetaError(f"{root} is not a directory")
    adapters: list[EngineAdapter] = (
        [get_adapter(engine)] if engine else [cls() for cls in available_adapters().values()]
    )
    scanned, owner = scan(root, adapters, max_depth=max_depth, exclude=exclude)
    project = Project(root=str(root))
    project.findings.extend(scanned.findings)
    project.load_errors.extend(scanned.load_errors)
    if not scanned.files:
        return project

    counts = Counter(owner[f.rel].name for f in scanned.files)
    logs = Counter(owner[f.rel].name for f in scanned.files if f.kind is FileKind.LOG)
    chosen_name = max(counts, key=lambda n: (logs.get(n, 0), counts[n]))
    adapter = next(a for a in adapters if a.name == chosen_name)
    files = [f for f in scanned.files if owner[f.rel] is adapter or adapter.claim_foreign(f)]
    # files of another engine that the chosen one does not use (adopted inputs don't count)
    # (engine-neutral formats such as PDB are never evidence of another engine)
    stray = Counter(
        owner[f.rel].name
        for f in scanned.files
        if owner[f.rel] is not adapter and f not in files and f.format not in GENERIC_FORMATS
    )
    if stray:
        project.add_finding(
            finding(
                "ENG001",
                f"files from several engines ({', '.join(sorted({chosen_name, *stray}))}); using {chosen_name}",
                data={chosen_name: counts[chosen_name], **dict(stray)},
            )
        )
    project.engine = adapter.name
    project.adapter = adapter
    project.parsed = {f.rel: f for f in files}

    steps, build_findings = adapter.build_steps(root, files)
    project.findings.extend(build_findings)
    project.steps = {s.id: s for s in steps}
    project.topologies = [f.ref() for f in files if f.kind is FileKind.TOPOLOGY and f.ok]
    project.coords_files = [f.ref() for f in files if f.kind is FileKind.COORDS and f.ok]
    project.build_logs = [f.ref() for f in files if f.kind is FileKind.BUILD_LOG and f.ok]

    _assign_topologies(project, files, adapter)
    for step in steps:
        if not step.role.is_known:  # keep roles the adapter knows from the engine itself
            step.role, role_findings = classify(step)
            project.findings.extend(role_findings)
        fam = family_index(step.name)
        if fam is not None:
            step.family, step.seq_index = fam
    for sid, role in roles.items():
        if sid not in project.steps:
            raise MddbmetaError(f"override steps.{sid}.role: no step '{sid}'")
        if role not in {r.value for r in PhaseRole}:
            raise MddbmetaError(f"override steps.{sid}.role: '{role}' is not a phase role")
        project.steps[sid].role = Mined(value=role, method="user", confidence=Confidence.EXACT)

    project.findings.extend(resolve_input_coords(project))
    declared = None if replica_glob else adapter.declare_replicas(steps)
    if declared is not None:
        reps = _declared_replicas(steps, declared)
    else:
        reps = infer_replicas(steps, replica_glob, edges=edges(project))
    project.findings.extend(reps.findings)
    for step in steps:
        step.replica = reps.assignment.get(step.id)
    project.replicas = reps.replicas
    project.findings.extend(check_forks(project))
    for rep in project.replicas:
        project.findings.extend(build_chains(project, rep))
        if rep.id in chains:
            unknown = [sid for sid in chains[rep.id] if sid not in rep.step_ids]
            if unknown:
                raise MddbmetaError(
                    f"override replicas.{rep.id}.production_chain: unknown steps {unknown}"
                )
            rep.production_chain = chains[rep.id]
            project.findings = [
                f for f in project.findings if not (f.code == "ORDER001" and f.subject == rep.id)
            ]
        rep.topology = _replica_topology(project, rep)
    missing_reps = set(chains) - {r.id for r in project.replicas}
    if missing_reps:
        raise MddbmetaError(
            f"override replicas.*.production_chain: unknown replica(s) {sorted(missing_reps)}"
        )

    system_top = _system_topology(project)
    if system_top is not None:
        project.system = adapter.system_info(project.parsed[system_top])
    project.findings.extend(run_all(project))
    project.findings.extend(mine(project))
    if reps.ambiguous:
        project.mined.clear()
        project.md_mined.clear()
    project.overrides = field_overrides
    project.findings = [
        f
        for f in project.findings
        if not (f.code in ("MINE001", "MINE002") and f.data.get("field") in field_overrides)
    ]
    _dedupe(project)
    return project


def load_project(target: str | Path, **options: Any) -> Project:
    """Discover a directory, or re-discover the project a manifest describes (with its overrides)."""
    path = Path(target).expanduser()
    options = {k: v for k, v in options.items() if v is not None}
    if path.is_dir():
        return discover(path, **options)
    if not path.exists():
        raise MddbmetaError(f"{path} does not exist")
    root, disc, overrides = read_manifest(path)
    opts = {
        k: v for k, v in disc.items() if k in ("engine", "replica_glob", "max_depth", "exclude")
    }
    opts.update(options)
    return discover(root, overrides=overrides, **opts)


def _assign_topologies(project: Project, files: list[ParsedFile], adapter: EngineAdapter) -> None:
    natoms: dict[str, int | None] = {}
    for f in files:
        if f.kind is FileKind.TOPOLOGY and f.ok:
            info = adapter.system_info(f)
            natoms[f.rel] = info.natom if info else None
    if not natoms:
        return
    for step in project.steps.values():
        if FileRole.TOPOLOGY in step.files:
            step.raw["topology_method"] = "engine log"
            continue
        natom = step.runtime.natom.value or step.trajectory.natom.value
        cands = [t for t, n in natoms.items() if natom is None or n == natom]
        if not cands:
            continue
        if len(cands) > 1:
            cands.sort(key=lambda t, d=step.dir: _closeness(d, t))
        pf = project.parsed[cands[0]]
        step.files[FileRole.TOPOLOGY] = pf.ref()
        step.raw["topology_method"] = (
            "only topology" if len(natoms) == 1 else "atom count / proximity"
        )


def _declared_replicas(steps: list[Step], declared: dict[str, str | None]) -> ReplicaResult:
    """Replicas an adapter knows from the engine itself (e.g. MELD ladder positions)."""
    ids = sorted({r for r in declared.values() if r}, key=natural_key)
    replicas = []
    for rid in ids:
        members = [s for s in steps if declared.get(s.id) == rid]
        dirs = {s.dir for s in members}
        replicas.append(Replica(id=rid, dir=dirs.pop() if len(dirs) == 1 else "",
                                step_ids=[s.id for s in members]))  # fmt: skip
    return ReplicaResult(replicas=replicas, assignment=dict(declared))


def _closeness(step_dir: str, topology: str) -> tuple[int, str]:
    """Sort key: the deepest topology directory that contains the step comes first."""
    tdir = posixpath.dirname(topology)
    inside = tdir == "" or step_dir == tdir or step_dir.startswith(tdir + "/")
    return (-(len(tdir.split("/")) if tdir else 0) if inside else 1, topology)


def _replica_topology(project: Project, rep) -> str | None:
    steps: list[Step] = [project.steps[s] for s in rep.production_chain] or [
        project.steps[s] for s in rep.step_ids
    ]
    counts = Counter(s.files[FileRole.TOPOLOGY].path for s in steps if FileRole.TOPOLOGY in s.files)
    return counts.most_common(1)[0][0] if counts else None


def _system_topology(project: Project) -> str | None:
    counts = Counter(rep.topology for rep in project.replicas if rep.topology)
    if counts:
        return counts.most_common(1)[0][0]
    if project.topologies:
        return project.topologies[0].path
    # no topology file (e.g. OpenMM from a PDB): the structure every production run read
    structs = Counter(
        project.steps[s].raw.get("structure")
        for rep in project.replicas
        for s in rep.production_chain
    )
    structs.pop(None, None)
    return structs.most_common(1)[0][0] if structs else None


def _dedupe(project: Project) -> None:
    seen = []
    for f in project.findings:
        if f not in seen:
            seen.append(f)
    project.findings = seen


def production_steps(project: Project) -> list[Step]:
    return [s for s in project.steps.values() if s.phase_role is PhaseRole.PRODUCTION]
