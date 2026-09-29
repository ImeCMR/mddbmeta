"""Mine MDDB inputs fields from a discovered Project.

Project-level fields go to ``project.mined``; per-replica values that differ
across replicas go to ``project.md_mined[replica]`` (MDDB-workflow reads them
as per-MD overrides) and the project-level value is left unset.
"""

from __future__ import annotations

from typing import Any

from mddbmeta.engines.base import BuildEvidence
from mddbmeta.findings import finding
from mddbmeta.mining.boxtype import mine_boxtype
from mddbmeta.mining.forcefield import mine_forcefields
from mddbmeta.mining.water import mine_water
from mddbmeta.model import FileKind, FileRole, Project, Replica, Step
from mddbmeta.provenance import Confidence, Finding, Mined, Source
from mddbmeta.units import clean_float, ps_to_fs, ps_to_ns

MAX_SOURCES = 6
FRAME_REL_TOL = 0.01

# fields that may differ per MD (MDDB-workflow METADATA_FIELDS)
PER_MD_FIELDS = ("framestep", "timestep", "temp", "ensemble")


def unify(values: list[Mined[Any]], what: str) -> tuple[Mined[Any], list[Any]]:
    """Merge agreeing values; returns (merged or missing, distinct values)."""
    known = [m for m in values if m.is_known]
    distinct: list[Any] = []
    for m in known:
        if m.value not in distinct:
            distinct.append(m.value)
    if not known:
        return Mined.missing(f"{what}: not stated by any production run"), []
    if len(distinct) > 1:
        return (
            Mined(
                value=None,
                method="conflict",
                confidence=Confidence.HEURISTIC,
                alternatives=tuple(distinct),
                note=f"{what}: production runs disagree ({', '.join(map(str, distinct))})",
            ),
            distinct,
        )
    sources: list[Source] = []
    for m in known:
        for s in m.sources:
            if s not in sources:
                sources.append(s)
    weakest = min((m.confidence for m in known), key=lambda c: c.rank)
    methods = sorted({m.method for m in known})
    note = None
    if len(sources) > MAX_SOURCES:
        note = f"{len(sources)} sources agree; first {MAX_SOURCES} listed"
        sources = sources[:MAX_SOURCES]
    return (
        Mined(
            value=distinct[0],
            unit=known[0].unit,
            sources=tuple(sources),
            method=" + ".join(methods),
            confidence=weakest,
            note=note,
        ),
        distinct,
    )


def _convert(m: Mined[Any], value: Any, unit: str, how: str) -> Mined[Any]:
    """A pure unit conversion keeps the input's confidence."""
    return m.with_value(value, unit=unit, method=f"{m.method}; {how}")


def _prod(project: Project, rep: Replica) -> list[Step]:
    return [s for s in project.production_steps(rep) if s.status in ("done", "orphan")]


def _replica_values(project: Project, rep: Replica) -> tuple[dict[str, Mined[Any]], list[Finding]]:
    findings: list[Finding] = []
    prod = _prod(project, rep)
    done = [s for s in prod if s.status == "done"]
    out: dict[str, Mined[Any]] = {}

    dt, _ = unify([s.settings.dt_ps for s in done], "timestep")
    out["timestep"] = _convert(dt, ps_to_fs(dt.value), "fs", "ps->fs") if dt.is_known else dt

    temp, _ = unify([s.settings.temp_K for s in done], "temperature")
    out["temp"] = temp
    ens, _ = unify([s.settings.ensemble for s in done], "ensemble")
    out["ensemble"] = ens

    # framestep: from the control files and, independently, from the trajectories themselves
    ctl_vals = []
    for s in done:
        if (
            s.settings.traj_interval_steps.is_known
            and s.settings.dt_ps.is_known
            and s.settings.traj_interval_steps.value
        ):
            ps = clean_float(
                float(s.settings.traj_interval_steps.value) * float(s.settings.dt_ps.value)
            )
            ctl_vals.append(
                Mined(
                    value=ps,
                    unit="ps",
                    sources=s.settings.traj_interval_steps.sources + s.settings.dt_ps.sources,
                    method="derived: ntwx*dt",
                    confidence=min(
                        (s.settings.traj_interval_steps.confidence, s.settings.dt_ps.confidence),
                        key=lambda c: c.rank,
                    ),
                )
            )
    ctl, _ = unify(ctl_vals, "frame interval (control files)")
    traj_vals = [
        s.trajectory.frame_interval_ps for s in prod if s.trajectory.frame_interval_ps.is_known
    ]
    trj = _unify_close(traj_vals, "frame interval (trajectories)")
    chosen = trj if trj.is_known else ctl
    if trj.is_known and ctl.is_known and abs(trj.value - ctl.value) > FRAME_REL_TOL * ctl.value:
        findings.append(
            finding(
                "CONS002",
                f"{rep.id}: control files say a frame every {ctl.value:g} ps, trajectories have "
                f"{trj.value:g} ps; using the trajectories",
                subject=rep.id,
            )
        )
    out["framestep"] = (
        _convert(chosen, ps_to_ns(chosen.value), "ns", "ps->ns") if chosen.is_known else chosen
    )
    return out, findings


def _unify_close(values: list[Mined[Any]], what: str) -> Mined[Any]:
    if not values:
        return Mined.missing(f"{what}: none")
    ref = values[0].value
    if all(abs(v.value - ref) <= FRAME_REL_TOL * abs(ref) for v in values):
        merged, _ = unify([v.with_value(ref) for v in values], what)
        return merged
    return Mined(value=None, method="conflict", confidence=Confidence.HEURISTIC,
                 alternatives=tuple(sorted({v.value for v in values})), note=f"{what}: disagree")  # fmt: skip


def mine(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    adapter = project.adapter
    mined: dict[str, Mined[Any]] = {}
    per_rep: dict[str, dict[str, Mined[Any]]] = {}
    for rep in project.replicas:
        values, f = _replica_values(project, rep)
        per_rep[rep.id] = values
        findings.extend(f)

    # fields that may be set per MD
    md_mined: dict[str, dict[str, Mined[Any]]] = {rep.id: {} for rep in project.replicas}
    for key in PER_MD_FIELDS:
        merged, distinct = unify([per_rep[r.id][key] for r in project.replicas], key)
        conflicted = any(per_rep[r.id][key].method == "conflict" for r in project.replicas)
        if len(distinct) > 1 and not conflicted:
            mined[key] = Mined(value=None, method="per-md", confidence=Confidence.HEURISTIC,
                               note=f"differs across replicas: {', '.join(map(str, distinct))}")  # fmt: skip
            for r in project.replicas:
                if per_rep[r.id][key].is_known:
                    md_mined[r.id][key] = per_rep[r.id][key]
            findings.append(
                finding(
                    "MINE003",
                    f"{key} differs across replicas ({', '.join(map(str, distinct))}); "
                    "written per MD",
                    data={"field": key},
                )  # fmt: skip
            )
        else:
            mined[key] = merged if not conflicted else Mined.missing(f"{key}: runs disagree")

    all_prod = [s for r in project.replicas for s in _prod(project, r)]
    done = [s for s in all_prod if s.status == "done"]

    # program / version
    programs = [s.runtime.program for s in done if s.runtime.program.is_known]
    if programs and adapter is not None:
        variants = sorted({s.runtime.variant.value for s in done if s.runtime.variant.is_known})
        mined["program"] = Mined(
            value=adapter.program,
            sources=tuple(dict.fromkeys(src for m in programs for src in m.sources))[:MAX_SOURCES],
            method=f"engine log header ({', '.join(sorted({m.value for m in programs}))})",
            confidence=Confidence.EXACT,
            note=f"executable: {', '.join(variants)}" if variants else None,
        )
    else:
        tprog = [s.trajectory.program for s in all_prod if s.trajectory.program.is_known]
        amberish = [m for m in tprog if str(m.value).lower().startswith(ENGINE_WRITERS)]
        if amberish and adapter is not None:
            mined["program"] = Mined(value=adapter.program, sources=amberish[0].sources,
                                     method="trajectory attribute program", confidence=Confidence.DERIVED)  # fmt: skip
        else:
            mined["program"] = Mined.missing("no engine log; the program is not recorded")
    version, distinct = unify([s.runtime.version for s in done], "version")
    if not version.is_known and not distinct:
        # only trust a trajectory's programVersion when the engine itself wrote the file
        # (converters such as cpptraj or MDAnalysis stamp their own name and version)
        tver = [
            s.trajectory.program_version.with_value(
                _clean_version(s.trajectory.program_version.value), confidence=Confidence.DERIVED
            )  # fmt: skip
            for s in all_prod
            if s.trajectory.program_version.is_known and _engine_wrote(s)
        ]
        version, _ = unify(tver, "version")
    mined["version"] = version

    # system-level fields
    system = project.system
    # the topology-side box (for GROMACS: a production structure) beats arbitrary coordinate files
    has_box = system is not None and system.box_complete
    mined["boxtype"] = mine_boxtype(system, None if has_box else _coords_box(project))
    evidence = _evidence(project)
    mined["wat"] = mine_water(system, evidence)
    tops = [t.path for t in project.topologies]
    chamber = bool(project.system and adapter and project.system.source in project.parsed
                   and adapter.is_chamber(project.parsed[project.system.source]))  # fmt: skip
    mined["ff"] = mine_forcefields(evidence, tops, chamber)

    with_traj = [
        r
        for r in project.replicas
        if any(FileRole.TRAJECTORY in s.files for s in _prod(project, r))
    ]
    ensembles = [s for s in all_prod if s.raw.get("mddb_type") == "ensemble"]
    if all_prod and len(ensembles) == len(all_prod):
        mined["type"] = Mined(value="ensemble", method="derived: frames are ensemble members, not a time series",
                              confidence=Confidence.DERIVED)  # fmt: skip
        mined["framestep"] = Mined(value=None, method="not applicable: ensemble", confidence=Confidence.DERIVED,
                                   note="MDDB: framestep may be None for an ensemble")  # fmt: skip
        for rep in project.replicas:
            md_mined[rep.id].pop("framestep", None)
    elif project.replicas and len(with_traj) == len(project.replicas):
        mined["type"] = Mined(value="trajectory", method="derived: every replica has a time-ordered trajectory",
                              confidence=Confidence.DERIVED)  # fmt: skip
    else:
        mined["type"] = Mined.missing("not every replica has a production trajectory")

    for key in ("method", "metadditions"):
        vals = [s.raw.get(f"mddb_{key}") for s in done if s.raw.get(f"mddb_{key}") is not None]
        if vals and all(v == vals[0] for v in vals) and len(vals) == len(done):
            src = done[0].files.get(FileRole.CONTROL) or done[0].files.get(FileRole.LOG)
            mined[key] = Mined(value=vals[0], sources=(Source(src.path, key),) if src else (),
                               method=f"derived: {project.engine} run files", confidence=Confidence.DERIVED)  # fmt: skip

    # a coordinate file the whole production shares (NAMD's COORDINATE PDB), for mwf's structure input
    structs = {s.raw.get("structure") for s in done if s.raw.get("structure")}
    if len(structs) == 1:
        path = structs.pop()
        mined["input_structure_filepath"] = Mined(value=path, sources=(Source(path, "coordinates"),),
                                                  method="structure every production run read",
                                                  confidence=Confidence.EXACT)  # fmt: skip

    # an OpenMM system built from a PDB + ForceField XML has no topology file MDDB-workflow can read
    if not project.topologies and "input_structure_filepath" in mined:
        mined["input_topology_filepath"] = Mined(
            value="no", method="no topology file (structure only)", confidence=Confidence.DERIVED,
            note="MDDB-workflow's 'no' flag: it skips charge-based analyses and guesses bonds",
        )  # fmt: skip
        findings.append(finding(
            "MINE004",
            "no topology file (the system was built from a structure and force-field files); MDDB-workflow "
            "will skip charge-based analyses and guess bonds -- export a prmtop/psf if you need them",
            data={"field": "input_topology_filepath"},
        ))  # fmt: skip

    # files and mds[]
    topo_by_rep = {r.id: r.topology for r in project.replicas}
    shared = len(set(topo_by_rep.values())) == 1 and None not in topo_by_rep.values()
    if shared and project.replicas:
        tp = next(iter(topo_by_rep.values()))
        mined["input_topology_filepath"] = Mined(value=tp, sources=(Source(tp, "topology"),),
                                                 method="shared by all replicas", confidence=Confidence.EXACT)  # fmt: skip
    elif project.replicas and None not in topo_by_rep.values():
        findings.append(finding("REP003", "replicas use different topologies; written per MD"))
    for i, rep in enumerate(project.replicas, start=1):
        md = md_mined[rep.id]
        md["name"] = Mined(value=rep.id, method="replica id", confidence=Confidence.DERIVED)
        shared_dir = sum(1 for r in project.replicas if r.dir == rep.dir) > 1
        md["mdir"] = Mined(
            value=(rep.id if shared_dir and rep.dir else rep.dir) or f"replica_{i}",
            method="replica directory"
            if rep.dir
            else "new directory (replica lives in the project root)",
            confidence=Confidence.DERIVED,
        )
        trajs = [
            s.files[FileRole.TRAJECTORY].path
            for s in _prod(project, rep)
            if FileRole.TRAJECTORY in s.files
        ]
        md["input_trajectory_filepaths"] = Mined(
            value=trajs,
            sources=tuple(Source(p, "production chain") for p in trajs)[:MAX_SOURCES],
            method="production chain order",
            confidence=Confidence.EXACT if trajs else Confidence.HEURISTIC,
        )
        if not shared and rep.topology:
            md["input_topology_filepath"] = Mined(value=rep.topology, sources=(Source(rep.topology, "topology"),),
                                                  method="replica topology", confidence=Confidence.EXACT)  # fmt: skip

    for key in (
        "program",
        "version",
        "timestep",
        "framestep",
        "temp",
        "ensemble",
        "boxtype",
        "wat",
        "ff",
        "type",
    ):
        m = mined.get(key)
        if m is None or m.method in (
            "per-md",
            "non-periodic",
            "no water residues",
            "not applicable: ensemble",
        ):
            continue
        if not m.is_known:
            findings.append(
                finding("MINE001", f"{key}: {m.note or 'could not be mined'}", data={"field": key})
            )
        elif m.confidence is Confidence.HEURISTIC:
            alts = f" (candidates: {', '.join(map(str, m.alternatives))})" if m.alternatives else ""
            findings.append(
                finding(
                    "MINE002", f"{key}: best guess {m.value}{alts}; confirm it", data={"field": key}
                )
            )

    project.mined = mined
    project.md_mined = md_mined
    return findings


ENGINE_WRITERS = ("pmemd", "sander")


def _engine_wrote(step: Step) -> bool:
    prog = step.trajectory.program.value
    return prog is not None and str(prog).lower().startswith(ENGINE_WRITERS)


def _clean_version(value: Any) -> str:
    text = str(value).strip()
    for prefix in ("Version ", "version "):
        if text.startswith(prefix):
            text = text[len(prefix) :]
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text


def _coords_box(project: Project):
    adapter = project.adapter
    if adapter is None:
        return None
    candidates = []
    if project.starting_structure is not None:
        candidates.append(project.starting_structure.path)
    candidates += [c.path for c in project.coords_files]
    for path in candidates:
        pf = project.parsed.get(path)
        if pf is None:
            continue
        box = adapter.coords_box(pf)
        if box is not None:
            return box[0], box[1], path
    return None


def _evidence(project: Project) -> list[BuildEvidence]:
    adapter = project.adapter
    out: list[BuildEvidence] = []
    if adapter is None:
        return out
    for ref in project.build_logs:
        pf = project.parsed.get(ref.path)
        if pf is None or pf.kind is not FileKind.BUILD_LOG:
            continue
        ev = adapter.build_evidence(pf)
        if ev is not None:
            out.append(ev)
    for ref in project.topologies:
        pf = project.parsed.get(ref.path)
        if pf is None:
            continue
        ev = adapter.topology_evidence(pf)
        if ev is not None and (ev.force_fields or ev.water_models):
            out.append(ev)
    prod = [project.steps[s] for r in project.replicas for s in r.production_chain]
    out.extend(ev for ev in adapter.run_evidence(prod) if ev.force_fields or ev.water_models)
    return out
