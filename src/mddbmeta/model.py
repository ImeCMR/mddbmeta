"""The engine-neutral model: Project -> Replica -> Phase -> Step.

* A **Step** is one engine invocation (one AMBER mdin/mdout pair, one GROMACS
  mdrun, ...). It owns its files, its settings (from the control file), its
  runtime facts (from the log), and a declared input-coordinate source.
* A **Phase** is a role-bearing grouping of steps (minimization, heating,
  equilibration, production). It owns no files.
* A **Replica** is one independent member of the project -- it becomes one
  ``mds[]`` entry in MDDB's inputs. Its production chain is the ordered list
  of production steps whose trajectories MDDB will merge.
* A **Project** is the whole directory: one engine, a topology pool, a
  starting structure, replicas, findings and mined MDDB fields.

All paths stored on the model are POSIX paths relative to ``Project.root``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from mddbmeta.provenance import FileLoadError, Finding, Mined


class FileKind(str, Enum):
    TOPOLOGY = "topology"
    CONTROL = "control"
    LOG = "log"
    COORDS = "coords"  # a single-frame structure / restart
    TRAJECTORY = "trajectory"
    BUILD_LOG = "build_log"  # e.g. tleap logs/scripts (force-field evidence)
    UNKNOWN = "unknown"


class FileRole(str, Enum):
    CONTROL = "control"
    LOG = "log"
    TOPOLOGY = "topology"
    INPUT_COORDS = "input_coords"
    OUTPUT_RESTART = "output_restart"
    TRAJECTORY = "trajectory"
    REF_COORDS = "ref_coords"


class PhaseRole(str, Enum):
    MINIMIZATION = "minimization"
    HEATING = "heating"
    EQUILIBRATION = "equilibration"
    PRODUCTION = "production"
    UNKNOWN = "unknown"


PHASE_ORDER = [
    PhaseRole.MINIMIZATION,
    PhaseRole.HEATING,
    PhaseRole.EQUILIBRATION,
    PhaseRole.PRODUCTION,
    PhaseRole.UNKNOWN,
]


def _missing() -> Mined[Any]:
    return Mined.missing()


@dataclass
class FileRef:
    path: str  # relative to the project root, POSIX
    kind: FileKind
    format: str | None = None  # e.g. "prmtop", "netcdf", "ascii_crd", "rst7"
    size: int | None = None
    sha256: str | None = None

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"path": self.path, "kind": self.kind.value}
        if self.format:
            out["format"] = self.format
        if self.sha256:
            out["sha256"] = self.sha256
        return out


@dataclass
class StepSettings:
    """What the control file asked for, in engine-neutral terms."""

    minimization: Mined[bool] = field(default_factory=_missing)
    nsteps: Mined[int] = field(default_factory=_missing)
    dt_ps: Mined[float] = field(default_factory=_missing)
    temp_K: Mined[float] = field(default_factory=_missing)
    temp_initial_K: Mined[float] = field(default_factory=_missing)
    thermostat: Mined[str] = field(default_factory=_missing)
    barostat: Mined[str] = field(default_factory=_missing)
    periodic: Mined[bool] = field(default_factory=_missing)
    constant_pressure: Mined[bool] = field(default_factory=_missing)
    restraints: Mined[bool] = field(default_factory=_missing)
    restart: Mined[bool] = field(default_factory=_missing)  # continues velocities from a restart
    traj_interval_steps: Mined[int] = field(default_factory=_missing)
    restart_interval_steps: Mined[int] = field(default_factory=_missing)
    ensemble: Mined[str] = field(default_factory=_missing)
    cutoff_A: Mined[float] = field(default_factory=_missing)
    temp_ramp: Mined[bool] = field(default_factory=_missing)
    implicit_solvent: Mined[bool] = field(default_factory=_missing)

    def to_dict(self) -> dict[str, Any]:
        return {k: v.to_dict() for k, v in self.__dict__.items() if v.is_known}


@dataclass
class StepRuntime:
    """What the log says actually happened."""

    program: Mined[str] = field(default_factory=_missing)
    version: Mined[str] = field(default_factory=_missing)
    variant: Mined[str] = field(default_factory=_missing)
    start_time_ps: Mined[float] = field(default_factory=_missing)
    end_time_ps: Mined[float] = field(default_factory=_missing)
    nsteps_completed: Mined[int] = field(default_factory=_missing)
    completed: Mined[bool] = field(default_factory=_missing)
    wall_time_s: Mined[float] = field(default_factory=_missing)
    natom: Mined[int] = field(default_factory=_missing)

    def to_dict(self) -> dict[str, Any]:
        return {k: v.to_dict() for k, v in self.__dict__.items() if v.is_known}


@dataclass
class TrajectoryInfo:
    """What the trajectory file itself says."""

    format: Mined[str] = field(default_factory=_missing)
    n_frames: Mined[int] = field(default_factory=_missing)
    natom: Mined[int] = field(default_factory=_missing)
    first_time_ps: Mined[float] = field(default_factory=_missing)
    last_time_ps: Mined[float] = field(default_factory=_missing)
    frame_interval_ps: Mined[float] = field(default_factory=_missing)
    program: Mined[str] = field(default_factory=_missing)
    program_version: Mined[str] = field(default_factory=_missing)
    # GROMACS writes the step-0 frame; AMBER's first frame is one interval in
    first_frame_at_start: bool = False

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            k: v.to_dict() for k, v in self.__dict__.items() if isinstance(v, Mined) and v.is_known
        }
        if out and self.first_frame_at_start:
            out["first_frame_at_start"] = True
        return out


class CoordSourceKind(str, Enum):
    STARTING_STRUCTURE = "starting_structure"
    STEP = "step"
    PATH = "path"
    UNKNOWN = "unknown"


@dataclass
class CoordSource:
    kind: CoordSourceKind = CoordSourceKind.UNKNOWN
    ref: str | None = None  # a Step.id when kind == STEP
    path: str | None = None  # the coordinate file actually read
    method: str | None = None  # how it was resolved (file_assignments, stem, user, ...)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"source": self.kind.value}
        if self.ref:
            out["ref"] = self.ref
        if self.path:
            out["path"] = self.path
        if self.method:
            out["method"] = self.method
        return out


@dataclass
class Step:
    id: str  # posix "dir/stem" relative to root; stable across scans
    name: str  # stem
    dir: str  # posix directory relative to root ("" for root)
    engine: str
    files: dict[FileRole, FileRef] = field(default_factory=dict)
    settings: StepSettings = field(default_factory=StepSettings)
    runtime: StepRuntime = field(default_factory=StepRuntime)
    trajectory: TrajectoryInfo = field(default_factory=TrajectoryInfo)
    role: Mined[str] = field(default_factory=_missing)
    input_coords: CoordSource = field(default_factory=CoordSource)
    family: str | None = None  # "prod" for prod_0004
    seq_index: int | None = None  # 4 for prod_0004
    replica: str | None = None
    status: str = "done"  # done | queued (control only) | orphan (trajectory only)
    load_errors: list[FileLoadError] = field(default_factory=list)
    raw: dict[str, Any] = field(default_factory=dict)  # engine-specific extras (not serialized)

    @property
    def phase_role(self) -> PhaseRole:
        try:
            return PhaseRole(self.role.value) if self.role.value else PhaseRole.UNKNOWN
        except ValueError:
            return PhaseRole.UNKNOWN

    def file(self, role: FileRole) -> FileRef | None:
        return self.files.get(role)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "id": self.id,
            "name": self.name,
            "dir": self.dir,
            "engine": self.engine,
            "status": self.status,
            "role": self.role.to_dict(),
            "files": {r.value: f.to_dict() for r, f in self.files.items()},
            "input_coords": self.input_coords.to_dict(),
        }
        if self.replica is not None:
            out["replica"] = self.replica
        if self.family is not None:
            out["family"] = self.family
            out["seq_index"] = self.seq_index
        settings = self.settings.to_dict()
        if settings:
            out["settings"] = settings
        runtime = self.runtime.to_dict()
        if runtime:
            out["runtime"] = runtime
        traj = self.trajectory.to_dict()
        if traj:
            out["trajectory"] = traj
        if self.load_errors:
            out["load_errors"] = [e.to_dict() for e in self.load_errors]
        return out


@dataclass
class SystemInfo:
    """Chemical-system facts, usually from the topology."""

    source: str  # relative path of the file these facts come from
    natom: int | None = None
    nres: int | None = None
    title: str | None = None
    residue_counts: dict[str, int] = field(default_factory=dict)
    residue_atoms: dict[str, int] = field(default_factory=dict)  # atoms per residue, by name
    residue_atom_names: dict[str, list[str]] = field(default_factory=dict)
    periodic_box_type: int | None = None  # AMBER IFBOX: 0 none, 1 rectangular, 2 trunc. oct.
    box_lengths: tuple[float, float, float] | None = None
    box_angles: tuple[float, float, float] | None = None
    box_complete: bool = False  # True when all three box vectors are known (not just AMBER's beta)
    n_extra_points: int | None = None

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"source": self.source}
        for key in ("natom", "nres", "title", "periodic_box_type", "n_extra_points"):
            value = getattr(self, key)
            if value is not None:
                out[key] = value
        if self.box_lengths:
            out["box_lengths"] = list(self.box_lengths)
        if self.box_angles:
            out["box_angles"] = list(self.box_angles)
        if self.residue_counts:
            out["residue_counts"] = dict(self.residue_counts)
        return out


@dataclass
class Replica:
    id: str
    dir: str  # posix relative to root; "" when the replica lives in the root
    step_ids: list[str] = field(default_factory=list)
    production_chain: list[str] = field(default_factory=list)
    lineages: list[list[str]] = field(default_factory=list)
    topology: str | None = None  # relative path

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "dir": self.dir,
            "steps": list(self.step_ids),
            "production_chain": list(self.production_chain),
            "lineages": [list(chain) for chain in self.lineages],
            "topology": self.topology,
        }


@dataclass
class Phase:
    role: PhaseRole
    step_ids: list[str]


@dataclass
class Project:
    root: str  # absolute
    engine: str | None = None
    steps: dict[str, Step] = field(default_factory=dict)
    replicas: list[Replica] = field(default_factory=list)
    topologies: list[FileRef] = field(default_factory=list)
    starting_structure: FileRef | None = None
    coords_files: list[FileRef] = field(default_factory=list)
    build_logs: list[FileRef] = field(default_factory=list)
    system: SystemInfo | None = None
    findings: list[Finding] = field(default_factory=list)
    load_errors: list[FileLoadError] = field(default_factory=list)
    mined: dict[str, Mined[Any]] = field(default_factory=dict)
    md_mined: dict[str, dict[str, Mined[Any]]] = field(default_factory=dict)
    overrides: dict[str, Any] = field(default_factory=dict)
    # discovery context (not serialized): parsed files by relative path, and the adapter
    parsed: dict[str, Any] = field(default_factory=dict, repr=False, compare=False)
    adapter: Any = field(default=None, repr=False, compare=False)

    def step(self, step_id: str) -> Step:
        return self.steps[step_id]

    def replica(self, replica_id: str) -> Replica:
        for rep in self.replicas:
            if rep.id == replica_id:
                return rep
        raise KeyError(replica_id)

    def phases(self, replica: Replica) -> list[Phase]:
        by_role: dict[PhaseRole, list[str]] = {}
        for sid in replica.step_ids:
            by_role.setdefault(self.steps[sid].phase_role, []).append(sid)
        return [Phase(role, by_role[role]) for role in PHASE_ORDER if role in by_role]

    def production_steps(self, replica: Replica) -> list[Step]:
        return [self.steps[sid] for sid in replica.production_chain]

    def add_finding(self, f: Finding) -> None:
        if f not in self.findings:
            self.findings.append(f)

    def to_dict(self) -> dict[str, Any]:
        return {
            "root": self.root,
            "engine": self.engine,
            "topologies": [t.to_dict() for t in self.topologies],
            "starting_structure": self.starting_structure.to_dict()
            if self.starting_structure
            else None,
            "build_logs": [b.to_dict() for b in self.build_logs],
            "system": self.system.to_dict() if self.system else None,
            "replicas": [r.to_dict() for r in self.replicas],
            "steps": [s.to_dict() for s in self.steps.values()],
            "mined": {k: v.to_dict() for k, v in self.mined.items()},
            "md_mined": {
                rid: {k: v.to_dict() for k, v in fields.items()}
                for rid, fields in self.md_mined.items()
            },
            "findings": [f.to_dict() for f in self.findings],
            "load_errors": [e.to_dict() for e in self.load_errors],
        }
