"""The engine adapter interface.

An adapter knows one MD engine's files. Discovery, replica inference, role
classification, lineage, validation and MDDB export are engine-neutral and
only talk to adapters through this interface, so adding GROMACS / NAMD /
OpenMM means writing one adapter, not touching the core.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

from mddbmeta.model import FileKind, FileRef, Step, SystemInfo, TrajectoryInfo
from mddbmeta.provenance import FileLoadError, Finding


@dataclass(frozen=True)
class SniffResult:
    kind: FileKind
    format: str  # adapter-specific format tag, e.g. "prmtop", "mdin", "netcdf_traj"
    score: float  # 0..1; content evidence scores higher than extension evidence


@dataclass
class ParsedFile:
    path: str  # absolute
    rel: str  # POSIX, relative to the project root
    kind: FileKind
    format: str
    result: Any = None  # adapter-specific parse result (None when parsing failed)
    errors: list[FileLoadError] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.result is not None

    def ref(self) -> FileRef:
        size = None
        try:
            size = Path(self.path).stat().st_size
        except OSError:
            pass
        return FileRef(path=self.rel, kind=self.kind, format=self.format, size=size)


@dataclass
class BuildEvidence:
    """What build logs (tleap, pdb2gmx, psfgen, ...) say about the system."""

    source: str  # relative path
    force_fields: list[str] = field(default_factory=list)
    water_models: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)  # topology files the log says it wrote
    linked: bool = False  # the evidence is part of the run itself (e.g. the OpenMM script that ran)


class EngineAdapter(ABC):
    name: ClassVar[str]
    program: ClassVar[str]  # the MDDB "program" value, e.g. "AMBER"

    @abstractmethod
    def sniff(self, path: Path, head: bytes) -> SniffResult | None:
        """Classify a file from its first bytes (content first, extension second)."""

    @abstractmethod
    def parse(self, path: Path, rel: str, sniffed: SniffResult) -> ParsedFile:
        """Parse one file. Must never raise: problems go in ParsedFile.errors."""

    @abstractmethod
    def build_steps(self, root: Path, files: list[ParsedFile]) -> tuple[list[Step], list[Finding]]:
        """Group parsed files into Steps and fill settings, runtime and trajectory facts.

        Each step's ``input_coords.path`` must be set when the engine records
        which coordinates the run read; the core resolves it into a lineage.
        """

    @abstractmethod
    def system_info(self, topology: ParsedFile) -> SystemInfo | None:
        """Chemical-system facts from a topology."""

    @abstractmethod
    def coords_natom_time(self, coords: ParsedFile) -> tuple[int | None, float | None]:
        """Atom count and simulation time (ps) stored in a coordinate/restart file."""

    def trajectory_info(self, traj: ParsedFile, natom: int | None) -> TrajectoryInfo:
        """Facts from a trajectory; ``natom`` helps formats that don't store it."""
        return TrajectoryInfo()

    def coords_box(
        self, coords: ParsedFile
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
        """(lengths, angles) of the periodic box stored in a coordinate file, if any."""
        return None

    def is_chamber(self, topology: ParsedFile) -> bool:
        """True for a topology converted from another force-field family (e.g. CHARMM)."""
        return False

    def build_evidence(self, build_log: ParsedFile) -> BuildEvidence | None:
        return None

    def declare_replicas(self, steps: list) -> dict[str, str | None] | None:
        """step id -> replica id, when the engine itself defines the replicas (e.g. REMD ladders)."""
        return None

    def run_evidence(self, steps: list) -> list[BuildEvidence]:
        """Force-field / water evidence recorded by the runs themselves (e.g. OpenMM scripts)."""
        return []

    def claim_foreign(self, pf: ParsedFile) -> bool:
        """Adopt a file another adapter sniffed (e.g. NAMD reading an AMBER prmtop)."""
        return False

    def topology_evidence(self, topology: ParsedFile) -> BuildEvidence | None:
        """Force-field / water evidence carried by a topology itself (e.g. GROMACS #include)."""
        return None

    def mddb_trajectory_format(self, fmt: str) -> str | None:
        """Map an adapter format tag to MDDB-workflow's trajectory format, or None."""
        return None

    def mddb_topology_format(self, fmt: str) -> str | None:
        return None
