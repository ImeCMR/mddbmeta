"""NAMD engine adapter.

One Step per NAMD log. The log's ``Info:`` lines name every file the run
used (structure, input coordinates, extended system, output and restart
names, dcd) and the resolved settings, including values a config computed
with Tcl. Continuation is exact: a run's ``BINARY COORDINATES`` is another
run's final output (``<outputName>.coor``) or restart (``<restartName>.coor``).

NAMD can also read AMBER (``amber on`` / ``parmfile``) topologies; those files
belong to the AMBER adapter's sniffer, so this adapter adopts them.
"""

from __future__ import annotations

import posixpath
from pathlib import Path

from mddbmeta.engines.amber import prmtop as amber_prmtop
from mddbmeta.engines.base import BuildEvidence, EngineAdapter, ParsedFile, SniffResult
from mddbmeta.engines.gromacs.structure import lengths_angles
from mddbmeta.engines.namd import config as cfgmod
from mddbmeta.engines.namd import files, log
from mddbmeta.engines.namd.semantics import settings_from_config, settings_from_log
from mddbmeta.findings import finding
from mddbmeta.model import (
    CoordSource,
    CoordSourceKind,
    FileKind,
    FileRole,
    Step,
    StepRuntime,
    SystemInfo,
    TrajectoryInfo,
)
from mddbmeta.provenance import (
    Confidence,
    FileLoadError,
    Finding,
    Mined,
    Source,
    classify_exception,
    stated,
)
from mddbmeta.units import clean_float

# CHARMM parameter-file stem -> (label, residue class that justifies listing it; None = always)
CHARMM_PARAMS = [
    ("par_all36m_prot", "CHARMM36m", None),
    ("par_all36_prot", "CHARMM36", None),
    ("par_all27_prot", "CHARMM27", None),
    ("par_all22_prot", "CHARMM22", None),
    ("par_all36_na", "CHARMM36 (nucleic acids)", "nucleic"),
    ("par_all36_lipid", "CHARMM36 (lipids)", "lipid"),
    ("par_all36_carb", "CHARMM36 (carbohydrates)", "carb"),
    ("par_all36_cgenff", "CGenFF", "other"),
]
PROTEIN = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HSD", "HSE", "HSP", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}  # fmt: skip
NUCLEIC = {"ADE", "CYT", "GUA", "THY", "URA", "DA", "DC", "DG", "DT", "A", "C", "G", "U"}
LIPID = {
    "POPC",
    "POPE",
    "POPS",
    "POPG",
    "DPPC",
    "DOPC",
    "DMPC",
    "DLPC",
    "CHL1",
    "PSM",
    "DOPE",
    "SAPI",
}
WATER_IONS = {
    "TIP3",
    "TIP4",
    "SPC",
    "SWM4",
    "SOD",
    "CLA",
    "POT",
    "CAL",
    "MG",
    "ZN2",
    "CES",
    "LIT",
    "RUB",
}
CARB_PREFIX = ("AGLC", "BGLC", "AMAN", "BMAN", "AGAL", "BGAL", "BGLCNA", "ANE5AC", "AFUC", "BXYL")


def _dir(rel: str) -> str:
    return posixpath.dirname(rel)


def _stem(rel: str) -> str:
    base = posixpath.basename(rel)
    return base.rsplit(".", 1)[0] if "." in base else base


class NamdAdapter(EngineAdapter):
    name = "namd"
    program = "NAMD"

    def __init__(self) -> None:
        self._by_rel: dict[str, ParsedFile] = {}

    # ------------------------------------------------------------------ sniff
    def sniff(self, path: Path, head: bytes) -> SniffResult | None:
        if files.is_dcd(head):
            return SniffResult(FileKind.TRAJECTORY, "dcd", 0.95)
        ext = path.suffix.lower()
        if ext in (".coor", ".vel"):
            try:
                size = path.stat().st_size
            except OSError:
                size = -1
            if files.binary_natom(head, size):
                return SniffResult(FileKind.COORDS, "namdvel" if ext == ".vel" else "namdbin", 0.9)
        if b"\x00" in head[:1024]:
            return None
        if files.looks_like_xsc(head):
            return SniffResult(FileKind.COORDS, "xsc", 0.9)
        if head.startswith(b"PSF"):
            return SniffResult(FileKind.TOPOLOGY, "psf", 0.95)
        if log.looks_like_log(head):
            return SniffResult(FileKind.LOG, "namdlog", 0.95)
        if cfgmod.looks_like_config(head, path.name):
            return SniffResult(FileKind.CONTROL, "namdconf", 0.8)
        if files.looks_like_pdb(head):
            return SniffResult(FileKind.COORDS, "pdb", 0.6)
        return None

    def claim_foreign(self, pf: ParsedFile) -> bool:
        """AMBER-format inputs of `amber on` runs (read with the AMBER parser)."""
        return pf.format in ("prmtop", "rst7") and pf.ok

    # ------------------------------------------------------------------ parse
    def parse(self, path: Path, rel: str, sniffed: SniffResult) -> ParsedFile:
        pf = ParsedFile(path=str(path), rel=rel, kind=sniffed.kind, format=sniffed.format)
        readers = {
            "dcd": files.read_dcd,
            "xsc": files.read_xsc,
            "psf": files.read_psf,
            "pdb": files.read_pdb,
            "namdlog": log.read_log,
            "namdconf": cfgmod.read_config,
            "namdbin": lambda p: files.NamdBinary(
                str(p), (Path(p).stat().st_size - 4) // 24, False
            ),
            "namdvel": lambda p: files.NamdBinary(str(p), (Path(p).stat().st_size - 4) // 24, True),
        }
        try:
            pf.result = readers[sniffed.format](path)
        except Exception as exc:
            pf.errors.append(
                FileLoadError(sniffed.kind.value, rel, classify_exception(exc), str(exc))
            )
            return pf
        if isinstance(pf.result, files.Dcd) and pf.result.truncated:
            pf.findings.append(
                finding("LOAD001", f"{rel}: last frame is truncated", subject=rel, paths=[rel])
            )
        return pf

    # ------------------------------------------------------------ build steps
    def build_steps(self, root: Path, parsed: list[ParsedFile]) -> tuple[list[Step], list[Finding]]:
        findings: list[Finding] = []
        by_rel = {f.rel: f for f in parsed}
        self._by_rel = by_rel
        claimed: set[str] = set()
        steps: list[Step] = []

        def find(run_dir: str, name: str | None) -> ParsedFile | None:
            if not name:
                return None
            for c in (posixpath.normpath(posixpath.join(run_dir, name)), posixpath.normpath(name)):
                if c in by_rel:
                    return by_rel[c]
            return None

        configs = [f for f in parsed if f.format == "namdconf" and f.ok]
        by_output: dict[tuple[str, str], ParsedFile] = {}
        for c in configs:
            out = c.result.get("outputname")
            if out:
                by_output[(_dir(c.rel), posixpath.normpath(posixpath.join(_dir(c.rel), out)))] = c

        outputs: dict[str, str] = {}  # coordinates file a run wrote -> step id
        for lf in sorted((f for f in parsed if f.kind is FileKind.LOG), key=lambda f: f.rel):
            run_dir, stem = _dir(lf.rel), _stem(lf.rel)
            step = Step(id=posixpath.join(run_dir, stem) if run_dir else stem, name=stem, dir=run_dir,
                        engine=self.name)  # fmt: skip
            step.files[FileRole.LOG] = lf.ref()
            claimed.add(lf.rel)
            step.load_errors.extend(lf.errors)
            nl: log.NamdLog | None = lf.result if lf.ok else None
            if nl is None:
                steps.append(step)
                continue
            out = nl.one("OUTPUT FILENAME")
            restart = nl.one("RESTART FILENAME")
            topo = find(
                run_dir,
                nl.one("STRUCTURE FILE") or nl.one("AMBER PARM FILE") or nl.one("PARM FILE"),
            )
            if topo is not None:
                step.files[FileRole.TOPOLOGY] = topo.ref()
            incoords = (
                nl.one("BINARY COORDINATES")
                or nl.one("COORDINATE PDB")
                or nl.one("AMBER COORDINATE FILE")
            )
            inc = find(run_dir, incoords)
            if inc is not None:
                step.files[FileRole.INPUT_COORDS] = inc.ref()
            elif incoords:
                step.raw["unresolved_assignments"] = {
                    "INPCRD": posixpath.join(run_dir, incoords) if run_dir else incoords
                }
            pdb = find(run_dir, nl.one("COORDINATE PDB"))
            if pdb is not None:
                step.raw["structure"] = pdb.rel
            traj = find(run_dir, nl.one("DCD FILENAME"))
            if traj is not None:
                step.files[FileRole.TRAJECTORY] = traj.ref()
                claimed.add(traj.rel)
            for base, is_restart in ((out, False), (restart, True)):
                if not base:
                    continue
                coor = posixpath.normpath(posixpath.join(run_dir, base + ".coor"))
                if is_restart:
                    outputs.setdefault(
                        coor, step.id
                    )  # a final output name wins over a restart name
                else:
                    outputs[coor] = step.id
            final = find(run_dir, out + ".coor") if out else None
            if final is not None:
                step.files[FileRole.OUTPUT_RESTART] = final.ref()
            ctl = (
                by_output.get((run_dir, posixpath.normpath(posixpath.join(run_dir, out))))
                if out
                else None
            )
            ctl = ctl or find(run_dir, f"{stem}.namd") or find(run_dir, f"{stem}.conf")
            if ctl is not None and ctl.format == "namdconf":
                step.files[FileRole.CONTROL] = ctl.ref()
                claimed.add(ctl.rel)
            step.raw["namd_xsc_out"] = (
                posixpath.normpath(posixpath.join(run_dir, out + ".xsc")) if out else None
            )
            self._fill_step(step, nl, lf.rel, by_rel)
            steps.append(step)

        # continuation: the coordinates a run read are another run's output or restart
        for step in steps:
            ref = step.files.get(FileRole.INPUT_COORDS)
            if ref is None:
                continue
            producer = outputs.get(ref.path)
            if producer and producer != step.id:
                kind = "restart" if ".restart." in ref.path else "final output"
                step.input_coords = CoordSource(CoordSourceKind.STEP, ref=producer, path=ref.path,
                                                method=f"binary coordinates = {kind} of the run")  # fmt: skip

        # configs that never ran
        for c in configs:
            if c.rel in claimed or not c.result.runs:
                continue
            run_dir, stem = _dir(c.rel), _stem(c.rel)
            step = Step(id=posixpath.join(run_dir, stem) if run_dir else stem, name=stem, dir=run_dir,
                        engine=self.name, status="queued")  # fmt: skip
            if any(s.id == step.id for s in steps):
                continue
            step.files[FileRole.CONTROL] = c.ref()
            step.settings = settings_from_config(c.result, c.rel)
            steps.append(step)
            findings.append(
                finding(
                    "STEP001",
                    f"{c.rel} has no NAMD log; treated as not run",
                    subject=step.id,
                    paths=[c.rel],
                )
            )

        for f in sorted(parsed, key=lambda f: f.rel):
            if f.kind is not FileKind.TRAJECTORY or f.rel in claimed:
                continue
            run_dir, stem = _dir(f.rel), _stem(f.rel)
            sid = posixpath.join(run_dir, stem) if run_dir else stem
            if any(s.id == sid for s in steps):
                sid = f.rel
            step = Step(id=sid, name=stem, dir=run_dir, engine=self.name, status="orphan")
            step.files[FileRole.TRAJECTORY] = f.ref()
            step.trajectory = self.trajectory_info(f, None)
            steps.append(step)
            findings.append(
                finding(
                    "STEP002",
                    f"{f.rel} has no NAMD log; its settings cannot be mined",
                    subject=sid,
                    paths=[f.rel],
                )
            )
        return steps, findings

    def _fill_step(
        self, step: Step, nl: log.NamdLog, rel: str, by_rel: dict[str, ParsedFile]
    ) -> None:
        step.settings = settings_from_log(nl, rel)
        rt = StepRuntime()
        rt.program = stated("NAMD", rel, "Info: NAMD")
        if nl.version:
            rt.version = stated(nl.version, rel, "Info: NAMD <version>")
        if nl.platform:
            rt.variant = stated(
                nl.platform + (" (GPU)" if nl.gpu else ""), rel, "Info: NAMD ... for <platform>"
            )
        if nl.natom is not None:
            rt.natom = stated(nl.natom, rel, "Info: N ATOMS")
        rt.completed = Mined(
            value=nl.completed,
            sources=(Source(rel, "FATAL ERROR / WallClock"),),
            method="stated" if nl.completed else "fatal error or no WallClock line",
            note="; ".join(nl.fatal[:2])
            or (None if nl.completed else "no WallClock line: stopped or still running"),
        )
        if nl.wallclock_s is not None:
            rt.wall_time_s = stated(nl.wallclock_s, rel, "WallClock", "s")
        dt = step.settings.dt_ps.value if step.settings.dt_ps.is_known else None
        first_ts = int(float(nl.one("FIRST TIMESTEP") or 0))
        if nl.last is not None:
            rt.nsteps_completed = Mined(value=nl.last.step - first_ts, sources=(Source(rel, "ENERGY rows"),),
                                        method="derived: last ENERGY step - FIRST TIMESTEP",
                                        confidence=Confidence.EXACT)  # fmt: skip
        if dt is not None and nl.dynamics_steps:
            rt.start_time_ps = Mined(value=clean_float(first_ts * dt), unit="ps",
                                     sources=(Source(rel, "Info: FIRST TIMESTEP x TIMESTEP"),),
                                     method="derived: FIRST TIMESTEP x TIMESTEP", confidence=Confidence.EXACT)  # fmt: skip
            if nl.last is not None:
                rt.end_time_ps = Mined(value=clean_float(nl.last.step * dt), unit="ps",
                                       sources=(Source(rel, "last ENERGY step x TIMESTEP"),),
                                       method="derived: last ENERGY step x TIMESTEP", confidence=Confidence.EXACT)  # fmt: skip
        step.runtime = rt
        traj = step.files.get(FileRole.TRAJECTORY)
        if traj is not None and traj.path in by_rel:
            step.trajectory = self.trajectory_info(by_rel[traj.path], None)
        step.raw["namd_parameters"] = nl.all("PARAMETERS")
        step.raw["namd_cell"] = [nl.one(f"PERIODIC CELL BASIS {i}") for i in (1, 2, 3)]

    # ---------------------------------------------------------------- queries
    def trajectory_info(self, traj: ParsedFile, natom: int | None) -> TrajectoryInfo:
        info = TrajectoryInfo(first_frame_at_start=False)
        info.format = Mined(value=traj.format, sources=(Source(traj.rel, "content"),))
        d = traj.result
        if not isinstance(d, files.Dcd):
            return info
        info.n_frames = stated(d.n_frames, traj.rel, "file size / frame size")
        info.natom = stated(d.natom, traj.rel, "header natoms")
        if d.timestep_fs and d.n_frames:
            info.first_time_ps = stated(d.time_ps(0), traj.rel, "ISTART x DELTA", "ps")
            info.last_time_ps = stated(
                d.time_ps(d.n_frames - 1), traj.rel, "last frame step x DELTA", "ps"
            )
            info.frame_interval_ps = Mined(value=clean_float(d.step_interval * d.timestep_fs / 1000.0), unit="ps",
                                           sources=(Source(traj.rel, "NSAVC x DELTA"),), method="derived: NSAVC x DELTA",
                                           confidence=Confidence.EXACT)  # fmt: skip
        return info

    def system_info(self, topology: ParsedFile) -> SystemInfo | None:
        r = topology.result
        if isinstance(r, files.Psf):
            info = SystemInfo(source=topology.rel, natom=r.natom, nres=r.nres,
                              residue_counts=dict(r.residue_counts), residue_atoms=dict(r.residue_atoms),
                              residue_atom_names={k: list(v) for k, v in r.residue_atom_names.items()})  # fmt: skip
        elif isinstance(r, amber_prmtop.Prmtop):
            atoms, names = r.residue_templates()
            info = SystemInfo(source=topology.rel, natom=r.natom, nres=r.nres, residue_counts=r.residue_counts(),
                              residue_atoms=atoms, residue_atom_names=names, n_extra_points=r.numextra,
                              periodic_box_type=r.ifbox)  # fmt: skip
        else:
            return None
        cell = self._cell()
        if cell is not None:
            info.box_lengths, info.box_angles = lengths_angles(cell)
            info.box_complete = True
            info.periodic_box_type = 1
        return info

    def _cell(self):
        for pf in self._by_rel.values():
            if pf.format == "namdlog" and pf.ok:
                vecs = [pf.result.one(f"PERIODIC CELL BASIS {i}") for i in (1, 2, 3)]
                if all(vecs):
                    try:
                        return tuple(tuple(float(x) for x in v.split()[:3]) for v in vecs)
                    except ValueError:
                        return None
        return None

    def coords_natom_time(self, coords: ParsedFile) -> tuple[int | None, float | None]:
        r = coords.result
        if isinstance(r, (files.NamdBinary, files.Pdb)):
            return r.natom, None
        if isinstance(r, files.Xsc):
            return None, None
        return None, None

    def coords_box(self, coords: ParsedFile):
        r = coords.result
        if isinstance(r, files.Xsc) and r.vectors_A:
            return lengths_angles(r.vectors_A)
        return None

    def topology_evidence(self, topology: ParsedFile) -> BuildEvidence | None:
        """Force fields from the CHARMM parameter files the runs loaded (from the logs)."""
        params: list[str] = []
        for pf in self._by_rel.values():
            if pf.format == "namdlog" and pf.ok:
                for p in pf.result.all("PARAMETERS"):
                    if p not in params:
                        params.append(p)
        if not params or not isinstance(topology.result, files.Psf):
            return None
        classes = residue_classes(topology.result.residue_counts)
        labels: list[str] = []
        waters: list[str] = []
        for p in params:
            stem = posixpath.basename(p).lower()
            for key, label, needs in CHARMM_PARAMS:
                if (
                    stem.startswith(key)
                    and (needs is None or needs in classes)
                    and label not in labels
                ):
                    labels.append(label)
                    break
            if (
                "water" in stem
                and "TIP3" in topology.result.residue_counts
                and "TIP3P" not in waters
            ):
                waters.append("TIP3P")  # CHARMM's modified TIP3P (with LJ on hydrogens)
        return BuildEvidence(source=topology.rel, force_fields=labels, water_models=waters,
                             outputs=[posixpath.basename(topology.rel)])  # fmt: skip

    def mddb_trajectory_format(self, fmt: str) -> str | None:
        return {"dcd": "dcd"}.get(fmt)

    def mddb_topology_format(self, fmt: str) -> str | None:
        return {"psf": "psf", "prmtop": "prmtop"}.get(fmt)


def residue_classes(counts: dict[str, int]) -> set[str]:
    classes = set()
    for name in counts:
        n = name.upper()
        if n in PROTEIN:
            classes.add("protein")
        elif n in NUCLEIC:
            classes.add("nucleic")
        elif n in LIPID:
            classes.add("lipid")
        elif n in WATER_IONS:
            continue
        elif n.startswith(CARB_PREFIX):
            classes.add("carb")
        else:
            classes.add("other")
    return classes
