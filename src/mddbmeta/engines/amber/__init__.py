"""AMBER (sander / pmemd) engine adapter."""

from __future__ import annotations

import posixpath
from pathlib import Path
from typing import Any

from mddbmeta.engines.amber import leaplog, mdin, mdout, netcdf, prmtop, rst7
from mddbmeta.engines.amber.semantics import CntrlResolver, settings_from_cntrl
from mddbmeta.engines.base import BuildEvidence, EngineAdapter, ParsedFile, SniffResult
from mddbmeta.errors import MissingDependencyError
from mddbmeta.findings import finding
from mddbmeta.model import (
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
    derived,
    stated,
)
from mddbmeta.units import clean_float

RESTART_EXTS = {".rst", ".rst7", ".restrt", ".ncrst", ".restart"}
NETCDF_TRAJ_EXTS = {".nc", ".ncdf", ".netcdf", ".cdf"}

# mdout File Assignments tag -> step file role
ASSIGNMENT_ROLES = {
    "MDIN": FileRole.CONTROL,
    "INPCRD": FileRole.INPUT_COORDS,
    "PARM": FileRole.TOPOLOGY,
    "RESTRT": FileRole.OUTPUT_RESTART,
    "MDCRD": FileRole.TRAJECTORY,
    "REFC": FileRole.REF_COORDS,
}
ROLE_KINDS = {
    FileRole.CONTROL: {FileKind.CONTROL},
    FileRole.INPUT_COORDS: {FileKind.COORDS},
    FileRole.TOPOLOGY: {FileKind.TOPOLOGY},
    FileRole.OUTPUT_RESTART: {FileKind.COORDS},
    FileRole.TRAJECTORY: {FileKind.TRAJECTORY},
    FileRole.REF_COORDS: {FileKind.COORDS},
}


def _stem(rel: str) -> str:
    base = posixpath.basename(rel)
    return base.rsplit(".", 1)[0] if "." in base else base


def _dir(rel: str) -> str:
    return posixpath.dirname(rel)


class AmberAdapter(EngineAdapter):
    name = "amber"
    program = "AMBER"

    # ------------------------------------------------------------------ sniff
    def sniff(self, path: Path, head: bytes) -> SniffResult | None:
        ext = path.suffix.lower()
        if netcdf.is_netcdf_bytes(head):
            kind = FileKind.COORDS if ext in RESTART_EXTS else FileKind.TRAJECTORY
            fmt = "netcdf_restart" if kind is FileKind.COORDS else "netcdf_traj"
            return SniffResult(kind, fmt, 0.7)
        if b"\x00" in head[:1024]:
            return None
        if prmtop.looks_like_prmtop(head):
            return SniffResult(FileKind.TOPOLOGY, "prmtop", 0.95)
        if mdout.looks_like_mdout(head):
            return SniffResult(FileKind.LOG, "mdout", 0.95)
        if mdin.looks_like_mdin(head):
            return SniffResult(FileKind.CONTROL, "mdin", 0.9)
        if leaplog.looks_like_leap(head, path.name):
            return SniffResult(FileKind.BUILD_LOG, "leap", 0.8)
        shape = rst7.sniff_ascii_coords(head)
        if shape == "rst7":
            return SniffResult(FileKind.COORDS, "rst7", 0.8)
        if shape == "mdcrd":
            return SniffResult(FileKind.TRAJECTORY, "mdcrd", 0.7)
        return None

    # ------------------------------------------------------------------ parse
    def parse(self, path: Path, rel: str, sniffed: SniffResult) -> ParsedFile:
        pf = ParsedFile(path=str(path), rel=rel, kind=sniffed.kind, format=sniffed.format)
        readers = {
            "prmtop": prmtop.read_prmtop,
            "mdout": mdout.read_mdout,
            "mdin": mdin.read_mdin,
            "leap": leaplog.read_leap,
            "rst7": rst7.read_rst7,
            "mdcrd": rst7.read_ascii_traj,
            "netcdf_traj": netcdf.read_netcdf,
            "netcdf_restart": netcdf.read_netcdf,
        }
        try:
            pf.result = readers[sniffed.format](path)
        except MissingDependencyError as exc:
            pf.errors.append(FileLoadError(sniffed.kind.value, rel, "unsupported", str(exc)))
            return pf
        except Exception as exc:  # parsers must not take discovery down
            pf.errors.append(
                FileLoadError(sniffed.kind.value, rel, classify_exception(exc), str(exc))
            )
            return pf
        # A NetCDF file's own Conventions attribute beats its extension.
        if sniffed.format.startswith("netcdf") and isinstance(pf.result, netcdf.NetcdfInfo):
            if pf.result.is_restart and pf.kind is FileKind.TRAJECTORY:
                pf.kind, pf.format = FileKind.COORDS, "netcdf_restart"
            elif not pf.result.is_restart and pf.kind is FileKind.COORDS:
                pf.kind, pf.format = FileKind.TRAJECTORY, "netcdf_traj"
        if isinstance(pf.result, mdin.Mdin) and pf.result.n_cntrl > 1:
            pf.findings.append(
                finding(
                    "ENG002",
                    f"{rel} holds {pf.result.n_cntrl} &cntrl namelists (groupfile/multi-run input); "
                    "only the first is used",
                    subject=rel,
                    paths=[rel],
                )
            )
        return pf

    # ------------------------------------------------------------ build steps
    def build_steps(self, root: Path, files: list[ParsedFile]) -> tuple[list[Step], list[Finding]]:
        findings: list[Finding] = []
        by_rel = {f.rel: f for f in files}
        by_dir_base: dict[tuple[str, str], ParsedFile] = {}
        by_dir_stem: dict[tuple[str, str], list[ParsedFile]] = {}
        for f in files:
            by_dir_base[(_dir(f.rel), posixpath.basename(f.rel))] = f
            by_dir_stem.setdefault((_dir(f.rel), _stem(f.rel)), []).append(f)

        claimed: set[str] = set()
        steps: list[Step] = []

        def resolve(name: str, run_dir: str) -> ParsedFile | None:
            cands: list[str] = []
            p = Path(name)
            if p.is_absolute():
                try:
                    cands.append(Path(p).resolve().relative_to(root).as_posix())
                except ValueError:
                    pass
            else:
                cands.append(posixpath.normpath(posixpath.join(run_dir, name)))
                cands.append(posixpath.normpath(name))
            for c in cands:
                if c in by_rel:
                    return by_rel[c]
            # runs are often executed elsewhere and copied: fall back to the basename
            return by_dir_base.get((run_dir, posixpath.basename(name)))

        def by_stem(run_dir: str, stem: str, kinds: set[FileKind]) -> ParsedFile | None:
            for f in by_dir_stem.get((run_dir, stem), []):
                if f.kind in kinds and f.rel not in claimed:
                    return f
            return None

        logs = sorted((f for f in files if f.kind is FileKind.LOG), key=lambda f: f.rel)
        for log in logs:
            run_dir, stem = _dir(log.rel), _stem(log.rel)
            step = Step(
                id=posixpath.join(run_dir, stem) if run_dir else stem,
                name=stem,
                dir=run_dir,
                engine=self.name,
            )
            step.files[FileRole.LOG] = log.ref()
            claimed.add(log.rel)
            step.load_errors.extend(log.errors)
            out: mdout.Mdout | None = log.result if log.ok else None
            assignments = out.assignments if out else {}
            unresolved: dict[str, str] = {}
            for tag, role in ASSIGNMENT_ROLES.items():
                name = assignments.get(tag)
                if not name:
                    continue
                target = resolve(name, run_dir)
                if target is not None and target.kind in ROLE_KINDS[role]:
                    step.files[role] = target.ref()
                    if role in (FileRole.CONTROL, FileRole.OUTPUT_RESTART, FileRole.TRAJECTORY):
                        claimed.add(target.rel)
                else:
                    unresolved[tag] = name
            # Stem fallbacks for files the log did not (resolvably) name.
            for role, kinds in (
                (FileRole.CONTROL, {FileKind.CONTROL}),
                (FileRole.TRAJECTORY, {FileKind.TRAJECTORY}),
                (FileRole.OUTPUT_RESTART, {FileKind.COORDS}),
            ):
                if role in step.files:
                    continue
                f = by_stem(run_dir, stem, kinds)
                if f is not None:
                    step.files[role] = f.ref()
                    claimed.add(f.rel)
                    step.raw.setdefault("stem_matched", []).append(role.value)
            step.raw["unresolved_assignments"] = unresolved
            self._fill_step(step, by_rel, out, log)
            steps.append(step)

        # Control files that no log claimed: prepared but not run (or log missing).
        for f in sorted(files, key=lambda f: f.rel):
            if f.kind is not FileKind.CONTROL or f.rel in claimed:
                continue
            run_dir, stem = _dir(f.rel), _stem(f.rel)
            step = Step(
                id=posixpath.join(run_dir, stem) if run_dir else stem,
                name=stem,
                dir=run_dir,
                engine=self.name,
                status="queued",
            )
            step.files[FileRole.CONTROL] = f.ref()
            claimed.add(f.rel)
            self._fill_step(step, by_rel, None, None)
            steps.append(step)
            findings.append(
                finding(
                    "STEP001",
                    f"{f.rel} has no engine log; treated as not run",
                    subject=step.id,
                    paths=[f.rel],
                )
            )

        # Trajectories that no log claimed.
        for f in sorted(files, key=lambda f: f.rel):
            if f.kind is not FileKind.TRAJECTORY or f.rel in claimed:
                continue
            run_dir, stem = _dir(f.rel), _stem(f.rel)
            step_id = posixpath.join(run_dir, stem) if run_dir else stem
            if any(s.id == step_id for s in steps):
                step_id = f.rel
            step = Step(id=step_id, name=stem, dir=run_dir, engine=self.name, status="orphan")
            step.files[FileRole.TRAJECTORY] = f.ref()
            claimed.add(f.rel)
            step.load_errors.extend(f.errors)
            self._fill_trajectory(step, by_rel, None)
            steps.append(step)
            findings.append(
                finding(
                    "STEP002",
                    f"{f.rel} has no engine log; its settings cannot be mined",
                    subject=step.id,
                    paths=[f.rel],
                )
            )

        # Ids must be unique even when two logs share a stem (prod.out + prod.log).
        seen: dict[str, int] = {}
        for s in steps:
            if s.id in seen:
                seen[s.id] += 1
                s.id = f"{s.id}~{seen[s.id]}"
            else:
                seen[s.id] = 0
        return steps, findings

    def _fill_step(
        self,
        step: Step,
        by_rel: dict[str, ParsedFile],
        out: mdout.Mdout | None,
        log: ParsedFile | None,
    ) -> None:
        control_ref = step.files.get(FileRole.CONTROL)
        control = by_rel.get(control_ref.path) if control_ref else None
        if control is not None:
            step.load_errors.extend(e for e in control.errors if e not in step.load_errors)
        stated_values: dict[str, Any] = {}
        stated_path: str | None = None
        prefix = "&cntrl"
        wt: list[dict[str, Any]] | None = None
        if control is not None and control.ok:
            parsed: mdin.Mdin = control.result
            stated_values, stated_path, wt = parsed.cntrl, control.rel, parsed.all("wt")
            step.raw["mdin_title"] = parsed.title
        elif out is not None and out.input_echo and log is not None:
            try:
                echo = mdin.parse_mdin_text(out.input_echo, path=log.rel)
                stated_values, stated_path, wt = echo.cntrl, log.rel, echo.all("wt")
                prefix = "input echo &cntrl"
            except mdin.MdinFormatError:
                pass
        resolver = CntrlResolver(
            stated_values,
            stated_path,
            stated_locator_prefix=prefix,
            reported=out.control if out else None,
            reported_path=log.rel if (out and log) else None,
        )
        if stated_path is not None or (out is not None and out.control):
            step.settings = settings_from_cntrl(resolver, wt)
            step.raw["cntrl"] = dict(stated_values)
        if out is not None and log is not None:
            step.runtime = self._runtime(step, out, log.rel, resolver)
        self._fill_trajectory(step, by_rel, out)

    def _runtime(self, step: Step, out: mdout.Mdout, rel: str, res: CntrlResolver) -> StepRuntime:
        rt = StepRuntime()
        if out.program:
            rt.program = stated(out.program, rel, "header")
        if out.version:
            rt.version = stated(out.version, rel, "header")
        if out.variant:
            rt.variant = stated(out.variant, rel, "Executable path")
        if out.natom is not None:
            rt.natom = stated(out.natom, rel, "RESOURCE USE:NATOM")
        rt.completed = Mined(
            value=out.completed,
            sources=(Source(rel, "TIMINGS"),),
            method="stated" if out.completed else "absent_final_timings",
            note=None
            if out.completed
            else "no final timing block: the run stopped or is still running",
        )
        if out.wall_time_s is not None:
            rt.wall_time_s = stated(out.wall_time_s, rel, "Total wall time", unit="s")
        if step.settings.minimization.value:
            if out.last_min_step is not None:
                rt.nsteps_completed = stated(out.last_min_step, rel, "last minimization record")
            return rt
        # start time: read from the input coordinates on a restart, else &cntrl t
        if out.begin_time_ps is not None:
            rt.start_time_ps = stated(
                out.begin_time_ps, rel, "begin time read from input coords", "ps"
            )
        else:
            restart = step.settings.restart.value
            if restart is False:
                t = res.get("t", unit="ps")
                rt.start_time_ps = t
        if out.last_record is not None:
            rt.nsteps_completed = stated(out.last_record.nstep, rel, "last energy record")
        nstlim, dt = step.settings.nsteps, step.settings.dt_ps
        if out.completed and rt.start_time_ps.is_known and nstlim.is_known and dt.is_known:
            end = clean_float(rt.start_time_ps.value + nstlim.value * dt.value)
            rt.end_time_ps = derived(
                end, "start + nstlim*dt", (rt.start_time_ps, nstlim, dt), unit="ps"
            )
            if out.last_record is not None and out.last_record.nstep < nstlim.value:
                rt.nsteps_completed = derived(
                    nstlim.value, "completed run: nstlim", (nstlim, rt.completed)
                )
        elif out.last_record is not None and out.last_record.time_ps is not None:
            rt.end_time_ps = stated(
                out.last_record.time_ps, rel, "last energy record TIME(PS)", "ps"
            )
        return rt

    def _fill_trajectory(
        self, step: Step, by_rel: dict[str, ParsedFile], out: mdout.Mdout | None
    ) -> None:
        ref = step.files.get(FileRole.TRAJECTORY)
        if ref is None:
            return
        traj = by_rel.get(ref.path)
        if traj is None:
            return
        natom = out.natom if out and out.natom else None
        step.trajectory = self.trajectory_info(traj, natom)
        step.load_errors.extend(e for e in traj.errors if e not in step.load_errors)

    # ---------------------------------------------------------------- queries
    def trajectory_info(self, traj: ParsedFile, natom: int | None) -> TrajectoryInfo:
        info = TrajectoryInfo()
        info.format = Mined(value=traj.format, sources=(Source(traj.rel, "content"),))
        r = traj.result
        if isinstance(r, netcdf.NetcdfInfo):
            if r.n_frames is not None:
                info.n_frames = stated(r.n_frames, traj.rel, "dimension frame")
            if r.natom is not None:
                info.natom = stated(r.natom, traj.rel, "dimension atom")
            # writers that don't track time (e.g. some converters) store 0 for every frame
            usable_times = r.first_time_ps is not None and not (
                (r.n_frames or 0) > 1 and (r.median_dt_ps is None or r.median_dt_ps <= 0)
            )
            if usable_times:
                info.first_time_ps = stated(r.first_time_ps, traj.rel, "time[0]", "ps")
                info.last_time_ps = stated(r.last_time_ps, traj.rel, "time[-1]", "ps")
            if usable_times and r.median_dt_ps is not None and r.median_dt_ps > 0:
                info.frame_interval_ps = Mined(
                    value=clean_float(r.median_dt_ps),
                    unit="ps",
                    sources=(Source(traj.rel, "median(diff(time))"),),
                    method="derived:median time delta",
                    confidence=Confidence.EXACT,
                )
            if r.program:
                info.program = stated(r.program, traj.rel, "attribute program")
            if r.program_version:
                info.program_version = stated(
                    r.program_version, traj.rel, "attribute programVersion"
                )
        elif isinstance(r, rst7.AsciiTraj) and natom:
            frames, _box = rst7.ascii_traj_frames(r, natom)
            if frames is not None:
                info.n_frames = Mined(
                    value=frames,
                    sources=(Source(traj.rel, "line count"),),
                    method="derived:line count / lines per frame",
                    confidence=Confidence.DERIVED,
                )
                info.natom = Mined(
                    value=natom, method="assumed from log/topology", confidence=Confidence.DERIVED
                )
        return info

    def system_info(self, topology: ParsedFile) -> SystemInfo | None:
        top = topology.result
        if not isinstance(top, prmtop.Prmtop):
            return None
        atoms, names = top.residue_templates()
        info = SystemInfo(
            source=topology.rel,
            natom=top.natom,
            nres=top.nres,
            title=top.title,
            residue_counts=top.residue_counts(),
            residue_atoms=atoms,
            residue_atom_names=names,
            periodic_box_type=top.ifbox,
            n_extra_points=top.numextra,
        )
        if top.box is not None:
            beta, a, b, c = top.box
            info.box_lengths = (a, b, c)
            info.box_angles = (beta, beta, beta) if top.ifbox == 2 else (90.0, beta, 90.0)
        if top.is_chamber:
            info.title = info.title or "CHAMBER topology"
        return info

    def coords_natom_time(self, coords: ParsedFile) -> tuple[int | None, float | None]:
        r = coords.result
        if isinstance(r, rst7.Rst7):
            return r.natom, r.time_ps
        if isinstance(r, netcdf.NetcdfInfo):
            return r.natom, r.first_time_ps
        return None, None

    def coords_box(
        self, coords: ParsedFile
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
        r = coords.result
        if isinstance(r, rst7.Rst7) and r.box:
            return tuple(r.box[:3]), tuple(r.box[3:])  # type: ignore[return-value]
        if isinstance(r, netcdf.NetcdfInfo) and r.cell_lengths and r.cell_angles:
            return r.cell_lengths, r.cell_angles
        return None

    def build_evidence(self, build_log: ParsedFile) -> BuildEvidence | None:
        r = build_log.result
        if not isinstance(r, leaplog.LeapEvidence):
            return None
        return BuildEvidence(
            source=build_log.rel,
            force_fields=list(r.force_fields),
            water_models=list(r.water_models),
            outputs=list(r.saved_topologies),
        )

    def is_chamber(self, topology: ParsedFile) -> bool:
        return isinstance(topology.result, prmtop.Prmtop) and topology.result.is_chamber

    def mddb_trajectory_format(self, fmt: str) -> str | None:
        return {"netcdf_traj": "nc", "mdcrd": "crd", "rst7": "rst7"}.get(fmt)

    def mddb_topology_format(self, fmt: str) -> str | None:
        return "prmtop" if fmt == "prmtop" else None
