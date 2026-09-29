"""GROMACS (gmx mdrun) engine adapter.

One Step per mdrun log. The log's command line names the files the run
used (``-deffnm``, ``-s``, ``-x``, ``-o``, ``-c``, ``-cpo``; ``-noappend``
adds ``.partNNNN`` to outputs). Settings come from the log's resolved
"Input Parameters", plus the matching .mdp for grompp-only options.

GROMACS does not record which coordinates grompp read, so continuations are
recognised from evidence the engine does record:

* ``Reading checkpoint file X`` (mdrun -cpi) with the checkpoint's time;
* energy fingerprints: a run continued with ``grompp -t prev.cpt`` starts in
  exactly the state the previous run ended in, so its step-0 energy terms
  equal the previous run's last energy record (announced as LIN003).
"""

from __future__ import annotations

import posixpath
import re
from pathlib import Path

from mddbmeta.engines.base import BuildEvidence, EngineAdapter, ParsedFile, SniffResult
from mddbmeta.engines.gromacs import binary, log, mdp, structure
from mddbmeta.engines.gromacs.semantics import GmxResolver, settings_from_params
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

_PART = re.compile(r"^(?P<base>.*)\.part(?P<n>\d{4})$")
_VERSION_NUM = re.compile(r"^(\d{4}(?:\.\d+)?|\d+\.\d+(?:\.\d+)?)")
MIN_COMMON_TERMS = 3
TIME_TOL_PS = 1e-3

# mdrun option -> (extension, default base name without -deffnm, gets .partNNNN with -noappend)
OUTPUTS = {
    "-s": (".tpr", "topol", False),
    "-x": (".xtc", "traj_comp", True),
    "-o": (".trr", "traj", True),
    "-c": (".gro", "confout", True),
    "-cpo": (".cpt", "state", False),
}


def _dir(rel: str) -> str:
    return posixpath.dirname(rel)


def _stem(rel: str) -> str:
    base = posixpath.basename(rel)
    return base.rsplit(".", 1)[0] if "." in base else base


def _last_digit_unit(text: str) -> float:
    """Value of one unit in the last printed digit, e.g. '4.57306e+03' -> 0.01."""
    mantissa, _, exp = text.lower().partition("e")
    decimals = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
    return 10.0 ** (int(exp or 0) - decimals)


def _same_printed(a: str, b: str) -> bool:
    """Equal as printed, allowing one unit of rounding in the last digit (pair-list noise)."""
    if a == b:
        return True
    try:
        fa, fb = float(a), float(b)
    except ValueError:
        return False
    return abs(fa - fb) <= 1.01 * max(_last_digit_unit(a), _last_digit_unit(b))


def clean_version(text: str | None) -> str | None:
    if not text:
        return None
    m = _VERSION_NUM.match(text.strip())
    return m.group(1) if m else text.strip()


class GromacsAdapter(EngineAdapter):
    name = "gromacs"
    program = "GROMACS"

    def __init__(self) -> None:
        self._gros: list[ParsedFile] = []
        self._tops: list[ParsedFile] = []

    # ------------------------------------------------------------------ sniff
    def sniff(self, path: Path, head: bytes) -> SniffResult | None:
        magic = binary.first_int(head)
        if magic == binary.XTC_MAGIC:
            return SniffResult(FileKind.TRAJECTORY, "xtc", 0.95)
        if magic == binary.TRR_MAGIC:
            return SniffResult(FileKind.TRAJECTORY, "trr", 0.95)
        if magic == binary.CPT_MAGIC:
            return SniffResult(FileKind.COORDS, "cpt", 0.9)
        if binary.looks_like_tpr(head):
            return SniffResult(FileKind.TOPOLOGY, "tpr", 0.95)
        if b"\x00" in head[:1024]:
            return None
        if log.looks_like_log(head):
            return SniffResult(FileKind.LOG, "gmxlog", 0.95)
        if structure.looks_like_gro(head):
            return SniffResult(FileKind.COORDS, "gro", 0.85)
        if structure.looks_like_top(head, path.name):
            return SniffResult(FileKind.TOPOLOGY, "top", 0.85)
        if mdp.looks_like_mdp(head, path.name):
            return SniffResult(FileKind.CONTROL, "mdp", 0.85)
        return None

    # ------------------------------------------------------------------ parse
    def parse(self, path: Path, rel: str, sniffed: SniffResult) -> ParsedFile:
        pf = ParsedFile(path=str(path), rel=rel, kind=sniffed.kind, format=sniffed.format)
        readers = {
            "xtc": binary.read_xtc,
            "trr": binary.read_trr,
            "tpr": binary.read_tpr_header,
            "gmxlog": log.read_log,
            "gro": structure.read_gro,
            "top": structure.read_top,
            "mdp": mdp.read_mdp,
            "cpt": lambda p: None,
        }
        try:
            pf.result = readers[sniffed.format](path)
        except Exception as exc:
            pf.errors.append(
                FileLoadError(sniffed.kind.value, rel, classify_exception(exc), str(exc))
            )
            return pf
        if sniffed.format == "cpt":
            pf.result = "checkpoint"
        if isinstance(pf.result, binary.Frames) and pf.result.truncated:
            pf.findings.append(
                finding("LOAD001", f"{rel}: last frame is truncated", subject=rel, paths=[rel])
            )
        if isinstance(pf.result, log.GmxLog) and not pf.result.sessions:
            pf.result = None
            pf.errors.append(FileLoadError("log", rel, "malformed", "no mdrun session found"))
        return pf

    # ------------------------------------------------------------ build steps
    def build_steps(self, root: Path, files: list[ParsedFile]) -> tuple[list[Step], list[Finding]]:
        findings: list[Finding] = []
        by_rel = {f.rel: f for f in files}
        self._gros = [f for f in files if f.format == "gro" and f.ok]
        self._tops = [f for f in files if f.format == "top" and f.ok]
        claimed: set[str] = set()
        steps: list[Step] = []

        def find(run_dir: str, name: str | None) -> ParsedFile | None:
            if not name:
                return None
            cands = [posixpath.normpath(posixpath.join(run_dir, name)), posixpath.normpath(name)]
            for c in cands:
                if c in by_rel:
                    return by_rel[c]
            base = posixpath.basename(name)
            return by_rel.get(posixpath.join(run_dir, base) if run_dir else base)

        for lf in sorted((f for f in files if f.kind is FileKind.LOG), key=lambda f: f.rel):
            run_dir, stem = _dir(lf.rel), _stem(lf.rel)
            step = Step(
                id=posixpath.join(run_dir, stem) if run_dir else stem,
                name=stem,
                dir=run_dir,
                engine=self.name,
            )
            step.files[FileRole.LOG] = lf.ref()
            claimed.add(lf.rel)
            step.load_errors.extend(lf.errors)
            glog: log.GmxLog | None = lf.result if lf.ok else None
            if glog is not None:
                names = self._output_names(glog, stem)
                step.raw["gmx_names"] = names
                for role, key, kinds in (
                    (FileRole.TOPOLOGY, "-s", {FileKind.TOPOLOGY}),
                    (FileRole.TRAJECTORY, "-x", {FileKind.TRAJECTORY}),
                    (FileRole.OUTPUT_RESTART, "-c", {FileKind.COORDS}),
                ):
                    f = find(run_dir, names.get(key))
                    if f is not None and f.kind in kinds:
                        step.files[role] = f.ref()
                        claimed.add(f.rel)
                if FileRole.TRAJECTORY not in step.files:
                    f = find(run_dir, names.get("-o"))
                    if f is not None and f.kind is FileKind.TRAJECTORY:
                        step.files[FileRole.TRAJECTORY] = f.ref()
                        claimed.add(f.rel)
                if FileRole.OUTPUT_RESTART not in step.files:
                    f = find(run_dir, names.get("-cpo"))
                    if f is not None and f.kind is FileKind.COORDS:
                        step.files[FileRole.OUTPUT_RESTART] = f.ref()
                # the .mdp that produced this run, when it shares the tpr/deffnm stem
                tpr_stem = _stem(names.get("-s", "")) if names.get("-s") else None
                for cand in (tpr_stem, names.get("deffnm"), _PART.sub(r"\g<base>", stem)):
                    f = find(run_dir, f"{cand}.mdp") if cand else None
                    if f is not None and f.kind is FileKind.CONTROL:
                        step.files[FileRole.CONTROL] = f.ref()
                        claimed.add(f.rel)
                        break
            self._fill_step(step, by_rel, glog, lf.rel)
            steps.append(step)

        # prepared but not run: a tpr no log used
        for f in sorted(files, key=lambda f: f.rel):
            if f.format != "tpr" or f.rel in claimed:
                continue
            run_dir, stem = _dir(f.rel), _stem(f.rel)
            step = Step(
                id=posixpath.join(run_dir, stem) if run_dir else stem,
                name=stem,
                dir=run_dir,
                engine=self.name,
                status="queued",
            )
            if any(s.id == step.id for s in steps):
                continue
            step.files[FileRole.TOPOLOGY] = f.ref()
            steps.append(step)
            findings.append(
                finding(
                    "STEP001",
                    f"{f.rel} has no mdrun log; treated as not run",
                    subject=step.id,
                    paths=[f.rel],
                )
            )

        for f in sorted(files, key=lambda f: f.rel):
            if f.kind is not FileKind.TRAJECTORY or f.rel in claimed:
                continue
            run_dir, stem = _dir(f.rel), _stem(f.rel)
            step_id = posixpath.join(run_dir, stem) if run_dir else stem
            if any(s.id == step_id for s in steps):
                step_id = f.rel
            step = Step(id=step_id, name=stem, dir=run_dir, engine=self.name, status="orphan")
            step.files[FileRole.TRAJECTORY] = f.ref()
            step.load_errors.extend(f.errors)
            step.trajectory = self.trajectory_info(f, None)
            steps.append(step)
            findings.append(
                finding(
                    "STEP002",
                    f"{f.rel} has no mdrun log; its settings cannot be mined",
                    subject=step.id,
                    paths=[f.rel],
                )
            )

        findings.extend(self._link(steps, by_rel))
        return steps, findings

    def _output_names(self, glog: log.GmxLog, log_stem: str) -> dict[str, str]:
        """File names mdrun used, from the first session's command line."""
        s = glog.first
        deffnm = s.option("-deffnm")
        part = _PART.match(log_stem)
        suffix = f".part{part.group('n')}" if part and s.has_flag("-noappend") else ""
        names: dict[str, str] = {}
        if deffnm:
            names["deffnm"] = deffnm
        for flag, (ext, default, parts) in OUTPUTS.items():
            explicit = s.option(flag)
            if explicit:
                if parts and suffix:
                    stem, dot, e = explicit.rpartition(".")
                    names[flag] = f"{stem}{suffix}.{e}" if dot else explicit + suffix
                else:
                    names[flag] = explicit
            else:
                base = deffnm or default
                names[flag] = f"{base}{suffix if parts else ''}{ext}"
        cpi = s.option("-cpi")
        if cpi or s.has_flag("-cpi"):
            names["-cpi"] = cpi or f"{deffnm or 'state'}.cpt"
        return names

    def _fill_step(
        self, step: Step, by_rel: dict[str, ParsedFile], glog: log.GmxLog | None, log_rel: str
    ) -> None:
        ctl_ref = step.files.get(FileRole.CONTROL)
        ctl = by_rel.get(ctl_ref.path) if ctl_ref else None
        stated_values = ctl.result if ctl is not None and ctl.ok else None
        if glog is None and stated_values is None:
            return
        res = GmxResolver(
            glog.params if glog else None,
            log_rel if glog else None,
            stated_values,
            ctl.rel if ctl is not None else None,
        )
        step.settings = settings_from_params(res)
        if glog is not None:
            step.runtime = self._runtime(glog, log_rel, step, by_rel)
            step.raw["gmx_sessions"] = len(glog.sessions)
            # appended sessions don't reprint Input Parameters, so after `convert-tpr -extend`
            # the logged nsteps is stale; the step records say how far the run really went
            _first, last = glog.records()
            init = glog.params.get("init-step", 0) or 0
            if (
                last is not None
                and step.settings.nsteps.is_known
                and isinstance(step.settings.nsteps.value, int)
                and last.step > step.settings.nsteps.value + init
            ):
                step.settings.nsteps = Mined(
                    value=last.step - init,
                    sources=(Source(log_rel, "last Step record"),),
                    method="derived: last Step record (run extended beyond the logged nsteps)",
                    confidence=Confidence.EXACT,
                )
        traj_ref = step.files.get(FileRole.TRAJECTORY)
        if traj_ref is not None and traj_ref.path in by_rel:
            tf = by_rel[traj_ref.path]
            step.trajectory = self.trajectory_info(tf, None)
            groups = glog.params.get("compressed-x-grps") if glog else None
            if tf.format == "xtc" and groups not in (None, "", "System", "system"):
                step.trajectory.natom = step.trajectory.natom.with_value(
                    step.trajectory.natom.value, method="stated (output group subset)"
                )

    def _runtime(
        self, glog: log.GmxLog, rel: str, step: Step, by_rel: dict[str, ParsedFile]
    ) -> StepRuntime:
        rt = StepRuntime()
        first_s, last_s = glog.first, glog.last
        rt.program = stated("GROMACS", rel, "banner")
        version = last_s.version or first_s.version
        if version:
            rt.version = Mined(
                value=clean_version(version),
                sources=(Source(rel, "GROMACS version"),),
                note=None
                if clean_version(version) == version
                else f"full version string: {version}",
            )
        exe = last_s.executable or first_s.executable
        if exe:
            rt.variant = stated(posixpath.basename(exe), rel, "Executable")
        rt.completed = Mined(
            value=last_s.finished,
            sources=(Source(rel, "Finished mdrun"),),
            method="stated" if last_s.finished else "absent: Finished mdrun",
            note=None
            if last_s.finished
            else "no 'Finished mdrun' line: the run stopped or is still running",
        )
        first, last = glog.records()
        if first is not None and last is not None:
            if not glog.minimization:
                rt.start_time_ps = stated(first.time_ps, rel, "first Step/Time record", "ps")
                rt.end_time_ps = stated(last.time_ps, rel, "last Step/Time record", "ps")
            rt.nsteps_completed = Mined(
                value=last.step - first.step,
                sources=(Source(rel, "Step records"),),
                method="derived: last step - first step",
                confidence=Confidence.EXACT,
            )
        tpr_ref = step.files.get(FileRole.TOPOLOGY)
        if tpr_ref is not None and tpr_ref.path in by_rel:
            hdr = by_rel[tpr_ref.path].result
            if isinstance(hdr, binary.TprHeader) and hdr.natom:
                rt.natom = stated(hdr.natom, tpr_ref.path, "tpr header natoms")
        step.raw["gmx_first_energies"] = first.energies if first else {}
        step.raw["gmx_last_energies"] = last.energies if last else {}
        cp = next((s for s in glog.sessions if s.checkpoint_file), None)
        if cp is not None:
            step.raw["gmx_checkpoint_read"] = (
                cp.checkpoint_file,
                cp.checkpoint_time_ps,
                cp.appending,
            )
        return rt

    # ---------------------------------------------------------- continuations
    def _link(self, steps: list[Step], by_rel: dict[str, ParsedFile]) -> list[Finding]:
        findings: list[Finding] = []
        runs = [s for s in steps if s.status == "done"]
        # (a) mdrun -cpi without appending: a new log continuing another run's checkpoint
        for s in runs:
            info = s.raw.get("gmx_checkpoint_read")
            if not info or info[2]:  # appending continues the same log/step
                continue
            cpt_name, t = posixpath.basename(info[0]), info[1]
            cands = [
                p
                for p in runs
                if p is not s
                and p.dir == s.dir
                and posixpath.basename(p.raw.get("gmx_names", {}).get("-cpo", "")) == cpt_name
                and (
                    t is None
                    or (
                        p.runtime.end_time_ps.is_known
                        and abs(p.runtime.end_time_ps.value - t) <= TIME_TOL_PS
                    )
                )
            ]
            if len(cands) == 1:
                s.input_coords = CoordSource(
                    CoordSourceKind.STEP, ref=cands[0].id, path=posixpath.join(s.dir, info[0]) if s.dir else info[0],
                    method="checkpoint read by mdrun -cpi",
                )  # fmt: skip
        # (b) energy fingerprints for everything else
        for s in runs:
            if s.input_coords.kind is CoordSourceKind.STEP:
                continue
            start = s.raw.get("gmx_first_energies") or {}
            if not start:
                continue
            matches = []
            for p in runs:
                if p is s:
                    continue
                end = p.raw.get("gmx_last_energies") or {}
                common = [t for t in log.STATE_TERMS if t in start and t in end]
                if len(common) >= MIN_COMMON_TERMS and all(
                    _same_printed(start[t], end[t]) for t in common
                ):
                    matches.append(p)
            if len(matches) == 1:
                p = matches[0]
                s.input_coords = CoordSource(
                    CoordSourceKind.STEP,
                    ref=p.id,
                    method="energy fingerprint (grompp -t continuation)",
                )
                findings.append(
                    finding(
                        "LIN003",
                        f"{s.id} continues {p.id}: its first energies equal {p.id}'s last energies",
                        subject=s.id,
                        data={"producer": p.id},
                    )
                )
            elif len(matches) > 1:
                findings.append(
                    finding(
                        "LIN001",
                        f"{s.id} starts in a state that several runs ended in: "
                        + ", ".join(m.id for m in matches),
                        subject=s.id,
                    )
                )
        return findings

    # ---------------------------------------------------------------- queries
    def trajectory_info(self, traj: ParsedFile, natom: int | None) -> TrajectoryInfo:
        info = TrajectoryInfo(first_frame_at_start=True)
        info.format = Mined(value=traj.format, sources=(Source(traj.rel, "content"),))
        r = traj.result
        if not isinstance(r, binary.Frames) or not r.n_frames:
            return info
        info.n_frames = stated(r.n_frames, traj.rel, "frame headers")
        if r.natom is not None:
            info.natom = stated(r.natom, traj.rel, "frame header natoms")
        usable = not (
            r.n_frames > 1 and (r.median_interval_ps is None or r.median_interval_ps <= 0)
        )
        if usable:
            info.first_time_ps = stated(r.first_time_ps, traj.rel, "first frame time", "ps")
            info.last_time_ps = stated(r.last_time_ps, traj.rel, "last frame time", "ps")
            if r.median_interval_ps:
                info.frame_interval_ps = Mined(
                    value=clean_float(r.median_interval_ps),
                    unit="ps",
                    sources=(Source(traj.rel, "median frame time delta"),),
                    method="derived:median time delta",
                    confidence=Confidence.EXACT,
                )
        return info

    def _gro_for(self, natom: int | None, near: str) -> ParsedFile | None:
        cands = [g for g in self._gros if natom is None or g.result.natom == natom]
        if not cands:
            return None
        near_dir = _dir(near)
        cands.sort(key=lambda g: (_dir(g.rel) != near_dir, g.rel))
        return cands[0]

    def system_info(self, topology: ParsedFile) -> SystemInfo | None:
        natom = None
        if isinstance(topology.result, binary.TprHeader):
            natom = topology.result.natom
        info = SystemInfo(source=topology.rel, natom=natom)
        gro = self._gro_for(natom, topology.rel)
        if gro is not None:
            g: structure.Gro = gro.result
            info.natom = info.natom or g.natom
            info.nres = g.nres
            info.title = g.title
            info.residue_counts = dict(g.residue_counts)
            info.residue_atoms = dict(g.residue_atoms)
            info.residue_atom_names = {k: list(v) for k, v in g.residue_atom_names.items()}
            la = g.box_lengths_angles()
            if la is not None:
                info.box_lengths, info.box_angles = la
                info.box_complete = True
                info.periodic_box_type = 1
            info.source = f"{topology.rel} + {gro.rel}"
        elif self._tops:
            top: structure.Top = self._tops[0].result
            info.residue_counts = {name: n for name, n in top.molecules}
        return info

    def coords_natom_time(self, coords: ParsedFile) -> tuple[int | None, float | None]:
        if isinstance(coords.result, structure.Gro):
            return coords.result.natom, coords.result.time_ps
        return None, None

    def coords_box(self, coords: ParsedFile):
        if isinstance(coords.result, structure.Gro):
            return coords.result.box_lengths_angles()
        return None

    def build_evidence(self, build_log: ParsedFile) -> BuildEvidence | None:
        return None

    def topology_evidence(self, topology: ParsedFile) -> BuildEvidence | None:
        if not isinstance(topology.result, structure.Top):
            return None
        ffs, waters = structure.evidence_from_includes(topology.result.includes)
        return BuildEvidence(
            source=topology.rel,
            force_fields=ffs,
            water_models=waters,
            outputs=[posixpath.basename(topology.rel)],
        )

    def mddb_trajectory_format(self, fmt: str) -> str | None:
        return {"xtc": "xtc", "trr": "trr"}.get(fmt)

    def mddb_topology_format(self, fmt: str) -> str | None:
        return {"tpr": "tpr", "top": "top"}.get(fmt)
