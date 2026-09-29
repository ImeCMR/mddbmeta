"""OpenMM engine adapter.

OpenMM has no input or log format of its own: a run is a Python script. So:

* each StateDataReporter log (CSV) is one Step -- it holds the engine's own
  record of steps and times;
* the script that wrote it is found by matching the log name against the
  scripts' StateDataReporter file templates (``f'prod_{i}.csv'``); the matched
  placeholder values (``i=2``) name the segment's other files (dcd/xtc,
  checkpoints, saved states) and the file it loaded (``prod_{i-1}.chk``);
* settings come from the script, statically (never executed), and from any
  XmlSerializer files it used (System, Integrator), which are exact;
* continuation: the state/checkpoint a segment loaded is another segment's
  saved state/checkpoint; a segment with no load continues the previous
  segment of the same script run (same Simulation, e.g. NVT -> NPT).

Topologies read with AmberPrmtopFile / CharmmPsfFile / GromacsTopFile belong
to other adapters' sniffers; this adapter adopts them.
"""

from __future__ import annotations

import posixpath
from pathlib import Path
from typing import Any

from mddbmeta.engines.amber import prmtop as amber_prmtop
from mddbmeta.engines.base import BuildEvidence, EngineAdapter, ParsedFile, SniffResult
from mddbmeta.engines.gromacs import binary as gmx_binary
from mddbmeta.engines.gromacs.structure import Gro, lengths_angles
from mddbmeta.engines.namd import files as namd_files
from mddbmeta.engines.openmm import files, script
from mddbmeta.findings import finding
from mddbmeta.model import (
    CoordSource,
    CoordSourceKind,
    FileKind,
    FileRole,
    Step,
    StepRuntime,
    StepSettings,
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

THERMOSTATTED = {"LangevinMiddleIntegrator", "LangevinIntegrator", "NoseHooverIntegrator", "BrownianIntegrator",
                 "VariableLangevinIntegrator", "DrudeLangevinIntegrator", "DrudeNoseHooverIntegrator"}  # fmt: skip

# ForceField XML (as passed to ForceField(...)) -> label; checked as prefixes of the path's stem
FORCE_FIELDS = [
    ("amber14-all", "Amber ff14SB"), ("amber14/protein.ff14sb", "Amber ff14SB"),
    ("amber14/protein.ff15ipq", "Amber ff15ipq"), ("amber19-all", "Amber ff19SB"),
    ("amber19/protein.ff19sb", "Amber ff19SB"), ("amber/protein.ff14sb", "Amber ff14SB"),
    ("amber/protein.ff19sb", "Amber ff19SB"), ("amber99sbildn", "Amber ff99SB-ILDN"),
    ("amber99sbnmr", "Amber ff99SB-NMR"), ("amber99sb", "Amber ff99SB"), ("amber03", "Amber ff03"),
    ("amber10", "Amber ff10"), ("amber96", "Amber ff96"), ("charmm36_2024", "CHARMM36 (2024)"),
    ("charmm36", "CHARMM36"), ("charmm/charmm36", "CHARMM36"), ("charmm_polar_2023", "CHARMM Drude 2023"),
    ("charmm_polar_2019", "CHARMM Drude 2019"), ("amoeba2018", "AMOEBA 2018"), ("amoeba2013", "AMOEBA 2013"),
    ("amoeba2009", "AMOEBA 2009"),
]  # fmt: skip
WATERS = {
    "tip3p": "TIP3P", "tip3pfb": "TIP3P-FB", "tip3p_standard": "TIP3P", "tip4pew": "TIP4P-Ew",
    "tip4pfb": "TIP4P-FB", "tip5p": "TIP5P", "spce": "SPC/E", "opc": "OPC", "opc3": "OPC3",
    "water": "TIP3P", "tip3p-pme-b": "TIP3P-PME-B", "tip3p-pme-f": "TIP3P-PME-F", "swm4ndp": "SWM4-NDP",
    "waters_ions_default": "TIP3P",
}  # fmt: skip


def _dir(rel: str) -> str:
    return posixpath.dirname(rel)


def _stem(rel: str) -> str:
    base = posixpath.basename(rel)
    return base.rsplit(".", 1)[0] if "." in base else base


def classify_forcefield_files(names: list[str]) -> tuple[list[str], list[str]]:
    ffs: list[str] = []
    waters: list[str] = []
    for name in names:
        low = name.lower().removesuffix(".xml")
        stem = low.rsplit("/", 1)[-1]
        if stem in WATERS:
            if WATERS[stem] not in waters:
                waters.append(WATERS[stem])
            continue
        if low.startswith("implicit/"):
            continue
        for key, label in FORCE_FIELDS:
            if low.startswith(key) or stem.startswith(key):
                if label not in ffs:
                    ffs.append(label)
                break
    return ffs, waters


def _q(v: Any, dim: str) -> float | None:
    if isinstance(v, script.Quantity):
        return v.to(dim)
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v)  # a bare number: OpenMM's default units (K, ps, bar, nm)
    return None


class OpenmmAdapter(EngineAdapter):
    name = "openmm"
    program = "OpenMM"

    def __init__(self) -> None:
        self._by_rel: dict[str, ParsedFile] = {}
        self._scripts: dict[str, script.ScriptInfo] = {}

    # ------------------------------------------------------------------ sniff
    def sniff(self, path: Path, head: bytes) -> SniffResult | None:
        if namd_files.is_dcd(head):
            return SniffResult(FileKind.TRAJECTORY, "dcd", 0.6)
        if gmx_binary.first_int(head) == gmx_binary.XTC_MAGIC:
            return SniffResult(FileKind.TRAJECTORY, "xtc", 0.6)
        if b"\x00" in head[:1024]:
            return None
        if script.looks_like_script(head, path.name):
            return SniffResult(FileKind.CONTROL, "openmm_script", 0.9)
        if files.looks_like_state_log(head):
            return SniffResult(FileKind.LOG, "statelog", 0.9)
        if files.looks_like_xml_state(head):
            return SniffResult(FileKind.COORDS, "openmm_xml", 0.9)
        if namd_files.looks_like_pdb(head):
            return SniffResult(FileKind.COORDS, "pdb", 0.55)
        return None

    def claim_foreign(self, pf: ParsedFile) -> bool:
        """Topologies and coordinates OpenMM scripts read with the AMBER/CHARMM/GROMACS readers."""
        return pf.ok and pf.format in (
            "prmtop",
            "rst7",
            "psf",
            "top",
            "gro",
            "pdb",
            "leap",
            "xtc",
            "dcd",
        )

    # ------------------------------------------------------------------ parse
    def parse(self, path: Path, rel: str, sniffed: SniffResult) -> ParsedFile:
        pf = ParsedFile(path=str(path), rel=rel, kind=sniffed.kind, format=sniffed.format)
        readers = {
            "dcd": namd_files.read_dcd,
            "xtc": gmx_binary.read_xtc,
            "openmm_script": script.read_script,
            "statelog": files.read_state_log,
            "openmm_xml": files.read_xml_state,
            "pdb": namd_files.read_pdb,
        }
        try:
            pf.result = readers[sniffed.format](path)
        except Exception as exc:
            pf.errors.append(
                FileLoadError(sniffed.kind.value, rel, classify_exception(exc), str(exc))
            )
            return pf
        if isinstance(pf.result, script.ScriptInfo) and not pf.result.is_openmm:
            pf.result = None
            pf.errors.append(
                FileLoadError("control", rel, "unsupported", "not an OpenMM run script")
            )
        return pf

    # ------------------------------------------------------------ build steps
    def build_steps(self, root: Path, parsed: list[ParsedFile]) -> tuple[list[Step], list[Finding]]:
        findings: list[Finding] = []
        by_rel = {f.rel: f for f in parsed}
        self._by_rel = by_rel
        scripts = [f for f in parsed if f.format == "openmm_script" and f.ok]
        self._scripts = {f.rel: f.result for f in scripts}
        steps: list[Step] = []
        claimed: set[str] = set()
        produced: dict[str, str] = {}  # file written (state/checkpoint) -> step id
        segment_step: dict[
            tuple[str, int, str], str
        ] = {}  # (script, segment index, placeholder key) -> step id
        pending_links: list[tuple[Step, str, int, dict[str, str]]] = []

        def resolve(
            run_dir: str, t: script.Template | None, values: dict[str, str]
        ) -> tuple[str | None, ParsedFile | None]:
            if t is None:
                return None, None
            name = t.fill(values)
            if name is None:
                return None, None
            rel = posixpath.normpath(posixpath.join(run_dir, name))
            return rel, by_rel.get(rel)

        logs = sorted(
            (f for f in parsed if f.kind is FileKind.LOG and f.format == "statelog"),
            key=lambda f: f.rel,
        )
        for lf in logs:
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
            match = self._find_segment(lf.rel, scripts)
            if match is not None:
                sf, idx, values = match
                info: script.ScriptInfo = sf.result
                seg = info.segments[idx]
                sdir = _dir(sf.rel)
                step.files[FileRole.CONTROL] = sf.ref()
                step.raw["openmm_values"] = values
                segment_step[(sf.rel, idx, repr(sorted(values.items())))] = step.id
                for _fmt, t, _interval in seg.trajectories:
                    _rel, tf = resolve(sdir, t, values)
                    if (
                        tf is not None
                        and tf.kind is FileKind.TRAJECTORY
                        and FileRole.TRAJECTORY not in step.files
                    ):
                        step.files[FileRole.TRAJECTORY] = tf.ref()
                        step.trajectory = self.trajectory_info(tf, None)
                        claimed.add(tf.rel)
                for t in seg.saves + seg.checkpoints_out:
                    rel, f = resolve(sdir, t, values)
                    if rel:
                        produced.setdefault(rel, step.id)
                        if (
                            f is not None
                            and FileRole.OUTPUT_RESTART not in step.files
                            and f.format == "openmm_xml"
                        ):
                            step.files[FileRole.OUTPUT_RESTART] = f.ref()
                for t in info.written_structures:
                    _rel, f = resolve(sdir, t, values)
                    if f is not None:
                        step.raw["structure"] = f.rel  # the modelled system as simulated
                        break
                for fmt, t in info.topology:
                    rel, f = resolve(sdir, t, values)
                    if f is None:
                        continue
                    if fmt in ("pdb", "pdbx") and info.modeller:
                        continue  # only the building input: Modeller changed the system
                    if fmt in ("pdb", "pdbx"):
                        step.raw["structure"] = f.rel
                        if FileRole.TOPOLOGY not in step.files:
                            step.raw["topology_candidate"] = f.rel
                    elif f.kind is FileKind.TOPOLOGY:
                        step.files[FileRole.TOPOLOGY] = f.ref()
                pending_links.append((step, sf.rel, idx, values))
                self._fill_from_script(step, sf, seg, values, by_rel)
            self._fill_from_log(step, lf)
            steps.append(step)

        # continuations, once every segment's outputs are known
        for step, srel, idx, values in pending_links:
            info = self._scripts[srel]
            seg = info.segments[idx]
            sdir = _dir(srel)
            chosen = None
            for prefer_loop in (True, False):
                for t, in_loop in seg.loads:
                    if in_loop is not prefer_loop:
                        continue
                    rel, f = resolve(sdir, t, values)
                    if rel and (rel in produced or f is not None) and produced.get(rel) != step.id:
                        chosen = rel
                        break
                if chosen:
                    break
            if chosen is not None:
                producer = produced.get(chosen)
                f = by_rel.get(chosen)
                if f is not None:
                    step.files[FileRole.INPUT_COORDS] = f.ref()
                    step.raw["loaded_state"] = chosen
                if producer:
                    step.input_coords = CoordSource(CoordSourceKind.STEP, ref=producer, path=chosen,
                                                    method=f"script loads {posixpath.basename(chosen)}")  # fmt: skip
                else:
                    step.input_coords = CoordSource(
                        CoordSourceKind.PATH, path=chosen, method="script load"
                    )
            elif seg.continues_previous and idx > 0:
                prev = segment_step.get((srel, idx - 1, repr(sorted(values.items()))))
                if prev:
                    step.input_coords = CoordSource(CoordSourceKind.STEP, ref=prev,
                                                    method="same script run continues (no state reload)")  # fmt: skip
            else:
                # fresh start: the coordinates the script read (inpcrd/crd/gro, else the PDB)
                if info.modeller:
                    starts = list(info.written_structures)
                else:
                    starts = list(info.coordinates) + [
                        t for fmt, t in info.topology if fmt in ("pdb", "pdbx")
                    ]
                for t in starts:
                    _rel, f = resolve(sdir, t, values)
                    if f is not None:
                        step.files[FileRole.INPUT_COORDS] = f.ref()  # the core resolves the source
                        break
            self._fill_times(step, by_rel)

        # scripts whose segments never logged anything
        logged_scripts = {
            s.files[FileRole.CONTROL].path for s in steps if FileRole.CONTROL in s.files
        }
        for sf in scripts:
            if sf.rel in logged_scripts or not sf.result.segments:
                continue
            run_dir, stem = _dir(sf.rel), _stem(sf.rel)
            step = Step(id=posixpath.join(run_dir, stem) if run_dir else stem, name=stem, dir=run_dir,
                        engine=self.name, status="queued")  # fmt: skip
            step.files[FileRole.CONTROL] = sf.ref()
            steps.append(step)
            findings.append(finding("STEP001", f"{sf.rel}: no StateDataReporter log found; treated as not run",
                                    subject=step.id, paths=[sf.rel]))  # fmt: skip

        for f in sorted(parsed, key=lambda f: f.rel):
            if (
                f.kind is not FileKind.TRAJECTORY
                or f.rel in claimed
                or f.format not in ("dcd", "xtc")
            ):
                continue
            run_dir, stem = _dir(f.rel), _stem(f.rel)
            sid = posixpath.join(run_dir, stem) if run_dir else stem
            if any(s.id == sid for s in steps):
                sid = f.rel
            step = Step(id=sid, name=stem, dir=run_dir, engine=self.name, status="orphan")
            step.files[FileRole.TRAJECTORY] = f.ref()
            step.trajectory = self.trajectory_info(f, None)
            steps.append(step)
            findings.append(finding("STEP002", f"{f.rel} has no StateDataReporter log or script",
                                    subject=sid, paths=[f.rel]))  # fmt: skip
        return steps, findings

    def _find_segment(self, log_rel: str, scripts: list[ParsedFile]):
        """(script, segment index, placeholder values) whose StateDataReporter wrote this log."""
        run_dir = _dir(log_rel)
        candidates = sorted(scripts, key=lambda s: (_dir(s.rel) != run_dir, s.rel))
        for sf in candidates:
            sdir = _dir(sf.rel)
            rel_to_script = posixpath.relpath(log_rel, sdir) if sdir else log_rel
            for idx, seg in enumerate(sf.result.segments):
                if seg.log is None:
                    continue
                values = seg.log.match(rel_to_script)
                if values is not None:
                    return sf, idx, values
        # a log captured from stdout: same stem as a script that reports to stdout
        for sf in candidates:
            if _dir(sf.rel) == run_dir and _stem(sf.rel) == _stem(log_rel):
                for idx, seg in enumerate(sf.result.segments):
                    if seg.log is None:
                        return sf, idx, {}
        return None

    # ------------------------------------------------------------- settings
    def _fill_from_script(self, step: Step, sf: ParsedFile, seg: script.Segment, values: dict[str, str],
                          by_rel: dict[str, ParsedFile]) -> None:  # fmt: skip
        info: script.ScriptInfo = sf.result
        rel = sf.rel
        s = StepSettings()

        def st(value: Any, locator: str, unit: str | None = None) -> Mined[Any]:
            return Mined(
                value=value,
                unit=unit,
                sources=(Source(rel, locator),),
                method="stated (script, static)",
            )

        system_xml = integrator_xml = None
        for t in info.deserialized:
            name = t.fill(values)
            f = by_rel.get(posixpath.normpath(posixpath.join(_dir(rel), name))) if name else None
            if f is not None and isinstance(f.result, files.XmlState):
                if f.result.kind == "System":
                    system_xml = f
                elif f.result.kind == "Integrator":
                    integrator_xml = f
        step.raw["openmm_system_xml"] = system_xml.rel if system_xml else None

        s.minimization = st(False, "step()")
        if seg.steps:
            s.nsteps = st(seg.steps, "step(N)")
        integ = dict(seg.integrator)
        if integrator_xml is not None:
            ia = integrator_xml.result.integrator
            integ = {"type": ia.get("type"), "temperature": _num(ia.get("temperature")),
                     "friction": _num(ia.get("friction")), "dt": _num(ia.get("stepSize"))}  # fmt: skip
        dt = _q(integ.get("dt"), "time")
        if dt is not None:
            s.dt_ps = st(round(dt, 9), f"{integ.get('type')}(stepSize)", "ps")
        itype = integ.get("type")
        thermo = None
        temp = None
        if itype in THERMOSTATTED:
            thermo = str(itype).removesuffix("Integrator").lower()
            temp = _q(integ.get("temperature"), "temperature")
        elif seg.thermostat:
            thermo = "andersen"
            temp = _q(seg.thermostat.get("temperature"), "temperature")
        s.thermostat = st(thermo or "none", f"{itype}")
        if seg.temperature_changes:
            temp = seg.temperature_changes[-1]
            s.temp_ramp = st(len(set(seg.temperature_changes)) > 1, "integrator.setTemperature()")
        else:
            s.temp_ramp = st(False, "integrator.setTemperature()")
        if thermo and temp is not None:
            where = (
                "AndersenThermostat(temperature)"
                if thermo == "andersen"
                else f"{itype}(temperature)"
            )
            s.temp_K = st(float(temp), where, "K")
        baro = seg.barostat
        baro_src = "addForce(*Barostat)"
        if system_xml is not None:
            xml_baro = next(
                (fo for fo in system_xml.result.forces if "Barostat" in fo.get("type", "")), None
            )
            if xml_baro is not None:
                baro = {"type": xml_baro.get("type")}
                baro_src = f"{system_xml.rel}:<Force type={xml_baro.get('type')}>"
        s.barostat = st(baro["type"] if baro else "none", baro_src)
        s.constant_pressure = st(bool(baro), baro_src)
        s.restart = st(not seg.new_velocities, "setVelocitiesToTemperature()")
        restr = seg.restraints
        if system_xml is not None and any(
            fo.get("type") == "CustomExternalForce" for fo in system_xml.result.forces
        ):
            restr = True
        s.restraints = st(restr, "CustomExternalForce")
        traj_interval = next((n for _f, _t, n in seg.trajectories if n), None)
        if traj_interval:
            s.traj_interval_steps = st(traj_interval, "Reporter(reportInterval)")
        cutoff = _q(info.system_options.get("nonbondedCutoff"), "length")
        if cutoff is None and system_xml is not None:
            nb = next(
                (fo for fo in system_xml.result.forces if fo.get("type") == "NonbondedForce"), None
            )
            cutoff = _num(nb.get("cutoff")) if nb else None
        if cutoff is not None:
            s.cutoff_A = st(cutoff * 10.0, "nonbondedCutoff", "A")
        if thermo:
            ens = "NPT" if baro else "NVT"
        else:
            ens = "NPH" if baro else "NVE"
        s.ensemble = Mined(value=ens, sources=(Source(rel, "integrator/barostat"),),
                           method="derived: integrator and barostat", confidence=Confidence.DERIVED)  # fmt: skip
        step.settings = s
        step.raw["openmm_report_interval"] = seg.report_interval
        step.raw["openmm_forcefields"] = list(info.force_fields) + list(info.charmm_parameters)
        step.raw["openmm_generators"] = list(info.template_generators)

    def _fill_from_log(self, step: Step, lf: ParsedFile) -> None:
        log: files.StateLog | None = lf.result if lf.ok else None
        rt = StepRuntime()
        rt.program = Mined(value="OpenMM", sources=(Source(lf.rel, "StateDataReporter format"),),
                           method="derived: StateDataReporter log", confidence=Confidence.DERIVED)  # fmt: skip
        step.runtime = rt
        if log is None or not log.first:
            return
        step.raw["openmm_log"] = log
        f, la = log.first, log.last
        if "Step" in f and "Time (ps)" in f and la["Step"] > f["Step"]:
            dt_ps = round((la["Time (ps)"] - f["Time (ps)"]) / (la["Step"] - f["Step"]), 9)
            step.raw["openmm_dt_from_log"] = dt_ps
            script_dt = step.settings.dt_ps
            if not script_dt.is_known or abs(float(script_dt.value) - dt_ps) > 1e-6 * max(
                dt_ps, 1e-9
            ):
                # the engine's own record beats a statically read script
                note = f"script states {script_dt.value} ps" if script_dt.is_known else None
                step.settings.dt_ps = Mined(value=dt_ps, unit="ps", sources=(Source(lf.rel, "Time/Step columns"),),
                                            method="engine_report: dTime/dStep", confidence=Confidence.EXACT,
                                            note=note)  # fmt: skip
        if not step.settings.traj_interval_steps.is_known and log.report_interval:
            step.raw["openmm_log_interval"] = log.report_interval

    def _fill_times(self, step: Step, by_rel: dict[str, ParsedFile]) -> None:
        log: files.StateLog | None = step.raw.get("openmm_log")
        rt = step.runtime
        loaded = step.raw.get("loaded_state")
        start_step = start_time = None
        lf = by_rel.get(loaded) if loaded else None
        if lf is not None and isinstance(lf.result, files.XmlState) and lf.result.kind == "State":
            start_step, start_time = lf.result.step, lf.result.time_ps
            if start_time is not None:
                rt.start_time_ps = stated(clean_float(start_time), lf.rel, "<State time>", "ps")
        if log is None or not log.first:
            return
        interval = log.report_interval or step.raw.get("openmm_report_interval")
        dt = step.raw.get("openmm_dt_from_log") or (
            step.settings.dt_ps.value if step.settings.dt_ps.is_known else None
        )
        first_step = log.first.get("Step")
        if start_step is None and first_step is not None and interval:
            start_step = int(first_step - interval)
        if not rt.start_time_ps.is_known and "Time (ps)" in log.first and interval and dt:
            rt.start_time_ps = Mined(value=clean_float(log.first["Time (ps)"] - interval * dt), unit="ps",
                                     sources=(Source(step.files[FileRole.LOG].path, "first row - interval"),),
                                     method="derived: first logged time - report interval", confidence=Confidence.DERIVED)  # fmt: skip
        if "Time (ps)" in log.last:
            rt.end_time_ps = Mined(value=clean_float(log.last["Time (ps)"]), unit="ps",
                                   sources=(Source(step.files[FileRole.LOG].path, "last row Time (ps)"),),
                                   method="engine_report", confidence=Confidence.EXACT)  # fmt: skip
        last_step = log.last.get("Step")
        if last_step is not None and start_step is not None:
            done = int(last_step - start_step)
            rt.nsteps_completed = Mined(value=done, sources=(Source(step.files[FileRole.LOG].path, "Step column"),),
                                        method="derived: last logged step - start step", confidence=Confidence.EXACT)  # fmt: skip
            planned = step.settings.nsteps.value if step.settings.nsteps.is_known else None
            if planned:
                rt.completed = Mined(value=done >= planned, sources=(Source(step.files[FileRole.LOG].path, "Step column"),),
                                     method="derived: logged steps vs step(N) in the script",
                                     confidence=Confidence.DERIVED,
                                     note=None if done >= planned else f"logged {done} of {planned} steps")  # fmt: skip
        # OpenMM's DCD/XTC reporters count time from their own creation, not the simulation clock:
        # shift the file's frame times onto the run clock.
        tr = step.trajectory
        if rt.start_time_ps.is_known and tr.first_time_ps.is_known and tr.last_time_ps.is_known:
            t0 = float(rt.start_time_ps.value)
            for attr in ("first_time_ps", "last_time_ps"):
                m = getattr(tr, attr)
                setattr(tr, attr, m.with_value(clean_float(t0 + float(m.value)),
                                               method="derived: reporter-relative time + run start",
                                               confidence=Confidence.DERIVED))  # fmt: skip
        version = self._version_for(step, by_rel)
        if version is not None:
            rt.version = version

    def _version_for(self, step: Step, by_rel: dict[str, ParsedFile]) -> Mined[Any] | None:
        ref = step.files.get(FileRole.OUTPUT_RESTART)
        f = by_rel.get(ref.path) if ref is not None else None
        if f is not None and isinstance(f.result, files.XmlState) and f.result.version:
            return stated(f.result.version, f.rel, "openmmVersion")  # written by this very run
        for rel in (step.raw.get("openmm_system_xml"), step.raw.get("loaded_state")):
            f = by_rel.get(rel) if rel else None
            if f is not None and isinstance(f.result, files.XmlState) and f.result.version:
                return Mined(value=f.result.version, sources=(Source(f.rel, "openmmVersion"),),
                             method="version of an XML file this run read", confidence=Confidence.DERIVED)  # fmt: skip
        return None

    # ---------------------------------------------------------------- queries
    def trajectory_info(self, traj: ParsedFile, natom: int | None) -> TrajectoryInfo:
        info = TrajectoryInfo(first_frame_at_start=False)  # reporters write after each interval
        info.format = Mined(value=traj.format, sources=(Source(traj.rel, "content"),))
        r = traj.result
        if isinstance(r, namd_files.Dcd) and r.n_frames:
            info.n_frames = stated(r.n_frames, traj.rel, "file size / frame size")
            info.natom = stated(r.natom, traj.rel, "header natoms")
            if r.timestep_fs:
                info.first_time_ps = stated(r.time_ps(0), traj.rel, "ISTART x DELTA", "ps")
                info.last_time_ps = stated(r.time_ps(r.n_frames - 1), traj.rel, "last frame", "ps")
                info.frame_interval_ps = Mined(value=clean_float(r.step_interval * r.timestep_fs / 1000.0), unit="ps",
                                               sources=(Source(traj.rel, "NSAVC x DELTA"),), method="derived: NSAVC x DELTA",
                                               confidence=Confidence.EXACT)  # fmt: skip
        elif isinstance(r, gmx_binary.Frames) and r.n_frames:
            info.n_frames = stated(r.n_frames, traj.rel, "frame headers")
            info.natom = stated(r.natom, traj.rel, "frame header natoms")
            if r.median_interval_ps and r.median_interval_ps > 0:
                info.first_time_ps = stated(r.first_time_ps, traj.rel, "first frame time", "ps")
                info.last_time_ps = stated(r.last_time_ps, traj.rel, "last frame time", "ps")
                info.frame_interval_ps = Mined(value=clean_float(r.median_interval_ps), unit="ps",
                                               sources=(Source(traj.rel, "median frame time delta"),),
                                               method="derived:median time delta", confidence=Confidence.EXACT)  # fmt: skip
        return info

    def system_info(self, topology: ParsedFile) -> SystemInfo | None:
        r = topology.result
        if isinstance(r, amber_prmtop.Prmtop):
            atoms, names = r.residue_templates()
            info = SystemInfo(source=topology.rel, natom=r.natom, nres=r.nres, residue_counts=r.residue_counts(),
                              residue_atoms=atoms, residue_atom_names=names, n_extra_points=r.numextra,
                              periodic_box_type=r.ifbox)  # fmt: skip
        elif isinstance(r, (namd_files.Psf, namd_files.Pdb, Gro)):
            info = SystemInfo(source=topology.rel, natom=r.natom, nres=r.nres, residue_counts=dict(r.residue_counts),
                              residue_atoms=dict(r.residue_atoms),
                              residue_atom_names={k: list(v) for k, v in r.residue_atom_names.items()})  # fmt: skip
        else:
            return None
        box = self._box()
        if box is not None:
            info.box_lengths, info.box_angles = lengths_angles(box)
            info.box_complete = True
            info.periodic_box_type = 1
        return info

    def _box(self):
        for pf in self._by_rel.values():
            if (
                pf.format == "openmm_xml"
                and pf.ok
                and pf.result.box_vectors_A
                and pf.result.kind == "State"
            ):
                return pf.result.box_vectors_A
        for pf in self._by_rel.values():
            if pf.format == "openmm_xml" and pf.ok and pf.result.box_vectors_A:
                return pf.result.box_vectors_A
        return None

    def coords_natom_time(self, coords: ParsedFile) -> tuple[int | None, float | None]:
        r = coords.result
        if isinstance(r, files.XmlState) and r.kind == "State":
            return r.natom, r.time_ps
        if isinstance(r, namd_files.Pdb):
            return r.natom, None
        return None, None

    def coords_box(self, coords: ParsedFile):
        r = coords.result
        if isinstance(r, files.XmlState) and r.box_vectors_A:
            return lengths_angles(r.box_vectors_A)
        return None

    def run_evidence(self, steps: list) -> list[BuildEvidence]:
        out: list[BuildEvidence] = []
        seen: set[str] = set()
        for s in steps:
            ctl = s.files.get(FileRole.CONTROL)
            if ctl is None or ctl.path in seen:
                continue
            seen.add(ctl.path)
            names = s.raw.get("openmm_forcefields") or []
            ffs, waters = classify_forcefield_files(names)
            for g in s.raw.get("openmm_generators") or []:
                label = _generator_label(g)
                if label and label not in ffs:
                    ffs.append(label)
            if ffs or waters:
                out.append(
                    BuildEvidence(
                        source=ctl.path, force_fields=ffs, water_models=waters, linked=True
                    )
                )
        return out

    def mddb_trajectory_format(self, fmt: str) -> str | None:
        return {"dcd": "dcd", "xtc": "xtc"}.get(fmt)

    def mddb_topology_format(self, fmt: str) -> str | None:
        return {"prmtop": "prmtop", "psf": "psf", "top": "top"}.get(fmt)


def _num(text: str | None) -> float | None:
    try:
        return float(text) if text is not None else None
    except ValueError:
        return None


def _generator_label(name: str) -> str | None:
    low = name.lower()
    if low.startswith("openff-"):
        return "OpenFF " + name.split("-", 1)[1].removesuffix(".offxml")
    if low.startswith("gaff-"):
        return "GAFF " + name.split("-", 1)[1]
    if low.startswith("espaloma"):
        return "espaloma " + name.split("-", 1)[1] if "-" in name else "espaloma"
    return None
