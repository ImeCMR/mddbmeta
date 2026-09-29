"""MELD (replica exchange on OpenMM) engine adapter.

A MELD run is one replica-exchange simulation of N ladder positions. After the
run, ``extract_trajectory extract_traj_dcd --replica i`` writes the frames of
ladder position i (fixed alpha slot: an *ensemble*, not a time series), and
``follow_dcd`` writes walkers (continuous, but crossing temperatures and
restraint strengths).

Each ladder DCD becomes one Replica with one production Step. Settings come
from the setup script (static), the leader log and the data store's alpha
history; the temperature of each ladder position follows from its alphas
through the script's temperature scaler. Walkers are recorded, not exported.
"""

from __future__ import annotations

import posixpath
import re
from pathlib import Path
from typing import Any

from mddbmeta.engines.amber import prmtop as amber_prmtop
from mddbmeta.engines.base import BuildEvidence, EngineAdapter, ParsedFile, SniffResult
from mddbmeta.engines.meld import files
from mddbmeta.engines.namd import files as namd_files
from mddbmeta.engines.openmm import classify_forcefield_files
from mddbmeta.errors import MissingDependencyError
from mddbmeta.findings import finding
from mddbmeta.model import (
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

_INDEX = re.compile(r"(\d+)\.dcd$")
TEMP_CHANGE_K = 1.0


class MeldAdapter(EngineAdapter):
    name = "meld"
    program = "MELD"

    def __init__(self) -> None:
        self._by_rel: dict[str, ParsedFile] = {}
        self.walkers: list[str] = []

    # ------------------------------------------------------------------ sniff
    def sniff(self, path: Path, head: bytes) -> SniffResult | None:
        if files.is_block(path, head):
            return SniffResult(FileKind.COORDS, "meld_block", 0.95)
        if b"\x00" in head[:1024]:
            return None
        if files.looks_like_setup(head, path.name):
            return SniffResult(FileKind.CONTROL, "meld_setup", 0.95)
        if files.looks_like_remd_log(head):
            return SniffResult(FileKind.LOG, "remd_log", 0.95)
        if files.looks_like_extract_script(head):
            return SniffResult(FileKind.BUILD_LOG, "meld_extract", 0.9)
        return None

    def claim_foreign(self, pf: ParsedFile) -> bool:
        return pf.ok and pf.format in ("dcd", "prmtop", "rst7", "pdb", "psf")

    # ------------------------------------------------------------------ parse
    def parse(self, path: Path, rel: str, sniffed: SniffResult) -> ParsedFile:
        pf = ParsedFile(path=str(path), rel=rel, kind=sniffed.kind, format=sniffed.format)
        readers = {
            "meld_setup": files.read_setup,
            "remd_log": files.read_remd_log,
            "meld_extract": files.read_extract_script,
            "meld_block": lambda p: str(p),  # read together, in build_steps
        }
        try:
            pf.result = readers[sniffed.format](path)
        except Exception as exc:
            pf.errors.append(
                FileLoadError(sniffed.kind.value, rel, classify_exception(exc), str(exc))
            )
        return pf

    # ------------------------------------------------------------ build steps
    def build_steps(self, root: Path, parsed: list[ParsedFile]) -> tuple[list[Step], list[Finding]]:
        findings: list[Finding] = []
        self._by_rel = {f.rel: f for f in parsed}
        steps: list[Step] = []
        setups = [f for f in parsed if f.format == "meld_setup" and f.ok]
        logs = [f for f in parsed if f.format == "remd_log" and f.ok]
        leader = next((f for f in logs if f.result.role == "leader"), None)
        setup_pf = setups[0] if setups else None
        setup: files.MeldSetup | None = setup_pf.result if setup_pf else None
        if len(setups) > 1:
            findings.append(finding("ENG002", f"several MELD setup scripts; using {setups[0].rel}",
                                    paths=[f.rel for f in setups]))  # fmt: skip

        # which DCDs are ladder positions and which are walkers
        roles: dict[str, str] = {}
        for f in parsed:
            if f.format == "meld_extract" and f.ok:
                for kind, target in f.result.targets:
                    d = posixpath.normpath(
                        posixpath.join(posixpath.dirname(f.rel), posixpath.dirname(target))
                    )
                    roles[d] = kind
        dcds = sorted((f for f in parsed if f.format == "dcd" and f.ok), key=lambda f: f.rel)
        ladder, walkers = [], []
        for f in dcds:
            kind = roles.get(posixpath.dirname(f.rel))
            base = posixpath.basename(f.rel)
            if kind is None:
                kind = (
                    "walker"
                    if base.startswith("follow")
                    else "ladder"
                    if base.startswith("trajectory")
                    else None
                )
            (walkers if kind == "walker" else ladder if kind == "ladder" else []).append(f)
        self.walkers = [f.rel for f in walkers]

        # the data store's alpha history: exact per-position alphas over the run
        blocks = [Path(f.path) for f in parsed if f.format == "meld_block"]
        history = None
        if blocks:
            try:
                history = files.read_alpha_history(blocks)
            except (MissingDependencyError, OSError, ValueError, KeyError) as exc:
                findings.append(
                    finding("LOAD001", f"MELD data store not read: {exc}", paths=["Data/Blocks"])
                )

        topology = next((f for f in parsed if f.format in ("prmtop", "psf") and f.ok), None)
        structure = self._structure(parsed, topology)
        scaler = setup.temperature_scaler if setup else None
        for f in ladder:
            m = _INDEX.search(posixpath.basename(f.rel))
            idx = int(m.group(1)) if m else None
            stem = posixpath.basename(f.rel).rsplit(".", 1)[0]
            rid = f"replica_{idx:02d}" if idx is not None else stem
            step = Step(id=posixpath.join(posixpath.dirname(f.rel), stem), name=stem,
                        dir=posixpath.dirname(f.rel), engine=self.name)  # fmt: skip
            step.replica = rid
            step.raw["meld_replica"] = rid
            step.raw["meld_index"] = idx
            step.files[FileRole.TRAJECTORY] = f.ref()
            if setup_pf is not None:
                step.files[FileRole.CONTROL] = setup_pf.ref()
            if leader is not None:
                step.files[FileRole.LOG] = leader.ref()
            if topology is not None:
                step.files[FileRole.TOPOLOGY] = topology.ref()
            if structure:
                step.raw["structure"] = structure
            step.raw["mddb_type"] = "ensemble"
            step.settings = self._settings(setup, setup_pf.rel if setup_pf else None)
            self._temperature(
                step, idx, history, scaler, setup_pf.rel if setup_pf else None, findings
            )
            step.runtime = self._runtime(leader, setup)
            step.trajectory = self.trajectory_info(f, None)
            step.role = Mined(
                value="production", method="MELD ladder position", confidence=Confidence.DERIVED
            )
            step.raw["mddb_method"] = "Replica exchange (MELD: temperature + Hamiltonian)"
            step.raw["mddb_metadditions"] = self._metadditions(setup, leader, history)
            steps.append(step)

        drifting = [s for s in steps if s.raw.get("meld_temp_range", 0) > TEMP_CHANGE_K]
        if drifting and history is not None:
            settled = max(history.settled_frame)
            findings.append(finding(
                "MELD001",
                f"the replica-exchange ladder adapted during the run: {len(drifting)} of {len(steps)} ladder "
                f"positions changed temperature by more than {TEMP_CHANGE_K:g} K (settled at frame {settled} of "
                f"{history.n_frames}); per-MD temp is the final value",
                data={"positions": [s.raw["meld_replica"] for s in drifting], "settled_frame": settled},
            ))  # fmt: skip
        if history is not None and ladder:
            extracted = {
                int(s.trajectory.n_frames.value) for s in steps if s.trajectory.n_frames.is_known
            }
            if extracted and extracted != {history.n_frames}:
                findings.append(finding(
                    "MELD002",
                    f"extracted DCDs hold {', '.join(map(str, sorted(extracted)))} frames; the data store holds "
                    f"{history.n_frames}",
                    data={"extracted": sorted(extracted), "stored": history.n_frames},
                ))  # fmt: skip
        if walkers:
            findings.append(finding(
                "MELD003",
                f"{len(walkers)} walker trajectories (follow_dcd) are recorded but not exported: each walker "
                "crosses temperatures and restraint strengths",
                paths=[w.rel for w in walkers[:5]],
            ))  # fmt: skip
        if leader is not None:
            restarts = sorted(
                {s for s in leader.result.sessions if s > 1 and s < (leader.result.max_steps or 0)}
            )
            if restarts:
                findings.append(finding("MELD004", f"the run was restarted at exchange(s) {', '.join(map(str, restarts))}",
                                        data={"restarts": restarts}))  # fmt: skip
        if not ladder:
            findings.append(finding("MINE001", "no extracted ladder DCDs (extract_trajectory extract_traj_dcd) found",
                                    data={"field": "input_trajectory_filepaths"}))  # fmt: skip
        return steps, findings

    def declare_replicas(self, steps: list[Step]) -> dict[str, str | None] | None:
        """Ladder positions are the replicas; directory/name inference cannot see them."""
        if not steps or not all(s.raw.get("meld_replica") for s in steps):
            return None
        return {s.id: s.raw["meld_replica"] for s in steps}

    def _structure(self, parsed: list[ParsedFile], topology: ParsedFile | None) -> str | None:
        natom = (
            topology.result.natom
            if topology is not None and hasattr(topology.result, "natom")
            else None
        )
        pdbs = [f for f in parsed if f.format == "pdb" and f.ok]
        pdbs.sort(key=lambda f: (posixpath.dirname(f.rel) != "", "min" not in f.rel, f.rel))
        for f in pdbs:
            if natom is None or f.result.natom == natom:
                return f.rel
        return None

    def _settings(self, setup: files.MeldSetup | None, rel: str | None) -> StepSettings:
        s = StepSettings()
        if setup is None or rel is None:
            return s

        def st(value: Any, locator: str, unit: str | None = None, conf: Confidence = Confidence.EXACT,
               method: str = "stated (setup script, static)") -> Mined[Any]:  # fmt: skip
            return Mined(
                value=value,
                unit=unit,
                sources=(Source(rel, locator),),
                method=method,
                confidence=conf,
            )

        s.minimization = st(False, "replica exchange")
        dt_fs, why = setup.timestep_fs
        s.dt_ps = st(dt_fs / 1000.0, f"{setup.builder or 'system'} options: {why}", "ps", Confidence.DERIVED,
                     "derived: MELD builder timestep rule (2.0 / 3.5 / 4.5 fs)")  # fmt: skip
        if setup.timesteps:
            # one stored frame per exchange; an ensemble has no trajectory write interval (see metadditions)
            if setup.n_steps:
                s.nsteps = st(setup.timesteps * setup.n_steps, "RunOptions.timesteps x N_STEPS")
        s.thermostat = st("langevin", "MELD builder integrator (Langevin, 1/ps)", conf=Confidence.DERIVED,
                          method="derived: MELD builders always use a Langevin integrator")  # fmt: skip
        s.barostat = st(
            "none", "replica exchange (constant volume)", conf=Confidence.DERIVED, method="derived"
        )
        s.constant_pressure = s.barostat.with_value(False)
        implicit = str(
            setup.builder_options.get("solvation_type", "")
        ).lower() == "implicit" or any(
            "implicit" in str(x) for x in setup.builder_options.get("base_forcefield_files") or []
        )
        s.implicit_solvent = st(implicit, "solvation_type / base_forcefield_files")
        s.periodic = st(not implicit, "solvation_type", conf=Confidence.DERIVED, method="derived")
        s.restraints = st(
            True, "MELD restraints (data-guided)", conf=Confidence.DERIVED, method="derived"
        )
        s.restart = st(True, "replica exchange continues each replica's state")
        s.ensemble = Mined(value="NVT", sources=(Source(rel, "Langevin, no barostat"),),
                           method="derived: Langevin integrator, no barostat", confidence=Confidence.DERIVED)  # fmt: skip
        return s

    def _temperature(self, step: Step, idx: int | None, history: files.AlphaHistory | None,
                     scaler: dict[str, Any] | None, rel: str | None, findings: list[Finding]) -> None:  # fmt: skip
        if idx is None or scaler is None:
            return
        if history is not None and idx < history.n_replicas:
            final = files.temperature_of(scaler, history.last[idx])
            lo = files.temperature_of(scaler, history.minimum[idx])
            hi = files.temperature_of(scaler, history.maximum[idx])
            if final is None:
                return
            rng = (hi - lo) if lo is not None and hi is not None else 0.0
            step.raw["meld_temp_range"] = rng
            note = None
            if rng > TEMP_CHANGE_K:
                note = (f"ladder adapted: {lo:.1f}-{hi:.1f} K over the run (alpha {history.minimum[idx]:.4f}-"
                        f"{history.maximum[idx]:.4f}); final {final:.2f} K from frame {history.settled_frame[idx]}")  # fmt: skip
            step.settings.temp_K = Mined(
                value=round(final, 2), unit="K",
                sources=(Source("Data/Blocks", f"alphas[{idx}] (last frame)"), Source(rel or "", "temperature scaler")),
                method="derived: final alpha through the temperature scaler", confidence=Confidence.DERIVED, note=note,
            )  # fmt: skip
            step.raw["meld_alpha_final"] = history.last[idx]
        else:
            findings.append(finding("MELD001", f"{step.id}: data store not read; temperature from the initial ladder",
                                    subject=step.id))  # fmt: skip

    def _runtime(self, leader: ParsedFile | None, setup: files.MeldSetup | None) -> StepRuntime:
        rt = StepRuntime()
        if leader is None:
            return rt
        log: files.RemdLog = leader.result
        rel = leader.rel
        rt.program = stated("MELD", rel, "Meld version is ...")
        if log.meld_version:
            rt.version = stated(log.meld_version, rel, "Meld version is ...")
        variant = f"OpenMM_Meld {log.openmm_meld_version}" if log.openmm_meld_version else None
        if variant:
            rt.variant = stated(
                variant + (f" ({log.platform})" if log.platform else ""),
                rel,
                "OpenMM_Meld version is ...",
            )
        done = log.finished_steps or log.last_step
        if done is not None:
            rt.nsteps_completed = Mined(value=done * (setup.timesteps if setup and setup.timesteps else 1),
                                        sources=(Source(rel, "replica exchange steps"),),
                                        method="derived: exchanges x MD steps per exchange", confidence=Confidence.DERIVED)  # fmt: skip
        rt.completed = Mined(value=log.finished_steps is not None and log.finished_steps >= (log.max_steps or 0),
                             sources=(Source(rel, "Finished N steps of replica exchange successfully"),),
                             method="stated" if log.finished_steps else "absent: Finished line",
                             note=None if log.finished_steps else f"stopped at exchange {log.last_step} of {log.max_steps}")  # fmt: skip
        return rt

    def _metadditions(self, setup: files.MeldSetup | None, leader: ParsedFile | None,
                      history: files.AlphaHistory | None) -> dict[str, Any]:  # fmt: skip
        md: dict[str, Any] = {"meld": {}}
        m = md["meld"]
        if leader is not None:
            m["openmm_meld_version"] = leader.result.openmm_meld_version
            m["platform"] = leader.result.platform
        if setup is not None:
            m["n_replicas"] = setup.n_replicas
            m["exchanges"] = setup.n_steps
            m["md_steps_per_exchange"] = setup.timesteps
            if setup.timesteps:
                m["exchange_interval_ps"] = clean_float(
                    setup.timesteps * setup.timestep_fs[0] / 1000.0
                )
            sc = setup.temperature_scaler or {}
            if sc:
                m["temperature_ladder"] = {k: v for k, v in sc.items()}
            m["adaptor"] = setup.adaptor
            opts = setup.builder_options
            m["builder"] = setup.builder
            if opts.get("grappa_model_tag"):
                m["grappa_model"] = opts["grappa_model_tag"]
            if opts.get("solvation_type"):
                m["solvation"] = opts["solvation_type"]
            implicit = [
                x for x in (opts.get("base_forcefield_files") or []) if "implicit" in str(x)
            ]
            if implicit:
                m["implicit_solvent"] = implicit[0]
        if history is not None:
            m["stored_frames"] = history.n_frames
            m["ladder_settled_frame"] = max(history.settled_frame)
        md["meld"] = {k: v for k, v in m.items() if v is not None}
        return md

    # ---------------------------------------------------------------- queries
    def trajectory_info(self, traj: ParsedFile, natom: int | None) -> TrajectoryInfo:
        info = TrajectoryInfo(first_frame_at_start=True)
        info.format = Mined(value=traj.format, sources=(Source(traj.rel, "content"),))
        r = traj.result
        if isinstance(r, namd_files.Dcd):
            # extract_trajectory writes DELTA=1, NSAVC=1: frame indices, not times
            info.n_frames = stated(r.n_frames, traj.rel, "file size / frame size")
            info.natom = stated(r.natom, traj.rel, "header natoms")
        return info

    def system_info(self, topology: ParsedFile) -> SystemInfo | None:
        r = topology.result
        if isinstance(r, amber_prmtop.Prmtop):
            atoms, names = r.residue_templates()
            return SystemInfo(source=topology.rel, natom=r.natom, nres=r.nres, residue_counts=r.residue_counts(),
                              residue_atoms=atoms, residue_atom_names=names, n_extra_points=r.numextra,
                              periodic_box_type=r.ifbox)  # fmt: skip
        if isinstance(r, namd_files.Psf):
            return SystemInfo(source=topology.rel, natom=r.natom, nres=r.nres, residue_counts=dict(r.residue_counts),
                              residue_atoms=dict(r.residue_atoms),
                              residue_atom_names={k: list(v) for k, v in r.residue_atom_names.items()},
                              periodic_box_type=0)  # fmt: skip
        return None

    def coords_natom_time(self, coords: ParsedFile) -> tuple[int | None, float | None]:
        r = coords.result
        return (r.natom, None) if isinstance(r, namd_files.Pdb) else (None, None)

    def run_evidence(self, steps: list) -> list[BuildEvidence]:
        for s in steps:
            ctl = s.files.get(FileRole.CONTROL)
            pf = self._by_rel.get(ctl.path) if ctl else None
            if pf is None or not isinstance(pf.result, files.MeldSetup):
                continue
            opts = pf.result.builder_options
            ffs: list[str] = []
            tag = opts.get("grappa_model_tag")
            if tag:
                ffs.append("Grappa " + str(tag).removeprefix("grappa-"))
            base, waters = classify_forcefield_files(
                [str(x) for x in opts.get("base_forcefield_files") or []]
            )
            ffs += [f for f in base if f not in ffs]
            return [
                BuildEvidence(source=pf.rel, force_fields=ffs, water_models=waters, linked=True)
            ]
        return []

    def mddb_trajectory_format(self, fmt: str) -> str | None:
        return {"dcd": "dcd"}.get(fmt)

    def mddb_topology_format(self, fmt: str) -> str | None:
        return {"prmtop": "prmtop", "psf": "psf"}.get(fmt)
