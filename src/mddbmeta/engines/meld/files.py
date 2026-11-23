"""MELD run files, read without MELD installed and without unpickling anything.

* setup script (``setup_*.py`` importing ``meld``): read statically with :mod:`ast`;
  run constants, ``RunOptions``, the system-builder options (Grappa/Amber), the
  temperature scaler and the replica-exchange adaptor. Both the current syntax
  (``AmberOptions(forcefield=..., use_big_timestep=...)``) and the older one
  (``SystemBuilder(forcefield=...)`` plus ``options = RunOptions(); options.X = ...``)
  are understood.
* ``Logs/remd_NNN.log`` (or the single ``remd.log`` of MELD <= 0.4, all ranks interleaved):
  MELD / OpenMM_Meld versions, ``launch_remd`` sessions
  (restarts), exchange progress and the "Finished ... successfully" line.
  Only the leader log (``... on leader``) carries progress; worker logs are
  read up to their banner.
* ``Data/Blocks/block_*.nc``: MELD's data store; only the ``alphas`` variable
  (each ladder position's alpha per stored frame) is read.
* ``extract_trajectory`` calls in shell/slurm scripts: which DCDs are ladder
  positions (``extract_traj_dcd``) and which are walkers (``follow_dcd``).

``Data/*.dat`` files are Python pickles and are never unpickled here: that can execute
arbitrary code (``Data/system.dat`` is read only by :mod:`.prmtop`'s inert loader).
"""

from __future__ import annotations

import ast
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta.engines.openmm.script import Quantity, _call_name, _Walker
from mddbmeta.errors import MissingDependencyError
from mddbmeta.io.safe import open_text

# timestep MELD's system builders use (meld/system/builders/{grappa,amber}: _create_integrator;
# meld/system/openmm_runner/runner.py in MELD <= 0.5). use_big_timestep is tested first.
MELD_TIMESTEP_FS = {"default": 2.0, "big": 3.5, "bigger": 4.5}
# the Amber builder's forcefield names -> the leaprc tleap sources (same table in MELD 0.4 to 0.6)
MELD_FORCEFIELDS = {
    "ff12sb": "leaprc.ff12SB",
    "ff14sb": "leaprc.protein.ff14SB",
    "ff14sbside": "leaprc.protein.ff14SBonlysc",
}
MELD_DEFAULT_FORCEFIELD = "ff14sbside"  # SystemBuilder / AmberOptions default
MELD_DEFAULT_IMPLICIT = "gbNeck2"  # RunOptions / AmberOptions default for implicit solvation
# MELD's implicit_solvent_model names -> the OpenMM GB model they select
MELD_GB_MODELS = {"obc": "OBC2", "gbneck": "GBn", "gbneck2": "GBn2"}


# ---------------------------------------------------------------- setup script


@dataclass
class MeldSetup:
    path: str
    constants: dict[str, Any] = field(default_factory=dict)
    n_replicas: int | None = None
    n_steps: int | None = None
    block_size: int | None = None
    timesteps: int | None = None  # MD steps per exchange (RunOptions.timesteps)
    minimize_steps: int | None = None
    builder: str | None = None  # "grappa" | "amber"
    builder_options: dict[str, Any] = field(default_factory=dict)
    run_options: dict[str, Any] = field(default_factory=dict)  # RunOptions(...) / options.X = ...
    legacy: bool = False  # the older SystemBuilder API (MELD <= 0.5)
    temperature_scaler: dict[str, Any] | None = None
    adaptor: str | None = None
    initial_alpha_rule: str | None = None  # e.g. "index / (N_REPLICAS - 1.0)"

    @property
    def options(self) -> dict[str, Any]:
        """Builder options over RunOptions: older MELD keeps timestep/solvent settings in RunOptions."""
        return {**self.run_options, **self.builder_options}

    @property
    def timestep_fs(self) -> tuple[float, str]:
        o = self.options
        if o.get("use_big_timestep") is True:  # MELD checks this one first
            return MELD_TIMESTEP_FS["big"], "use_big_timestep=True"
        if o.get("use_bigger_timestep") is True:
            return MELD_TIMESTEP_FS["bigger"], "use_bigger_timestep=True"
        return MELD_TIMESTEP_FS[
            "default"
        ], "use_big_timestep/use_bigger_timestep False (MELD default)"

    @property
    def explicit_solvent(self) -> bool:
        o = self.options
        if self.builder == "grappa":
            return str(o.get("solvation_type", "")).lower() == "explicit"
        return (
            o.get("explicit_solvent") is True or str(o.get("solvation", "")).lower() == "explicit"
        )

    @property
    def implicit_solvent_model(self) -> tuple[str | None, bool]:
        """(MELD's implicit_solvent_model, stated in the script) for an Amber-built system;
        (None, ...) for explicit solvent, vacuum or a Grappa build."""
        if self.builder != "amber" or self.explicit_solvent:
            return None, True
        model = self.options.get("implicit_solvent_model")
        stated = model is not None
        model = model if stated else MELD_DEFAULT_IMPLICIT
        return (None if str(model).lower() == "vacuum" else str(model)), stated

    @property
    def forcefield(self) -> tuple[str | None, bool]:
        """(Amber builder forcefield name, stated in the script); MELD's default when not given."""
        if self.builder != "amber":
            return None, True
        ff = self.builder_options.get("forcefield")
        return (str(ff), True) if ff is not None else (MELD_DEFAULT_FORCEFIELD, False)

    @property
    def leaprcs(self) -> list[str]:
        """The leaprc files MELD's Amber builder sourced in tleap for this setup."""
        ff, _ = self.forcefield
        out = [MELD_FORCEFIELDS.get(ff.lower(), f"leaprc.{ff}")] if ff else []
        if self.builder == "amber" and self.explicit_solvent:
            out.append(f"leaprc.water.{self.options.get('solvent_forcefield') or 'tip3p'}")
        return out


class _MeldWalker(_Walker):
    """The OpenMM walker's value/template machinery, plus MELD calls and function bodies."""

    def __init__(self, path: str):
        super().__init__(path)
        self.setup = MeldSetup(path=path)
        self._run_option_vars: set[str] = set()  # names bound to RunOptions(...)

    def walk(self, body: list[ast.stmt]) -> None:
        for stmt in body:
            if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)):
                self.walk(stmt.body)  # MELD setup scripts do their work inside setup_system()
                continue
            if isinstance(stmt, ast.Assign) and isinstance(
                stmt.value, (ast.Constant, ast.BinOp, ast.UnaryOp)
            ):
                for t in stmt.targets:
                    if isinstance(t, ast.Name) and t.id.isupper():
                        v = self.value(stmt.value)
                        if v is not None:
                            self.setup.constants[t.id] = v
            if isinstance(stmt, ast.Assign):
                if isinstance(stmt.value, ast.Call) and _call_name(stmt.value) == "RunOptions":
                    self._run_option_vars |= {t.id for t in stmt.targets if isinstance(t, ast.Name)}
                for t in stmt.targets:
                    if isinstance(t, ast.Attribute) and t.attr == "alpha":
                        self.setup.initial_alpha_rule = ast.unparse(stmt.value)
                    elif (isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
                          and t.value.id in self._run_option_vars):  # fmt: skip
                        self.setup.run_options[t.attr] = self.value(stmt.value)
            super().walk([stmt])

    def call(self, node: ast.Call) -> None:
        name = _call_name(node)
        kw = {k.arg: k.value for k in node.keywords if k.arg}
        args = node.args
        s = self.setup
        if name == "RunOptions":
            s.timesteps = _int(self.value(kw.get("timesteps"))) or s.timesteps
            s.minimize_steps = _int(self.value(kw.get("minimize_steps"))) or s.minimize_steps
            s.run_options.update(self._kwargs(kw))
            if args:  # older MELD: RunOptions(solvation="implicit"|"explicit")
                s.run_options.setdefault("solvation", self.value(args[0]))
        elif name in ("GrappaOptions", "AmberOptions"):
            s.builder = "grappa" if name == "GrappaOptions" else "amber"
            s.builder_options = self._kwargs(kw)
        elif name == "SystemBuilder":  # older MELD: SystemBuilder(forcefield="ff14sbside", ...)
            s.builder, s.legacy = "amber", True
            s.builder_options = self._kwargs(kw)
            if args:
                s.builder_options.setdefault("forcefield", self.value(args[0]))
        elif name in (
            "GeometricTemperatureScaler",
            "LinearTemperatureScaler",
            "ConstantTemperatureScaler",
        ):
            vals = [self.value(a) for a in args]
            kind = name.removesuffix("TemperatureScaler").lower()
            if kind == "constant":
                s.temperature_scaler = {
                    "kind": kind,
                    "temperature": _temp(vals[0] if vals else self.value(kw.get("temperature"))),
                }
            elif len(vals) >= 4:
                s.temperature_scaler = {"kind": kind, "alpha_min": _float(vals[0]), "alpha_max": _float(vals[1]),
                                        "t_min": _temp(vals[2]), "t_max": _temp(vals[3])}  # fmt: skip
        elif name in (
            "EqualAcceptanceAdaptor",
            "SwitchingCompositeAdaptor",
            "NullAdaptor",
            "FluxAdaptor",
        ):
            s.adaptor = name
        elif name == "LeaderReplicaExchangeRunner":
            if args:
                s.n_replicas = s.n_replicas or _int(self.value(args[0]))
            s.n_steps = s.n_steps or _int(self.value(kw.get("max_steps")))
        elif name == "DataStore":
            if len(args) >= 2:
                s.n_replicas = s.n_replicas or _int(self.value(args[1]))
            s.block_size = s.block_size or _int(self.value(kw.get("block_size")))
        super().call(node)

    def _kwargs(self, kw: dict[str, ast.expr]) -> dict[str, Any]:
        opts: dict[str, Any] = {}
        for key, v in kw.items():
            if isinstance(v, (ast.List, ast.Tuple)):
                opts[key] = [self.value(e) for e in v.elts]
            else:
                opts[key] = self.value(v)
        return opts


def _int(v: Any) -> int | None:
    return int(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _float(v: Any) -> float | None:
    if isinstance(v, Quantity):
        return v.value
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _temp(v: Any) -> float | None:
    if isinstance(v, Quantity):
        return v.to("temperature") if v.dimension == "temperature" else v.value
    return _float(v)


def read_setup(path: str | Path) -> MeldSetup:
    with open_text(path) as fh:
        tree = ast.parse(fh.read(), filename=str(path))
    w = _MeldWalker(str(path))
    w.walk(tree.body)
    s = w.setup
    c = s.constants
    s.timesteps = s.timesteps or _int(s.run_options.get("timesteps"))
    s.minimize_steps = s.minimize_steps or _int(s.run_options.get("minimize_steps"))
    s.n_replicas = s.n_replicas or _int(c.get("N_REPLICAS"))
    s.n_steps = s.n_steps or _int(c.get("N_STEPS"))
    s.block_size = s.block_size or _int(c.get("BLOCK_SIZE"))
    return s


def looks_like_setup(head: bytes, name: str) -> bool:
    if not name.lower().endswith(".py"):
        return False
    text = head.decode("latin-1", errors="replace")
    return bool(re.search(r"^\s*(import meld\b|from meld\b)", text, re.MULTILINE))


def temperature_of(scaler: dict[str, Any] | None, alpha: float) -> float | None:
    """MELD's temperature for an alpha (meld/system/temperature.py semantics)."""
    if not scaler:
        return None
    kind = scaler.get("kind")
    if kind == "constant":
        return scaler.get("temperature")
    a0, a1, t0, t1 = (scaler.get(k) for k in ("alpha_min", "alpha_max", "t_min", "t_max"))
    if None in (a0, a1, t0, t1):
        return None
    if alpha <= a0:
        return t0
    if alpha > a1:
        return t1
    frac = (alpha - a0) / (a1 - a0)
    if kind == "linear":
        return t0 + frac * (t1 - t0)
    return math.exp((math.log(t1) - math.log(t0)) * frac + math.log(t0))


# ---------------------------------------------------------------- logs

_LAUNCH = re.compile(r"^Launching replica exchange on (leader|worker|follower)")
_VERSION = re.compile(r"^Meld version is (\S+)")
_OMM = re.compile(r"^OpenMM_Meld version is (\S+)")
_STEP = re.compile(r"^Running replica exchange step (\d+) of (\d+)\.")
_FINISHED = re.compile(r"^Finished (\d+) steps of replica exchange successfully")
_PLATFORM = re.compile(r"^Using (\w+) platform\.")
# MELD <= 0.4 writes one remd.log for all ranks, each line prefixed with
# "host:rank date time LEVEL logger: ", e.g. "c44a-s9.ufhpc:000 2021-12-22 11:00:54 INFO meld.remd.launch: "
_PREFIX = re.compile(r"^\S+:\d+ \d{4}-\d\d-\d\d \d\d:\d\d:\d\d \w+ [\w.]+: ")
# the only lines read; everything else (MELD's DEBUG timing lines, the bulk of a log) is skipped unparsed
_KEYS = (
    b"Launching replica",
    b"version is",
    b"Running replica exchange step",
    b"Finished",
    b"platform",
)


@dataclass
class RemdLog:
    path: str
    role: str | None = None  # leader | worker
    meld_version: str | None = None
    openmm_meld_version: str | None = None
    platform: str | None = None
    sessions: list[int] = field(default_factory=list)  # first exchange step of each launch
    last_step: int | None = None
    max_steps: int | None = None
    finished_steps: int | None = None


def read_remd_log(path: str | Path) -> RemdLog:
    """Read a MELD log: ``Logs/remd_NNN.log`` (one file per rank, MELD >= 0.5) or the single
    ``remd.log`` of MELD <= 0.4, which interleaves every rank behind a ``host:rank date ...`` prefix.
    A combined log counts as the leader's when any line of it comes from the leader."""
    log = RemdLog(path=str(path))
    new_session = True  # the first exchange step seen opens a session, launch line or not
    combined = False
    with open(path, "rb") as fh:
        for raw in fh:
            if not any(k in raw for k in _KEYS):
                continue
            line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
            m = _PREFIX.match(line)
            if m:
                combined = True
                line = line[m.end() :]
            m = _LAUNCH.match(line)
            if m:
                role = "worker" if m.group(1) == "follower" else m.group(1)
                if role == "leader":
                    log.role = "leader"
                    new_session = True
                elif log.role is None:
                    log.role = "worker"
                if role == "worker" and log.meld_version and not combined:
                    break  # a worker's own log carries no progress: stop after the banner
                continue
            if log.meld_version is None:
                m = _VERSION.match(line)
                if m:
                    log.meld_version = m.group(1)
                    continue
            if log.openmm_meld_version is None:
                m = _OMM.match(line)
                if m:
                    log.openmm_meld_version = m.group(1)
                    if log.role == "worker" and not combined:
                        break
                    continue
            if line.startswith("Running replica exchange step"):
                m = _STEP.match(line)
                if m:
                    log.role = "leader"  # only the leader logs exchange steps
                    step = int(m.group(1))
                    if new_session:
                        log.sessions.append(step)
                        new_session = False
                    log.last_step, log.max_steps = step, int(m.group(2))
                continue
            if line.startswith("Finished"):
                m = _FINISHED.match(line)
                if m:
                    log.finished_steps = int(m.group(1))
            elif log.platform is None and line.startswith("Using "):
                m = _PLATFORM.match(line)
                if m:
                    log.platform = m.group(1)
    return log


def looks_like_remd_log(head: bytes) -> bool:
    text = head.decode("latin-1", errors="replace")
    return "Launching replica exchange on" in text and "Meld version" in text


# ---------------------------------------------------------------- data store


_BLOCK = re.compile(r"^block_(\d+)\.nc$")


def is_block(path: Path, head: bytes) -> bool:
    return bool(_BLOCK.match(path.name)) and head[:4] in (b"CDF\x01", b"CDF\x02", b"\x89HDF")


@dataclass
class AlphaHistory:
    n_replicas: int
    n_frames: int
    first: list[float]
    last: list[float]
    minimum: list[float]
    maximum: list[float]
    blocks: int
    settled_frame: list[int]  # last frame at which each position's alpha still changed (> tol)


def _block_alphas(ds: Any, start: int):
    """Columns ``start:`` of a block's ``alphas`` as floats (NaN for frames never written)."""
    import numpy as np

    v = ds.variables["alphas"]
    return np.ma.filled(np.ma.asarray(v[:, start:], dtype=float), np.nan)


def read_alpha_history(blocks: list[Path], tol: float = 1e-6) -> AlphaHistory:
    """Scan ``alphas[n_replicas, timesteps]`` of every block, in block order; nothing else is read.

    MELD's blocks share one growing time axis: block k's ``alphas`` runs from frame 0 to the end
    of block k, and only the columns after the previous block's end hold its frames (the rest are
    empty). Each block is read from there, so the cost per block stays constant; reading the whole
    variable grows with the run and, with one-column HDF5 chunks, took most of an hour on long runs."""
    try:
        import netCDF4
        import numpy as np
    except ImportError as exc:  # MELD blocks are NetCDF4/HDF5: scipy cannot read them
        raise MissingDependencyError("Reading MELD data-store blocks", "netcdf") from exc
    ordered = sorted(blocks, key=lambda p: int(_BLOCK.match(p.name).group(1)))  # type: ignore[union-attr]
    arrays = []
    previous = 0  # length of the previous block's time axis
    for path in ordered:
        with netCDF4.Dataset(str(path), "r") as ds:
            length = ds.variables["alphas"].shape[1]
            start = previous if 0 < previous < length else 0  # axis not growing: read it all
            arrays.append(_block_alphas(ds, start))
        previous = length
    arrays = [a for a in arrays if a.ndim == 2 and a.shape[1]]
    if not arrays:
        raise ValueError("no frames in the MELD data store")
    a = np.concatenate(arrays, axis=1)  # [n_replicas, all stored frames]
    a = a[:, ~np.isnan(a).any(axis=0)]  # frames the store allocated but never wrote
    if a.shape[1] == 0:
        raise ValueError("no frames in the MELD data store")
    n, frames = a.shape
    # settled_frame[i]: the last frame at which position i's alpha changed (by more than tol)
    changed = np.abs(np.diff(a, axis=1)) > tol  # True at k: frame k + 1 differs from frame k
    settled = [int(np.flatnonzero(row)[-1]) + 1 if row.any() else 0 for row in changed]
    return AlphaHistory(n_replicas=n, n_frames=frames, first=a[:, 0].tolist(), last=a[:, -1].tolist(),
                        minimum=a.min(axis=1).tolist(), maximum=a.max(axis=1).tolist(), blocks=len(blocks),
                        settled_frame=settled)  # fmt: skip


# ---------------------------------------------------------------- extraction scripts

_EXTRACT = re.compile(r"extract_trajectory\s+(extract_traj_dcd|follow_dcd)\b([^\n]*)")


@dataclass
class ExtractScript:
    path: str
    targets: list[tuple[str, str]] = field(default_factory=list)  # (kind, output path pattern)


def read_extract_script(path: str | Path) -> ExtractScript:
    out = ExtractScript(path=str(path))
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                continue
            m = _EXTRACT.search(s)
            if m:
                tokens = [
                    t for t in m.group(2).split() if not t.startswith("-") and not t.startswith("$")
                ]
                target = next((t for t in tokens if t.endswith(".dcd")), "")
                out.targets.append(
                    ("ladder" if m.group(1) == "extract_traj_dcd" else "walker", target)
                )
    return out


def looks_like_extract_script(head: bytes) -> bool:
    return b"extract_trajectory" in head
