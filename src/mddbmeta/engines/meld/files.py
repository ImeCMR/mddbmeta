"""MELD run files, read without MELD installed and without unpickling anything.

* setup script (``setup_*.py`` importing ``meld``): read statically with :mod:`ast`;
  run constants, ``RunOptions``, the system-builder options (Grappa/Amber), the
  temperature scaler and the replica-exchange adaptor.
* ``Logs/remd_NNN.log``: MELD / OpenMM_Meld versions, ``launch_remd`` sessions
  (restarts), exchange progress and the "Finished ... successfully" line.
  Only the leader log (``... on leader``) carries progress; worker logs are
  read up to their banner.
* ``Data/Blocks/block_*.nc``: MELD's data store; only the ``alphas`` variable
  (each ladder position's alpha per stored frame) is read.
* ``extract_trajectory`` calls in shell/slurm scripts: which DCDs are ladder
  positions (``extract_traj_dcd``) and which are walkers (``follow_dcd``).

``Data/*.dat`` files are Python pickles and are deliberately never opened:
unpickling can execute arbitrary code.
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

# timestep MELD's system builders use (meld/system/builders/{grappa,amber}: _create_integrator)
MELD_TIMESTEP_FS = {"default": 2.0, "big": 3.5, "bigger": 4.5}


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
    temperature_scaler: dict[str, Any] | None = None
    adaptor: str | None = None
    initial_alpha_rule: str | None = None  # e.g. "index / (N_REPLICAS - 1.0)"

    @property
    def timestep_fs(self) -> tuple[float, str]:
        o = self.builder_options
        if o.get("use_bigger_timestep") is True:
            return MELD_TIMESTEP_FS["bigger"], "use_bigger_timestep=True"
        if o.get("use_big_timestep") is True:
            return MELD_TIMESTEP_FS["big"], "use_big_timestep=True"
        return MELD_TIMESTEP_FS[
            "default"
        ], "use_big_timestep/use_bigger_timestep False (MELD default)"


class _MeldWalker(_Walker):
    """The OpenMM walker's value/template machinery, plus MELD calls and function bodies."""

    def __init__(self, path: str):
        super().__init__(path)
        self.setup = MeldSetup(path=path)

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
                for t in stmt.targets:
                    if isinstance(t, ast.Attribute) and t.attr == "alpha":
                        self.setup.initial_alpha_rule = ast.unparse(stmt.value)
            super().walk([stmt])

    def call(self, node: ast.Call) -> None:
        name = _call_name(node)
        kw = {k.arg: k.value for k in node.keywords if k.arg}
        args = node.args
        s = self.setup
        if name == "RunOptions":
            s.timesteps = _int(self.value(kw.get("timesteps"))) or s.timesteps
            s.minimize_steps = _int(self.value(kw.get("minimize_steps"))) or s.minimize_steps
        elif name in ("GrappaOptions", "AmberOptions"):
            s.builder = "grappa" if name == "GrappaOptions" else "amber"
            opts: dict[str, Any] = {}
            for key, v in kw.items():
                if isinstance(v, (ast.List, ast.Tuple)):
                    opts[key] = [self.value(e) for e in v.elts]
                else:
                    opts[key] = self.value(v)
            s.builder_options = opts
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

_LAUNCH = re.compile(r"^Launching replica exchange on (leader|worker)")
_VERSION = re.compile(r"^Meld version is (\S+)")
_OMM = re.compile(r"^OpenMM_Meld version is (\S+)")
_STEP = re.compile(r"^Running replica exchange step (\d+) of (\d+)\.")
_FINISHED = re.compile(r"^Finished (\d+) steps of replica exchange successfully")
_PLATFORM = re.compile(r"^Using (\w+) platform\.")


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
    log = RemdLog(path=str(path))
    new_session = False
    with open_text(path) as fh:
        for line in fh:
            m = _LAUNCH.match(line)
            if m:
                log.role = log.role or m.group(1)
                new_session = True
                if m.group(1) == "worker" and log.meld_version:
                    break  # workers carry no progress: stop after the banner
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
                    if log.role == "worker":
                        break
                    continue
            if line.startswith("Running replica exchange step"):
                m = _STEP.match(line)
                if m:
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


def read_alpha_history(blocks: list[Path], tol: float = 1e-6) -> AlphaHistory:
    """Scan ``alphas[n_replicas, timesteps]`` of every block, in block order; nothing else is read."""
    try:
        import netCDF4
    except ImportError as exc:  # MELD blocks are NetCDF4/HDF5: scipy cannot read them
        raise MissingDependencyError("Reading MELD data-store blocks", "netcdf") from exc
    first = last = mn = mx = None
    settled: list[int] = []
    frames = 0
    for path in sorted(blocks, key=lambda p: int(_BLOCK.match(p.name).group(1))):  # type: ignore[union-attr]
        with netCDF4.Dataset(str(path), "r") as ds:
            a = ds.variables["alphas"][:].tolist()  # [n_replicas][timesteps]
        n, t = len(a), len(a[0]) if a else 0
        for j in range(t):
            raw = [a[i][j] for i in range(n)]
            if any(v is None for v in raw):
                continue  # masked: a frame the store allocated but never wrote (end of the last block)
            col = [float(v) for v in raw]
            if first is None:
                first, last, mn, mx = list(col), list(col), list(col), list(col)
                settled = [0] * n
            else:
                for i in range(n):
                    if abs(col[i] - last[i]) > tol:
                        settled[i] = frames
                    mn[i] = min(mn[i], col[i])
                    mx[i] = max(mx[i], col[i])
                last = col
            frames += 1
    if first is None:
        raise ValueError("no frames in the MELD data store")
    return AlphaHistory(n_replicas=len(first), n_frames=frames, first=first, last=last, minimum=mn,
                        maximum=mx, blocks=len(blocks), settled_frame=settled)  # fmt: skip


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
