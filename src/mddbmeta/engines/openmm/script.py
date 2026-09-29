"""Static analysis of OpenMM run scripts (Python), without executing them.

The script is parsed with :mod:`ast`. Only literals, simple arithmetic and
OpenMM unit expressions (``300*kelvin``, ``0.002*picoseconds``,
``1/picosecond``) are evaluated. The walk records, in statement order:

* inputs: ``PDBFile``/``AmberPrmtopFile``/``CharmmPsfFile``/... and
  ``ForceField(...)`` / ``CharmmParameterSet(...)`` files;
* the system options (``createSystem`` keywords), integrators, thermostats and
  barostats, and ``setTemperature`` calls (temperature ramps);
* reporters (StateDataReporter, DCD/XTC/PDB, Checkpoint), ``step(N)``,
  ``minimizeEnergy``, ``load*`` / ``save*`` calls.

A *segment* is what one StateDataReporter logged: the step() calls issued
after it was created and before the next one. File names may be f-strings
(``f'prod_{i}.dcd'``); they are kept as templates and matched to real files.
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta.io.safe import open_text

# unit name -> (dimension, factor to the canonical unit: ps, K, bar, nm, amu, 1/ps)
UNITS: dict[str, tuple[str, float]] = {
    "femtosecond": ("time", 1e-3), "femtoseconds": ("time", 1e-3), "fs": ("time", 1e-3),
    "picosecond": ("time", 1.0), "picoseconds": ("time", 1.0), "ps": ("time", 1.0),
    "nanosecond": ("time", 1e3), "nanoseconds": ("time", 1e3), "ns": ("time", 1e3),
    "kelvin": ("temperature", 1.0), "K": ("temperature", 1.0),
    "bar": ("pressure", 1.0), "bars": ("pressure", 1.0), "atmosphere": ("pressure", 1.01325),
    "atmospheres": ("pressure", 1.01325), "atm": ("pressure", 1.01325),
    "nanometer": ("length", 1.0), "nanometers": ("length", 1.0), "nm": ("length", 1.0),
    "angstrom": ("length", 0.1), "angstroms": ("length", 0.1),
    "amu": ("mass", 1.0), "dalton": ("mass", 1.0), "daltons": ("mass", 1.0),
    "molar": ("concentration", 1.0),
}  # fmt: skip

# OpenMM names used as option values (nonbondedMethod=PME, constraints=HBonds, ...)
CONSTANTS = {"PME", "LJPME", "Ewald", "NoCutoff", "CutoffPeriodic", "CutoffNonPeriodic", "HBonds",
             "AllBonds", "HAngles", "OBC1", "OBC2", "GBn", "GBn2", "HCT"}  # fmt: skip

INTEGRATORS = {
    # name: positional argument roles
    "LangevinMiddleIntegrator": ("temperature", "friction", "dt"),
    "LangevinIntegrator": ("temperature", "friction", "dt"),
    "NoseHooverIntegrator": ("temperature", "collision", "dt"),
    "BrownianIntegrator": ("temperature", "friction", "dt"),
    "VerletIntegrator": ("dt",),
    "VariableLangevinIntegrator": ("temperature", "friction", "error_tolerance"),
    "VariableVerletIntegrator": ("error_tolerance",),
    "DrudeLangevinIntegrator": (
        "temperature",
        "friction",
        "drude_temperature",
        "drude_friction",
        "dt",
    ),
    "DrudeNoseHooverIntegrator": (
        "temperature",
        "collision",
        "drude_temperature",
        "drude_collision",
        "dt",
    ),
}
BAROSTATS = {"MonteCarloBarostat", "MonteCarloAnisotropicBarostat", "MonteCarloMembraneBarostat",
             "MonteCarloFlexibleBarostat"}  # fmt: skip
TOPOLOGY_READERS = {"AmberPrmtopFile": "prmtop", "CharmmPsfFile": "psf", "GromacsTopFile": "top",
                    "PDBFile": "pdb", "PDBxFile": "pdbx"}  # fmt: skip
COORD_READERS = {"AmberInpcrdFile", "CharmmCrdFile", "GromacsGroFile"}
TRAJ_REPORTERS = {
    "DCDReporter": "dcd",
    "XTCReporter": "xtc",
    "PDBReporter": "pdb",
    "PDBxReporter": "pdbx",
}


@dataclass(frozen=True)
class Quantity:
    value: float
    dimension: str | None  # None: plain number

    def to(self, dimension: str) -> float | None:
        return self.value if self.dimension == dimension else None


@dataclass
class Template:
    """A file name, possibly with f-string placeholders (``prod_{i}.dcd``)."""

    text: str
    placeholders: list[str] = field(default_factory=list)  # source of each placeholder expression

    @property
    def regex(self) -> re.Pattern[str]:
        parts = re.split(r"(\{\d+\})", self.text)
        out = []
        for p in parts:
            m = re.fullmatch(r"\{(\d+)\}", p)
            out.append(f"(?P<g{m.group(1)}>.+?)" if m else re.escape(p))
        return re.compile("^" + "".join(out) + "$")

    def match(self, name: str) -> dict[str, str] | None:
        """Placeholder expression -> value, if ``name`` fits this template."""
        m = self.regex.match(name)
        if m is None:
            return None
        return {self.placeholders[int(k[1:])]: v for k, v in m.groupdict().items()}

    def fill(self, values: dict[str, str]) -> str | None:
        """Substitute known placeholder values; expressions like ``i-1`` are evaluated."""
        text = self.text
        for i, expr in enumerate(self.placeholders):
            val = _eval_expr(expr, values)
            if val is None:
                return None
            text = text.replace("{" + str(i) + "}", str(val))
        return text

    @property
    def is_literal(self) -> bool:
        return not self.placeholders


def _eval_expr(expr: str, values: dict[str, str]) -> Any:
    if expr in values:
        return values[expr]
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError:
        return None
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    env: dict[str, Any] = {}
    for n in names:
        if n not in values:
            return None
        try:
            env[n] = int(values[n])
        except ValueError:
            env[n] = values[n]
    try:
        result = _arith(tree.body, env)
    except (ValueError, TypeError, ZeroDivisionError):
        return None
    if isinstance(result, float) and result.is_integer():
        result = int(result)
    return result


def _arith(node: ast.AST, env: dict[str, Any]) -> Any:
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return env[node.id]
    if isinstance(node, ast.BinOp):
        a, b = _arith(node.left, env), _arith(node.right, env)
        ops = {ast.Add: lambda: a + b, ast.Sub: lambda: a - b, ast.Mult: lambda: a * b,
               ast.Div: lambda: a / b, ast.FloorDiv: lambda: a // b, ast.Mod: lambda: a % b}  # fmt: skip
        if type(node.op) in ops:
            return ops[type(node.op)]()
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_arith(node.operand, env)
    raise ValueError("unsupported expression")


@dataclass(frozen=True)
class _OpenFile:
    template: Template


@dataclass
class Event:
    index: int  # statement order
    kind: str
    data: dict[str, Any] = field(default_factory=dict)
    loop: tuple[str, int | None] | None = None  # (loop variable, iterations) when inside a loop


@dataclass
class Segment:
    """One StateDataReporter's worth of simulation."""

    log: Template | None  # None: reported to stdout
    report_interval: int | None
    steps: int = 0
    in_loop: bool = False
    trajectories: list[tuple[str, Template, int | None]] = field(
        default_factory=list
    )  # (fmt, template, interval)
    checkpoints_out: list[Template] = field(default_factory=list)
    saves: list[Template] = field(default_factory=list)
    loads: list[tuple[Template, bool]] = field(
        default_factory=list
    )  # (file, loaded inside the loop)
    continues_previous: bool = False  # no load between the previous segment and this one
    integrator: dict[str, Any] = field(default_factory=dict)
    barostat: dict[str, Any] | None = None
    thermostat: dict[str, Any] | None = None
    temperature_changes: list[float] = field(default_factory=list)
    new_velocities: bool = False
    minimized: bool = False
    restraints: bool = False


@dataclass
class ScriptInfo:
    path: str
    force_fields: list[str] = field(default_factory=list)
    charmm_parameters: list[str] = field(default_factory=list)
    template_generators: list[str] = field(default_factory=list)  # e.g. "openff-2.1.0", "gaff-2.11"
    topology: list[tuple[str, Template]] = field(default_factory=list)  # (format, file)
    coordinates: list[Template] = field(default_factory=list)
    deserialized: list[Template] = field(
        default_factory=list
    )  # XmlSerializer.deserialize(open(f).read())
    system_options: dict[str, Any] = field(default_factory=dict)
    platform: str | None = None
    modeller: bool = (
        False  # the system was built/changed with Modeller (solvent, hydrogens, membrane)
    )
    written_structures: list[Template] = field(
        default_factory=list
    )  # PDBFile/PDBxFile.writeFile targets
    segments: list[Segment] = field(default_factory=list)
    events: list[Event] = field(default_factory=list)
    is_openmm: bool = False


class _Walker:
    def __init__(self, path: str):
        self.info = ScriptInfo(path=path)
        self.env: dict[str, Any] = {}
        self.index = 0
        self.loop: tuple[str, int | None] | None = None

    # -------------------------------------------------------------- values
    def value(self, node: ast.AST | None) -> Any:
        if node is None:
            return None
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id in self.env:
                return self.env[node.id]
            if node.id in UNITS:
                return Quantity(UNITS[node.id][1], UNITS[node.id][0])
            return node.id if node.id in CONSTANTS else None
        if isinstance(node, ast.Attribute):
            if node.attr in UNITS:
                return Quantity(UNITS[node.attr][1], UNITS[node.attr][0])
            return node.attr if node.attr in CONSTANTS else None
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            v = self.value(node.operand)
            return -v if isinstance(v, (int, float)) else None
        if isinstance(node, ast.BinOp):
            a, b = self.value(node.left), self.value(node.right)
            return _combine(a, b, node.op)
        if isinstance(node, ast.Call):
            name = _call_name(node)
            if name == "Quantity" and len(node.args) == 2:
                return _combine(self.value(node.args[0]), self.value(node.args[1]), ast.Mult())
            if name in ("int", "float") and node.args:
                v = self.value(node.args[0])
                return v if isinstance(v, (int, float)) else None
        return None

    def template(self, node: ast.AST | None) -> Template | None:
        if node is None:
            return None
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return Template(node.value)
        if isinstance(node, ast.Name) and isinstance(self.env.get(node.id), (str, Template)):
            v = self.env[node.id]
            return v if isinstance(v, Template) else Template(v)
        if isinstance(node, ast.JoinedStr):
            text, ph = "", []
            for part in node.values:
                if isinstance(part, ast.Constant):
                    text += str(part.value)
                elif isinstance(part, ast.FormattedValue):
                    v = self.value(part.value)
                    if isinstance(v, (str, int)) and part.format_spec is None:
                        text += str(v)
                    else:
                        text += "{" + str(len(ph)) + "}"
                        ph.append(ast.unparse(part.value))
            return Template(text, ph)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            a, b = self.template(node.left), self.template(node.right)
            if a is None or b is None:
                return None
            shift = len(a.placeholders)
            btext = re.sub(r"\{(\d+)\}", lambda m: "{" + str(int(m.group(1)) + shift) + "}", b.text)
            return Template(a.text + btext, a.placeholders + b.placeholders)
        if isinstance(node, ast.Call) and _call_name(node) == "join" and len(node.args) >= 2:
            parts = [self.template(a) for a in node.args]
            if all(parts):
                out = parts[0]
                for p in parts[1:]:
                    shift = len(out.placeholders)
                    ptext = re.sub(
                        r"\{(\d+)\}",
                        lambda m, s=shift: "{" + str(int(m.group(1)) + s) + "}",
                        p.text,
                    )
                    out = Template(
                        out.text.rstrip("/") + "/" + ptext, out.placeholders + p.placeholders
                    )
                return out
        return None

    # -------------------------------------------------------------- walking
    def event(self, kind: str, **data: Any) -> None:
        self.info.events.append(Event(self.index, kind, data, self.loop))

    def walk(self, body: list[ast.stmt]) -> None:
        for stmt in body:
            self.index += 1
            if isinstance(stmt, ast.Assign):
                self.calls(stmt.value)
                val = self.value(stmt.value)
                tmpl = self.template(stmt.value)
                for target in stmt.targets:
                    if isinstance(target, ast.Name):
                        if tmpl is not None and not isinstance(val, (int, float, Quantity)):
                            self.env[target.id] = tmpl if tmpl.placeholders else tmpl.text
                        elif val is not None:
                            self.env[target.id] = val
                        else:
                            self.env.pop(target.id, None)
                continue
            if isinstance(stmt, (ast.For, ast.AsyncFor)):
                var = stmt.target.id if isinstance(stmt.target, ast.Name) else None
                iters = self._iterations(stmt.iter)
                outer = self.loop
                self.loop = (var or "?", iters)
                if var:
                    self.env.pop(var, None)
                self.walk(stmt.body)
                self.loop = outer
                continue
            if isinstance(stmt, ast.While):
                outer = self.loop
                self.loop = ("<while>", None)
                self.walk(stmt.body)
                self.loop = outer
                continue
            if isinstance(stmt, ast.If):
                self.walk(stmt.body)
                self.walk(stmt.orelse)
                continue
            if isinstance(stmt, (ast.With, ast.Try)):
                if isinstance(stmt, ast.With):
                    for item in stmt.items:
                        self.calls(item.context_expr)
                        ce = item.context_expr
                        if (isinstance(ce, ast.Call) and _call_name(ce) == "open" and ce.args
                                and isinstance(item.optional_vars, ast.Name)):  # fmt: skip
                            t = self.template(ce.args[0])
                            if t is not None:
                                self.env[item.optional_vars.id] = _OpenFile(t)
                self.walk(stmt.body)
                if isinstance(stmt, ast.Try):
                    for h in stmt.handlers:
                        self.walk(h.body)
                    self.walk(stmt.orelse)
                    self.walk(stmt.finalbody)
                continue
            if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue  # not executed at import; simulate() helpers are rare in run scripts
            if isinstance(stmt, (ast.Import, ast.ImportFrom)):
                mods = [a.name for a in stmt.names] + (
                    [stmt.module] if isinstance(stmt, ast.ImportFrom) and stmt.module else []
                )
                if any(m and (m.startswith("openmm") or m.startswith("simtk")) for m in mods):
                    self.info.is_openmm = True
                continue
            for node in ast.iter_child_nodes(stmt):
                self.calls(node)

    def _iterations(self, node: ast.AST) -> int | None:
        if isinstance(node, ast.Call) and _call_name(node) == "range":
            args = [self.value(a) for a in node.args]
            if all(isinstance(a, int) for a in args) and args:
                start, stop, step = (
                    (0, args[0], 1)
                    if len(args) == 1
                    else (args[0], args[1], args[2] if len(args) > 2 else 1)
                )
                return max(0, (stop - start + (step - 1 if step > 0 else step + 1)) // step)
        if isinstance(node, (ast.List, ast.Tuple)):
            return len(node.elts)
        return None

    def calls(self, node: ast.AST) -> None:
        if node is None:
            return
        for sub in ast.walk(node):
            if isinstance(sub, ast.Call):
                self.call(sub)

    def call(self, node: ast.Call) -> None:
        name = _call_name(node)
        kw = {k.arg: k.value for k in node.keywords if k.arg}
        args = node.args
        if name == "ForceField":
            self.info.force_fields += [t.text for t in (self.template(a) for a in args) if t]
        elif name == "CharmmParameterSet":
            self.info.charmm_parameters += [t.text for t in (self.template(a) for a in args) if t]
        elif name in (
            "SMIRNOFFTemplateGenerator",
            "GAFFTemplateGenerator",
            "EspalomaTemplateGenerator",
        ):
            ff = self.value(kw.get("forcefield"))
            self.info.template_generators.append(str(ff) if ff else name)
        elif name in TOPOLOGY_READERS and args:
            t = self.template(args[0])
            if t:
                self.info.topology.append((TOPOLOGY_READERS[name], t))
        elif name in COORD_READERS and args:
            t = self.template(args[0])
            if t:
                self.info.coordinates.append(t)
        elif name == "deserialize" and args:
            inner = args[0]
            if (
                isinstance(inner, ast.Call)
                and _call_name(inner) == "read"
                and isinstance(inner.func, ast.Attribute)
            ):
                target = inner.func.value
                if isinstance(target, ast.Call) and _call_name(target) == "open" and target.args:
                    t = self.template(target.args[0])
                    if t:
                        self.info.deserialized.append(t)
                elif isinstance(target, ast.Name) and isinstance(
                    self.env.get(target.id), _OpenFile
                ):
                    self.info.deserialized.append(self.env[target.id].template)
        elif (
            name in ("addSolvent", "addHydrogens", "addMembrane", "addExtraParticles")
            or name == "Modeller"
        ):
            if name != "Modeller":
                self.info.modeller = True
        elif name == "writeFile" and len(args) >= 3:
            target = args[2]
            t = None
            if isinstance(target, ast.Name) and isinstance(self.env.get(target.id), _OpenFile):
                t = self.env[target.id].template
            elif isinstance(target, ast.Call) and _call_name(target) == "open" and target.args:
                t = self.template(target.args[0])
            if t is not None:
                self.info.written_structures.append(t)
        elif name == "createSystem":
            for key in ("nonbondedMethod", "nonbondedCutoff", "constraints", "hydrogenMass", "rigidWater",
                        "implicitSolvent", "ewaldErrorTolerance"):  # fmt: skip
                if key in kw:
                    self.info.system_options[key] = self.value(kw[key])
            self.info.is_openmm = True
        elif name in INTEGRATORS:
            roles = INTEGRATORS[name]
            data: dict[str, Any] = {"type": name}
            for role, a in zip(roles, args, strict=False):
                data[role] = self.value(a)
            for k, v in kw.items():
                data[k] = self.value(v)
            self.event("integrator", **data)
        elif name == "CustomExternalForce" or (name == "addForce" and args and isinstance(args[0], ast.Name)
                                                and "restr" in args[0].id.lower()):  # fmt: skip
            self.event("restraint")
        elif name == "AndersenThermostat":
            self.event("thermostat", type=name, temperature=self.value(args[0]) if args else None)
        elif name in BAROSTATS:
            self.event("barostat", type=name, pressure=self.value(args[0]) if args else None,
                       temperature=self.value(args[1]) if len(args) > 1 else None)  # fmt: skip
        elif name == "setTemperature" and args:
            self.event("set_temperature", temperature=self.value(args[0]))
        elif name == "StateDataReporter":
            target = args[0] if args else kw.get("file")
            t = None
            if not (
                isinstance(target, (ast.Attribute, ast.Name))
                and ast.unparse(target) in ("stdout", "sys.stdout")
            ):
                t = self.template(target)
            interval = self.value(args[1] if len(args) > 1 else kw.get("reportInterval"))
            self.event("log", file=t, interval=interval if isinstance(interval, int) else None,
                       stdout=t is None)  # fmt: skip
        elif name in TRAJ_REPORTERS and args:
            interval = self.value(args[1] if len(args) > 1 else kw.get("reportInterval"))
            self.event("trajectory", fmt=TRAJ_REPORTERS[name], file=self.template(args[0]),
                       interval=interval if isinstance(interval, int) else None)  # fmt: skip
        elif name == "CheckpointReporter" and args:
            self.event("checkpoint_out", file=self.template(args[0]))
        elif name == "step" and isinstance(node.func, ast.Attribute) and args:
            n = self.value(args[0])
            self.event("step", n=n if isinstance(n, int) else None)
        elif name == "minimizeEnergy":
            self.event("minimize")
        elif name in ("loadState", "loadCheckpoint") and args:
            self.event("load", file=self.template(args[0]), how=name)
        elif name in ("saveState", "saveCheckpoint") and args:
            self.event("save", file=self.template(args[0]), how=name)
        elif name == "setVelocitiesToTemperature":
            self.event("velocities")
        elif name == "getPlatformByName" and args:
            p = self.value(args[0])
            self.info.platform = p if isinstance(p, str) else None
        elif name == "Simulation":
            self.info.is_openmm = True


def _call_name(node: ast.Call) -> str | None:
    f = node.func
    if isinstance(f, ast.Name):
        return f.id
    if isinstance(f, ast.Attribute):
        return f.attr
    return None


def _combine(a: Any, b: Any, op: ast.operator) -> Any:
    def q(x: Any) -> Quantity | None:
        if isinstance(x, Quantity):
            return x
        if isinstance(x, (int, float)) and not isinstance(x, bool):
            return Quantity(float(x), None)
        return None

    qa, qb = q(a), q(b)
    if qa is None or qb is None:
        return None
    if isinstance(op, ast.Mult):
        dim = qa.dimension or qb.dimension
        if qa.dimension is None and qb.dimension is None:
            return a * b
        return Quantity(qa.value * qb.value, dim)
    if isinstance(op, ast.Div):
        if qb.dimension == "time" and qa.dimension is None:
            return Quantity(qa.value / qb.value, "rate")  # 1/picosecond
        if qb.dimension is None:
            if qa.dimension is None:
                return a / b
            return Quantity(qa.value / qb.value, qa.dimension)
        return None
    if isinstance(op, (ast.Add, ast.Sub)) and qa.dimension == qb.dimension:
        if qa.dimension is None:
            return a + b if isinstance(op, ast.Add) else a - b
        v = qa.value + qb.value if isinstance(op, ast.Add) else qa.value - qb.value
        return Quantity(v, qa.dimension)
    return None


def _segments(info: ScriptInfo) -> list[Segment]:
    """Group the event stream into StateDataReporter segments."""
    segments: list[Segment] = []
    integrator: dict[str, Any] = {}
    barostat: dict[str, Any] | None = None
    thermostat: dict[str, Any] | None = None
    restraints = False
    pending_traj: list[tuple[str, Template, int | None]] = []
    pending_ckpt: list[Template] = []
    pending_loads: list[tuple[Template, bool]] = []
    pending_velocities = False
    pending_min = False
    loaded_since_last_step = False
    cur: Segment | None = None
    stepped_since_log = False
    for ev in info.events:
        d = ev.data
        if ev.kind == "integrator":
            integrator = dict(d)
        elif ev.kind == "barostat":
            barostat = dict(d)
        elif ev.kind == "thermostat":
            thermostat = dict(d)
        elif ev.kind == "restraint":
            restraints = True
        elif ev.kind == "trajectory" and d.get("file") is not None:
            if cur is not None and not stepped_since_log:
                cur.trajectories.append((d["fmt"], d["file"], d["interval"]))
            else:
                pending_traj.append((d["fmt"], d["file"], d["interval"]))
        elif ev.kind == "checkpoint_out" and d.get("file") is not None:
            if cur is not None and not stepped_since_log:
                cur.checkpoints_out.append(d["file"])
            else:
                pending_ckpt.append(d["file"])
        elif ev.kind == "log":
            cur = Segment(
                log=d.get("file"), report_interval=d.get("interval"), in_loop=ev.loop is not None
            )
            cur.trajectories, pending_traj = pending_traj, []
            cur.checkpoints_out, pending_ckpt = pending_ckpt, []
            cur.loads, pending_loads = pending_loads, []
            cur.continues_previous = bool(segments) and not loaded_since_last_step
            cur.new_velocities, pending_velocities = pending_velocities, False
            cur.minimized, pending_min = pending_min, False
            segments.append(cur)
            stepped_since_log = False
            loaded_since_last_step = False
        elif ev.kind == "load" and d.get("file") is not None:
            pending_loads.append((d["file"], ev.loop is not None))
            loaded_since_last_step = True
            if cur is not None and not stepped_since_log:
                cur.loads.append((d["file"], ev.loop is not None))
        elif ev.kind == "velocities":
            pending_velocities = True
            if cur is not None and not stepped_since_log:
                cur.new_velocities = True
        elif ev.kind == "minimize":
            pending_min = True
        elif ev.kind == "set_temperature":
            t = d.get("temperature")
            if cur is not None and isinstance(t, Quantity):
                cur.temperature_changes.append(t.value)
        elif ev.kind == "step" and cur is not None:
            n = d.get("n") or 0
            loop_iters = ev.loop[1] if ev.loop else 1
            cur.steps += n * (loop_iters or 1) if (ev.loop and not cur.in_loop) else n
            cur.integrator = dict(integrator)
            cur.barostat = dict(barostat) if barostat else None
            cur.thermostat = dict(thermostat) if thermostat else None
            cur.restraints = restraints
            stepped_since_log = True
            loaded_since_last_step = False
        elif ev.kind == "save" and cur is not None and d.get("file") is not None:
            cur.saves.append(d["file"])
    return segments


def read_script(path: str | Path) -> ScriptInfo:
    with open_text(path) as fh:
        source = fh.read()
    tree = ast.parse(source, filename=str(path))
    walker = _Walker(str(path))
    walker.walk(tree.body)
    info = walker.info
    info.segments = _segments(info)
    return info


def looks_like_script(head: bytes, name: str) -> bool:
    if not name.lower().endswith(".py"):
        return False
    text = head.decode("latin-1", errors="replace")
    return ("openmm" in text or "simtk" in text) and (
        "Simulation" in text or "createSystem" in text or "Integrator" in text or "Reporter" in text
    )
