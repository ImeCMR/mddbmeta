"""Write an AMBER prmtop from a MELD run's ``Data/system.dat`` -- without MELD.

``system.dat`` is a pickle of MELD's ``System``. It is read with :class:`InertLoader`, which never
imports a module the file names: every class becomes an inert placeholder that only holds the saved
values, so no code from MELD, OpenMM or the file runs. Only numpy's array reconstruction is real
(it rebuilds plain number arrays). From the placeholders we take

* the OpenMM ``Topology`` (chains, residues, atoms, elements, bonds), rebuilt as a real one, and
* the OpenMM ``System``, which pickles itself as its XML serialization.

The system is then written by ParmEd. Grappa splits the bonded terms over several forces (two
``HarmonicBondForce``, a zeroed ff14SB ``PeriodicTorsionForce`` plus Grappa propers and impropers);
ParmEd would read a second bond force as CHARMM Urey-Bradley terms, so the forces are merged first.
Implicit solvent (``CustomGBForce``/``GBSAOBCForce``) has no prmtop form and is left out.

Every conversion is checked: the bond, angle, torsion and nonbonded energies of the written prmtop
must match the MELD system's at the saved coordinates, or nothing is written.

The saved coordinates (MELD's template coordinates, the run's starting structure) can also be
written as a PDB with the same topology, so its atoms match the prmtop one to one.
"""

from __future__ import annotations

import collections
import math
import os
import pickle
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta.errors import MddbmetaError

# forces a prmtop cannot hold but that are not part of the molecule's force field
DROPPABLE = {
    "CMMotionRemover": "centre-of-mass motion remover",
    "CustomGBForce": "implicit solvent (GB)",
    "GBSAOBCForce": "implicit solvent (GB)",
    "MonteCarloBarostat": "barostat",
    "MonteCarloAnisotropicBarostat": "barostat",
    "MonteCarloMembraneBarostat": "barostat",
    "AndersenThermostat": "thermostat",
}
DEFAULT_NAME = (
    "meld_system.prmtop"  # not topology.*: MDDB-workflow writes its own topology.* outputs
)
TITLE = "mddbmeta meld-prmtop from Data/system.dat"
PDB_NAME = "meld_system.pdb"  # the starting structure, for runs that have no PDB of their own
PDB_MARK = "REMARK   1 mddbmeta: starting structure from Data/system.dat"
ENERGY_TOL = 1e-3  # kJ/mol, absolute ...
ENERGY_RTOL = 1e-6  # ... plus relative


class _Inert:
    """Placeholder for any class a pickle names; it only stores what the pickle hands it."""

    _qualname = "?"

    def __new__(cls, *args: Any, **kwargs: Any) -> _Inert:
        obj = object.__new__(cls)
        obj.__dict__["_args"] = args
        return obj

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        if args:
            self.__dict__["_args"] = args

    def __setstate__(self, state: Any) -> None:
        if isinstance(state, tuple) and len(state) == 2 and isinstance(state[1], dict):
            state = {**(state[0] or {}), **state[1]}  # (dict, slots) form
        if isinstance(state, dict):
            self.__dict__.update(state)
        else:
            self.__dict__["_state"] = state

    def __call__(self, *args: Any, **kwargs: Any) -> _Inert:  # a placeholder function
        return _placeholder(f"{self._qualname}()")(*args)

    def __repr__(self) -> str:
        return f"<inert {self._qualname}>"


def _placeholder(qualname: str) -> type:
    return type(qualname.rsplit(".", 1)[-1], (_Inert,), {"_qualname": qualname})


def _numpy_allowlist() -> dict[tuple[str, str], Any]:
    import numpy as np

    try:
        from numpy._core.multiarray import _reconstruct, scalar
    except ImportError:  # numpy < 2
        from numpy.core.multiarray import _reconstruct, scalar  # type: ignore[no-redef]
    allowed: dict[tuple[str, str], Any] = {}
    for mod in ("numpy.core.multiarray", "numpy._core.multiarray"):
        allowed[(mod, "_reconstruct")] = _reconstruct
        allowed[(mod, "scalar")] = scalar
    allowed[("numpy", "ndarray")] = np.ndarray
    allowed[("numpy", "dtype")] = np.dtype
    return allowed


class InertLoader(pickle.Unpickler):
    """An unpickler that imports nothing: unknown classes become :class:`_Inert` placeholders."""

    def __init__(self, fh: Any) -> None:
        super().__init__(fh)
        self._allowed = _numpy_allowlist()
        self._allowed[("collections", "OrderedDict")] = collections.OrderedDict

    def find_class(self, module: str, name: str) -> Any:
        real = self._allowed.get((module, name))
        return real if real is not None else _placeholder(f"{module}.{name}")


def load_inert(path: str | Path) -> Any:
    with open(path, "rb") as fh:
        try:
            return InertLoader(fh).load()
        except Exception as exc:  # corrupt or not a pickle
            raise MddbmetaError(f"{path}: not a readable MELD system pickle ({exc})") from exc


@dataclass
class MeldSystem:
    topology: Any  # openmm.app.Topology
    system: Any  # openmm.System
    positions: Any | None  # numpy array (nm) or None
    source: str


def _qual(obj: Any) -> str:
    return getattr(type(obj), "_qualname", "")


def _find(spec: Any, want: str) -> Any | None:
    """The first attribute of the MELD System whose class name ends with ``want``."""
    for value in vars(spec).values() if hasattr(spec, "__dict__") else ():
        if _qual(value).endswith(want):
            return value
    return None


def _system_xml(obj: Any) -> str | None:
    state = vars(obj).get("_state") if hasattr(obj, "__dict__") else None
    if isinstance(state, str) and "<System" in state[:500]:
        return state
    for value in vars(obj).values() if hasattr(obj, "__dict__") else ():
        if isinstance(value, str) and "<System" in value[:500]:
            return value
    return None


def read_meld_system(path: str | Path) -> MeldSystem:
    """Rebuild the OpenMM Topology and System saved in a MELD ``system.dat``."""
    import numpy as np
    import openmm
    import openmm.app as app

    spec = load_inert(path)
    if not _qual(spec).endswith(".System") or "meld" not in _qual(spec):
        raise MddbmetaError(f"{path}: not a MELD System pickle (found {_qual(spec) or type(spec)})")
    otop = _find(spec, "openmm.app.topology.Topology")
    osys = next((v for v in vars(spec).values() if _system_xml(v)), None) if otop else None
    if otop is None or osys is None:
        raise MddbmetaError(
            f"{path}: no OpenMM topology/system inside (this MELD version stores them differently)"
        )
    system = openmm.XmlSerializer.deserialize(_system_xml(osys))

    top = app.Topology()
    atoms: dict[int, Any] = {}
    try:
        for c in otop._chains:
            chain = top.addChain(str(getattr(c, "id", "A")))
            for r in c._residues:
                res = top.addResidue(r.name, chain, str(getattr(r, "id", "")))
                for a in r._atoms:
                    el = getattr(a, "element", None)
                    sym = el._args[0] if isinstance(el, _Inert) and el._args else None
                    element = app.Element.getBySymbol(sym) if sym else None
                    atoms[a.index] = top.addAtom(a.name, element, res, str(getattr(a, "id", "")))
        for b in otop._bonds:
            a1, a2 = (b._args[0], b._args[1]) if b._args else (b.atom1, b.atom2)
            top.addBond(atoms[a1.index], atoms[a2.index])
    except (AttributeError, KeyError, IndexError, TypeError) as exc:
        raise MddbmetaError(f"{path}: unexpected OpenMM topology layout ({exc})") from exc
    box = getattr(otop, "_periodicBoxVectors", None)
    if box is not None and system.usesPeriodicBoundaryConditions():
        top.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
    if top.getNumAtoms() != system.getNumParticles():
        raise MddbmetaError(
            f"{path}: topology has {top.getNumAtoms()} atoms but the system {system.getNumParticles()}"
        )

    positions = None
    for key in ("_template_coordinates", "template_coordinates", "_coordinates"):
        value = vars(spec).get(key)
        if isinstance(value, np.ndarray) and value.shape == (system.getNumParticles(), 3):
            positions = value
            break
    return MeldSystem(topology=top, system=system, positions=positions, source=str(path))


def merge_forces(system: Any) -> tuple[Any, list[str]]:
    """One force per bonded type (as AMBER has), zeroed torsions dropped; returns (system, dropped)."""
    import openmm
    from openmm import unit

    merged = openmm.System()
    for i in range(system.getNumParticles()):
        merged.addParticle(system.getParticleMass(i))
    for i in range(system.getNumConstraints()):
        merged.addConstraint(*system.getConstraintParameters(i))
    if system.usesPeriodicBoundaryConditions():
        merged.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    bonds, angles, torsions = (
        openmm.HarmonicBondForce(),
        openmm.HarmonicAngleForce(),
        openmm.PeriodicTorsionForce(),
    )
    dropped: list[str] = []
    for f in system.getForces():
        name = type(f).__name__
        if isinstance(f, openmm.HarmonicBondForce):
            for j in range(f.getNumBonds()):
                bonds.addBond(*f.getBondParameters(j))
        elif isinstance(f, openmm.HarmonicAngleForce):
            for j in range(f.getNumAngles()):
                angles.addAngle(*f.getAngleParameters(j))
        elif isinstance(f, openmm.PeriodicTorsionForce):
            for j in range(f.getNumTorsions()):
                p = f.getTorsionParameters(j)
                if p[-1].value_in_unit(unit.kilojoule_per_mole) != 0.0:
                    torsions.addTorsion(*p)
        elif isinstance(f, (openmm.NonbondedForce, openmm.CMAPTorsionForce)):
            merged.addForce(_copy(f))
        elif name in DROPPABLE:
            dropped.append(name)
        else:
            raise MddbmetaError(
                f"the MELD system has a {name}, which a prmtop cannot hold; "
                "the prmtop would not describe the simulated system"
            )
    for f in (bonds, angles, torsions):
        merged.addForce(f)
    return merged, dropped


def _copy(force: Any) -> Any:
    import openmm

    return openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(force))


TERMS = ("bond", "angle", "torsion", "cmap", "nonbonded")


def _term(force: Any) -> str | None:
    import openmm

    return {
        openmm.HarmonicBondForce: "bond",
        openmm.HarmonicAngleForce: "angle",
        openmm.PeriodicTorsionForce: "torsion",
        openmm.CMAPTorsionForce: "cmap",
        openmm.NonbondedForce: "nonbonded",
    }.get(type(force))


def term_energies(system: Any, positions: Any) -> dict[str, float]:
    """Energy (kJ/mol) per force-field term, on the CPU reference platform."""
    import openmm
    from openmm import unit

    system = _copy_system(system)
    groups: dict[str, int] = {}
    for f in system.getForces():
        term = _term(f)
        if term is None:
            f.setForceGroup(31)
            continue
        f.setForceGroup(groups.setdefault(term, len(groups)))
    ctx = openmm.Context(
        system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference")
    )
    ctx.setPositions(positions)
    out = {}
    for term, g in groups.items():
        state = ctx.getState(getEnergy=True, groups={g})
        out[term] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    del ctx
    return out


def _copy_system(system: Any) -> Any:
    import openmm

    return openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))


def _prmtop_system(path: str, like: Any) -> Any:
    """The prmtop as an OpenMM System with the MELD system's nonbonded settings."""
    import openmm
    import openmm.app as app

    nb = next((f for f in like.getForces() if isinstance(f, openmm.NonbondedForce)), None)
    methods = {
        openmm.NonbondedForce.NoCutoff: app.NoCutoff,
        openmm.NonbondedForce.CutoffNonPeriodic: app.CutoffNonPeriodic,
        openmm.NonbondedForce.CutoffPeriodic: app.CutoffPeriodic,
        openmm.NonbondedForce.Ewald: app.Ewald,
        openmm.NonbondedForce.PME: app.PME,
        openmm.NonbondedForce.LJPME: app.LJPME,
    }
    kwargs: dict[str, Any] = {"removeCMMotion": False}
    if nb is not None:
        kwargs["nonbondedMethod"] = methods[nb.getNonbondedMethod()]
        if nb.getNonbondedMethod() != openmm.NonbondedForce.NoCutoff:
            kwargs["nonbondedCutoff"] = nb.getCutoffDistance()
    prm = app.AmberPrmtopFile(path)
    if like.usesPeriodicBoundaryConditions():
        prm.topology.setPeriodicBoxVectors(like.getDefaultPeriodicBoxVectors())
    system = prm.createSystem(**kwargs)
    if nb is not None:
        mine = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        mine.setReactionFieldDielectric(nb.getReactionFieldDielectric())
        mine.setUseDispersionCorrection(nb.getUseDispersionCorrection())
        mine.setUseSwitchingFunction(nb.getUseSwitchingFunction())
        if nb.getUseSwitchingFunction():
            mine.setSwitchingDistance(nb.getSwitchingDistance())
        mine.setEwaldErrorTolerance(nb.getEwaldErrorTolerance())
        mine.setPMEParameters(*nb.getPMEParameters())
    if like.usesPeriodicBoundaryConditions():
        system.setDefaultPeriodicBoxVectors(*like.getDefaultPeriodicBoxVectors())
    return system


@dataclass
class PrmtopReport:
    output: str
    source: str
    natom: int
    nres: int
    nbond: int
    dropped: list[str]
    energies: dict[str, tuple[float, float]] = field(default_factory=dict)  # term -> (meld, prmtop)
    checked: bool = True


def find_system_dat(target: Path) -> Path:
    if target.is_file():
        return target
    for cand in (target / "Data" / "system.dat", target / "system.dat"):
        if cand.is_file():
            return cand
    raise MddbmetaError(f"no Data/system.dat under {target} (is it a MELD run directory?)")


def write_prmtop(target: str | Path, output: str | Path | None = None, force: bool = False
                 ) -> PrmtopReport:  # fmt: skip
    """Convert ``target`` (a MELD run directory or its ``system.dat``) to an AMBER prmtop."""
    import parmed

    src = find_system_dat(Path(target).expanduser())
    run_dir = src.parent.parent if src.parent.name == "Data" else src.parent
    out = Path(output).expanduser() if output else run_dir / DEFAULT_NAME
    if out.exists() and not force:
        raise MddbmetaError(f"{out} exists; pass --force to overwrite it")
    if not out.parent.is_dir():
        raise MddbmetaError(f"output directory {out.parent} does not exist")

    meld = read_meld_system(src)
    merged, dropped = merge_forces(meld.system)
    struct = parmed.openmm.load_topology(meld.topology, merged)
    if struct.urey_bradleys or struct.impropers or struct.rb_torsions:
        raise MddbmetaError(
            "the MELD system has CHARMM-style terms; it cannot be written as an AMBER prmtop"
        )

    fd, tmp = tempfile.mkstemp(prefix=".mddbmeta-", suffix=".prmtop", dir=out.parent)
    os.close(fd)
    try:
        parm = parmed.amber.AmberParm.from_structure(struct)
        parm.parm_data["TITLE"] = [TITLE]
        parm.write_parm(tmp)
        report = PrmtopReport(output=str(out), source=str(src), natom=len(struct.atoms),
                              nres=len(struct.residues), nbond=len(struct.bonds), dropped=dropped)  # fmt: skip
        if meld.positions is None:
            report.checked = False
        else:
            want = term_energies(merged, meld.positions)
            got = term_energies(_prmtop_system(tmp, merged), meld.positions)
            bad = []
            for term in TERMS:
                if term not in want and term not in got:
                    continue
                a, b = want.get(term, 0.0), got.get(term, 0.0)
                report.energies[term] = (a, b)
                if not math.isclose(a, b, rel_tol=ENERGY_RTOL, abs_tol=ENERGY_TOL):
                    bad.append(f"{term} {a:.3f} vs {b:.3f} kJ/mol")
            if bad:
                raise MddbmetaError(
                    "the prmtop does not reproduce the MELD system ("
                    + "; ".join(bad)
                    + "); not written"
                )
        _publish(tmp, out)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return report


def _publish(tmp: str, out: Path) -> None:
    """Move a finished temp file into place, with the permissions a new file would get."""
    umask = os.umask(0)
    os.umask(umask)
    os.chmod(tmp, 0o666 & ~umask)  # mkstemp makes 0600
    os.replace(tmp, out)


def write_structure(target: str | Path, output: str | Path | None = None, force: bool = False
                    ) -> tuple[Path, int]:  # fmt: skip
    """Write the starting structure saved in a MELD ``system.dat`` as a PDB; returns (path, natom)."""
    import io

    import openmm.app as app
    import openmm.unit as unit

    src = find_system_dat(Path(target).expanduser())
    run_dir = src.parent.parent if src.parent.name == "Data" else src.parent
    out = Path(output).expanduser() if output else run_dir / PDB_NAME
    if out.exists() and not force:
        raise MddbmetaError(f"{out} exists; pass --force to overwrite it")
    if not out.parent.is_dir():
        raise MddbmetaError(f"output directory {out.parent} does not exist")
    meld = read_meld_system(src)
    if meld.positions is None:
        raise MddbmetaError(f"{src} holds no starting coordinates")
    buf = io.StringIO()
    app.PDBFile.writeFile(meld.topology, meld.positions * unit.nanometer, buf)
    fd, tmp = tempfile.mkstemp(prefix=".mddbmeta-", suffix=".pdb", dir=out.parent)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(PDB_MARK + "\n" + buf.getvalue())
        _publish(tmp, out)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return out, meld.topology.getNumAtoms()


def is_generated(path: str | Path, marker: str = TITLE) -> bool:
    """True for a file this module wrote (the prmtop TITLE / PDB REMARK carries the marker)."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(400)
    except OSError:
        return False
    return marker.encode() in head


def ensure_structure(system_dat: Path) -> tuple[Path, bool]:
    """The run's generated starting-structure PDB: reused when newer than ``system_dat``, else
    (re)written. Returns (path, written). Raises MddbmetaError when it cannot be written."""
    run_dir = system_dat.parent.parent if system_dat.parent.name == "Data" else system_dat.parent
    out = run_dir / PDB_NAME
    if out.exists():
        if not is_generated(out, PDB_MARK):
            raise MddbmetaError(f"{out} exists and was not written by mddbmeta; not overwriting it")
        if out.stat().st_mtime >= system_dat.stat().st_mtime:
            return out, False
    try:
        write_structure(system_dat, out, force=True)
    except OSError as exc:
        raise MddbmetaError(f"cannot write {out}: {exc.strerror or exc}") from exc
    return out, True


def ensure_prmtop(system_dat: Path) -> tuple[Path, PrmtopReport | None]:
    """The run's generated prmtop: reused when it is newer than ``system_dat``, else (re)written.

    Returns (path, report); report is None when an up-to-date file was reused.
    Raises MddbmetaError when it cannot be written (a foreign file in the way, a failed check, ...).
    """
    run_dir = system_dat.parent.parent if system_dat.parent.name == "Data" else system_dat.parent
    out = run_dir / DEFAULT_NAME
    if out.exists():
        if not is_generated(out):
            raise MddbmetaError(f"{out} exists and was not written by mddbmeta; not overwriting it")
        if out.stat().st_mtime >= system_dat.stat().st_mtime:
            return out, None
    try:
        return out, write_prmtop(system_dat, out, force=True)
    except OSError as exc:
        raise MddbmetaError(f"cannot write {out}: {exc.strerror or exc}") from exc
