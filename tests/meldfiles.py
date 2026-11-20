"""Write a MELD-like Data/system.dat: a Grappa-like OpenMM system pickled as meld.system.meld_system.System."""

from __future__ import annotations

import contextlib
import pickle
import sys
import types
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import openmm
import openmm.app as app

PDB = Path(__file__).resolve().parent / "data" / "openmm_real" / "prep" / "ala_ala.pdb"
AMBER = (
    Path(__file__).resolve().parent / "data" / "openmm_amber"
)  # ala_ala.prmtop/.inpcrd, 33 atoms
NATOM = 33


def grappa_like(system: openmm.System) -> openmm.System:
    """Split the bonds over two forces and replace the torsions by a zeroed set plus a copy,
    the way MELD's Grappa builder leaves them."""
    out = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))
    bonds = next(f for f in out.getForces() if isinstance(f, openmm.HarmonicBondForce))
    second = openmm.HarmonicBondForce()
    for j in range(bonds.getNumBonds()):
        if j % 2:
            a, b, r, k = bonds.getBondParameters(j)
            second.addBond(a, b, r, k)
            bonds.setBondParameters(j, a, b, r, 0.0)
    out.addForce(second)
    tors = next(f for f in out.getForces() if isinstance(f, openmm.PeriodicTorsionForce))
    copy = openmm.PeriodicTorsionForce()
    for j in range(tors.getNumTorsions()):
        *atoms, per, phase, k = tors.getTorsionParameters(j)
        copy.addTorsion(*atoms, per, phase, k)
        tors.setTorsionParameters(j, *atoms, per, phase, 0.0)
    out.addForce(copy)
    return out


@contextlib.contextmanager
def _fake_meld_modules(module: str = "meld.system.meld_system") -> Iterator[type]:
    """A module tree named like MELD's, so an object pickles as ``<module>.System``."""
    mod = types.ModuleType(module)

    class System:
        pass

    System.__module__ = mod.__name__
    System.__qualname__ = "System"
    mod.System = System  # type: ignore[attr-defined]
    names = ("meld", "meld.system", mod.__name__)
    saved = {n: sys.modules.get(n) for n in names}
    sys.modules.update({"meld": types.ModuleType("meld"), "meld.system": types.ModuleType("meld.system"),
                        mod.__name__: mod})  # fmt: skip
    try:
        yield System
    finally:
        for n, m in saved.items():
            if m is None:
                sys.modules.pop(n, None)
            else:
                sys.modules[n] = m


def write_system_dat(
    run: Path, system: openmm.System | None = None, with_positions: bool = True
) -> Path:
    """Write ``run/Data/system.dat`` for the ala-ala test molecule; returns the path."""
    pdb = app.PDBFile(str(PDB))
    if system is None:
        ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
        system = grappa_like(ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic))
    path = run / "Data" / "system.dat"
    path.parent.mkdir(parents=True, exist_ok=True)
    with _fake_meld_modules() as cls:
        spec = cls()
        spec._openmm_topology = pdb.topology
        spec._openmm_system = system
        if with_positions:
            spec._template_coordinates = np.array(pdb.getPositions(asNumpy=True)._value)
        spec.builder_info = {"builder": "grappa"}
        with open(path, "wb") as fh:
            pickle.dump(spec, fh)
    return path


def write_legacy_system_dat(run: Path, with_coordinates: bool = True) -> Path:
    """``run/Data/system.dat`` as older MELD saved it: ``meld.system.system.System`` holding the
    AMBER prmtop and inpcrd as text (and the coordinates in angstrom), no OpenMM system."""
    path = run / "Data" / "system.dat"
    path.parent.mkdir(parents=True, exist_ok=True)
    crd = (AMBER / "ala_ala.inpcrd").read_text()
    with _fake_meld_modules("meld.system.system") as cls:
        spec = cls()
        spec._top_string = (AMBER / "ala_ala.prmtop").read_text()
        spec._mdcrd_string = crd if with_coordinates else None
        spec._n_atoms = NATOM
        with open(path, "wb") as fh:
            pickle.dump(spec, fh)
    return path
