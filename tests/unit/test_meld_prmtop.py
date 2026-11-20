"""MELD system.dat -> prmtop: inert loading, Grappa force merging, the energy check, the CLI."""

from __future__ import annotations

import os
import pickle
from pathlib import Path

import openmm
import openmm.app as app
import pytest

from mddbmeta.cli.main import main
from mddbmeta.engines.meld.prmtop import load_inert, read_meld_system, write_prmtop
from mddbmeta.errors import MddbmetaError
from meldfiles import PDB, write_system_dat

PROTEIN_L = os.environ.get("MDDBMETA_MELD_RUN")


def _run_dir(
    tmp_path: Path, system: openmm.System | None = None, with_positions: bool = True
) -> Path:
    run = tmp_path / "run"
    write_system_dat(run, system, with_positions)
    return run


def test_grappa_like_system_round_trips_exactly(tmp_path: Path) -> None:
    run = _run_dir(tmp_path)
    report = write_prmtop(run)
    assert Path(report.output) == run / "meld_system.prmtop"
    assert report.checked and report.natom == 33
    assert {"bond", "angle", "torsion", "nonbonded"} <= set(report.energies)
    for meld_e, prm_e in report.energies.values():
        assert meld_e == pytest.approx(prm_e, abs=1e-3)
    assert "CustomGBForce" in report.dropped
    prm = app.AmberPrmtopFile(report.output)
    assert prm.topology.getNumAtoms() == 33
    source = app.PDBFile(str(PDB)).topology  # same atom order and names as the MELD topology
    assert [a.name for a in prm.topology.atoms()] == [a.name for a in source.atoms()]


def test_refuses_to_overwrite(tmp_path: Path) -> None:
    run = _run_dir(tmp_path)
    write_prmtop(run)
    with pytest.raises(MddbmetaError, match="--force"):
        write_prmtop(run)
    write_prmtop(run, force=True)


def test_unsupported_force_is_refused(tmp_path: Path) -> None:
    pdb = app.PDBFile(str(PDB))
    system = app.ForceField("amber14/protein.ff14SB.xml").createSystem(pdb.topology)
    system.addForce(openmm.CustomExternalForce("x"))
    run = _run_dir(tmp_path, system=system)
    with pytest.raises(MddbmetaError, match="CustomExternalForce"):
        write_prmtop(run)
    assert not (run / "meld_system.prmtop").exists()


def test_without_positions_the_check_is_skipped(tmp_path: Path) -> None:
    report = write_prmtop(_run_dir(tmp_path, with_positions=False))
    assert not report.checked and not report.energies


def test_loader_never_runs_code(tmp_path: Path) -> None:
    marker = tmp_path / "pwned"

    class Evil:
        def __reduce__(self):
            return (os.system, (f"touch {marker}",))

    path = tmp_path / "evil.dat"
    path.write_bytes(pickle.dumps(Evil()))
    obj = load_inert(path)
    assert not marker.exists()
    assert "system" in repr(obj)
    with pytest.raises(MddbmetaError, match="not a MELD System"):
        read_meld_system(path)


def test_not_a_pickle(tmp_path: Path) -> None:
    path = tmp_path / "system.dat"
    path.write_text("hello")
    with pytest.raises(MddbmetaError, match="not a readable"):
        write_prmtop(path, tmp_path / "out.prmtop")


def test_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    run = _run_dir(tmp_path)
    assert main(["meld-prmtop", str(run), "-o", str(tmp_path / "t.prmtop")]) == 0
    out = capsys.readouterr().out
    assert "wrote" in out and "nonbonded" in out and "implicit solvent" in out
    assert main(["meld-prmtop", str(tmp_path / "nowhere")]) == 1
    assert "no Data/system.dat" in capsys.readouterr().err


@pytest.mark.skipif(not PROTEIN_L, reason="set MDDBMETA_MELD_RUN to a real MELD run directory")
def test_real_meld_run(tmp_path: Path) -> None:
    report = write_prmtop(PROTEIN_L, tmp_path / "topology.prmtop")
    assert report.checked
    for meld_e, prm_e in report.energies.values():
        assert meld_e == pytest.approx(prm_e, abs=1e-3)
