from __future__ import annotations

import struct
from pathlib import Path

import pytest

from mddbmeta.engines.namd import NamdAdapter, residue_classes
from mddbmeta.engines.namd.config import read_config
from mddbmeta.engines.namd.files import binary_natom, read_dcd, read_pdb, read_psf, read_xsc
from mddbmeta.engines.namd.log import read_log
from mddbmeta.engines.namd.semantics import settings_from_config, settings_from_log
from mddbmeta.model import FileKind

REAL = Path(__file__).resolve().parents[1] / "data" / "namd_real"


def test_config_variables_sources_and_tcl() -> None:
    cfg = read_config(REAL / "rep2" / "prod.namd")
    assert cfg.get("bincoordinates") == "../prep/eq.coor"  # $inputname substituted
    assert cfg.get("langevintemp") == "300" and cfg.get("outputname") == "prod"
    assert len(cfg.parameters) == 6 and cfg.sources  # via `source ../common.inc`
    assert cfg.runs == [("run", 2500)]
    tcl = read_config(REAL / "rep2" / "prod_2.namd")
    assert "firsttimestep" in tcl.unresolved and tcl.has_tcl_blocks  # proc + [get_first_ts ...]
    heat = read_config(REAL / "prep" / "heat.namd")
    assert heat.has_tcl_blocks and heat.runs == [
        ("run", 100)
    ]  # the loop's runs are left to the log


def test_config_settings_for_unrun_input(tmp_path: Path) -> None:
    p = tmp_path / "next.namd"
    p.write_text("structure a.psf\ncoordinates a.pdb\ntimestep 2.0\nlangevin on\nlangevinTemp 310\n"
                 "langevinPiston on\ndcdfreq 5000\nrun 500000\n")  # fmt: skip
    s = settings_from_config(read_config(p), "next.namd")
    assert (s.dt_ps.value, s.temp_K.value, s.ensemble.value) == (0.002, 310.0, "NPT")
    assert s.nsteps.value == 500000 and s.traj_interval_steps.value == 5000


def test_log_info_tcl_and_energies() -> None:
    log = read_log(REAL / "rep1" / "md_2.log")
    assert (log.version, log.platform, log.natom) == ("2.14b1", "Linux-x86_64-verbs-smp", 1032)
    assert log.one("FIRST TIMESTEP") == "1500" and log.one("BINARY COORDINATES") == "md_1.coor"
    assert log.runs == [("Running", 1500)] and log.completed
    # the first ENERGY rows come before any ETITLE header and must not be lost
    assert log.first.step == 1500 and "POTENTIAL" in log.first.values and log.last.step == 3000


def test_heating_ramp_and_restraints_from_log() -> None:
    heat = settings_from_log(read_log(REAL / "prep" / "heat.log"), "heat.log")
    assert heat.temp_ramp.value is True and heat.temp_K.value == 300.0  # last langevinTemp set
    assert (
        heat.temp_initial_K.value == 0.0
        and heat.nsteps.value == 600
        and heat.ensemble.value == "NVT"
    )
    eq = settings_from_log(read_log(REAL / "prep" / "eq.log"), "eq.log")
    assert (
        eq.restraints.value is True
        and eq.barostat.value == "langevin-piston"
        and eq.ensemble.value == "NPT"
    )
    assert eq.cutoff_A.value == 10.0
    mini = settings_from_log(read_log(REAL / "prep" / "min.log"), "min.log")
    assert mini.minimization.value is True and mini.nsteps.value == 200


def test_failed_run_is_not_complete(tmp_path: Path) -> None:
    text = (REAL / "rep1" / "md_1.log").read_text().split("\nENERGY:", 1)[0] + "\n"
    p = tmp_path / "bad.log"
    p.write_text(text + "FATAL ERROR: number of steps must be a multiple of stepsPerCycle\n"
                 "WallClock: 1.0  CPUTime: 1.0  Memory: 1 MB\n[Partition 0][Node 0] End of program\n")  # fmt: skip
    log = read_log(p)
    assert log.fatal and not log.completed  # NAMD prints WallClock even when it dies


def test_dcd_header_and_times() -> None:
    d = read_dcd(REAL / "rep1" / "md_2.dcd")
    assert (d.natom, d.n_frames, d.first_step, d.step_interval, d.timestep_fs) == (
        1032,
        6,
        1750,
        250,
        2.0,
    )
    assert d.time_ps(0) == 3.5 and d.time_ps(5) == 6.0 and d.has_cell and not d.truncated


def test_dcd_truncation(tmp_path: Path) -> None:
    data = (REAL / "rep1" / "md_2.dcd").read_bytes()
    p = tmp_path / "cut.dcd"
    p.write_bytes(data[:-500])
    d = read_dcd(p)
    assert d.truncated and d.n_frames == 5


def test_xsc_psf_pdb_and_binary() -> None:
    x = read_xsc(REAL / "rep1" / "md_1.xsc")
    assert x.step == 1500 and x.vectors_A[0][0] == pytest.approx(24.3243677112)
    psf = read_psf(REAL / "system" / "solv.psf")
    assert psf.natom == 1032 and psf.residue_counts == {"ALA": 3, "TIP3": 333}
    assert psf.residue_atom_names["TIP3"] == ["OH2", "H1", "H2"]
    assert read_pdb(REAL / "system" / "solv.pdb").natom == 1032
    coor = (REAL / "rep1" / "md_1.coor").read_bytes()
    assert binary_natom(coor[:4], len(coor)) == 1032
    assert binary_natom(struct.pack("<i", 5) + b"\0" * 10, 14) is None


def test_residue_classes_gate_force_field_labels() -> None:
    assert residue_classes({"ALA": 3, "TIP3": 10, "SOD": 1}) == {"protein"}
    assert residue_classes({"POPC": 100, "CHL1": 10}) == {"lipid"}
    assert residue_classes({"LIG": 1}) == {"other"}


def test_sniff_real_files() -> None:
    ad = NamdAdapter()
    expected = {
        "prep/eq.log": FileKind.LOG, "prep/eq.namd": FileKind.CONTROL, "system/solv.psf": FileKind.TOPOLOGY,
        "system/solv.pdb": FileKind.COORDS, "rep1/md_1.dcd": FileKind.TRAJECTORY, "rep1/md_1.xsc": FileKind.COORDS,
        "rep1/md_1.coor": FileKind.COORDS,
    }  # fmt: skip
    for rel, kind in expected.items():
        path = REAL / rel
        res = ad.sniff(path, path.read_bytes()[:4096])
        assert res is not None and res.kind is kind, rel
