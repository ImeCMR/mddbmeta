"""MELD replica exchange: synthetic trees (unit + scenario) and the real run (env-gated)."""

from __future__ import annotations

import os
import shutil
import struct
from pathlib import Path

import numpy as np
import pytest

from amberfiles import protein_residues, write_prmtop
from mddbmeta.engines.meld.files import read_remd_log, read_setup, temperature_of
from mddbmeta.export.mddb import build_export
from mddbmeta.model import FileRole
from mddbmeta.pipeline import discover
from mddbmeta.validate.schema_lite import check_inputs
from meldfiles import NATOM, PDB, write_legacy_system_dat, write_system_dat

SETUP = """
import meld
from meld import comm, vault, system
from meld.system.builders.grappa import GrappaOptions, GrappaSystemBuilder
from openmm import unit as u

N_REPLICAS = 3
N_STEPS = 10
BLOCK_SIZE = 5

def gen_state(s, index):
    state = s.get_state_template()
    state.alpha = index / (N_REPLICAS - 1.0)
    return state

def setup_system():
    options = GrappaOptions(solvation_type="implicit", grappa_model_tag="grappa-1.4.0",
                            base_forcefield_files=['amber14/protein.ff14SB.xml', 'implicit/gbn2.xml'],
                            default_temperature=300.0 * u.kelvin, use_big_timestep={big})
    s.temperature_scaler = meld.system.temperature.GeometricTemperatureScaler(0, 0.5, 300.*u.kelvin, 400.*u.kelvin)
    options = meld.RunOptions(timesteps=5000, minimize_steps=100)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)
    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)

setup_system()
"""


def _log(role: str, sessions: list[int], finished: bool) -> str:
    out = []
    for start in sessions:
        out += [f"Launching replica exchange on {role}", "Meld version is 0.6.1",
                "OpenMM_Meld version is 8.0.0.dev-a780005", "Loading system", "Using CUDA platform."]  # fmt: skip
        if role == "leader":
            out += [f"Running replica exchange step {i} of 10." for i in range(start, 11)]
    if finished and role == "leader":
        out.append("Finished 10 steps of replica exchange successfully.")
    return "\n".join(out) + "\n"


def _dcd(path: Path, natom: int, frames: int) -> None:
    def rec(b: bytes) -> bytes:
        return struct.pack("<i", len(b)) + b + struct.pack("<i", len(b))

    icntrl = [frames, 0, 1, frames] + [0] * 5
    head = (
        b"CORD"
        + struct.pack("<9i", *icntrl)
        + struct.pack("<f", 1.0)
        + struct.pack("<10i", *([0] * 9 + [24]))
    )
    title = struct.pack("<i", 1) + b"REMARKS extracted by MELD".ljust(80)
    body = b"".join(rec(struct.pack(f"<{natom}f", *([0.0] * natom))) * 3 for _ in range(frames))
    path.write_bytes(rec(head) + rec(title) + rec(struct.pack("<i", natom)) + body)


def _blocks(d: Path, alphas: list[list[float]], growing: bool = True) -> None:
    """MELD's data store: 5 frames per block. As MELD writes it (``growing``), every block's time
    axis runs from frame 0 to its own end and only its own columns hold data."""
    netCDF4 = pytest.importorskip("netCDF4")
    d.mkdir(parents=True)
    n = len(alphas[0])
    for b in range(0, len(alphas), 5):
        chunk = alphas[b : b + 5]
        with netCDF4.Dataset(str(d / f"block_{b // 5:06d}.nc"), "w") as ds:
            ds.createDimension("n_replicas", n)
            ds.createDimension("timesteps", None)
            v = ds.createVariable("alphas", "f8", ("n_replicas", "timesteps"))
            for j, col in enumerate(chunk):
                v[:, (b + j) if growing else j] = col
            if growing and len(chunk) < 5:  # the last block is allocated in full, partly written
                v[:, b + 4] = np.ma.masked


def make_meld_tree(root: Path, big: bool = False, finished: bool = True) -> Path:
    residues = protein_residues(4)
    write_prmtop(root / "system.prmtop", residues, ifbox=0)
    natom = sum(len(a) for _, a in residues)
    (root / "setup_meld.py").write_text(SETUP.format(big=big))
    (root / "Logs").mkdir()
    (root / "Logs" / "remd_000.log").write_text(_log("leader", [1, 7], finished))
    for i in (1, 2):
        (root / "Logs" / f"remd_00{i}.log").write_text(_log("worker", [1, 7], finished))
    # 11 stored frames; the middle position adapts from alpha 0.5 to 0.25 at frame 6
    alphas = [[0.0, 0.5 if f < 6 else 0.25, 1.0] for f in range(11)]
    _blocks(root / "Data" / "Blocks", alphas)
    (root / "extract_trajs.slurm").write_text(
        "extract_trajectory extract_traj_dcd trajectories/trajectory.$b.dcd --replica 0\n"
    )
    (root / "extract_walkers.sh").write_text(
        "extract_trajectory follow_dcd --replica 0 walkers/follow.$b.dcd\n"
    )
    for sub, stem in (("trajectories", "trajectory"), ("walkers", "follow")):
        (root / sub).mkdir()
        for i in range(3):
            _dcd(root / sub / f"{stem}.{i:02d}.dcd", natom, 10)
    return root


def test_setup_timestep_rule_and_scaler(tmp_path: Path) -> None:
    (tmp_path / "a.py").write_text(SETUP.format(big=False))
    (tmp_path / "b.py").write_text(SETUP.format(big=True))
    a, b = read_setup(tmp_path / "a.py"), read_setup(tmp_path / "b.py")
    assert (a.n_replicas, a.n_steps, a.block_size, a.timesteps) == (3, 10, 5, 5000)
    assert a.timestep_fs[0] == 2.0 and b.timestep_fs[0] == 3.5
    assert (
        a.adaptor == "EqualAcceptanceAdaptor"
        and a.initial_alpha_rule == "index / (N_REPLICAS - 1.0)"
    )
    sc = a.temperature_scaler
    assert temperature_of(sc, 0.0) == 300.0 and temperature_of(sc, 0.9) == 400.0
    assert temperature_of(sc, 0.25) == pytest.approx(300.0 * (400.0 / 300.0) ** 0.5)
    lin = {"kind": "linear", "alpha_min": 0.0, "alpha_max": 1.0, "t_min": 300.0, "t_max": 400.0}
    assert temperature_of(lin, 0.5) == 350.0


# the older MELD API (<= 0.5), as in real setup scripts: SystemBuilder + RunOptions attributes
OLD_SETUP = """
import numpy as np
from meld import system
import meld.system.montecarlo as mc

N_REPLICAS = 3
N_STEPS = 10
BLOCK_SIZE = 5

def setup_system():
    p = system.ProteinMoleculeFromSequence(sequence)
    b = system.SystemBuilder({builder_args})
    s = b.build_system_from_molecules([p])
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 0.5, 300., 400.)
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.use_big_timestep = {big}
    options.use_bigger_timestep = True
    options.cutoff = 1.8
    options.timesteps = 11111
    options.minimize_steps = 20000

setup_system()
"""

AMBER_SETUP = """
import meld
from meld.system.builders.amber import AmberOptions, AmberSystemBuilder
def setup_system():
    options = AmberOptions({args})
    builder = AmberSystemBuilder(options)
"""


def test_older_setup_syntax(tmp_path: Path) -> None:
    (tmp_path / "a.py").write_text(
        OLD_SETUP.format(builder_args='forcefield="ff14sbside"', big=False)
    )
    a = read_setup(tmp_path / "a.py")
    assert (a.builder, a.legacy) == ("amber", True)
    assert a.timestep_fs == (4.5, "use_bigger_timestep=True")
    assert (a.timesteps, a.minimize_steps, a.n_replicas) == (11111, 20000, 3)
    assert a.forcefield == ("ff14sbside", True) and a.leaprcs == ["leaprc.protein.ff14SBonlysc"]
    assert a.implicit_solvent_model == ("gbNeck2", True)
    # MELD tests use_big_timestep first
    (tmp_path / "b.py").write_text(OLD_SETUP.format(builder_args="", big=True))
    b = read_setup(tmp_path / "b.py")
    assert b.timestep_fs[0] == 3.5
    assert b.forcefield == ("ff14sbside", False)  # MELD's default, not stated
    (tmp_path / "c.py").write_text(OLD_SETUP.format(
        builder_args='"ff14sb", explicit_solvent=True, solvent_forcefield="opc"', big=False))  # fmt: skip
    c = read_setup(tmp_path / "c.py")
    assert c.leaprcs == ["leaprc.protein.ff14SB", "leaprc.water.opc"]
    assert c.implicit_solvent_model == (None, True)


def test_amber_options_forcefield(tmp_path: Path) -> None:
    (tmp_path / "a.py").write_text(
        AMBER_SETUP.format(args='forcefield="ff12sb", use_big_timestep=True')
    )
    a = read_setup(tmp_path / "a.py")
    assert (a.builder, a.legacy, a.timestep_fs[0]) == ("amber", False, 3.5)
    assert a.leaprcs == ["leaprc.ff12SB"] and a.implicit_solvent_model == ("gbNeck2", False)
    (tmp_path / "b.py").write_text(AMBER_SETUP.format(args='solvation="explicit"'))
    b = read_setup(tmp_path / "b.py")
    assert b.leaprcs == ["leaprc.protein.ff14SBonlysc", "leaprc.water.tip3p"]


def test_older_setup_in_a_run(tmp_path: Path) -> None:
    root = make_meld_tree(tmp_path)
    (root / "setup_meld.py").write_text(
        OLD_SETUP.format(builder_args='forcefield="ff14sbside"', big=False)
    )
    p = discover(root)
    m = p.mined
    assert m["ff"].value == ["Amber ff14SBonlysc"]
    assert m["timestep"].value == 4.5
    ex = build_export(p)
    meld = ex.inputs["metadditions"]["meld"]
    assert (
        meld["forcefield"] == "ff14sbside" and meld["implicit_solvent"] == "gbNeck2 (OpenMM GBn2)"
    )
    assert meld["exchange_interval_ps"] == pytest.approx(11111 * 4.5 / 1000)


def _combined_log() -> str:
    """MELD <= 0.4: one remd.log, every rank interleaved behind a host:rank date time prefix."""

    def ln(rank: str, level: str, logger: str, msg: str) -> str:
        return f"c45a-s21.ufhpc:{rank} 2021-12-19 18:12:32 {level} {logger}: {msg}"

    out = []
    for r in ("022", "014"):
        out += [ln(r, "INFO", "meld.remd.launch", "Launching replica exchange on follower"),
                ln(r, "INFO", "meld.remd.launch", "Meld version is 0.4.14"),
                ln(r, "INFO", "meld.remd.launch", "OpenMM_Meld version is 7.5.0.dev-Unknown")]  # fmt: skip
    out += [
        ln("000", "INFO", "meld.remd.leader", f"Running replica exchange step {i} of 40.")
        for i in (5, 6, 7)
    ]
    out += [ln("022", "DEBUG", "meld.comm", "receive_alpha_from_leader took 293.820 ms")] * 50
    out += [ln("000", "INFO", "meld.remd.launch", "Launching replica exchange on leader"),
            ln("000", "INFO", "meld.remd.launch", "Meld version is 0.4.14")]  # fmt: skip
    out += [
        ln("000", "INFO", "meld.remd.leader", f"Running replica exchange step {i} of 40.")
        for i in (7, 8)
    ]
    return "\n".join(out) + "\n"


def test_combined_remd_log(tmp_path: Path) -> None:
    p = tmp_path / "remd.log"
    p.write_text(_combined_log())
    log = read_remd_log(p)
    assert (log.role, log.meld_version, log.openmm_meld_version) == (
        "leader",
        "0.4.14",
        "7.5.0.dev-Unknown",
    )
    assert log.sessions == [5, 7] and (log.last_step, log.max_steps) == (8, 40)
    assert log.finished_steps is None


def test_combined_remd_log_in_a_run(tmp_path: Path) -> None:
    root = make_meld_tree(tmp_path)
    shutil.rmtree(root / "Logs")
    (root / "remd.log").write_text(_combined_log())
    p = discover(root)
    assert (p.mined["program"].value, p.mined["version"].value) == ("MELD", "0.4.14")


def test_uneven_ladder_dcds_are_flagged(tmp_path: Path) -> None:
    root = make_meld_tree(tmp_path)
    natom = sum(len(a) for _, a in protein_residues(4))
    _dcd(root / "trajectories" / "trajectory.01.dcd", natom, 4)  # an extraction that stopped early
    p = discover(root)
    msg = next(f.message for f in p.findings if f.code == "MELD009")
    assert "1 of 3" in msg and "down to 4" in msg


def test_data_store_layouts_read_the_same(tmp_path: Path) -> None:
    from mddbmeta.engines.meld.files import read_alpha_history

    alphas = [[0.0, 0.5 if f < 6 else 0.25, 1.0] for f in range(11)]
    _blocks(tmp_path / "grow", alphas, growing=True)
    _blocks(tmp_path / "flat", alphas, growing=False)
    g = read_alpha_history(sorted((tmp_path / "grow").iterdir()))
    f = read_alpha_history(sorted((tmp_path / "flat").iterdir()))
    assert (g.n_frames, g.last, g.settled_frame) == (f.n_frames, f.last, f.settled_frame)
    assert (g.n_frames, g.last, g.settled_frame) == (11, [0.0, 0.25, 1.0], [0, 6, 0])
    assert (g.minimum[1], g.maximum[1]) == (0.25, 0.5)


def test_remd_log_sessions(tmp_path: Path) -> None:
    p = tmp_path / "remd_000.log"
    p.write_text(_log("leader", [1, 7], True))
    log = read_remd_log(p)
    assert (log.role, log.meld_version, log.platform) == ("leader", "0.6.1", "CUDA")
    assert log.sessions == [1, 7] and log.finished_steps == 10
    w = tmp_path / "remd_001.log"
    w.write_text(_log("worker", [1], True))
    assert read_remd_log(w).role == "worker"


def test_meld_tree(tmp_path: Path) -> None:
    p = discover(make_meld_tree(tmp_path))
    assert p.engine == "meld"
    assert [r.id for r in p.replicas] == ["replica_00", "replica_01", "replica_02"]
    m = p.mined
    assert (m["program"].value, m["version"].value, m["type"].value) == (
        "MELD",
        "0.6.1",
        "ensemble",
    )
    assert (
        m["timestep"].value == 2.0 and m["ensemble"].value == "NVT" and not m["framestep"].is_known
    )
    assert m["ff"].value == ["Grappa 1.4.0", "Amber ff14SB"]
    temps = {rid: f["temp"].value for rid, f in p.md_mined.items()}
    assert temps["replica_00"] == 300.0 and temps["replica_02"] == 400.0
    assert temps["replica_01"] == pytest.approx(round(300.0 * (400.0 / 300.0) ** 0.5, 2))
    codes = {f.code for f in p.findings}
    assert {
        "MELD001",
        "MELD002",
        "MELD003",
        "MELD004",
    } <= codes  # adapted ladder, 10 vs 11 frames, walkers, restart
    ex = build_export(p)
    assert ex.ok and check_inputs(ex.inputs) == []
    assert [md["input_trajectory_filepaths"] for md in ex.inputs["mds"]] == [
        [f"trajectories/trajectory.{i:02d}.dcd"] for i in range(3)
    ]
    assert ex.inputs["metadditions"]["meld"]["exchange_interval_ps"] == 10.0  # 5000 x 2 fs
    assert "walkers" not in str(ex.inputs["mds"])


def test_unfinished_run(tmp_path: Path) -> None:
    p = discover(make_meld_tree(tmp_path, finished=False))
    step = next(iter(p.steps.values()))
    assert step.runtime.completed.value is False
    # one finding for the run, not one per ladder position
    from mddbmeta.validate import run_all

    cmp = [f for f in run_all(p) if f.code == "CMP002"]
    assert len(cmp) == 1 and "3 runs logged in Logs/remd_000.log" in cmp[0].message
    assert "stopped at exchange 10 of 10" in cmp[0].message


def test_unreadable_data_store_leaves_temps_blank(tmp_path: Path) -> None:
    root = make_meld_tree(tmp_path)
    for b in (root / "Data" / "Blocks").iterdir():
        b.write_bytes(b"CDF\x01" + b"\0" * 8)  # still looks like a block, but cannot be read
    p = discover(root)
    load = [f for f in p.findings if f.code == "LOAD001"]
    assert len(load) == 1 and "left blank" in load[0].message
    assert not [f for f in p.findings if f.code == "MELD001"]
    assert all("temp" not in f or not f["temp"].is_known for f in p.md_mined.values())


def test_missing_topology_points_to_meld_prmtop(tmp_path: Path) -> None:
    root = make_meld_tree(tmp_path)
    (root / "system.prmtop").unlink()
    p = discover(root)
    hints = [f for f in p.findings if f.code == "MINE004"]
    assert hints and "mddbmeta meld-prmtop" in hints[0].message


def _meld_with_system_dat(root: Path) -> Path:
    """A MELD tree whose run also saved Data/system.dat (33 atoms; DCDs to match)."""
    make_meld_tree(root)
    for sub, stem in (("trajectories", "trajectory"), ("walkers", "follow")):
        for i in range(3):
            _dcd(root / sub / f"{stem}.{i:02d}.dcd", NATOM, 10)
    (root / "system.prmtop").unlink()
    write_system_dat(root)
    return root


def test_topology_is_written_from_system_dat_then_used(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    p = discover(root)
    assert (root / "topology.prmtop").is_file()
    umask = os.umask(0)
    os.umask(umask)
    assert (root / "topology.prmtop").stat().st_mode & 0o777 == 0o666 & ~umask
    assert {s.files[FileRole.TOPOLOGY].path for s in p.steps.values()} == {"topology.prmtop"}
    assert p.mined["input_topology_filepath"].value == "topology.prmtop"
    assert p.system is not None and p.system.natom == NATOM
    msg = next(f.message for f in p.findings if f.code == "MELD005")
    assert "written" in msg and "energies match" in msg
    assert not [f for f in p.findings if f.code in ("MINE004", "MELD006")]

    # a second discovery reuses the file instead of rewriting it
    stamp = (root / "topology.prmtop").stat().st_mtime_ns
    p = discover(root)
    assert (root / "topology.prmtop").stat().st_mtime_ns == stamp
    assert "reused" in next(f.message for f in p.findings if f.code == "MELD005")


def test_written_topology_wins_over_another_prmtop(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    write_prmtop(
        root / "system.prmtop", [("ALA", ["N", "CA"])] * 1, ifbox=0
    )  # a stand-in from tleap
    p = discover(root)
    assert p.mined["input_topology_filepath"].value == "topology.prmtop"
    assert "not used: system.prmtop" in next(f.message for f in p.findings if f.code == "MELD005")


def test_starting_structure_is_written_when_the_run_has_no_pdb(tmp_path: Path) -> None:
    import numpy as np
    import openmm.app as app
    import openmm.unit as unit

    root = _meld_with_system_dat(tmp_path)
    p = discover(root)
    pdb = root / "topology.pdb"
    assert pdb.is_file()
    assert p.mined["input_structure_filepath"].value == "topology.pdb"
    assert "written" in next(f.message for f in p.findings if f.code == "MELD007")
    # same atoms, in the same order, as the written topology; the saved starting coordinates
    got = app.PDBFile(str(pdb))
    prm = app.AmberPrmtopFile(str(root / "topology.prmtop"))
    assert [a.name for a in got.topology.atoms()] == [a.name for a in prm.topology.atoms()]
    want = app.PDBFile(str(PDB)).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    assert np.allclose(
        got.getPositions(asNumpy=True).value_in_unit(unit.nanometer), want, atol=1e-3
    )

    stamp = pdb.stat().st_mtime_ns
    p = discover(root)
    assert pdb.stat().st_mtime_ns == stamp
    assert "reused" in next(f.message for f in p.findings if f.code == "MELD007")


def test_starting_structure_beats_analysis_pdbs_in_subfolders(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    (root / "cluster").mkdir()
    shutil.copy(PDB, root / "cluster" / "best_rmsd.pdb")  # a clustering result, same atom count
    p = discover(root)
    assert p.mined["input_structure_filepath"].value == "topology.pdb"


def test_a_pdb_mddbmeta_wrote_under_an_old_name_is_ignored(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    discover(root)
    shutil.copy(
        root / "topology.pdb", root / "meld_system.pdb"
    )  # left over from an earlier version
    p = discover(root)
    assert p.mined["input_structure_filepath"].value == "topology.pdb"


def test_the_runs_own_pdb_is_preferred(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    shutil.copy(PDB, root / "protein_min.pdb")
    p = discover(root)
    assert p.mined["input_structure_filepath"].value == "protein_min.pdb"
    assert not [f for f in p.findings if f.code == "MELD007"]


def test_foreign_structure_is_not_overwritten(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    (root / "topology.pdb").write_text("REMARK someone else's file\nEND\n")
    p = discover(root)
    assert (root / "topology.pdb").read_text() == "REMARK someone else's file\nEND\n"
    assert any(f.code == "MELD008" and "not overwriting" in f.message for f in p.findings)
    assert p.mined["input_topology_filepath"].value == "topology.prmtop"


def test_older_meld_layout_gives_topology_and_structure(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    write_legacy_system_dat(root)  # replaces the newer system.dat
    p = discover(root)
    assert p.mined["input_topology_filepath"].value == "topology.prmtop"
    assert p.mined["input_structure_filepath"].value == "topology.pdb"
    assert p.system is not None and p.system.natom == NATOM
    assert "written out unchanged" in next(f.message for f in p.findings if f.code == "MELD005")
    assert not [f for f in p.findings if f.code in ("MINE004", "MELD006", "MELD008")]


def test_read_only_discovery_writes_nothing(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    p = discover(root, write_topology=False)
    assert not (root / "topology.prmtop").exists()
    assert not (root / "topology.pdb").exists()
    assert "mddbmeta meld-prmtop" in next(f.message for f in p.findings if f.code == "MINE004")


def test_foreign_file_is_not_overwritten(tmp_path: Path) -> None:
    root = _meld_with_system_dat(tmp_path)
    (root / "topology.prmtop").write_text("someone else's file\n")
    p = discover(root)
    assert (root / "topology.prmtop").read_text() == "someone else's file\n"
    assert any(f.code == "MELD006" and "not overwriting" in f.message for f in p.findings)


@pytest.mark.external_data
def test_real_meld_run() -> None:
    """MDDBMETA_EXTERNAL_DATA=/path/to/protein_G_3GB1 pytest -m external_data"""
    root = next((Path(d) for d in os.environ.get("MDDBMETA_EXTERNAL_DATA", "").split(":")
                 if (Path(d) / "Data" / "Blocks").is_dir()), None)  # fmt: skip
    if root is None:
        pytest.skip("no MELD directory in MDDBMETA_EXTERNAL_DATA")
    p = discover(root)
    assert p.engine == "meld" and len(p.replicas) == 30
    m = p.mined
    assert (m["version"].value, m["timestep"].value, m["type"].value) == ("0.6.1", 2.0, "ensemble")
    assert (
        p.md_mined["replica_00"]["temp"].value == 300.0
        and p.md_mined["replica_29"]["temp"].value == 550.0
    )
    assert ("MELD004", [10400, 13850]) in [(f.code, f.data.get("restarts")) for f in p.findings]
    assert check_inputs(build_export(p).inputs) == []
