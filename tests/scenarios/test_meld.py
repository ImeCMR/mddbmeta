"""MELD replica exchange: synthetic trees (unit + scenario) and the real run (env-gated)."""

from __future__ import annotations

import os
import struct
from pathlib import Path

import pytest

from amberfiles import protein_residues, write_prmtop
from mddbmeta.engines.meld.files import read_remd_log, read_setup, temperature_of
from mddbmeta.export.mddb import build_export
from mddbmeta.pipeline import discover
from mddbmeta.validate.schema_lite import check_inputs

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


def _blocks(d: Path, alphas: list[list[float]]) -> None:
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
                v[:, j] = col


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
