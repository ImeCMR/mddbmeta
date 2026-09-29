"""OpenMM discovery on real OpenMM 8.6.1 output."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from mddbmeta.export.mddb import build_export
from mddbmeta.model import CoordSourceKind
from mddbmeta.pipeline import discover
from mddbmeta.provenance import Confidence

DATA = Path(__file__).resolve().parents[1] / "data"
REAL = DATA / "openmm_real"


@pytest.fixture(scope="module")
def real():
    return discover(REAL)


def codes(project, code):
    return [f for f in project.findings if f.code == code]


def test_layout_and_lineage(real) -> None:
    assert real.engine == "openmm"
    assert [r.id for r in real.replicas] == ["rep1", "rep2"]
    assert {s.id for s in real.steps.values() if s.replica is None} == {"prep/nvt", "prep/npt"}
    src = {s.id: s.input_coords for s in real.steps.values()}
    assert src["prep/npt"].ref == "prep/nvt" and "same script" in src["prep/npt"].method
    assert (
        src["rep1/prod_1"].ref == "prep/npt"
    )  # loadState('../prep/eq_state.xml') outside the loop
    assert src["rep1/prod_3"].ref == "rep1/prod_2"  # loadCheckpoint(f'prod_{i-1}.chk') inside it
    assert (
        src["rep2/prod2"].ref == "rep2/prod"
    )  # restart.py loads prod.chk from prod.py's CheckpointReporter
    assert (
        real.starting_structure.path == "prep/solvated.pdb"
    )  # the Modeller output, not ala_ala.pdb
    assert real.replicas[0].production_chain == ["rep1/prod_1", "rep1/prod_2", "rep1/prod_3"]


def test_settings_and_times(real) -> None:
    npt, nvt = real.steps["prep/npt"].settings, real.steps["prep/nvt"].settings
    assert (
        nvt.ensemble.value == "NVT" and npt.ensemble.value == "NPT"
    )  # barostat added between them
    assert real.steps["rep2/prod"].settings.restart.value is False  # setVelocitiesToTemperature
    p1 = real.steps["rep1/prod_1"]
    assert p1.runtime.start_time_ps.value == pytest.approx(
        10.0
    ) and p1.runtime.end_time_ps.value == pytest.approx(15.0)
    assert p1.runtime.completed.value is True
    # DCD/XTC reporters count from their own creation; shifted onto the run clock
    assert p1.trajectory.first_time_ps.value == pytest.approx(11.0)
    assert {f.subject for f in codes(real, "CONT000")} >= {
        "rep1/prod_2",
        "rep1/prod_3",
        "rep2/prod2",
    }
    assert not codes(real, "CMP003") and not codes(real, "CMP004") and not codes(real, "CONS004")


def test_mined_fields_and_export(real) -> None:
    m = real.mined
    assert (m["program"].value, m["version"].value) == ("OpenMM", "8.6.1")
    assert (m["timestep"].value, m["temp"].value, m["ensemble"].value) == (2.0, 300.0, "NPT")
    assert m["framestep"].value == pytest.approx(0.001)
    assert m["ff"].value == ["Amber ff14SB"] and m["ff"].confidence is Confidence.DERIVED
    assert m["wat"].value == "TIP3P-FB" and m["boxtype"].value == "Cubic"
    ex = build_export(real)
    assert ex.ok and ex.inputs["input_topology_filepath"] == "no"
    assert ex.inputs["input_structure_filepath"] == "prep/solvated.pdb"
    assert [f.code for f in ex.findings] == []
    assert codes(real, "MINE004")


def test_amber_topology_script_with_stdout_log() -> None:
    p = discover(DATA / "openmm_amber")
    assert p.engine == "openmm" and not codes(p, "ENG001")
    (step,) = [s for s in p.steps.values() if s.status == "done"]
    assert step.files and step.input_coords.kind is CoordSourceKind.STARTING_STRUCTURE
    m = p.mined
    assert (m["timestep"].value, m["temp"].value, m["ensemble"].value) == (1.0, 310.0, "NVT")
    assert m["boxtype"].method == "non-periodic"
    assert build_export(p).inputs["input_topology_filepath"] == "ala_ala.prmtop"


def test_incomplete_segment(tmp_path: Path) -> None:
    shutil.copytree(REAL, tmp_path / "p")
    log = tmp_path / "p" / "rep2" / "prod2.log"
    log.write_text(
        "".join(log.read_text().splitlines(keepends=True)[:3])
    )  # the run was killed early
    p = discover(tmp_path / "p")
    rt = p.steps["rep2/prod2"].runtime
    assert rt.completed.value is False and "1000 of 2500" in rt.completed.note
    assert codes(p, "CMP002")


def test_log_timestep_beats_script(tmp_path: Path) -> None:
    """A script edited after the run (or misread statically) must not override the engine's record."""
    shutil.copytree(REAL, tmp_path / "p")
    run = tmp_path / "p" / "rep1" / "run.py"
    run.write_text(run.read_text().replace("0.002*picoseconds", "0.004*picoseconds"))
    dt = discover(tmp_path / "p").steps["rep1/prod_2"].settings.dt_ps
    assert (
        dt.value == pytest.approx(0.002)
        and dt.method.startswith("engine_report")
        and "0.004" in dt.note
    )
