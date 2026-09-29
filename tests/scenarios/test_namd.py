"""NAMD discovery on real NAMD 2.14 output."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from mddbmeta.export.mddb import build_export
from mddbmeta.model import CoordSourceKind
from mddbmeta.pipeline import discover
from mddbmeta.provenance import Confidence

REAL = Path(__file__).resolve().parents[1] / "data" / "namd_real"


@pytest.fixture(scope="module")
def real():
    return discover(REAL)


def codes(project, code):
    return [f for f in project.findings if f.code == code]


def test_layout_roles_and_lineage(real) -> None:
    assert real.engine == "namd"
    assert [r.id for r in real.replicas] == ["rep1", "rep2"]
    assert {s.id for s in real.steps.values() if s.replica is None} == {
        "prep/min",
        "prep/heat",
        "prep/eq",
    }
    roles = {s.id: s.role.value for s in real.steps.values()}
    assert roles["prep/min"] == "minimization" and roles["prep/heat"] == "heating"
    assert roles["prep/eq"] == "equilibration" and roles["rep2/prod_2"] == "production"
    src = {s.id: s.input_coords for s in real.steps.values()}
    assert src["prep/min"].kind is CoordSourceKind.STARTING_STRUCTURE
    assert real.starting_structure.path == "system/solv.pdb"
    assert (src["rep1/md_3"].ref, src["rep2/prod"].ref) == ("rep1/md_2", "prep/eq")
    assert src["rep2/prod_2"].ref == "rep2/prod" and "restart" in src["rep2/prod_2"].method
    assert real.replicas[1].production_chain == ["rep2/prod", "rep2/prod_2"]


def test_continuity(real) -> None:
    healthy = {f.subject for f in codes(real, "CONT000")}
    assert {"rep1/md_2", "rep1/md_3", "rep2/prod_2"} <= healthy  # firsttimestep carried the clock
    resets = codes(real, "CONT005")
    assert resets and all(f.severity.value == "info" for f in resets)  # only prep -> production
    assert not codes(real, "CMP003") and not codes(real, "ORDER001")


def test_mined_fields_and_export(real) -> None:
    m = real.mined
    assert (m["program"].value, m["version"].value) == ("NAMD", "2.14b1")
    assert m["timestep"].value == 2.0 and m["framestep"].value == pytest.approx(0.0005)
    assert (m["temp"].value, m["ensemble"].value, m["boxtype"].value) == (
        300.0,
        "NPT",
        "Rectangular",
    )
    assert m["ff"].value == [
        "CHARMM36"
    ]  # lipid/NA/carb/CGenFF files loaded but not used by this system
    assert m["wat"].value == "TIP3P" and m["wat"].confidence is Confidence.DERIVED
    ex = build_export(real)
    assert ex.ok and not [f for f in ex.findings if f.code == "MRG001"]
    assert ex.inputs["input_topology_filepath"] == "system/solv.psf"
    assert ex.inputs["input_structure_filepath"] == "system/solv.pdb"
    assert ex.inputs["mds"][0]["input_trajectory_filepaths"] == [
        "rep1/md_1.dcd",
        "rep1/md_2.dcd",
        "rep1/md_3.dcd",
    ]


def test_missing_restart_file_is_reported(tmp_path: Path) -> None:
    shutil.copytree(REAL, tmp_path / "p")
    (tmp_path / "p" / "rep1" / "md_1.coor").unlink()
    p = discover(tmp_path / "p")
    assert p.steps["rep1/md_2"].input_coords.kind is CoordSourceKind.UNKNOWN
    assert any(f.subject == "rep1/md_2" for f in codes(p, "CONT002"))


def test_amber_topology_run() -> None:
    """`amber on`: the prmtop/inpcrd belong to the AMBER sniffer; NAMD adopts them."""
    p = discover(REAL.parent / "namd_amber")
    assert p.engine == "namd" and not codes(p, "ENG001")
    assert [t.path for t in p.topologies] == ["ala_ala.prmtop"]
    assert p.starting_structure.path == "ala_ala.inpcrd"
    m = p.mined
    assert (m["timestep"].value, m["temp"].value, m["ensemble"].value) == (1.0, 300.0, "NVT")
    assert m["framestep"].value == pytest.approx(0.0001)
    assert m["boxtype"].method == "non-periodic" and m["wat"].method == "no water residues"
    assert build_export(p).inputs["input_topology_filepath"] == "ala_ala.prmtop"
