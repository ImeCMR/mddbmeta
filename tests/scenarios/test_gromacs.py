"""GROMACS discovery on real GROMACS 2022.3 output and on MDDB-workflow's own test data."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from mddbmeta.export.mddb import build_export
from mddbmeta.model import CoordSourceKind
from mddbmeta.pipeline import discover
from mddbmeta.provenance import Confidence

REAL = Path(__file__).resolve().parents[1] / "data" / "gromacs_real"
MDDB = Path(__file__).resolve().parents[3] / "MDDB-workflow" / "test" / "data" / "input"


@pytest.fixture(scope="module")
def real():
    return discover(REAL)


def codes(project, code):
    return [f for f in project.findings if f.code == code]


def test_layout_roles_and_replicas(real) -> None:
    assert real.engine == "gromacs"
    assert [r.id for r in real.replicas] == [
        "rep1",
        "rep2",
    ]  # disjoint run names, independent lineages
    shared = {s.id for s in real.steps.values() if s.replica is None}
    assert shared == {"prep/em", "prep/nvt", "prep/npt"}
    roles = {s.id: s.role.value for s in real.steps.values()}
    assert (
        roles["prep/em"] == "minimization"
        and real.steps["prep/em"].role.confidence is Confidence.EXACT
    )
    assert roles["prep/nvt"] == roles["prep/npt"] == "equilibration"
    assert roles["rep2/prod"] == roles["rep1/md_2"] == "production"


def test_continuations_from_engine_evidence(real) -> None:
    ref = {s.id: (s.input_coords.ref, s.input_coords.method) for s in real.steps.values()
           if s.input_coords.kind is CoordSourceKind.STEP}  # fmt: skip
    assert ref["rep1/md_2"][0] == "rep1/md_1" and "energy fingerprint" in ref["rep1/md_2"][1]
    assert ref["rep1/md_3"][0] == "rep1/md_2"
    assert ref["rep1/md_1"][0] == ref["rep2/prod"][0] == "prep/npt"
    assert ref["rep2/prod.part0003"] == ("rep2/prod", "checkpoint read by mdrun -cpi")
    assert len(codes(real, "LIN003")) == 5  # every fingerprint inference is announced
    assert real.replicas[0].production_chain == ["rep1/md_1", "rep1/md_2", "rep1/md_3"]
    assert real.replicas[1].production_chain == ["rep2/prod", "rep2/prod.part0003"]
    assert not codes(real, "ORDER001")


def test_clock_resets_and_healthy_checkpoint_continuation(real) -> None:
    resets = codes(real, "CONT005")
    prod_resets = {f.subject for f in resets if f.severity.value == "warning"}
    assert prod_resets == {"rep1/md_2", "rep1/md_3"}  # prep-stage resets are info only
    healthy = [f for f in codes(real, "CONT000") if f.subject == "rep2/prod.part0003"]
    assert healthy and healthy[0].data["gap_ps"] == 0.0
    assert not codes(real, "CMP003")  # GROMACS writes the step-0 frame: 2500/500 + 1 = 6 frames


def test_mined_fields(real) -> None:
    m = real.mined
    assert (m["program"].value, m["version"].value) == ("GROMACS", "2022.3")
    assert m["timestep"].value == 2.0 and m["timestep"].method.startswith("engine_report")
    assert m["framestep"].value == pytest.approx(0.001) and m["temp"].value == 300.0
    assert m["ensemble"].value == "NPT" and m["boxtype"].value == "Dodecahedron"
    assert m["ff"].value == ["Amber ff99SB-ILDN"] and m["ff"].confidence is Confidence.DERIVED
    assert m["wat"].value == "TIP3P" and m["wat"].confidence is Confidence.DERIVED


def test_export_reports_duplicated_join_frames(real) -> None:
    ex = build_export(real)
    assert ex.ok  # nothing here stops MDDB-workflow
    (mrg,) = [f for f in ex.findings if f.code == "MRG001"]
    assert mrg.subject == "rep2" and mrg.data["joins"] == ["rep2/prod -> rep2/prod.part0003"]
    assert mrg.severity.value == "warning"
    assert [md["input_topology_filepath"] for md in ex.inputs["mds"]] == [
        "rep1/md_1.tpr",
        "rep2/prod.tpr",
    ]


def test_mddb_dummy_gromacs_project() -> None:
    d = MDDB / "dummy" / "gromacs"
    if not d.exists():
        pytest.skip("MDDB-workflow test data not available")
    p = discover(d)
    assert p.engine == "gromacs"
    orphans = [s for s in p.steps.values() if s.status == "orphan"]
    assert {s.id for s in orphans} == {"raw_trajectory", "rmsd_jump"}
    assert not p.mined["framestep"].is_known  # the xtc stores t=0 for every frame
    assert p.mined["ff"].value == ["Amber ff99SB-ILDN"]
    assert p.mined["wat"].method == "no water residues"  # vacuum dipeptide


def test_mddb_raw_project_matches_its_inputs(tmp_path: Path) -> None:
    d = MDDB / "raw_project"
    if not d.exists():
        pytest.skip("MDDB-workflow test data not available")
    work = tmp_path / "raw_project"
    work.mkdir()
    for name in ("raw_topology.tpr", "raw_trajectory.xtc", "2_frames.mdcrd"):
        shutil.copy(d / name, work / name)
    p = discover(work)
    assert p.engine == "gromacs" and codes(p, "ENG001")  # the stray AMBER mdcrd is noticed
    # MDDB-workflow's hand-written inputs.yaml for this project says framestep: 0.01 (ns)
    assert p.mined["framestep"].value == pytest.approx(0.01)
    ex = build_export(p)
    assert ex.inputs["mds"][0]["input_trajectory_filepaths"] == ["raw_trajectory.xtc"]
    assert ex.inputs["input_topology_filepath"] == "raw_topology.tpr"


def test_extended_run_nsteps_comes_from_records(real) -> None:
    prod = real.steps["rep2/prod"]
    assert prod.settings.nsteps.value == 5000 and "extended" in prod.settings.nsteps.method
    assert prod.runtime.nsteps_completed.value == 5000 and prod.raw["gmx_sessions"] == 2


def test_linked_chunk_dirs_are_one_replica(tmp_path: Path) -> None:
    """Directories named like replicas but chained by restarts are consecutive chunks."""
    shutil.copytree(REAL / "prep", tmp_path / "prep")
    for i in (1, 2, 3):
        d = tmp_path / f"run{i}"
        d.mkdir()
        for ext in ("log", "tpr", "xtc"):
            shutil.copy(REAL / "rep1" / f"md_{i}.{ext}", d / f"md_{i}.{ext}")
    p = discover(tmp_path)
    assert len(p.replicas) == 1
    assert p.replicas[0].production_chain == ["run1/md_1", "run2/md_2", "run3/md_3"]


def test_independent_labelled_dirs_are_replicas(tmp_path: Path) -> None:
    shutil.copytree(REAL / "prep", tmp_path / "prep")
    shutil.copytree(REAL / "rep1", tmp_path / "rep1")
    shutil.copytree(REAL / "rep2", tmp_path / "rep2")
    p = discover(tmp_path)
    assert [r.id for r in p.replicas] == ["rep1", "rep2"]
    # with no restart link between them, disjoint run names don't make them one replica
    assert not any(s.replica is None for s in p.steps.values() if s.dir.startswith("rep"))
