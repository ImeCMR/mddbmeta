"""End-to-end discovery scenarios on synthetic AMBER trees."""

from __future__ import annotations

from pathlib import Path

import pytest

from amberfiles import (
    EQUIL,
    HEAT,
    MIN,
    PROD,
    WATER_OPC,
    RunSpec,
    TreeBuilder,
    full_protocol,
    prod_chain,
    protein_residues,
    write_mdcrd,
    write_nc_traj,
    write_run,
)
from mddbmeta.errors import MddbmetaError
from mddbmeta.export.mddb import build_export
from mddbmeta.model import CoordSourceKind, FileRole
from mddbmeta.pipeline import discover
from mddbmeta.provenance import Confidence

pytestmark = pytest.mark.netcdf


def codes(project, code: str | None = None):
    found = [f.code for f in project.findings]
    return found if code is None else [f for f in project.findings if f.code == code]


def errors(project):
    return [f for f in project.findings if f.severity.value == "error"]


def test_linear_protocol(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", full_protocol(3), st, top)
    p = discover(tree.root)
    roles = {s.id: s.role.value for s in p.steps.values()}
    assert roles == {
        "min": "minimization", "heat": "heating", "equil": "equilibration",
        "prod_0001": "production", "prod_0002": "production", "prod_0003": "production",
    }  # fmt: skip
    assert p.steps["min"].role.confidence is Confidence.EXACT
    assert p.steps["prod_0002"].input_coords.kind is CoordSourceKind.STEP
    assert p.steps["prod_0002"].input_coords.ref == "prod_0001"
    assert p.starting_structure.path == "system.inpcrd"
    (rep,) = p.replicas
    assert rep.production_chain == ["prod_0001", "prod_0002", "prod_0003"]
    assert rep.lineages == [["min", "heat", "equil", "prod_0001", "prod_0002", "prod_0003"]]
    assert len(codes(p, "CONT000")) == 4
    assert not errors(p)
    m = p.mined
    assert m["timestep"].value == 2.0 and m["timestep"].unit == "fs"
    assert m["framestep"].value == pytest.approx(0.001) and m["framestep"].unit == "ns"
    assert m["temp"].value == 300.0 and m["ensemble"].value == "NPT"
    assert m["program"].value == "AMBER" and m["version"].value == "22"
    assert m["boxtype"].value == "Cubic" and m["type"].value == "trajectory"
    assert m["wat"].confidence is Confidence.HEURISTIC  # no build log: only a guess
    assert not m["ff"].is_known


def test_sequence_hole_and_missing_restart(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(5), st, top)
    for f in tree.root.glob("prod_0003.*"):
        f.unlink()
    p = discover(tree.root)
    (hole,) = codes(p, "SEQ001")
    assert hole.data["missing"] == [3]
    assert p.steps["prod_0004"].input_coords.kind is CoordSourceKind.UNKNOWN
    assert any(f.subject == "prod_0004" for f in codes(p, "CONT002"))
    assert p.replicas[0].production_chain == ["prod_0001", "prod_0002", "prod_0004", "prod_0005"]
    # the trajectory join where frames are missing is reported for the merge
    assert any(f.subject == "prod_0004" and "missing" in f.message for f in codes(p, "CONT001"))


def test_replicas_share_topology(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    end = tree.chain("prep", [("min", MIN), ("heat", HEAT), ("equil", EQUIL)], st, top)
    for r in ("rep1", "rep2", "rep3"):
        tree.chain(r, prod_chain(2), "prep/equil.rst7", top, start_time=end)
    p = discover(tree.root)
    assert [r.id for r in p.replicas] == ["rep1", "rep2", "rep3"]
    assert {s.id for s in p.steps.values() if s.replica is None} == {
        "prep/min",
        "prep/heat",
        "prep/equil",
    }
    assert p.replicas[1].production_chain == ["rep2/prod_0001", "rep2/prod_0002"]
    assert "LIN001" not in codes(p)  # a fork across replicas is how replicas start
    ex = build_export(p)
    assert ex.inputs["input_topology_filepath"] == "system.prmtop"
    assert [md["mdir"] for md in ex.inputs["mds"]] == ["rep1", "rep2", "rep3"]
    assert ex.inputs["mds"][2]["input_trajectory_filepaths"] == [
        "rep3/prod_0001.nc",
        "rep3/prod_0002.nc",
    ]


def test_crashed_replica_is_reported(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("rep1", prod_chain(3), st, top)
    tree.chain("rep2", prod_chain(1), st, top)
    p = discover(tree.root)
    assert [r.id for r in p.replicas] == ["rep1", "rep2"]
    (f,) = codes(p, "REP002")
    assert f.data["missing"] == {"rep2": ["prod_0002", "prod_0003"]}


def test_nested_layout_is_refused_then_declared(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    for cond in ("300K", "310K"):
        for r in ("rep1", "rep2"):
            tree.chain(f"{cond}/{r}", prod_chain(2), st, top)
    p = discover(tree.root)
    assert codes(p, "REP001")
    with pytest.raises(MddbmetaError, match="ambiguous"):
        build_export(p)
    p2 = discover(tree.root, replica_glob="*/rep*")
    assert [r.id for r in p2.replicas] == ["300K/rep1", "300K/rep2", "310K/rep1", "310K/rep2"]
    assert not codes(p2, "REP001")


def test_truncated_runs(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(3), st, top, overrides={"prod_0003": {"completed": False}})
    p = discover(tree.root)
    (f,) = codes(p, "CMP002")
    assert f.subject == "prod_0003" and f.severity.value == "warning"


def test_truncated_run_mid_chain_is_an_error(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(3), st, top, overrides={"prod_0002": {"completed": False}})
    # a crashed run can leave a checkpoint restart (ntwr) behind, and the next run read it
    (tree.root / "prod_0002.rst7").write_text((tree.root / "prod_0001.rst7").read_text())
    p = discover(tree.root)
    f = [x for x in codes(p, "CMP002") if x.subject == "prod_0002"]
    assert f and f[0].severity.value == "error"


def test_temperature_differs_per_replica(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("rep1", prod_chain(2), st, top)
    tree.chain("rep2", prod_chain(2, cntrl={**PROD, "temp0": 310.0}), st, top)
    p = discover(tree.root)
    assert p.mined["temp"].value is None
    assert p.md_mined["rep1"]["temp"].value == 300.0 and p.md_mined["rep2"]["temp"].value == 310.0
    assert codes(p, "MINE003")
    ex = build_export(p)
    assert "temp" not in ex.inputs
    assert [md.get("temp") for md in ex.inputs["mds"]] == [300.0, 310.0]


def test_order_follows_lineage_not_names(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(11, prefix="md", width=1), st, top)
    p = discover(tree.root)
    chain = p.replicas[0].production_chain
    assert chain == [f"md{i}" for i in range(1, 12)]  # md10 after md9, not after md1
    assert "ORDER001" not in codes(p)


def test_orphan_trajectories_ordered_by_time(tmp_path: Path) -> None:
    tree = TreeBuilder(tmp_path)
    top = tree.topology()
    # no logs at all: trajectories named so that name order is wrong
    write_nc_traj(tmp_path / "b_part.nc", tree.natom, [11.0, 12.0, 13.0])
    write_nc_traj(tmp_path / "a_part.nc", tree.natom, [14.0, 15.0, 16.0])
    write_nc_traj(tmp_path / "c_part.nc", tree.natom, [8.0, 9.0, 10.0])
    p = discover(tmp_path)
    assert p.replicas[0].production_chain == ["c_part", "b_part", "a_part"]
    assert "ORDER001" not in codes(p)
    assert len(codes(p, "STEP002")) == 3
    assert p.mined["framestep"].value == pytest.approx(0.001)
    assert p.mined["program"].value == "AMBER"  # from the NetCDF program attribute
    assert p.mined["program"].confidence is Confidence.DERIVED
    assert build_export(p).inputs["input_topology_filepath"] == top


def test_orphan_trajectories_without_times_fall_back_to_names(tmp_path: Path) -> None:
    tree = TreeBuilder(tmp_path)
    tree.topology()
    for n in (10, 2, 1):
        write_mdcrd(tmp_path / f"part{n}.mdcrd", tree.natom, 3)
    p = discover(tmp_path)
    assert p.replicas[0].production_chain == ["part1", "part2", "part10"]
    assert codes(p, "ORDER001")


def test_mixed_trajectory_formats_block_export(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(2), st, top, overrides={"prod_0002": {"traj_format": "mdcrd"}})
    p = discover(tree.root)
    ex = build_export(p)
    assert any(f.code == "FMT001" for f in ex.findings)
    assert not ex.ok


def test_replicas_from_name_tokens(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", [("rep1_prod_001", PROD), ("rep1_prod_002", PROD)], st, top)
    tree.chain("", [("rep2_prod_001", PROD), ("rep2_prod_002", PROD)], st, top)
    p = discover(tree.root)
    assert [r.id for r in p.replicas] == ["rep1", "rep2"]
    assert codes(p, "REP004")
    ex = build_export(p)
    assert [md["mdir"] for md in ex.inputs["mds"]] == ["replica_1", "replica_2"]


def test_frame_count_mismatch(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(2), st, top, overrides={"prod_0002": {"extra_frames": -2}})
    p = discover(tree.root)
    (f,) = codes(p, "CMP003")
    assert f.data == {"frames": 8, "expected": 10}


def test_time_gap_between_runs(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(1), st, top)
    # prod_0002 claims to start 50 ps after prod_0001 ended (e.g. a restart from a later run)
    write_run(
        tree.root,
        "prod_0002",
        RunSpec(cntrl=dict(PROD), start_time_ps=60.0),
        "prod_0001.rst7",
        top,
        tree.natom,
    )
    p = discover(tree.root)
    gaps = [f for f in codes(p, "CONT001") if "prod_0001 -> prod_0002" in f.message]
    assert gaps and gaps[0].data["gap_ps"] == pytest.approx(50.0)


def test_atom_count_mismatch(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    write_run(tree.root, "prod_0001", RunSpec(cntrl=dict(PROD)), st, top, tree.natom + 3)
    p = discover(tree.root)
    assert codes(p, "CONS004")


def test_build_log_and_extra_point_water(tmp_path: Path) -> None:
    tree = TreeBuilder(tmp_path, residues=protein_residues(3) + [WATER_OPC] * 5)
    top, st = tree.topology(numextra=5), tree.starting()
    tree.chain("", prod_chain(1), st, top)
    p = discover(tmp_path)
    assert p.mined["wat"].value == "OPC" and p.mined["wat"].confidence is Confidence.HEURISTIC
    (tmp_path / "leap.log").write_text(
        "Welcome to LEaP!\n----- Source: /amber/dat/leap/cmd/leaprc.protein.ff19SB\n"
        "----- Source: /amber/dat/leap/cmd/leaprc.water.opc\n> saveamberparm m system.prmtop system.inpcrd\n"
    )
    p = discover(tmp_path)
    assert p.mined["wat"].value == "OPC" and p.mined["wat"].confidence is Confidence.DERIVED
    assert (
        p.mined["ff"].value == ["Amber ff19SB"] and p.mined["ff"].confidence is Confidence.DERIVED
    )


def test_truncated_octahedron_and_nonperiodic(tmp_path: Path) -> None:
    tree = TreeBuilder(tmp_path)
    tree.topology(ifbox=2, box=(109.4712206, 60.0, 60.0, 60.0))
    write_nc_traj(tmp_path / "traj.nc", tree.natom, [1.0, 2.0], cell=None)
    assert discover(tmp_path).mined["boxtype"].value == "Truncated Octahedron"
    tree.topology(ifbox=0)
    m = discover(tmp_path).mined["boxtype"]
    assert m.value is None and m.method == "non-periodic"


def test_queued_control_file(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(2), st, top)
    (tree.root / "prod_0003.mdin").write_text((tree.root / "prod_0002.mdin").read_text())
    p = discover(tree.root)
    assert p.steps["prod_0003"].status == "queued"
    assert codes(p, "STEP001")
    assert p.replicas[0].production_chain == ["prod_0001", "prod_0002"]


def test_mwf_output_files_are_ignored(tree: TreeBuilder) -> None:
    top, st = tree.topology(), tree.starting()
    tree.chain("", prod_chain(1), st, top)
    (tree.root / ".mwf_cache.json").write_text("{}")
    tree.topology("topology.prmtop")  # the processed copy mwf writes
    p = discover(tree.root)
    assert [t.path for t in p.topologies] == ["system.prmtop"]


def test_relative_root_is_normalized(tree: TreeBuilder, monkeypatch: pytest.MonkeyPatch) -> None:
    top = tree.topology()
    tree.starting()
    tree.chain("runs", prod_chain(2), "system.inpcrd", top)
    monkeypatch.chdir(tree.root.parent)
    p = discover(tree.root.name)
    assert Path(p.root).is_absolute()
    assert p.replicas[0].production_chain == ["runs/prod_0001", "runs/prod_0002"]
    assert p.steps["runs/prod_0002"].files[FileRole.TRAJECTORY].path == "runs/prod_0002.nc"


def test_zero_time_trajectory_and_converter_version(tmp_path: Path) -> None:
    """Converters (MDAnalysis, cpptraj) may store time=0 everywhere and stamp their own version."""
    tree = TreeBuilder(tmp_path)
    tree.topology()
    write_nc_traj(tmp_path / "raw.nc", tree.natom, [0.0, 0.0, 0.0],
                  program="MDAnalysis.coordinates.TRJ.NCDFWriter", version="2.9.0")  # fmt: skip
    p = discover(tmp_path)
    assert not p.mined["framestep"].is_known  # never 0
    assert not p.mined["version"].is_known and not p.mined["program"].is_known
    from mddbmeta.validate.schema_lite import check_inputs

    assert check_inputs(build_export(p).inputs) == []


def test_real_mddb_dummy_project(mddb_dummy_amber: Path) -> None:
    p = discover(mddb_dummy_amber)
    (rep,) = p.replicas
    assert rep.production_chain == ["raw_trajectory"]
    ex = build_export(p)
    assert ex.inputs["input_topology_filepath"] == "ala_ala.prmtop"
    assert ex.inputs["mds"][0]["input_trajectory_filepaths"] == ["raw_trajectory.nc"]
    assert "boxtype" not in ex.inputs  # non-periodic dipeptide
