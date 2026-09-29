"""Exporter, fill/reconcile properties, manifest round-trip, CLI exit codes, golden files."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
import yaml

from amberfiles import EQUIL, HEAT, MIN, PROD, TreeBuilder, full_protocol, prod_chain
from mddbmeta.cli.main import main
from mddbmeta.export.manifest import write_manifest
from mddbmeta.export.mddb import build_export
from mddbmeta.export.reconcile import fill, reconcile
from mddbmeta.pipeline import discover, load_project
from mddbmeta.validate.schema_lite import check_inputs

pytestmark = pytest.mark.netcdf
GOLDEN = Path(__file__).resolve().parents[1] / "golden"


def replica_tree(tree: TreeBuilder) -> Path:
    top, st = tree.topology(), tree.starting()
    end = tree.chain("prep", [("min", MIN), ("heat", HEAT), ("equil", EQUIL)], st, top)
    for r in ("rep1", "rep2"):
        tree.chain(r, prod_chain(3), "prep/equil.rst7", top, start_time=end)
    (tree.root / "tleap.in").write_text(
        "source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\nsaveamberparm m system.prmtop system.inpcrd\n"
    )
    return tree.root


def _normalize(obj, root: str):
    text = json.dumps(obj, indent=2, sort_keys=True, default=str).replace(root, "<ROOT>")
    return re.sub(r'"generated_at": "[^"]+"', '"generated_at": "<TIME>"', text)


@pytest.mark.parametrize("scenario", ["linear", "replicas"])
def test_golden(tmp_path: Path, scenario: str, update_golden: bool) -> None:
    tree = TreeBuilder(tmp_path)
    if scenario == "linear":
        tree.chain("", full_protocol(3), tree.starting(), tree.topology())
    else:
        replica_tree(tree)
    p = discover(tmp_path)
    ex = build_export(p)
    outputs = {
        "discover.json": _normalize(p.to_dict(), str(tmp_path)),
        "inputs.json": _normalize(ex.inputs, str(tmp_path)),
        "sidecar.json": _normalize(ex.sidecar, str(tmp_path)),
    }
    gdir = GOLDEN / scenario
    for name, text in outputs.items():
        path = gdir / name
        if update_golden or not path.exists():
            gdir.mkdir(parents=True, exist_ok=True)
            path.write_text(text + "\n")
        assert text + "\n" == path.read_text(), (
            f"{scenario}/{name} changed; rerun with --update-golden"
        )


def test_export_passes_schema_and_has_template_order(tree: TreeBuilder) -> None:
    ex = build_export(discover(replica_tree(tree)))
    assert check_inputs(ex.inputs) == []
    keys = list(ex.inputs)
    assert keys.index("program") < keys.index("framestep") < keys.index("mds") < keys.index("mdref")
    assert ex.inputs["ff"] == ["Amber ff14SB"] and ex.inputs["wat"] == "TIP3P"


def test_heuristics_need_opt_in(tree: TreeBuilder) -> None:
    tree.chain("", prod_chain(2), tree.starting(), tree.topology())
    p = discover(tree.root)
    assert "wat" not in build_export(p).inputs
    assert build_export(p, accept_heuristic=True).inputs["wat"] == "TIP3P"


def test_unstated_values_never_claim_to_be_stated(tree: TreeBuilder) -> None:
    # temp0 and dt absent from the mdin: values come from the engine's own report or defaults
    cntrl = {k: v for k, v in PROD.items() if k not in ("temp0",)}
    tree.chain("", prod_chain(2, cntrl=cntrl), tree.starting(), tree.topology())
    for f in tree.root.glob("*.mdin"):  # also drop the echo's source of truth
        f.write_text(f.read_text())
    p = discover(tree.root)
    for step in p.steps.values():
        temp = step.settings.temp_K
        assert temp.method != "stated", temp
        assert temp.method in ("engine_report", "engine_default")


def test_fill_never_overwrites_and_is_idempotent(tree: TreeBuilder) -> None:
    ex = build_export(discover(replica_tree(tree)))
    user = {"name": "mine", "temp": 310, "wat": None, "mds": [{"name": "rep1", "mdir": "rep1",
            "input_trajectory_filepaths": ["rep1/prod_0002.nc", "rep1/prod_0001.nc"]}]}  # fmt: skip
    first = fill(user, ex)
    assert first.inputs["temp"] == 310 and first.inputs["name"] == "mine"
    assert first.inputs["wat"] == "TIP3P" and "wat" in first.filled
    assert first.inputs["mds"][0]["input_trajectory_filepaths"] == [
        "rep1/prod_0002.nc",
        "rep1/prod_0001.nc",
    ]
    assert first.skipped_mds == ["rep2"]
    second = fill(first.inputs, ex)
    assert second.inputs == first.inputs and second.filled == []
    added = fill(user, ex, add_mds=True)
    assert [m["name"] for m in added.inputs["mds"]] == ["rep1", "rep2"]


def test_reconcile_detects_order_and_value_conflicts(tree: TreeBuilder) -> None:
    root = replica_tree(tree)
    p = discover(root)
    ex = build_export(p, accept_heuristic=True)
    user = {"temp": 300, "timestep": 2, "ff": "Amber ff14SB", "mds": [
        {"name": "rep1", "mdir": "rep1", "input_trajectory_filepaths": ["prod_0003.nc", "prod_0001.nc", "prod_0002.nc"]},
        {"name": "rep2", "mdir": "rep2", "input_trajectory_filepaths": ["rep2/prod_0001.nc", "rep2/prod_0002.nc", "rep2/prod_0003.nc"]},
    ]}  # fmt: skip
    items, findings = reconcile(user, p, ex)
    status = {i.field: i.status for i in items}
    assert status["temp"] == "agree" and status["timestep"] == "agree" and status["ff"] == "agree"
    assert status["mds.rep1.input_trajectory_filepaths"] == "conflict"
    assert status["mds.rep2.input_trajectory_filepaths"] == "agree"
    assert [f.code for f in findings] == ["ORDER002"]


def test_manifest_round_trip_with_overrides(tree: TreeBuilder, tmp_path: Path) -> None:
    root = replica_tree(tree)
    p = discover(root)
    manifest = root / "project.yaml"
    write_manifest(p, manifest)
    data = yaml.safe_load(manifest.read_text())
    assert data["discovery"]["root"] == "."
    assert build_export(load_project(manifest)).inputs == build_export(p).inputs
    data["overrides"] = {
        "wat": "OPC",
        "mds.rep2.temp": 310,
        "steps.rep1/prod_0003.role": "equilibration",
    }
    manifest.write_text(yaml.safe_dump(data))
    q = load_project(manifest)
    ex = build_export(q)
    assert ex.inputs["wat"] == "OPC"
    assert ex.sidecar["fields"]["wat"]["method"] == "user"
    assert q.replicas[0].production_chain == ["rep1/prod_0001", "rep1/prod_0002"]
    assert ex.inputs["mds"][0]["input_trajectory_filepaths"] == [
        "rep1/prod_0001.nc",
        "rep1/prod_0002.nc",
    ]
    assert ex.inputs["mds"][1]["temp"] == 310


def test_cli_exit_codes(tree: TreeBuilder, capsys: pytest.CaptureFixture[str]) -> None:
    root = replica_tree(tree)
    assert main(["validate", str(root)]) == 0
    assert main(["discover", str(root), "--write", str(root / "project.yaml")]) == 0
    out = root / "inputs.yaml"
    assert main(["export", str(root / "project.yaml"), "--mddb", str(out)]) == 0
    assert out.exists() and (root / "inputs.yaml.mddbmeta.json").exists()
    exported = yaml.safe_load(out.read_text())
    assert exported["mds"][1]["name"] == "rep2"
    # user file with wrong merge order -> conflict exit code
    bad = {"mds": [{"name": "rep1", "mdir": "rep1", "input_trajectory_filepaths": "prod_*.nc"}]}
    (root / "user.yaml").write_text(yaml.safe_dump(bad))
    code = main(["reconcile", str(root / "user.yaml")])
    assert code in (0, 3)  # glob order is filesystem dependent: agree (with a warning) or conflict
    for f in root.glob("rep1/prod_0003.*"):
        f.unlink()
    for f in root.glob("rep1/prod_0002.*"):
        f.unlink()
    for f in root.glob("rep1/prod_0001.mdout"):
        f.unlink()
    assert main(["validate", str(root), "--strict"]) == 2
    assert main(["info", str(root / "system.prmtop")]) == 0
    assert main(["info", str(root / "nope")]) == 1
    capsys.readouterr()


def test_cli_fill_in_place(tree: TreeBuilder, capsys: pytest.CaptureFixture[str]) -> None:
    root = replica_tree(tree)
    user = root / "inputs.yaml"
    user.write_text("name: my run\ntemp: 305\n")
    assert main(["export", str(root), "--fill", str(user), "--in-place", "--add-mds"]) == 0
    data = yaml.safe_load(user.read_text())
    assert data["name"] == "my run" and data["temp"] == 305
    assert data["timestep"] == 2.0 and len(data["mds"]) == 2
    assert (root / "inputs.yaml.bak").read_text() == "name: my run\ntemp: 305\n"
    assert "Kept your value for: temp" in capsys.readouterr().out
