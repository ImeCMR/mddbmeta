"""Validate exports against MDDB-workflow's real pydantic schema, when it can be imported.

MDDB-workflow's package import needs its full environment (rdkit, gromacs, pytraj),
so this test only runs inside that environment (e.g. `conda activate mwf_env`).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from amberfiles import EQUIL, HEAT, MIN, TreeBuilder, prod_chain
from mddbmeta.export.mddb import build_export
from mddbmeta.pipeline import discover

pytestmark = [pytest.mark.mwf, pytest.mark.netcdf]


def _schema():
    try:
        from mddb_workflow.core.inputs_schema import validate_inputs
    except Exception as exc:  # ImportError, or RuntimeError when gromacs is missing
        pytest.skip(f"mddb_workflow not importable here: {exc}")
    return validate_inputs


def test_export_is_accepted_by_mwf_schema(tmp_path: Path) -> None:
    validate_inputs = _schema()
    tree = TreeBuilder(tmp_path)
    top, st = tree.topology(), tree.starting()
    end = tree.chain("prep", [("min", MIN), ("heat", HEAT), ("equil", EQUIL)], st, top)
    for r in ("rep1", "rep2"):
        tree.chain(r, prod_chain(2), "prep/equil.rst7", top, start_time=end)
    inputs = build_export(discover(tmp_path), accept_heuristic=True).inputs
    validate_inputs(dict(inputs), strict_unknown=True)


def _schema_from_source():
    """Load MDDB-workflow's inputs_schema.py source with its three light imports stubbed.

    The schema module itself only needs pydantic plus InputError, warn and
    PDB_ID_FORMAT; stubbing those avoids importing the whole mwf package.
    """
    import importlib.util
    import re
    import sys
    import types

    pytest.importorskip("pydantic")
    src = (
        Path(__file__).resolve().parents[3]
        / "MDDB-workflow"
        / "mddb_workflow"
        / "core"
        / "inputs_schema.py"
    )
    if not src.exists():
        pytest.skip("MDDB-workflow source not next to mddbmeta")
    constants = (src.parents[1] / "utils" / "constants.py").read_text()
    pdb_format = re.search(r"PDB_ID_FORMAT = r'([^']+)'", constants).group(1)

    class InputError(Exception):
        pass

    stubs = {
        "mddb_workflow": types.ModuleType("mddb_workflow"),
        "mddb_workflow.utils": types.ModuleType("mddb_workflow.utils"),
        "mddb_workflow.utils.auxiliar": types.ModuleType("mddb_workflow.utils.auxiliar"),
        "mddb_workflow.utils.constants": types.ModuleType("mddb_workflow.utils.constants"),
    }
    stubs["mddb_workflow.utils.auxiliar"].InputError = InputError
    stubs["mddb_workflow.utils.auxiliar"].warn = lambda *a, **k: None
    stubs["mddb_workflow.utils.constants"].PDB_ID_FORMAT = pdb_format
    saved = {k: sys.modules.get(k) for k in stubs}
    sys.modules.update(stubs)
    try:
        spec = importlib.util.spec_from_file_location("mwf_inputs_schema", src)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    return module.validate_inputs, InputError


@pytest.mark.parametrize("accept", [False, True])
def test_export_is_accepted_by_mwf_schema_source(tmp_path: Path, accept: bool) -> None:
    validate_inputs, _err = _schema_from_source()
    tree = TreeBuilder(tmp_path)
    top, st = tree.topology(), tree.starting()
    end = tree.chain("prep", [("min", MIN), ("heat", HEAT), ("equil", EQUIL)], st, top)
    for r in ("rep1", "rep2"):
        tree.chain(r, prod_chain(2, cntrl={**__import__("amberfiles").PROD, "temp0": 290.0 + len(r)}),
                   "prep/equil.rst7", top, start_time=end)  # fmt: skip
    inputs = build_export(discover(tmp_path), accept_heuristic=accept).inputs
    validate_inputs(dict(inputs), strict_unknown=True)


def test_schema_source_rejects_what_schema_lite_rejects() -> None:
    validate_inputs, InputError = _schema_from_source()
    from mddbmeta.validate.schema_lite import check_inputs

    for bad in (
        {"framestep": 0},
        {"type": "movie"},
        {"orientation": [1, 2]},
        {"pdb_ids": ["XXXXX"]},
    ):
        assert check_inputs(bad), bad
        with pytest.raises(InputError):
            validate_inputs(dict(bad))


def test_gromacs_export_is_accepted_by_mwf_schema_source() -> None:
    validate_inputs, _err = _schema_from_source()
    real = Path(__file__).resolve().parents[1] / "data" / "gromacs_real"
    inputs = build_export(discover(real)).inputs
    validate_inputs(dict(inputs), strict_unknown=True)


def test_namd_export_is_accepted_by_mwf_schema_source() -> None:
    validate_inputs, _err = _schema_from_source()
    real = Path(__file__).resolve().parents[1] / "data" / "namd_real"
    validate_inputs(dict(build_export(discover(real)).inputs), strict_unknown=True)


@pytest.mark.parametrize("fixture", ["openmm_real", "openmm_amber"])
def test_openmm_export_is_accepted_by_mwf_schema_source(fixture: str) -> None:
    validate_inputs, _err = _schema_from_source()
    data = Path(__file__).resolve().parents[1] / "data" / fixture
    validate_inputs(dict(build_export(discover(data)).inputs), strict_unknown=True)
