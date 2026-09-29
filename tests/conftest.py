from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from amberfiles import TreeBuilder, netcdf_available  # noqa: E402

MDDB_WORKFLOW = Path(__file__).resolve().parents[2] / "MDDB-workflow"


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--update-golden", action="store_true", help="rewrite golden files")


@pytest.fixture
def update_golden(request: pytest.FixtureRequest) -> bool:
    return bool(request.config.getoption("--update-golden"))


@pytest.fixture
def tree(tmp_path: Path) -> TreeBuilder:
    return TreeBuilder(tmp_path)


@pytest.fixture
def mddb_dummy_amber() -> Path:
    path = MDDB_WORKFLOW / "test" / "data" / "input" / "dummy" / "amber"
    if not path.is_dir():
        pytest.skip("MDDB-workflow test data not available")
    return path


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    skip_nc = pytest.mark.skip(reason="no NetCDF backend (netCDF4 or scipy)")
    skip_ext = pytest.mark.skip(reason="set MDDBMETA_EXTERNAL_DATA to run")
    for item in items:
        if "netcdf" in item.keywords and not netcdf_available():
            item.add_marker(skip_nc)
        if "external_data" in item.keywords and not os.environ.get("MDDBMETA_EXTERNAL_DATA"):
            item.add_marker(skip_ext)
