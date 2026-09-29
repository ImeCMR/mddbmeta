"""Smoke-test real run directories without vendoring them.

MDDBMETA_EXTERNAL_DATA=/path/to/run_dir[:/another] pytest -m external_data
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from mddbmeta.pipeline import discover

pytestmark = pytest.mark.external_data


def _dirs() -> list[Path]:
    return [Path(p) for p in os.environ.get("MDDBMETA_EXTERNAL_DATA", "").split(":") if p]


@pytest.mark.parametrize("index", range(8))
def test_real_directory_discovers_cleanly(index: int) -> None:
    dirs = _dirs()
    if index >= len(dirs):
        pytest.skip("no directory at this index")
    project = discover(dirs[index])
    assert project.steps, f"no runs found in {dirs[index]}"
    assert not project.load_errors, project.load_errors
    assert project.mined.get("program") is not None
