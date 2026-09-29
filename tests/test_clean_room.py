"""mddbmeta is a clean-room implementation: it must never import or embed AmberMeta code."""

from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "mddbmeta"


def test_no_ambermeta_imports() -> None:
    offenders = []
    for path in SRC.rglob("*.py"):
        for n, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith(("import ", "from ")) and "ambermeta" in stripped:
                offenders.append(f"{path.relative_to(SRC)}:{n}: {stripped}")
    assert not offenders, offenders


def test_no_ambermeta_dependency() -> None:
    pyproject = (SRC.parents[1] / "pyproject.toml").read_text().lower()
    assert "ambermeta" not in pyproject
