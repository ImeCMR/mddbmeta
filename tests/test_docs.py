from pathlib import Path

from mddbmeta.findings import CATALOG

DOC = Path(__file__).resolve().parents[1] / "docs" / "findings.md"


def test_every_finding_code_is_documented() -> None:
    text = DOC.read_text()
    missing = [code for code in CATALOG if f"`{code}`" not in text]
    assert not missing, f"document these codes in docs/findings.md: {missing}"
