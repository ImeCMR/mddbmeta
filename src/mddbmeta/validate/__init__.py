"""Validation: continuity, sequence holes, completeness and consistency.

Every check reports findings; none raises. ``run_all`` is what discovery calls.
"""

from __future__ import annotations

from mddbmeta.model import Project
from mddbmeta.provenance import Finding
from mddbmeta.validate.completeness import check_completeness
from mddbmeta.validate.consistency import check_consistency
from mddbmeta.validate.continuity import check_continuity, check_sequences


def run_all(project: Project) -> list[Finding]:
    findings: list[Finding] = []
    findings.extend(check_continuity(project))
    findings.extend(check_sequences(project))
    findings.extend(check_completeness(project))
    findings.extend(check_consistency(project))
    return findings
