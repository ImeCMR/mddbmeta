"""Catalogue of stable finding codes.

Tests and downstream tools match on ``code``; messages may change freely.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from mddbmeta.provenance import Finding, Severity

S = Severity


@dataclass(frozen=True)
class Code:
    code: str
    severity: Severity
    title: str


CATALOG: dict[str, Code] = {
    c.code: c
    for c in [
        # Continuity
        Code("CONT000", S.INFO, "Healthy hand-off between chained steps"),
        Code("CONT001", S.ERROR, "Time gap between a step and the step it continues from"),
        Code("CONT003", S.WARNING, "Step overlaps the step it continues from"),
        Code(
            "CONT002", S.WARNING, "Restart without a resolvable producer, or time reset mid-chain"
        ),
        Code("CONT004", S.INFO, "Continuity not measured (times unavailable)"),
        Code("CONT005", S.WARNING, "Step continues the previous state but restarts the clock"),
        # Sequences
        Code("SEQ001", S.ERROR, "Numbered run family is missing members"),
        Code("SEQ002", S.ERROR, "Numbered run family has duplicate members"),
        # Completeness
        Code("CMP001", S.WARNING, "Step writes trajectory frames but no trajectory file was found"),
        Code("CMP002", S.WARNING, "Engine log is incomplete (run did not finish)"),
        Code("CMP003", S.WARNING, "Trajectory frame count disagrees with the log"),
        Code("CMP004", S.WARNING, "Trajectory time span disagrees with the log"),
        # Consistency
        Code("CONS001", S.WARNING, "Production steps disagree on a parameter"),
        Code("CONS002", S.WARNING, "Frame interval from control file disagrees with trajectory"),
        Code("CONS003", S.WARNING, "Production steps disagree on engine or version"),
        Code("CONS004", S.ERROR, "Atom counts disagree between files"),
        # Replicas / layout
        Code("REP001", S.ERROR, "Replica layout is ambiguous"),
        Code("REP002", S.INFO, "Replicas contain different step sets"),
        Code("REP003", S.INFO, "Replicas use different topologies"),
        Code("REP004", S.SUGGESTION, "Replicas inferred from the directory layout"),
        # Steps / grouping / roles / ordering
        Code("STEP001", S.WARNING, "Control file without an engine log (not run yet?)"),
        Code("STEP002", S.WARNING, "Trajectory without an engine log"),
        Code("ROLE000", S.SUGGESTION, "Step role inferred from name or content"),
        Code("ROLE001", S.WARNING, "Step name and content suggest different roles"),
        Code("ORDER001", S.WARNING, "Production order inferred without lineage evidence"),
        Code("ORDER002", S.ERROR, "User trajectory order differs from the continuity chain"),
        Code("LIN001", S.WARNING, "One restart feeds several steps (fork)"),
        Code("LIN002", S.ERROR, "Lineage cycle"),
        Code("LIN003", S.INFO, "Continuation inferred from matching energies"),
        # Mining
        Code("MINE001", S.SUGGESTION, "Field could not be mined; fill it by hand"),
        Code("MINE002", S.SUGGESTION, "Field has heuristic candidates only"),
        Code("MINE003", S.WARNING, "Field differs across replicas; written per MD"),
        Code("MINE004", S.SUGGESTION, "No topology file; MDDB-workflow runs without one"),
        # MELD / replica exchange
        Code("MELD001", S.WARNING, "Replica-exchange ladder temperatures changed during the run"),
        Code(
            "MELD002",
            S.INFO,
            "Extracted trajectories hold a different number of frames than the data store",
        ),
        Code("MELD003", S.INFO, "Walker trajectories recorded but not exported"),
        Code("MELD004", S.INFO, "The replica-exchange run was restarted"),
        # Formats / schema / loading / engines
        Code("FMT001", S.ERROR, "Trajectory parts of one MD have different formats"),
        Code("FMT002", S.WARNING, "File format not supported by MDDB-workflow"),
        Code(
            "MRG001", S.WARNING, "MDDB-workflow's trajectory merge will duplicate frames at joins"
        ),
        Code("SCH001", S.ERROR, "Exported inputs violate the MDDB inputs schema"),
        Code("LOAD001", S.WARNING, "A file could not be read or parsed"),
        Code("ENG001", S.WARNING, "Files from more than one MD engine"),
        Code("ENG002", S.WARNING, "Unsupported engine feature"),
        # Reconcile
        Code("REC001", S.WARNING, "User value conflicts with the mined value"),
    ]
}


def finding(
    code: str,
    message: str,
    subject: str | None = None,
    paths: tuple[str, ...] | list[str] = (),
    severity: Severity | None = None,
    data: dict[str, Any] | None = None,
    **extra: Any,
) -> Finding:
    """Build a Finding with the catalogue's default severity; ``data`` and ``extra`` merge."""
    entry = CATALOG[code]
    return Finding(
        code=code,
        severity=severity or entry.severity,
        message=message,
        subject=subject,
        paths=tuple(paths),
        data={**(data or {}), **extra},
    )


def worst(findings: list[Finding]) -> Severity | None:
    if not findings:
        return None
    return max((f.severity for f in findings), key=lambda s: s.rank)
