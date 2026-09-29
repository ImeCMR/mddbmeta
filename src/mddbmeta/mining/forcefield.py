"""Force-field labels from build logs (tleap and friends).

A build log that says it saved the project's own topology is linked
evidence (DERIVED); an unlinked log in the directory is only HEURISTIC --
it may describe a different system.
"""

from __future__ import annotations

import posixpath

from mddbmeta.engines.base import BuildEvidence
from mddbmeta.provenance import Confidence, Mined, Source


def mine_forcefields(
    evidence: list[BuildEvidence], topology_paths: list[str], chamber: bool = False
) -> Mined[list[str]]:
    top_names = {posixpath.basename(p) for p in topology_paths}
    linked = [
        ev
        for ev in evidence
        if ev.linked or top_names & {posixpath.basename(o) for o in ev.outputs}
    ]
    pool = linked or evidence
    labels: list[str] = []
    for ev in pool:
        for ff in ev.force_fields:
            if ff not in labels:
                labels.append(ff)
    if labels:
        return Mined(
            value=labels,
            sources=tuple(
                Source(ev.source, "force-field declaration") for ev in pool if ev.force_fields
            ),
            method="declared (linked to the topology)"
            if linked
            else "declared (not linked to the topology)",
            confidence=Confidence.DERIVED if linked else Confidence.HEURISTIC,
        )
    if chamber:
        return Mined(
            value=["CHARMM"],
            sources=tuple(Source(p, "CHAMBER topology") for p in topology_paths[:1]),
            method="heuristic: CHAMBER-converted topology",
            confidence=Confidence.HEURISTIC,
            note="topology was converted from CHARMM files; the exact CHARMM version is not recorded",
        )
    return Mined.missing(
        "no build log (tleap) found; force fields are not recorded in AMBER run files"
    )
