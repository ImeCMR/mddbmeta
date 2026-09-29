"""Water-model detection: build-log evidence first, then residue-shape heuristics.

A residue name and site count can narrow the candidates (3-site vs 4-site vs
5-site) but cannot tell TIP3P from SPC/E; those stay HEURISTIC and are never
exported unless the user accepts heuristics.
"""

from __future__ import annotations

from mddbmeta.engines.base import BuildEvidence
from mddbmeta.model import SystemInfo
from mddbmeta.provenance import Confidence, Mined, Source

WATER_RESNAMES = {"WAT", "HOH", "TIP3", "TP3", "T3P", "SOL", "H2O", "OPC", "OPC3", "T4P", "TP4", "T4E",
                  "TIP4", "T5P", "TP5", "TIP5", "SPC", "SPCE", "SPE", "FB3", "FB4", "OP3"}  # fmt: skip

# residue-name hints (strong) and site-count candidates (weak), best first
NAME_HINTS = {
    "OPC": "OPC", "OPC3": "OPC3", "OP3": "OPC3", "T4E": "TIP4P-Ew", "T4P": "TIP4P", "TP4": "TIP4P",
    "TIP4": "TIP4P", "T5P": "TIP5P", "TP5": "TIP5P", "TIP5": "TIP5P", "SPC": "SPC/E", "SPCE": "SPC/E",
    "SPE": "SPC/E", "TIP3": "TIP3P", "TP3": "TIP3P", "T3P": "TIP3P", "FB3": "TIP3P-FB", "FB4": "TIP4P-FB",
}  # fmt: skip
BY_SITES = {
    3: ["TIP3P", "SPC/E", "OPC3", "TIP3P-FB"],
    4: ["OPC", "TIP4P-Ew", "TIP4P", "TIP4P-D", "TIP4P-FB"],
    5: ["TIP5P"],
}
SITES = {m: n for n, models in BY_SITES.items() for m in models}
SITES.update({"SPC/Eb": 3, "OPC3-pol": 3, "TIP4P/2005": 4, "SPC": 3, "TIPS3P": 3})


def water_residues(system: SystemInfo) -> dict[str, int]:
    return {name: n for name, n in system.residue_counts.items() if name.upper() in WATER_RESNAMES}


def mine_water(system: SystemInfo | None, evidence: list[BuildEvidence]) -> Mined[str]:
    if system is None:
        return Mined.missing("no topology to inspect")
    waters = water_residues(system)
    if not waters:
        return Mined(
            value=None,
            method="no water residues",
            confidence=Confidence.EXACT,
            sources=(Source(system.source, "residue labels"),),
            note="no water residues: implicit solvent or vacuum",
        )
    resname = max(waters, key=waters.get)  # the dominant water residue
    sites = system.residue_atoms.get(resname)
    atom_names = [a.upper() for a in system.residue_atom_names.get(resname, [])]
    extra_points = any(a.startswith(("EP", "MW", "LP")) or a in ("M", "OM") for a in atom_names)
    src = Source(system.source, f"residue {resname} ({sites} sites)")

    declared = [w for ev in evidence for w in ev.water_models]
    if declared:
        model = declared[0]
        linked = [ev for ev in evidence if model in ev.water_models]
        consistent = sites is None or SITES.get(model) in (None, sites)
        if consistent:
            return Mined(
                value=model,
                sources=tuple(Source(ev.source, "declared water model") for ev in linked) + (src,),
                method="declared in build log / topology",
                confidence=Confidence.DERIVED,
                note=None
                if len(set(declared)) == 1
                else f"several water models sourced: {sorted(set(declared))}",
            )
    if resname.upper() in NAME_HINTS:
        hinted = NAME_HINTS[resname.upper()]
        if sites is None or SITES.get(hinted) == sites:
            candidates = [hinted] + [m for m in BY_SITES.get(sites or 0, []) if m != hinted]
            return Mined(
                value=hinted,
                sources=(src,),
                method="heuristic: residue name",
                confidence=Confidence.HEURISTIC,
                alternatives=tuple(candidates),
                note=f"residue name {resname} suggests {hinted}",
            )
    candidates = list(BY_SITES.get(sites or 0, []))
    if sites == 4 and not extra_points:
        candidates = []
    if not candidates:
        return Mined.missing(f"water residue {resname} with {sites} sites matches no known model")
    return Mined(
        value=candidates[0],
        sources=(src,),
        method=f"heuristic: {sites}-site water",
        confidence=Confidence.HEURISTIC,
        alternatives=tuple(candidates),
        note=f"{sites}-site water; candidates: {', '.join(candidates)}",
    )
