"""tleap scripts and logs: force-field and water-model evidence.

Looks for ``source leaprc.<family>.<name>`` commands (scripts), ``----- Source:
<path>/leaprc.<...>`` lines (leap.log) and ``saveamberparm <unit> <prmtop>
<inpcrd>`` (to link the log to the topology it produced).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from mddbmeta.io.safe import open_text

_SOURCE_CMD = re.compile(r"^\s*source\s+(\S*leaprc\S*)", re.IGNORECASE | re.MULTILINE)
_SOURCE_LOG = re.compile(r"-----\s*Source:\s*(\S*leaprc\S*)", re.IGNORECASE)
_SAVE = re.compile(
    r"^\s*(?:>\s*)?saveamberparm\s+\S+\s+(\S+)\s+(\S+)", re.IGNORECASE | re.MULTILINE
)

# leaprc stem -> MDDB-style force-field label (unknown stems fall back to a generic label)
FORCE_FIELDS = {
    "ff14sb": "Amber ff14SB",
    "ff19sb": "Amber ff19SB",
    "ff99sb": "Amber ff99SB",
    "ff99sbildn": "Amber ff99SB-ILDN",
    "ff15ipq": "Amber ff15ipq",
    "fb15": "Amber ff15FB",
    "ff03": "Amber ff03",
    "ff03.r1": "Amber ff03",
    "ff10": "Amber ff10",
    "ff12sb": "Amber ff12SB",
    "ol15": "Amber OL15",
    "ol21": "Amber OL21",
    "bsc1": "Amber parmbsc1",
    "bsc0": "Amber parmbsc0",
    "ol3": "Amber OL3",
    "yil": "Amber RNA.YIL",
    "gaff": "GAFF",
    "gaff2": "GAFF2",
    "glycam_06j-1": "GLYCAM-06j",
    "glycam_06j": "GLYCAM-06j",
    "glycam_06ewh": "GLYCAM-06EWH",
    "lipid21": "Lipid21",
    "lipid17": "Lipid17",
    "lipid14": "Lipid14",
    "phosaa19sb": "phosaa19SB",
    "phosaa14sb": "phosaa14SB",
    "modrna08": "modrna08",
}
WATER_MODELS = {
    "tip3p": "TIP3P",
    "tip3pfb": "TIP3P-FB",
    "tip4pew": "TIP4P-Ew",
    "tip4pd": "TIP4P-D",
    "tip4pfb": "TIP4P-FB",
    "tip5p": "TIP5P",
    "opc": "OPC",
    "opc3": "OPC3",
    "opc3pol": "OPC3-pol",
    "spce": "SPC/E",
    "spceb": "SPC/Eb",
    "fb3": "TIP3P-FB",
    "fb4": "TIP4P-FB",
}


@dataclass
class LeapEvidence:
    path: str
    leaprcs: list[str] = field(default_factory=list)
    force_fields: list[str] = field(default_factory=list)
    water_models: list[str] = field(default_factory=list)
    saved_topologies: list[str] = field(default_factory=list)


def classify_leaprc(name: str) -> tuple[str, str] | None:
    """('ff'|'water', label) for a leaprc file name, or None."""
    base = Path(name).name
    if not base.lower().startswith("leaprc"):
        return None
    parts = base.split(".", 1)
    rest = parts[1] if len(parts) > 1 else ""
    low = rest.lower()
    if low.startswith("water."):
        key = low[len("water.") :]
        return "water", WATER_MODELS.get(key, key.upper())
    for prefix in (
        "protein.",
        "dna.",
        "rna.",
        "lipid",
        "gaff",
        "glycam",
        "phosaa",
        "modrna",
        "ff",
        "constph",
        "conste",
    ):
        if low.startswith(prefix):
            key = low.split(".", 1)[1] if prefix.endswith(".") else low
            if key in FORCE_FIELDS:
                return "ff", FORCE_FIELDS[key]
            if low in FORCE_FIELDS:
                return "ff", FORCE_FIELDS[low]
            return "ff", rest
    if low in FORCE_FIELDS:
        return "ff", FORCE_FIELDS[low]
    return ("ff", rest) if rest else None


def parse_leap_text(text: str, path: str) -> LeapEvidence:
    ev = LeapEvidence(path=path)
    names = _SOURCE_CMD.findall(text) + _SOURCE_LOG.findall(text)
    for name in names:
        base = Path(name).name
        if base in ev.leaprcs:
            continue
        ev.leaprcs.append(base)
        kind = classify_leaprc(base)
        if kind is None:
            continue
        target = ev.water_models if kind[0] == "water" else ev.force_fields
        if kind[1] not in target:
            target.append(kind[1])
    for prmtop, _crd in _SAVE.findall(text):
        name = Path(prmtop).name
        if name not in ev.saved_topologies:
            ev.saved_topologies.append(name)
    return ev


def read_leap(path: str | Path) -> LeapEvidence:
    with open_text(path) as fh:
        return parse_leap_text(fh.read(), str(path))


def looks_like_leap(head: bytes, name: str) -> bool:
    text = head.decode("latin-1", errors="replace")
    low = name.lower()
    if "Welcome to LEaP" in text or "-----  Source:" in text or "----- Source:" in text:
        return True
    if re.search(r"^\s*source\s+\S*leaprc", text, re.IGNORECASE | re.MULTILINE):
        return True
    return "leap" in low and ("saveamberparm" in text.lower() or "loadpdb" in text.lower())
