"""A stdlib mirror of MDDB-workflow's inputs schema (mddb_workflow/core/inputs_schema.py).

The real schema can't be imported without MDDB-workflow's heavy environment
(rdkit, gromacs, pytraj), so the checks that matter for generated files are
reproduced here. An optional integration test validates exports against the
real schema when MDDB-workflow is importable.
"""

from __future__ import annotations

import re
from typing import Any

from mddbmeta.findings import finding
from mddbmeta.provenance import Finding

PDB_ID_FORMAT = re.compile(r"^[1-9][a-zA-Z0-9]{3}$")
POSITIVE = ("framestep", "timestep", "temp")
STR_OR_LIST = (
    "authors",
    "groups",
    "ff",
    "boxtype",
    "pdb_ids",
    "multimeric",
    "collections",
    "input_trajectory_filepaths",
)
TYPES = ("trajectory", "ensemble")

# mddb_workflow.utils.constants
TOPOLOGY_SUPPORTED_FORMATS = {"tpr", "top", "prmtop", "psf"}
TRAJECTORY_SUPPORTED_FORMATS = {"xtc", "trr", "nc", "dcd", "crd", "pdb", "rst7"}
EXTENSION_FORMATS = {
    "tpr": "tpr", "top": "top", "psf": "psf", "prmtop": "prmtop", "parm7": "prmtop", "prm7": "prmtop",
    "pdb": "pdb", "gro": "gro", "cif": "cif",
    "xtc": "xtc", "trr": "trr", "dcd": "dcd", "nc": "nc", "cdf": "nc", "netcdf": "nc",
    "crd": "crd", "mdcrd": "crd", "trj": "crd", "rst7": "rst7",
}  # fmt: skip


def mwf_format(path: str) -> str | None:
    """The format MDDB-workflow will assign to a file, from its extension (as mwf does)."""
    ext = path.rsplit(".", 1)[-1].lower() if "." in path else ""
    return EXTENSION_FORMATS.get(ext)


def _positive(where: str, key: str, value: Any, out: list[Finding]) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        out.append(
            finding("SCH001", f"{where}{key} must be a number (got {value!r})", data={"field": key})
        )
    elif value <= 0:
        out.append(
            finding("SCH001", f"{where}{key} must be positive (got {value})", data={"field": key})
        )


def check_inputs(inputs: dict[str, Any]) -> list[Finding]:
    out: list[Finding] = []
    for key in POSITIVE:
        _positive("", key, inputs.get(key), out)
    if inputs.get("type") is not None and inputs["type"] not in TYPES:
        out.append(finding("SCH001", f"type must be one of {TYPES} (got {inputs['type']!r})"))
    for key in STR_OR_LIST:
        v = inputs.get(key)
        if v is not None and not isinstance(v, (str, list)):
            out.append(
                finding("SCH001", f"{key} must be a string or a list (got {type(v).__name__})")
            )
    ver = inputs.get("version")
    if ver is not None and (isinstance(ver, bool) or not isinstance(ver, (str, int, float))):
        out.append(finding("SCH001", f"version must be a string or number (got {ver!r})"))
    orient = inputs.get("orientation")
    if orient is not None and (not isinstance(orient, list) or len(orient) != 16):
        out.append(finding("SCH001", "orientation must be a list of 16 numbers"))
    pdb_ids = inputs.get("pdb_ids")
    for pid in [pdb_ids] if isinstance(pdb_ids, str) else (pdb_ids or []):
        if not PDB_ID_FORMAT.match(str(pid)):
            out.append(finding("SCH001", f"{pid!r} does not look like a PDB id"))

    mds = inputs.get("mds")
    if mds is not None:
        if not isinstance(mds, list):
            out.append(finding("SCH001", "mds must be a list"))
            mds = []
        names, dirs = [], []
        for i, md in enumerate(mds):
            where = f"mds[{i}]."
            if not isinstance(md, dict):
                out.append(finding("SCH001", f"mds[{i}] must be a mapping"))
                continue
            for key in POSITIVE:
                _positive(where, key, md.get(key), out)
            if md.get("name") is not None:
                names.append(md["name"])
            if md.get("mdir") is not None:
                dirs.append(str(md["mdir"]).rstrip("/"))
                if str(md["mdir"]).strip("/") in ("", "."):
                    out.append(finding("SCH001", f"{where}mdir must not be the project directory"))
            trajs = md.get("input_trajectory_filepaths")
            paths = [trajs] if isinstance(trajs, str) else (trajs or [])
            if paths:
                absolute = {str(p).startswith("/") for p in paths}
                if len(absolute) > 1:
                    out.append(
                        finding(
                            "SCH001",
                            f"{where}input_trajectory_filepaths mixes absolute and relative paths",
                        )
                    )
        for label, values in (("name", names), ("mdir", dirs)):
            dupes = sorted({v for v in values if values.count(v) > 1})
            if dupes:
                out.append(
                    finding("SCH001", f"duplicate MD {label}(s): {', '.join(map(str, dupes))}")
                )
        mdref = inputs.get("mdref")
        if mdref is not None and (not isinstance(mdref, int) or not 0 <= mdref < max(len(mds), 1)):
            out.append(finding("SCH001", f"mdref {mdref!r} is not a valid index into mds"))
    return out
