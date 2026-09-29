"""Generator for MDDB-workflow's `mwf dataset inputs -ig <this file>`.

MDDB-workflow loads this file, calls ``inputs_generator(project_dir)`` and
renders the returned keys into a Jinja2 template. Jinja would print a Python
``None`` as the text "None" and a list as its Python repr, so every value
returned here is already a YAML literal string:

* ``<field>_yaml`` for each mined field (``framestep_yaml``, ``temp_yaml``, ...),
  ``"null"`` when unknown;
* ``mds_yaml``: the mds list as one flow-style YAML value;
* ``mddbmeta_block``: all mined keys plus mds/mdref as block YAML, ready to
  paste at the top level of a template;
* ``mddbmeta_ok`` (``"true"``/``"false"``) and ``mddbmeta_findings`` (comment lines).

Environment knobs: ``MDDBMETA_REPLICA_GLOB``, ``MDDBMETA_ACCEPT_HEURISTIC=1``,
``MDDBMETA_NO_SIDECAR=1``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from mddbmeta.errors import MddbmetaError
from mddbmeta.export.mddb import MINED_KEYS, build_export
from mddbmeta.io import yamlio
from mddbmeta.pipeline import discover
from mddbmeta.provenance import Severity

SIDECAR_NAME = "inputs.yaml.mddbmeta.json"


def _comment(lines: list[str]) -> str:
    return "\n".join(f"# {line}" for line in lines) if lines else "# mddbmeta: no findings"


def inputs_generator(project_dir: str) -> dict[str, str]:
    accept = os.environ.get("MDDBMETA_ACCEPT_HEURISTIC", "") not in ("", "0", "false")
    out: dict[str, str] = {f"{key}_yaml": "null" for key in MINED_KEYS}
    out["mds_yaml"] = "[]"
    try:
        project = discover(
            project_dir, replica_glob=os.environ.get("MDDBMETA_REPLICA_GLOB") or None
        )
        export = build_export(project, accept_heuristic=accept)
    except MddbmetaError as exc:
        out["mddbmeta_ok"] = "false"
        out["mddbmeta_block"] = f"# mddbmeta could not mine this project: {exc}"
        out["mddbmeta_findings"] = _comment([str(exc)])
        return out
    for key in MINED_KEYS:
        out[f"{key}_yaml"] = yamlio.dump_scalar(export.inputs.get(key))
    out["mds_yaml"] = json.dumps(export.inputs.get("mds", []), ensure_ascii=False)
    out["mddbmeta_block"] = yamlio.dumps_yaml(export.inputs).rstrip("\n")
    problems = [
        f for f in project.findings + export.findings if f.severity.rank >= Severity.WARNING.rank
    ]
    out["mddbmeta_ok"] = "false" if any(f.severity is Severity.ERROR for f in problems) else "true"
    out["mddbmeta_findings"] = _comment(
        [f"[{f.severity.value}] {f.code} {f.message}" for f in problems]
    )
    if os.environ.get("MDDBMETA_NO_SIDECAR", "") in ("", "0", "false"):
        yamlio.dump(export.sidecar, Path(project_dir) / SIDECAR_NAME)
    return out
