"""Build MDDB-workflow inputs (and the provenance sidecar) from a Project."""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mddbmeta import __version__
from mddbmeta.errors import MddbmetaError
from mddbmeta.findings import finding
from mddbmeta.io import yamlio
from mddbmeta.model import FileRole, Project
from mddbmeta.provenance import Confidence, Finding, Mined, Severity
from mddbmeta.validate.schema_lite import (
    TOPOLOGY_SUPPORTED_FORMATS,
    TRAJECTORY_SUPPORTED_FORMATS,
    check_inputs,
    mwf_format,
)

SIDECAR_SCHEMA_VERSION = 1

# Key order of MDDB-workflow's inputs template (resources/inputs_file_template.yml)
TEMPLATE_ORDER = [
    "name", "description", "authors", "groups", "contact", "program", "version", "type", "method",
    "license", "linkcense", "citation", "thanks", "accession", "links", "pdb_ids", "forced_references",
    "ligands", "framestep", "timestep", "temp", "ensemble", "ff", "wat", "boxtype", "interactions",
    "pbc_selection", "cg_selection", "dummy_selection", "forced_class_selections", "chainnames",
    "customs", "orientation", "multimeric", "collections", "cv19_unit", "cv19_startconf", "cv19_abs",
    "cv19_nanobs", "input_topology_filepath", "input_structure_filepath", "input_trajectory_filepaths",
    "mds", "mdref", "dataset_path", "metadditions",
]  # fmt: skip
MINED_KEYS = ["program", "version", "type", "framestep", "timestep", "temp", "ensemble", "ff", "wat",
              "boxtype", "method", "input_topology_filepath", "input_structure_filepath", "metadditions"]  # fmt: skip
MD_KEYS = ["name", "mdir", "input_topology_filepath", "input_trajectory_filepaths",
           "framestep", "timestep", "temp", "ensemble"]  # fmt: skip


@dataclass
class Export:
    inputs: dict[str, Any]
    sidecar: dict[str, Any]
    findings: list[Finding] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)  # field -> why it was not written

    @property
    def ok(self) -> bool:
        return not any(f.severity.value == "error" for f in self.findings)


def ordered(inputs: dict[str, Any]) -> dict[str, Any]:
    rank = {k: i for i, k in enumerate(TEMPLATE_ORDER)}
    return dict(sorted(inputs.items(), key=lambda kv: rank.get(kv[0], len(rank))))


def build_export(project: Project, accept_heuristic: bool = False) -> Export:
    if any(f.code == "REP001" for f in project.findings):
        raise MddbmetaError(
            "the replica layout is ambiguous (REP001); declare replicas with --replica-glob "
            "or in the manifest before exporting"
        )
    if not project.replicas:
        raise MddbmetaError(f"no MD runs found under {project.root}")
    threshold = Confidence.HEURISTIC if accept_heuristic else Confidence.DERIVED
    inputs: dict[str, Any] = {}
    fields: dict[str, Any] = {}
    skipped: dict[str, str] = {}
    findings: list[Finding] = []

    for key in MINED_KEYS:
        m = project.mined.get(key)
        if m is None:
            continue
        m = _override(project, key, m)
        fields[key] = m.to_dict()
        if m.meets(threshold):
            inputs[key] = m.value
        elif m.is_known:
            skipped[key] = (
                f"only a {m.confidence.value} guess ({m.value}); use --accept-heuristic or set it"
            )
        elif m.method not in ("per-md",):
            skipped[key] = m.note or "not mined"

    mds = []
    md_fields: dict[str, dict[str, Any]] = {}
    for rep in project.replicas:
        mined = project.md_mined.get(rep.id, {})
        entry: dict[str, Any] = {}
        md_fields[rep.id] = {}
        for key in MD_KEYS:
            m = mined.get(key)
            if m is None and f"mds.{rep.id}.{key}" in project.overrides:
                m = Mined.missing()
            if m is None:
                continue
            m = _override(project, f"mds.{rep.id}.{key}", m)
            md_fields[rep.id][key] = m.to_dict()
            if m.is_known and (
                m.meets(threshold) or key in ("name", "mdir", "input_trajectory_filepaths")
            ):
                entry[key] = m.value
        if not entry.get("input_trajectory_filepaths"):
            findings.append(
                finding(
                    "MINE001",
                    f"{rep.id}: no production trajectory found; set input_trajectory_filepaths by hand",
                    subject=rep.id,
                    data={"field": "input_trajectory_filepaths"},
                )
            )
        findings.extend(_format_checks(project, rep.id, entry))
        findings.extend(_merge_check(project, rep))
        mds.append(entry)
    inputs["mds"] = mds
    inputs["mdref"] = 0
    topo = inputs.get("input_topology_filepath")
    if (
        topo
        and str(topo).lower() not in ("no", "not", "na")
        and mwf_format(topo) not in TOPOLOGY_SUPPORTED_FORMATS
    ):
        findings.append(
            finding(
                "FMT002",
                f"MDDB-workflow will not recognize the topology {topo} by its extension; "
                "symlink it as e.g. topology.prmtop",
                paths=[topo],
            )
        )
    findings.extend(check_inputs(inputs))
    sidecar = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "mddbmeta_version": __version__,
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "root": project.root,
        "engine": project.engine,
        "accept_heuristic": accept_heuristic,
        "fields": fields,
        "mds": md_fields,
        "skipped": skipped,
        "findings": [f.to_dict() for f in project.findings + findings],
        "load_errors": [e.to_dict() for e in project.load_errors],
    }
    return Export(inputs=ordered(inputs), sidecar=sidecar, findings=findings, skipped=skipped)


def _override(project: Project, key: str, m: Mined[Any]) -> Mined[Any]:
    if key in project.overrides:
        return Mined(value=project.overrides[key], method="user", confidence=Confidence.EXACT,
                     note="set in the manifest overrides")  # fmt: skip
    return m


def _format_checks(project: Project, rep_id: str, entry: dict[str, Any]) -> list[Finding]:
    findings = []
    paths = entry.get("input_trajectory_filepaths") or []
    formats = {p: mwf_format(p) for p in paths}
    unsupported = [p for p, f in formats.items() if f not in TRAJECTORY_SUPPORTED_FORMATS]
    if unsupported:
        findings.append(
            finding(
                "FMT002",
                f"{rep_id}: MDDB-workflow will not recognize {', '.join(unsupported)} by extension; "
                "symlink with a supported extension (.nc, .crd, .xtc, ...)",
                subject=rep_id,
                paths=unsupported,
            )
        )
    if len({f for f in formats.values() if f}) > 1:
        findings.append(
            finding(
                "FMT001",
                f"{rep_id}: trajectory parts have different formats ({', '.join(sorted({str(f) for f in formats.values()}))}); "
                "MDDB-workflow can only merge parts of one format",
                subject=rep_id,
                paths=paths,
            )
        )
    return findings


def duplicated_join_frames(parts: list[tuple[float | None, float | None]]) -> list[int]:
    """Indices of joins where the next part re-writes the previous part's last frame.

    MDDB-workflow merges parts with MDTraj's ``mdconvert -o out part1 part2 ...``,
    which concatenates every frame in list order (verified on GROMACS xtc and
    NAMD dcd parts). Nothing is dropped, so a GROMACS ``-noappend`` continuation,
    whose first frame equals the previous part's last frame, leaves that frame
    twice in the merged trajectory. ``parts`` are (first frame time, last frame time).
    """
    dupes = []
    for i in range(1, len(parts)):
        prev_last, first = parts[i - 1][1], parts[i][0]
        if prev_last is not None and first is not None and abs(first - prev_last) <= 1e-6:
            dupes.append(i)
    return dupes


def _merge_check(project: Project, rep) -> list[Finding]:
    steps = [s for s in project.production_steps(rep) if FileRole.TRAJECTORY in s.files]
    if len(steps) < 2:
        return []
    parts = [
        (
            float(s.trajectory.first_time_ps.value)
            if s.trajectory.first_time_ps.is_known
            else None,
            float(s.trajectory.last_time_ps.value) if s.trajectory.last_time_ps.is_known else None,
        )
        for s in steps
    ]
    dupes = duplicated_join_frames(parts)
    if not dupes:
        return []
    joins = [f"{steps[i - 1].id} -> {steps[i].id}" for i in dupes]
    return [
        finding(
            "MRG001",
            f"{rep.id}: MDDB-workflow concatenates the parts frame by frame (mdconvert), so the "
            f"boundary frame shared at {len(dupes)} join(s) will appear twice: {', '.join(joins)}",
            subject=rep.id,
            data={"duplicated_frames": len(dupes), "joins": joins},
        )
    ]


HEADER = (
    "Generated by mddbmeta {version} from {root}\n"
    "Every value below was mined from the run files; provenance: {sidecar}\n"
    "Fill in the project metadata (name, authors, ...) by hand; see `mwf inputs` for the template."
)


def write_export(export: Export, out: Path, sidecar: bool = True) -> Path | None:
    out = Path(out)
    side_path = out.with_name(out.name + ".mddbmeta.json") if sidecar else None
    header = HEADER.format(version=__version__, root=export.sidecar["root"],
                           sidecar=side_path.name if side_path else "(no sidecar)")  # fmt: skip
    yamlio.dump(export.inputs, out, header=header)
    if side_path is not None:
        yamlio.dump(export.sidecar, side_path)
    return side_path


def relative_root_warning(project: Project, out: Path) -> Finding | None:
    """Paths are relative to the project root; warn when the file is written elsewhere."""
    if Path(out).resolve().parent == Path(project.root):
        return None
    return finding(
        "SCH001",
        f"{out} is outside the project directory; its relative paths assume {project.root}",
        severity=Severity.WARNING,
    )
