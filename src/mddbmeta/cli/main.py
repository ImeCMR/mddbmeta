"""mddbmeta command line.

Exit codes: 0 ok; 1 usage or expected failure; 2 validation errors
(warnings too under --strict); 3 reconcile conflicts.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import shutil
import sys
from pathlib import Path
from typing import Any

from mddbmeta import __version__
from mddbmeta.cli.render import findings_block, render_project, render_reconcile
from mddbmeta.engines.registry import available_adapters, sniff_all
from mddbmeta.errors import MddbmetaError
from mddbmeta.export.manifest import write_manifest
from mddbmeta.export.mddb import build_export, relative_root_warning, write_export
from mddbmeta.export.reconcile import fill, has_conflicts, reconcile
from mddbmeta.io import yamlio
from mddbmeta.io.safe import read_head
from mddbmeta.pipeline import load_project
from mddbmeta.provenance import FileLoadError, Finding, Severity

EXIT_OK, EXIT_FAIL, EXIT_INVALID, EXIT_CONFLICT = 0, 1, 2, 3
PACKAGE_DIR = Path(__file__).resolve().parents[1]


def _discovery_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--replica-glob", help="declare replica directories, e.g. 'rep*' or '*/rep*'")
    p.add_argument("--engine", help="restrict to one engine adapter (default: all)")
    p.add_argument(
        "--exclude", action="append", default=None, metavar="GLOB", help="skip matching paths"
    )
    p.add_argument(
        "--max-depth", type=int, default=None, help="directory depth to scan (default 4)"
    )


def _options(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "replica_glob": args.replica_glob,
        "engine": args.engine,
        "exclude": args.exclude,
        "max_depth": args.max_depth,
    }


def _exit_for(findings: list[Finding], strict: bool) -> int:
    worst = max((f.severity.rank for f in findings), default=-1)
    if worst >= Severity.ERROR.rank or (strict and worst >= Severity.WARNING.rank):
        return EXIT_INVALID
    return EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mddbmeta",
        description="Mine MDDB inputs metadata from raw MD run files, with provenance.",
    )
    parser.add_argument("--version", action="version", version=f"mddbmeta {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="also show info findings")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-v", "--verbose", action="store_true", default=argparse.SUPPRESS,
                        help="also show info findings")  # fmt: skip
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("info", parents=[common], help="parse one file and print what it says")
    p.add_argument("file")
    p.add_argument("--json", action="store_true")

    p = sub.add_parser("discover", parents=[common], help="reconstruct the project in a directory")
    p.add_argument("target", help="project directory or project.yaml manifest")
    p.add_argument(
        "--write", metavar="PROJECT_YAML", help="write a manifest (snapshot + overrides)"
    )
    p.add_argument("--json", action="store_true")
    _discovery_args(p)

    p = sub.add_parser(
        "validate", parents=[common], help="report continuity/completeness/consistency findings"
    )
    p.add_argument("target", help="project directory or project.yaml manifest")
    p.add_argument("--strict", action="store_true", help="fail on warnings too")
    p.add_argument("--json", action="store_true")
    _discovery_args(p)

    p = sub.add_parser("export", parents=[common], help="write an MDDB-workflow inputs file")
    p.add_argument("target", help="project directory or project.yaml manifest")
    p.add_argument("--mddb", metavar="OUT_YAML", help="inputs file to write (e.g. inputs.yaml)")
    p.add_argument(
        "--fill", metavar="EXISTING", help="fill only the empty fields of an existing inputs file"
    )
    p.add_argument(
        "--in-place", action="store_true", help="with --fill: rewrite EXISTING (keeps a .bak)"
    )
    p.add_argument("--add-mds", action="store_true", help="with --fill: append MDs the file lacks")
    p.add_argument(
        "--accept-heuristic", action="store_true", help="also write best guesses (wat, ff, ...)"
    )
    p.add_argument("--no-sidecar", action="store_true", help="do not write the provenance sidecar")
    p.add_argument(
        "--force", action="store_true", help="write even if the result violates the schema"
    )
    p.add_argument("--strict", action="store_true", help="exit 2 on warnings too")
    _discovery_args(p)

    p = sub.add_parser(
        "reconcile", parents=[common], help="compare an existing inputs file with the run files"
    )
    p.add_argument("inputs", help="existing inputs.yaml")
    p.add_argument("--project", help="project directory (default: the inputs file's directory)")
    p.add_argument("--json", action="store_true")
    _discovery_args(p)

    sub.add_parser(
        "template-path", help="print the Jinja2 inputs template for `mwf dataset inputs -it`"
    )
    sub.add_parser("generator-path", help="print the generator for `mwf dataset inputs -ig`")
    return parser


def cmd_info(args: argparse.Namespace) -> int:
    path = Path(args.file)
    if not path.is_file():
        raise MddbmetaError(f"{path} is not a file")
    head = read_head(path)
    if isinstance(head, FileLoadError):
        raise MddbmetaError(f"{path}: {head.error_type}: {head.detail}")
    adapters = [cls() for cls in available_adapters().values()]
    claim = sniff_all(path, head, adapters)
    if claim is None:
        raise MddbmetaError(f"{path}: not a file any engine adapter recognizes")
    adapter, sniffed = claim
    pf = adapter.parse(path, path.name, sniffed)
    payload: dict[str, Any] = {
        "file": str(path),
        "engine": adapter.name,
        "kind": pf.kind.value,
        "format": pf.format,
        "errors": [e.to_dict() for e in pf.errors],
    }
    if pf.ok:
        payload["content"] = _summarize(pf.result)
    if args.json:
        print(json.dumps(payload, indent=2, default=str))
    else:
        print(f"{path}  [{adapter.name} {pf.kind.value}/{pf.format}]")
        for err in pf.errors:
            print(f"  error: {err.error_type}: {err.detail}")
        for key, value in (payload.get("content") or {}).items():
            print(f"  {key}: {value}")
    return EXIT_OK if pf.ok else EXIT_FAIL


def _summarize(obj: Any) -> dict[str, Any]:
    data = dataclasses.asdict(obj) if dataclasses.is_dataclass(obj) else dict(vars(obj))
    out: dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, (list, tuple)) and len(value) > 12:
            out[key] = f"[{len(value)} items: {', '.join(map(str, value[:8]))}, ...]"
        elif isinstance(value, dict) and len(value) > 12:
            items = list(value.items())[:8]
            out[key] = f"{{{len(value)} entries: {', '.join(f'{k}={v}' for k, v in items)}, ...}}"
        elif isinstance(value, str) and len(value) > 200:
            out[key] = value[:200] + "..."
        else:
            out[key] = value
    return out


def cmd_discover(args: argparse.Namespace) -> int:
    project = load_project(args.target, **_options(args))
    if args.write:
        disc = {k: v for k, v in _options(args).items() if v is not None and k != "engine"}
        write_manifest(project, args.write, discovery=disc)
    if args.json:
        print(json.dumps(project.to_dict(), indent=2, default=str))
    else:
        print(render_project(project, args.verbose))
        if args.write:
            print(f"\nWrote manifest: {args.write}")
    return EXIT_OK


def cmd_validate(args: argparse.Namespace) -> int:
    project = load_project(args.target, **_options(args))
    if args.json:
        print(json.dumps([f.to_dict() for f in project.findings], indent=2, default=str))
    else:
        print("\n".join(findings_block(project.findings, args.verbose)))
    return _exit_for(project.findings, args.strict)


def cmd_export(args: argparse.Namespace) -> int:
    if not args.mddb and not (args.fill and args.in_place):
        raise MddbmetaError("give --mddb OUT_YAML (or --fill EXISTING --in-place)")
    project = load_project(args.target, **_options(args))
    export = build_export(project, accept_heuristic=args.accept_heuristic)
    out = Path(args.fill) if (args.fill and args.in_place) else Path(args.mddb)
    extra: list[Finding] = []
    w = relative_root_warning(project, out)
    if w is not None:
        extra.append(w)
    schema_errors = [
        f
        for f in export.findings
        if f.code in ("SCH001", "FMT001") and f.severity is Severity.ERROR
    ]
    if schema_errors and not args.force:
        print("\n".join(findings_block(project.findings + export.findings, args.verbose)))
        raise MddbmetaError(
            "not writing an inputs file that MDDB-workflow would reject (use --force to override)"
        )
    if args.fill:
        existing = yamlio.load(args.fill) or {}
        if not isinstance(existing, dict):
            raise MddbmetaError(f"{args.fill} is not a mapping")
        res = fill(existing, export, add_mds=args.add_mds)
        export.inputs = res.inputs
        if args.in_place:
            shutil.copy2(args.fill, args.fill + ".bak")
        sidecar = write_export(export, out, sidecar=not args.no_sidecar)
        print(f"Filled {len(res.filled)} field(s): {', '.join(res.filled) or '-'}")
        if res.kept:
            print(f"Kept your value for: {', '.join(res.kept)}  (see `mddbmeta reconcile`)")
        if res.added_mds:
            print(f"Added MDs: {', '.join(res.added_mds)}")
        if res.skipped_mds:
            print(f"MDs not in the file (use --add-mds): {', '.join(res.skipped_mds)}")
    else:
        sidecar = write_export(export, out, sidecar=not args.no_sidecar)
    print(f"Wrote {out}" + (f" (provenance: {sidecar})" if sidecar else ""))
    for key, why in export.skipped.items():
        print(f"  not written: {key}: {why}")
    all_findings = project.findings + export.findings + extra
    print("\n".join(findings_block(all_findings, args.verbose)))
    return _exit_for(all_findings, args.strict)


def cmd_reconcile(args: argparse.Namespace) -> int:
    inputs_path = Path(args.inputs)
    existing = yamlio.load(inputs_path) or {}
    if not isinstance(existing, dict):
        raise MddbmetaError(f"{inputs_path} is not a mapping")
    target = args.project or str(inputs_path.resolve().parent)
    project = load_project(target, **_options(args))
    export = build_export(project, accept_heuristic=True)
    items, findings = reconcile(existing, project, export)
    if args.json:
        print(json.dumps({"items": [i.to_dict() for i in items], "findings": [f.to_dict() for f in findings]},
                         indent=2, default=str))  # fmt: skip
    else:
        print(render_reconcile(items, findings))
    return EXIT_CONFLICT if has_conflicts(items, findings) else EXIT_OK


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "template-path":
            print(PACKAGE_DIR / "integrations" / "templates" / "inputs.mddbmeta.yml.j2")
            return EXIT_OK
        if args.command == "generator-path":
            print(PACKAGE_DIR / "integrations" / "mwf_generator.py")
            return EXIT_OK
        handler = {
            "info": cmd_info,
            "discover": cmd_discover,
            "validate": cmd_validate,
            "export": cmd_export,
            "reconcile": cmd_reconcile,
        }[args.command]
        return handler(args)
    except MddbmetaError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return EXIT_FAIL


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
