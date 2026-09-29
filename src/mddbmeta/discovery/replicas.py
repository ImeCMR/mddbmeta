"""Replica inference from the directory layout (or an explicit glob).

A replica is one independent member of the project (one ``mds[]`` entry).
The inference looks for one directory level whose sibling directories hold
"the same runs": either their names look like replica labels (``rep1``,
``replica_2``, ``run03``, ``seed4``, ``01``) and share at least one run, or
their contents are nearly identical (Jaccard >= 0.8). Directories at that
level that don't qualify (``common/``, ``prep/``) hold shared steps.

Anything it cannot reconcile -- nested replica levels (``300K/rep1``),
several parents -- is refused with REP001 rather than guessed; the user then
passes ``--replica-glob`` or edits the manifest.
"""

from __future__ import annotations

import fnmatch
import posixpath
import re
from dataclasses import dataclass, field
from itertools import combinations

from mddbmeta.findings import finding
from mddbmeta.model import Replica, Step
from mddbmeta.provenance import Finding

REPLICA_DIR = re.compile(
    r"^(?:rep|replica|repl|r|run|seed|sim|copy|traj|md)?[_.\-]?\d+$", re.IGNORECASE
)
REPLICA_IN_NAME = re.compile(
    r"(?:^|[_.\-])(?P<tag>(?:rep|replica|seed|run)[_.\-]?\d+)(?=$|[_.\-])", re.IGNORECASE
)
SIMILAR = 0.8
SINGLE_ID = "replica_1"


@dataclass
class ReplicaResult:
    replicas: list[Replica]
    assignment: dict[str, str | None]  # step id -> replica id (None = shared)
    findings: list[Finding] = field(default_factory=list)
    ambiguous: bool = False


def _segs(step: Step) -> list[str]:
    return step.dir.split("/") if step.dir else []


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def _detect_level(
    steps: list[Step], start: int = 0, edges: dict[str, str] | None = None
) -> tuple[int, dict[str, list[Step]], str] | None:
    """First directory level (>= start) whose siblings are replicas.

    Returns (level, {replica_dir: steps}, how) or None.
    """
    depth = max((len(_segs(s)) for s in steps), default=0)
    for k in range(start, depth):
        deep = [s for s in steps if len(_segs(s)) > k]
        prefixes = {tuple(_segs(s)[:k]) for s in deep}
        if len(prefixes) != 1:
            continue
        groups: dict[str, list[Step]] = {}
        for s in deep:
            groups.setdefault("/".join(_segs(s)[: k + 1]), []).append(s)
        if len(groups) < 2:
            continue
        contents = {
            g: {posixpath.join(*(_segs(s)[k + 1 :] or [""]), s.name) for s in members}
            for g, members in groups.items()
        }
        labelled = [g for g in groups if REPLICA_DIR.match(g.rsplit("/", 1)[-1])]
        if len(labelled) >= 2:
            connected = all(
                any(contents[g] & contents[h] for h in labelled if h != g) for g in labelled
            )
            # disjoint run names are still replicas when no run continues another member's run
            # (sequential chunk directories like run1/ -> run2/ are linked by their restarts)
            independent = edges is not None and not _crossing(
                edges, {g: groups[g] for g in labelled}
            )
            if connected or independent:
                return k, {g: groups[g] for g in labelled}, "labelled"
        names = list(groups)
        if all(_jaccard(contents[a], contents[b]) >= SIMILAR for a, b in combinations(names, 2)):
            return k, groups, "similar"
    return None


def _crossing(edges: dict[str, str], groups: dict[str, list[Step]]) -> bool:
    member = {s.id: g for g, ms in groups.items() for s in ms}
    return any(
        member.get(consumer) and member.get(producer) and member[consumer] != member[producer]
        for consumer, producer in edges.items()
    )


def infer_replicas(
    steps: list[Step], replica_glob: str | None = None, edges: dict[str, str] | None = None
) -> ReplicaResult:
    if not steps:
        return ReplicaResult(replicas=[], assignment={})
    if replica_glob:
        return _from_glob(steps, replica_glob)

    found = _detect_level(steps, edges=edges)
    if found is not None:
        level, groups, how = found
        # nested replica levels are ambiguous: which level is "the" replica?
        for members in groups.values():
            if _detect_level(members, start=level + 1, edges=edges) is not None:
                return ReplicaResult(
                    replicas=[],
                    assignment={s.id: None for s in steps},
                    findings=[
                        finding(
                            "REP001",
                            "replica directories are nested (e.g. condition/replica); "
                            "declare them with --replica-glob (e.g. '*/rep*') or in the manifest",
                            data={"level_dirs": sorted(groups)},
                        )
                    ],
                    ambiguous=True,
                )
        return _build(steps, groups, how)

    by_tag: dict[str, list[Step]] = {}
    for s in steps:
        m = REPLICA_IN_NAME.search(s.name)
        if m:
            by_tag.setdefault(m.group("tag").lower(), []).append(s)
    if len(by_tag) >= 2:
        stripped = {
            tag: {REPLICA_IN_NAME.sub("", s.name) for s in members}
            for tag, members in by_tag.items()
        }
        tags = list(by_tag)
        if all(stripped[a] & stripped[b] for a, b in combinations(tags, 2)):
            assignment: dict[str, str | None] = {s.id: None for s in steps}
            replicas = []
            common_dir = _common_dir(steps)
            for tag in sorted(by_tag, key=_natural):
                replicas.append(Replica(id=tag, dir=common_dir))
                for s in by_tag[tag]:
                    assignment[s.id] = tag
            _fill(replicas, steps, assignment)
            return ReplicaResult(
                replicas=replicas,
                assignment=assignment,
                findings=[
                    finding(
                        "REP004",
                        f"replicas inferred from run names: {', '.join(r.id for r in replicas)}",
                        data={"method": "name token"},
                    )
                ],
            )

    rep = Replica(id=SINGLE_ID, dir=_common_dir(steps))
    assignment = {s.id: SINGLE_ID for s in steps}
    _fill([rep], steps, assignment)
    return ReplicaResult(replicas=[rep], assignment=assignment)


def _build(steps: list[Step], groups: dict[str, list[Step]], how: str) -> ReplicaResult:
    assignment: dict[str, str | None] = {s.id: None for s in steps}
    replicas = []
    for gdir in sorted(groups, key=_natural):
        rid = gdir.rsplit("/", 1)[-1]
        replicas.append(Replica(id=rid, dir=gdir))
        for s in groups[gdir]:
            assignment[s.id] = rid
    if len({r.id for r in replicas}) != len(replicas):  # same label under different parents
        for r in replicas:
            r.id = r.dir
        for gdir, members in groups.items():
            for s in members:
                assignment[s.id] = gdir
    _fill(replicas, steps, assignment)
    findings = []
    if how == "similar":
        findings.append(
            finding(
                "REP004",
                f"replicas inferred from matching directory contents: {', '.join(r.dir for r in replicas)}",
                data={"method": "similar contents"},
            )
        )
    sets = {r.id: {posixpath.relpath(sid, r.dir) for sid in r.step_ids} for r in replicas}
    if len({frozenset(v) for v in sets.values()}) > 1:
        union = set().union(*sets.values())
        missing = {rid: sorted(union - v) for rid, v in sets.items() if union - v}
        findings.append(
            finding(
                "REP002",
                "replicas do not hold the same runs: "
                + "; ".join(f"{rid} lacks {', '.join(m)}" for rid, m in missing.items()),
                data={"missing": missing},
            )
        )
    return ReplicaResult(replicas=replicas, assignment=assignment, findings=findings)


def _from_glob(steps: list[Step], pattern: str) -> ReplicaResult:
    dirs = set()
    for s in steps:
        segs = _segs(s)
        for i in range(1, len(segs) + 1):
            d = "/".join(segs[:i])
            if fnmatch.fnmatch(d, pattern):
                dirs.add(d)
    # keep outermost matches only
    roots = sorted(
        (d for d in dirs if not any(d != o and d.startswith(o + "/") for o in dirs)), key=_natural
    )
    assignment: dict[str, str | None] = {s.id: None for s in steps}
    replicas = [Replica(id=d.rsplit("/", 1)[-1], dir=d) for d in roots]
    if len({r.id for r in replicas}) != len(replicas):
        for r in replicas:
            r.id = r.dir
    for s in steps:
        for r in replicas:
            if s.dir == r.dir or s.dir.startswith(r.dir + "/"):
                assignment[s.id] = r.id
    _fill(replicas, steps, assignment)
    findings = []
    if not replicas:
        findings.append(finding("REP001", f"--replica-glob '{pattern}' matched no run directory"))
    return ReplicaResult(
        replicas=replicas, assignment=assignment, findings=findings, ambiguous=not replicas
    )


def _fill(replicas: list[Replica], steps: list[Step], assignment: dict[str, str | None]) -> None:
    for r in replicas:
        r.step_ids = [s.id for s in steps if assignment.get(s.id) == r.id]


def steps_by_id(steps: list[Step]) -> dict[str, Step]:
    return {s.id: s for s in steps}


def _common_dir(steps: list[Step]) -> str:
    dirs = [s.dir for s in steps]
    if not dirs or any(d == "" for d in dirs):
        return ""
    common = posixpath.commonpath(dirs)
    return "" if common == "." else common


def _natural(text: str) -> list:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


natural_key = _natural
