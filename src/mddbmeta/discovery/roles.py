"""The one role classifier (engine-neutral).

Precedence:
1. authoritative content -- a minimization flag in the control file;
2. name tokens (step name, then its directories, nearest first), matched as
   whole tokens so ``minor`` never matches ``min`` and ``premin`` is not
   ``min``;
3. content heuristics -- temperature ramp/rise -> heating, positional
   restraints -> equilibration, otherwise production.

Every non-content role is announced (ROLE000); a name that contradicts a
strong content signal is flagged (ROLE001).
"""

from __future__ import annotations

import re

from mddbmeta.findings import finding
from mddbmeta.model import FileRole, PhaseRole, Step
from mddbmeta.provenance import Confidence, Finding, Mined, Source

NAME_TOKENS: dict[str, PhaseRole] = {}
for _role, _tokens in {
    PhaseRole.MINIMIZATION: ["min", "mini", "minim", "minimization", "minimisation", "em", "opt"],
    PhaseRole.HEATING: ["heat", "heating", "warm", "warmup", "therm", "thermalization", "anneal"],
    PhaseRole.EQUILIBRATION: [
        "eq", "equi", "equil", "equilibration", "equilibrate", "relax", "nvt", "npt", "density",
    ],
    PhaseRole.PRODUCTION: ["prod", "production", "production_run"],
}.items():  # fmt: skip
    for _t in _tokens:
        NAME_TOKENS[_t] = _role

_SPLIT = re.compile(r"[^A-Za-z0-9]+|(?<=[A-Za-z])(?=\d)|(?<=\d)(?=[A-Za-z])")

# A run this long, unrestrained and not ramping, is production-like (heuristic only).
LONG_RUN_PS = 1000.0


def tokens(text: str) -> list[str]:
    return [t.lower() for t in _SPLIT.split(text) if t]


def role_from_name(step: Step) -> tuple[PhaseRole, str] | None:
    candidates = [step.name] + list(reversed(step.dir.split("/"))) if step.dir else [step.name]
    for i, text in enumerate(candidates):
        for tok in tokens(text):
            role = NAME_TOKENS.get(tok)
            if role is not None:
                return role, ("name" if i == 0 else "directory")
    return None


def _content_role(step: Step) -> tuple[PhaseRole | None, str | None, bool]:
    """(role, reason, strong) from the settings; strong means hard to argue with."""
    s = step.settings
    if s.temp_ramp.value:
        return PhaseRole.HEATING, "temperature ramp (&wt / annealing)", True
    if (
        s.temp_initial_K.is_known
        and s.temp_K.is_known
        and s.restart.value is False
        and s.temp_initial_K.value < s.temp_K.value - 1.0
    ):
        return PhaseRole.HEATING, "starts colder than its target temperature", False
    if s.restraints.value:
        return PhaseRole.EQUILIBRATION, "positional restraints active", True
    if s.nsteps.is_known and s.dt_ps.is_known:
        length_ps = s.nsteps.value * s.dt_ps.value
        if length_ps >= LONG_RUN_PS:
            return PhaseRole.PRODUCTION, f"unrestrained {length_ps:g} ps run", False
    if s.minimization.is_known and s.minimization.value is False:
        return PhaseRole.PRODUCTION, "unrestrained MD", False
    return None, None, False


def _source_of(step: Step) -> tuple[Source, ...]:
    for role in (FileRole.CONTROL, FileRole.LOG, FileRole.TRAJECTORY):
        ref = step.files.get(role)
        if ref is not None:
            return (Source(ref.path, "name/content"),)
    return ()


def classify(step: Step) -> tuple[Mined[str], list[Finding]]:
    findings: list[Finding] = []
    s = step.settings
    if s.minimization.value is True:
        return (
            Mined(
                value=PhaseRole.MINIMIZATION.value,
                sources=s.minimization.sources,
                method="content:minimization flag",
                confidence=Confidence.EXACT,
            ),
            findings,
        )
    by_name = role_from_name(step)
    content_role, reason, strong = _content_role(step)
    sources = _source_of(step)

    if by_name is not None:
        role, where = by_name
        if role is PhaseRole.MINIMIZATION and s.minimization.value is False:
            findings.append(
                finding(
                    "ROLE001",
                    f"{step.id}: named like a minimization but runs MD; using the content",
                    subject=step.id,
                )
            )
            role = content_role or PhaseRole.UNKNOWN
            where = "content"
        elif content_role is not None and strong and content_role is not role:
            findings.append(
                finding(
                    "ROLE001",
                    f"{step.id}: {where} says {role.value} but content says {content_role.value} "
                    f"({reason}); keeping {role.value}",
                    subject=step.id,
                )
            )
        return (
            Mined(
                value=role.value,
                sources=sources,
                method=f"{where} token",
                confidence=Confidence.HEURISTIC,
            ),
            findings,
        )
    if content_role is not None:
        findings.append(
            finding(
                "ROLE000",
                f"{step.id}: role {content_role.value} inferred from {reason}",
                subject=step.id,
            )
        )
        return (
            Mined(
                value=content_role.value,
                sources=sources,
                method=f"content heuristic: {reason}",
                confidence=Confidence.HEURISTIC,
            ),
            findings,
        )
    if step.status == "orphan":
        findings.append(
            finding(
                "ROLE000",
                f"{step.id}: trajectory without a log, assumed to be production",
                subject=step.id,
            )
        )
        return (
            Mined(
                value=PhaseRole.PRODUCTION.value,
                sources=sources,
                method="assumed: orphan trajectory",
                confidence=Confidence.HEURISTIC,
            ),
            findings,
        )
    return Mined(
        value=PhaseRole.UNKNOWN.value, method="unclassified", confidence=Confidence.HEURISTIC
    ), findings


_FAMILY = re.compile(r"^(?P<base>.*?)[_.\-]?(?P<idx>\d+)$")


def family_index(name: str) -> tuple[str, int] | None:
    """('prod', 4) for 'prod_0004'; ('ntp_prod', 1) for 'ntp_prod_0001'; None otherwise."""
    m = _FAMILY.match(name)
    if m is None:
        return None
    return m.group("base"), int(m.group("idx"))
