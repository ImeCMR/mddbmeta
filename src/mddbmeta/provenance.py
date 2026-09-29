"""Provenance primitives: every mined value knows where it came from.

The central rule of mddbmeta: a value that a file does not state is ``None``,
never a parser default. When an engine documents a default for an unstated
parameter, the value is recorded with ``method="engine_default"`` so the
difference between "the user wrote this" and "the engine assumed this" is
never lost.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Generic, TypeVar

T = TypeVar("T")


class Confidence(str, Enum):
    """How much a mined value can be trusted."""

    EXACT = "exact"  # read verbatim from a file (or from an explicit user override)
    DERIVED = "derived"  # computed deterministically from exact values or engine defaults
    HEURISTIC = "heuristic"  # a best guess (names, patterns); never auto-filled by default

    @property
    def rank(self) -> int:
        return {"exact": 2, "derived": 1, "heuristic": 0}[self.value]


class Severity(str, Enum):
    INFO = "info"
    SUGGESTION = "suggestion"
    WARNING = "warning"
    ERROR = "error"

    @property
    def rank(self) -> int:
        return {"info": 0, "suggestion": 1, "warning": 2, "error": 3}[self.value]


@dataclass(frozen=True)
class Source:
    """A pointer into a file: which file, and where inside it."""

    path: str
    locator: str | None = None  # e.g. "&cntrl:dt", "%FLAG POINTERS[IFBOX]", "header"

    def to_dict(self) -> dict[str, Any]:
        return {"path": self.path, "locator": self.locator}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Source:
        return cls(path=data["path"], locator=data.get("locator"))


@dataclass(frozen=True)
class Mined(Generic[T]):
    """A value plus its provenance."""

    value: T | None
    unit: str | None = None
    sources: tuple[Source, ...] = ()
    method: str = "stated"
    confidence: Confidence = Confidence.EXACT
    note: str | None = None
    alternatives: tuple[Any, ...] = ()

    @classmethod
    def missing(cls, note: str | None = None, alternatives: tuple[Any, ...] = ()) -> Mined[Any]:
        """A value nobody stated and nothing could derive."""
        return cls(
            value=None,
            method="unstated",
            confidence=Confidence.HEURISTIC,
            note=note,
            alternatives=alternatives,
        )

    @property
    def is_known(self) -> bool:
        return self.value is not None

    def meets(self, threshold: Confidence) -> bool:
        return self.is_known and self.confidence.rank >= threshold.rank

    def with_value(self, value: Any, **changes: Any) -> Mined[Any]:
        data = {
            "value": value,
            "unit": self.unit,
            "sources": self.sources,
            "method": self.method,
            "confidence": self.confidence,
            "note": self.note,
            "alternatives": self.alternatives,
        }
        data.update(changes)
        return Mined(**data)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "value": _jsonable(self.value),
            "unit": self.unit,
            "method": self.method,
            "confidence": self.confidence.value,
            "sources": [s.to_dict() for s in self.sources],
        }
        if self.note:
            out["note"] = self.note
        if self.alternatives:
            out["alternatives"] = [_jsonable(a) for a in self.alternatives]
        return out

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Mined[Any]:
        return cls(
            value=data.get("value"),
            unit=data.get("unit"),
            sources=tuple(Source.from_dict(s) for s in data.get("sources", ())),
            method=data.get("method", "stated"),
            confidence=Confidence(data.get("confidence", "exact")),
            note=data.get("note"),
            alternatives=tuple(data.get("alternatives", ())),
        )


def stated(
    value: Any, path: str, locator: str | None = None, unit: str | None = None
) -> Mined[Any]:
    """Shorthand for a value read verbatim from a file."""
    return Mined(value=value, unit=unit, sources=(Source(path, locator),))


def engine_default(
    value: Any,
    path: str,
    locator: str | None = None,
    unit: str | None = None,
    note: str | None = None,
) -> Mined[Any]:
    """Shorthand for a documented engine default applied to an unstated parameter."""
    return Mined(
        value=value,
        unit=unit,
        sources=(Source(path, locator),),
        method="engine_default",
        confidence=Confidence.DERIVED,
        note=note or "not stated in the file; engine default applied",
    )


def derived(
    value: Any,
    method: str,
    inputs: tuple[Mined[Any], ...] = (),
    unit: str | None = None,
    note: str | None = None,
    confidence: Confidence | None = None,
) -> Mined[Any]:
    """A value computed from other mined values; inherits their sources.

    Confidence is the weakest of the inputs, capped at DERIVED.
    """
    sources: list[Source] = []
    for m in inputs:
        for s in m.sources:
            if s not in sources:
                sources.append(s)
    if confidence is None:
        confidence = Confidence.DERIVED
        for m in inputs:
            if m.confidence.rank < confidence.rank:
                confidence = m.confidence
    return Mined(
        value=value,
        unit=unit,
        sources=tuple(sources),
        method=method,
        confidence=confidence,
        note=note,
    )


@dataclass(frozen=True)
class Finding:
    """Something discovery or validation wants the user to know."""

    code: str
    severity: Severity
    message: str
    subject: str | None = None  # a step id, replica id, or field name
    paths: tuple[str, ...] = ()
    data: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "code": self.code,
            "severity": self.severity.value,
            "message": self.message,
        }
        if self.subject is not None:
            out["subject"] = self.subject
        if self.paths:
            out["paths"] = list(self.paths)
        if self.data:
            out["data"] = _jsonable(self.data)
        return out

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Finding:
        return cls(
            code=data["code"],
            severity=Severity(data["severity"]),
            message=data["message"],
            subject=data.get("subject"),
            paths=tuple(data.get("paths", ())),
            data=dict(data.get("data", {})),
        )


ERROR_TYPES = ("missing", "permission", "decode", "malformed", "unsupported")


@dataclass(frozen=True)
class FileLoadError:
    """A single file that could not be read or parsed. Failures are data."""

    kind: str
    path: str
    error_type: str  # one of ERROR_TYPES
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "path": self.path,
            "error_type": self.error_type,
            "detail": self.detail,
        }


def classify_exception(exc: BaseException) -> str:
    if isinstance(exc, FileNotFoundError):
        return "missing"
    if isinstance(exc, PermissionError):
        return "permission"
    if isinstance(exc, UnicodeDecodeError):
        return "decode"
    return "malformed"


def _jsonable(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_jsonable(v) for v in value]
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    return value
