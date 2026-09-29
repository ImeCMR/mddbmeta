"""Exceptions the CLI turns into a clean message and a non-zero exit code."""


class MddbmetaError(Exception):
    """Base class for expected, user-facing failures (no traceback)."""


class StrictModeError(MddbmetaError):
    """Raised when --strict promotes a recorded problem to a hard failure."""


class MissingDependencyError(MddbmetaError):
    """An optional extra is needed for the requested operation."""

    def __init__(self, feature: str, extra: str):
        super().__init__(
            f"{feature} needs the optional '{extra}' extra: pip install 'mddbmeta[{extra}]'"
        )
        self.feature = feature
        self.extra = extra
