"""AMBER NetCDF trajectories and restarts (optional backends: netCDF4, then scipy).

Conventions: global attribute ``Conventions`` is ``AMBER`` (trajectory) or
``AMBERRESTART`` (restart); dimensions ``frame`` (trajectory only) and
``atom``; variables ``time`` (ps), ``cell_lengths`` and ``cell_angles``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any

from mddbmeta.errors import MissingDependencyError

NETCDF_MAGICS = (b"CDF\x01", b"CDF\x02", b"CDF\x05", b"\x89HDF")


def is_netcdf_bytes(head: bytes) -> bool:
    return head[:4] in NETCDF_MAGICS


def backend() -> str | None:
    try:
        import netCDF4  # noqa: F401

        return "netCDF4"
    except ImportError:
        pass
    try:
        import scipy.io  # noqa: F401

        return "scipy"
    except ImportError:
        return None


@dataclass
class NetcdfInfo:
    path: str
    conventions: str | None
    is_restart: bool
    program: str | None
    program_version: str | None
    title: str | None
    n_frames: int | None
    natom: int | None
    first_time_ps: float | None
    last_time_ps: float | None
    median_dt_ps: float | None
    cell_lengths: tuple[float, float, float] | None
    cell_angles: tuple[float, float, float] | None


def _attr(ds: Any, name: str) -> str | None:
    value = getattr(ds, name, None)
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("latin-1")
    return str(value).strip() or None


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    s = sorted(values)
    mid = len(s) // 2
    return s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2


def _open(path: str) -> tuple[Any, dict[str, int | None], dict[str, Any]]:
    kind = backend()
    if kind is None:
        raise MissingDependencyError("Reading NetCDF files", "netcdf")
    if kind == "netCDF4":
        import netCDF4

        ds = netCDF4.Dataset(path, "r")
        dims: dict[str, int | None] = {n: len(d) for n, d in ds.dimensions.items()}
        return ds, dims, dict(ds.variables)
    from scipy.io import netcdf_file

    ds = netcdf_file(path, "r", mmap=False)
    dims = dict(ds.dimensions)
    return ds, dims, dict(ds.variables)


def _values(var: Any, limit: int) -> list[float]:
    shape = getattr(var, "shape", ())
    if not shape:  # scalar (restart time)
        data = var.getValue() if hasattr(var, "getValue") else var[...]
        return [float(data)]
    n = min(shape[0], limit)
    return [float(x) for x in var[:n]]


def _row(var: Any) -> tuple[float, ...]:
    shape = getattr(var, "shape", ())
    data = var[0] if len(shape) == 2 else var[:]
    return tuple(float(x) for x in data)


def read_netcdf(path: str | Path, max_time_values: int = 200000) -> NetcdfInfo:
    ds, dims, variables = _open(str(path))
    try:
        conventions = _attr(ds, "Conventions")
        is_restart = bool(conventions and "RESTART" in conventions.upper()) or "frame" not in dims
        times = _values(variables["time"], max_time_values) if "time" in variables else []
        n_frames = None
        if not is_restart:
            n_frames = dims.get("frame")
            if n_frames is None and "time" in variables:  # unlimited dimension (scipy)
                n_frames = variables["time"].shape[0]
        cell_l = _row(variables["cell_lengths"]) if "cell_lengths" in variables else None
        cell_a = _row(variables["cell_angles"]) if "cell_angles" in variables else None
        deltas = [b - a for a, b in pairwise(times)]
        return NetcdfInfo(
            path=str(path),
            conventions=conventions,
            is_restart=is_restart,
            program=_attr(ds, "program"),
            program_version=_attr(ds, "programVersion"),
            title=_attr(ds, "title"),
            n_frames=n_frames,
            natom=dims.get("atom"),
            first_time_ps=times[0] if times else None,
            last_time_ps=times[-1] if times else None,
            median_dt_ps=_median(deltas),
            cell_lengths=cell_l,  # type: ignore[arg-type]
            cell_angles=cell_a,  # type: ignore[arg-type]
        )
    finally:
        ds.close()
