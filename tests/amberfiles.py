"""Writers for small, valid, synthetic AMBER files (test fixtures only).

Everything the tests need is generated here, so no third-party run data is
vendored into the repository.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------- prmtop


def _fmt_block(values: list[Any], per_line: int, width: int, kind: str) -> str:
    lines = []
    for i in range(0, max(len(values), 1), per_line):
        chunk = values[i : i + per_line]
        if kind == "a":
            lines.append("".join(f"{str(v):<{width}}"[:width] for v in chunk))
        elif kind == "i":
            lines.append("".join(f"{int(v):>{width}d}" for v in chunk))
        else:
            lines.append("".join(f"{float(v):>{width}.8E}" for v in chunk))
    return "\n".join(lines) + "\n"


WATER_TIP3P = ("WAT", ["O", "H1", "H2"])
WATER_OPC = ("WAT", ["O", "H1", "H2", "EPW"])


def protein_residues(n: int) -> list[tuple[str, list[str]]]:
    names = ["ALA", "GLY", "SER", "LEU"]
    return [(names[i % 4], ["N", "H", "CA", "HA", "C", "O"]) for i in range(n)]


def write_prmtop(
    path: Path,
    residues: list[tuple[str, list[str]]],
    ifbox: int = 1,
    box: tuple[float, float, float, float] = (90.0, 40.0, 40.0, 40.0),
    title: str = "synthetic",
    numextra: int = 0,
    chamber: bool = False,
) -> Path:
    atom_names: list[str] = []
    labels: list[str] = []
    pointers_res: list[int] = []
    for name, atoms in residues:
        pointers_res.append(len(atom_names) + 1)
        labels.append(name)
        atom_names.extend(atoms)
    natom = len(atom_names)
    pointers = [0] * 32
    pointers[0] = natom
    pointers[11] = len(residues)
    pointers[27] = ifbox
    pointers[30] = numextra
    out = ["%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/26  00:00:00\n"]

    def section(flag: str, fmt: str, body: str) -> None:
        out.append(f"%FLAG {flag:<74}\n%FORMAT({fmt})\n{body}")

    if chamber:
        section("CTITLE", "a80", f"{title}\n")
        section("FORCE_FIELD_TYPE", "i2,a78", " 1 CHARMM  36\n")
    else:
        section("TITLE", "20a4", f"{title}\n")
    section("POINTERS", "10I8", _fmt_block(pointers, 10, 8, "i"))
    section("ATOM_NAME", "20a4", _fmt_block(atom_names, 20, 4, "a"))
    section("MASS", "5E16.8", _fmt_block([1.008] * natom, 5, 16, "e"))
    section("RESIDUE_LABEL", "20a4", _fmt_block(labels, 20, 4, "a"))
    section("RESIDUE_POINTER", "10I8", _fmt_block(pointers_res, 10, 8, "i"))
    if ifbox:
        n_solute = sum(1 for n, _ in residues if n != "WAT")
        section(
            "SOLVENT_POINTERS",
            "3I8",
            _fmt_block([n_solute, len(residues), n_solute + 1], 3, 8, "i"),
        )
        section("BOX_DIMENSIONS", "5E16.8", _fmt_block(list(box), 5, 16, "e"))
    path.write_text("".join(out))
    return path


def natom_of(residues: list[tuple[str, list[str]]]) -> int:
    return sum(len(a) for _, a in residues)


# ---------------------------------------------------------------- mdin / mdout


def mdin_text(
    cntrl: dict[str, Any], title: str = "run", wt: list[dict[str, Any]] | None = None
) -> str:
    def val(v: Any) -> str:
        if isinstance(v, str):
            return f"'{v}'"
        return str(v)

    items = ", ".join(f"{k}={val(v)}" for k, v in cntrl.items())
    text = f"{title}\n &cntrl\n  {items},\n /\n"
    for w in wt or []:
        text += " &wt " + ", ".join(f"{k}={val(v)}" for k, v in w.items()) + " /\n"
    return text


@dataclass
class RunSpec:
    cntrl: dict[str, Any]
    start_time_ps: float = 0.0
    completed: bool = True
    program: str = "PMEMD"
    version: str = "22"
    executable: str = "/opt/amber22/bin/pmemd.cuda"
    ntpr: int | None = None
    write_restart: bool = True
    write_traj: bool = True
    traj_format: str = "nc"  # "nc" | "mdcrd"
    extra_frames: int = 0  # frames to add (or drop, if negative) vs what the log implies
    wt: list[dict[str, Any]] = field(default_factory=list)
    mdin_on_disk: bool = True
    begin_time_in_log: bool | None = None  # default: when irest=1


def mdout_text(
    spec: RunSpec,
    names: dict[str, str],
    natom: int,
) -> str:
    c = spec.cntrl
    imin = int(c.get("imin", 0))
    dt = float(c.get("dt", 0.001))
    nstlim = int(c.get("nstlim", 1))
    ntpr = spec.ntpr or int(c.get("ntpr", max(1, nstlim // 5 if nstlim >= 5 else 1)))
    banner = f"Amber {spec.version} {spec.program}"
    lines = [
        "",
        "          -------------------------------------------------------",
        f"          {banner:<46}{2000 + int(spec.version) if spec.version.isdigit() else ''}",
        "          -------------------------------------------------------",
        "",
        f"| {spec.program} implementation of SANDER, Release {spec.version}",
        "",
        "|  Compiled date/time: Thu Apr 14 2022",
        f"|   Executable path: {spec.executable}",
        "",
        "File Assignments:",
    ]
    for tag in (
        "MDIN",
        "MDOUT",
        "INPCRD",
        "PARM",
        "RESTRT",
        "REFC",
        "MDVEL",
        "MDEN",
        "MDCRD",
        "MDINFO",
    ):
        lines.append(f"|{tag:>7}: {names.get(tag, tag.lower())}")
    lines += ["", "", " Here is the input file:", ""]
    lines += mdin_text(c, wt=spec.wt).splitlines()
    lines += [
        "",
        "",
        "|--------------------- INFORMATION ----------------------",
        "| GPU (CUDA) Version of PMEMD in use: NVIDIA GPU IN USE.",
        "",
        "--------------------------------------------------------------------------------",
        "   1.  RESOURCE   USE: ",
        "--------------------------------------------------------------------------------",
        "",
        f" NATOM  = {natom:>7d} NTYPES =       2 NBONH =       0 MBONA  =       0",
        "",
        "--------------------------------------------------------------------------------",
        "   2.  CONTROL  DATA  FOR  THE  RUN",
        "--------------------------------------------------------------------------------",
        "",
        "General flags:",
        f"     imin    = {imin:>7d}, nmropt  = {int(c.get('nmropt', 0)):>7d}",
        "",
        "Nature and format of input:",
        f"     ntx     = {int(c.get('ntx', 1)):>7d}, irest   = {int(c.get('irest', 0)):>7d}, ntrx    =       1",
        "",
        "Nature and format of output:",
        f"     ntxo    =       2, ntpr    = {ntpr:>7d}, ntrx    =       1, ntwr    = {int(c.get('ntwr', nstlim)):>7d}",
        f"     iwrap   =       0, ntwx    = {int(c.get('ntwx', 0)):>7d}, ntwv    =       0, ntwe    =       0",
        "",
        "Potential function:",
        f"     ntf     = {int(c.get('ntf', 1)):>7d}, ntb     = {int(c.get('ntb', 2 if c.get('ntp') else 1)):>7d}, igb     = {int(c.get('igb', 0)):>7d}, nsnb    =      25",
        f"     cut     = {float(c.get('cut', 8.0)):>9.5f}, intdiel =       0",
        "",
    ]
    if not imin:
        lines += [
            "Molecular dynamics:",
            f"     nstlim  = {nstlim:>7d}, nscm    =    1000, nrespa  =       1",
            f"     t       = {float(c.get('t', 0.0)):>9.5f}, dt      = {dt:>9.5f}, vlimit  =  -1.00000",
            "",
            "Langevin dynamics temperature regulation:",
            "     ig      =  -1",
            f"     temp0   = {float(c.get('temp0', 300.0)):>9.5f}, tempi   = {float(c.get('tempi', 0.0)):>9.5f}, gamma_ln=   1.00000",
            "",
        ]
    else:
        lines += [
            "Energy minimization:",
            f"     maxcyc  = {int(c.get('maxcyc', 1)):>7d}, ncyc    =      10, ntmin   =       1",
            "",
        ]
    lines += [
        "--------------------------------------------------------------------------------",
        "   3.  ATOMIC COORDINATES AND VELOCITIES",
        "--------------------------------------------------------------------------------",
        "",
        "default_name",
    ]
    begin = spec.begin_time_in_log
    if begin is None:
        begin = int(c.get("irest", 0)) == 1
    if begin and not imin:
        lines.append(f" begin time read from input coords = {spec.start_time_ps:10.3f} ps")
    lines += [
        "",
        "--------------------------------------------------------------------------------",
        "   4.  RESULTS",
        "--------------------------------------------------------------------------------",
        "",
    ]
    steps_done = nstlim if spec.completed else max(ntpr, nstlim // 2)
    if imin:
        lines += [
            "   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER",
            f"{int(c.get('maxcyc', 1)):>7d}      -1.2345E+03     1.2345E+00     5.0000E+00     CA         12",
            "",
        ]
    else:
        start = spec.start_time_ps
        for n in range(ntpr, steps_done + 1, ntpr):
            t = start + n * dt
            lines += [
                f" NSTEP = {n:>8d}   TIME(PS) = {t:>12.3f}  TEMP(K) = {float(c.get('temp0', 300.0)):>8.2f}  PRESS =     0.0",
                " Etot   =    -12345.6789  EKtot   =      1234.5678  EPtot      =    -13580.2467",
                " ------------------------------------------------------------------------------",
                "",
            ]
        if spec.completed:
            t = start + nstlim * dt
            lines += [
                "      A V E R A G E S   O V E R       5 S T E P S",
                "",
                f" NSTEP = {nstlim:>8d}   TIME(PS) = {t:>12.3f}  TEMP(K) = {float(c.get('temp0', 300.0)):>8.2f}  PRESS =     0.0",
                " ------------------------------------------------------------------------------",
                "",
                "      R M S  F L U C T U A T I O N S",
                "",
                f" NSTEP = {nstlim:>8d}   TIME(PS) = {t:>12.3f}  TEMP(K) =     1.00  PRESS =     0.0",
                "",
            ]
    if spec.completed:
        lines += [
            "--------------------------------------------------------------------------------",
            "   5.  TIMINGS",
            "--------------------------------------------------------------------------------",
            "",
            "|  Final Performance Info:",
            "|     ns/day =     100.00   seconds/ns =     864.00",
            "|  Total wall time:          42    seconds     0.01 hours",
        ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- coordinates


def write_rst7(
    path: Path,
    natom: int,
    time_ps: float | None,
    box: tuple[float, ...] | None = (40.0, 40.0, 40.0, 90.0, 90.0, 90.0),
    velocities: bool = True,
) -> Path:
    lines = ["default_name"]
    lines.append(f"{natom:6d}" + (f"{time_ps:15.7E}" if time_ps is not None else ""))
    coords = [1.0 + 0.001 * i for i in range(3 * natom)]
    for block in [coords] + ([coords] if velocities else []):
        for i in range(0, len(block), 6):
            lines.append("".join(f"{v:12.7f}" for v in block[i : i + 6]))
    if box is not None:
        lines.append("".join(f"{v:12.7f}" for v in box))
    path.write_text("\n".join(lines) + "\n")
    return path


def write_mdcrd(path: Path, natom: int, n_frames: int, box: bool = True) -> Path:
    lines = ["Cpptraj Generated trajectory"]
    vals = [1.0 + 0.01 * i for i in range(3 * natom)]
    for _ in range(n_frames):
        for i in range(0, len(vals), 10):
            lines.append("".join(f"{v:8.3f}" for v in vals[i : i + 10]))
        if box:
            lines.append("".join(f"{v:8.3f}" for v in (40.0, 40.0, 40.0)))
    path.write_text("\n".join(lines) + "\n")
    return path


def netcdf_available() -> bool:
    try:
        import scipy.io  # noqa: F401

        return True
    except ImportError:
        return False


def write_nc_traj(
    path: Path,
    natom: int,
    times: list[float],
    program: str = "pmemd.cuda",
    version: str = "22.0",
    cell: tuple[float, float, float] | None = (40.0, 40.0, 40.0),
) -> Path:
    from scipy.io import netcdf_file

    with netcdf_file(str(path), "w", version=2) as nc:
        nc.Conventions = "AMBER"
        nc.ConventionVersion = "1.0"
        nc.program = program
        nc.programVersion = version
        nc.title = "default_name"
        nc.createDimension("frame", None)
        nc.createDimension("spatial", 3)
        nc.createDimension("atom", natom)
        t = nc.createVariable("time", "f", ("frame",))
        t.units = "picosecond"
        xyz = nc.createVariable("coordinates", "f", ("frame", "atom", "spatial"))
        xyz.units = "angstrom"
        if cell is not None:
            nc.createDimension("cell_spatial", 3)
            nc.createDimension("cell_angular", 3)
            cl = nc.createVariable("cell_lengths", "d", ("frame", "cell_spatial"))
            ca = nc.createVariable("cell_angles", "d", ("frame", "cell_angular"))
        for i, time in enumerate(times):
            t[i] = time
            xyz[i] = [[0.0, 0.0, 0.0]] * natom
            if cell is not None:
                cl[i] = list(cell)
                ca[i] = [90.0, 90.0, 90.0]
    return path


def write_nc_restart(path: Path, natom: int, time_ps: float) -> Path:
    from scipy.io import netcdf_file

    with netcdf_file(str(path), "w", version=2) as nc:
        nc.Conventions = "AMBERRESTART"
        nc.ConventionVersion = "1.0"
        nc.program = "pmemd.cuda"
        nc.programVersion = "22.0"
        nc.createDimension("spatial", 3)
        nc.createDimension("atom", natom)
        t = nc.createVariable("time", "d", ())
        t[...] = time_ps
        xyz = nc.createVariable("coordinates", "d", ("atom", "spatial"))
        xyz[:] = [[0.0, 0.0, 0.0]] * natom
    return path


# ---------------------------------------------------------------- whole runs


def write_run(
    run_dir: Path,
    stem: str,
    spec: RunSpec,
    inpcrd: str,
    prmtop: str,
    natom: int,
) -> dict[str, Path]:
    """Write one AMBER run (mdin, mdout, restart, trajectory) as pmemd would."""
    run_dir.mkdir(parents=True, exist_ok=True)
    c = spec.cntrl
    dt = float(c.get("dt", 0.001))
    nstlim = int(c.get("nstlim", 1))
    imin = int(c.get("imin", 0))
    ntwx = int(c.get("ntwx", 0))
    written: dict[str, Path] = {}
    names = {
        "MDIN": f"{stem}.mdin",
        "MDOUT": f"{stem}.mdout",
        "INPCRD": inpcrd,
        "PARM": prmtop,
        "RESTRT": f"{stem}.rst7",
    }
    if spec.write_traj and ntwx > 0 and not imin:
        names["MDCRD"] = f"{stem}.{spec.traj_format}"
    if spec.mdin_on_disk:
        written["mdin"] = run_dir / names["MDIN"]
        written["mdin"].write_text(mdin_text(c, title=stem, wt=spec.wt))
    written["mdout"] = run_dir / names["MDOUT"]
    written["mdout"].write_text(mdout_text(spec, names, natom))
    steps_done = nstlim if spec.completed else nstlim // 2
    end_time = spec.start_time_ps + (0 if imin else steps_done * dt)
    if spec.write_restart and spec.completed:
        written["restart"] = write_rst7(
            run_dir / names["RESTRT"], natom, None if imin else round(end_time, 6)
        )
    if "MDCRD" in names:
        n_frames = steps_done // ntwx + spec.extra_frames
        times = [
            round(spec.start_time_ps + (i + 1) * ntwx * dt, 6) for i in range(max(n_frames, 0))
        ]
        target = run_dir / names["MDCRD"]
        if spec.traj_format == "nc":
            written["traj"] = write_nc_traj(target, natom, times)
        else:
            written["traj"] = write_mdcrd(target, natom, len(times))
    return written


PROD = {"imin": 0, "irest": 1, "ntx": 5, "nstlim": 5000, "dt": 0.002, "ntt": 3, "temp0": 300.0,
        "ntb": 2, "ntp": 1, "barostat": 2, "ntwx": 500, "ntpr": 1000, "cut": 9.0}  # fmt: skip
MIN = {"imin": 1, "maxcyc": 500, "ntb": 1, "ntr": 1, "cut": 9.0}
HEAT = {"imin": 0, "irest": 0, "ntx": 1, "nstlim": 2000, "dt": 0.002, "ntt": 3, "tempi": 0.0,
        "temp0": 300.0, "ntb": 1, "ntr": 1, "ntwx": 500, "ntpr": 500, "cut": 9.0}  # fmt: skip
EQUIL = {"imin": 0, "irest": 1, "ntx": 5, "nstlim": 2000, "dt": 0.002, "ntt": 3, "temp0": 300.0,
         "ntb": 2, "ntp": 1, "ntr": 1, "ntwx": 500, "ntpr": 500, "cut": 9.0}  # fmt: skip


class TreeBuilder:
    """Build realistic AMBER project trees in a temp directory."""

    def __init__(self, root: Path, residues: list[tuple[str, list[str]]] | None = None):
        self.root = root
        self.residues = residues or (protein_residues(4) + [WATER_TIP3P] * 6)
        self.natom = natom_of(self.residues)

    def topology(self, rel: str = "system.prmtop", **kw: Any) -> str:
        path = self.root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        write_prmtop(path, self.residues, **kw)
        return rel

    def starting(self, rel: str = "system.inpcrd") -> str:
        path = self.root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        write_rst7(path, self.natom, None, velocities=False)
        return rel

    def chain(
        self,
        run_dir: str,
        stems: list[tuple[str, dict[str, Any]]],
        inpcrd: str,
        prmtop: str,
        start_time: float = 0.0,
        overrides: dict[str, dict[str, Any]] | None = None,
    ) -> float:
        """Write consecutive runs, each reading the previous restart. Returns the end time."""
        d = self.root / run_dir if run_dir else self.root
        prev = inpcrd
        t = start_time
        for stem, cntrl in stems:
            opts = dict((overrides or {}).get(stem, {}))
            cntrl = {**cntrl, **opts.pop("cntrl", {})}
            if int(cntrl.get("irest", 0)) == 0 and not cntrl.get("imin"):
                t = float(cntrl.get("t", 0.0))
            spec = RunSpec(cntrl=cntrl, start_time_ps=t, **opts)
            write_run(d, stem, spec, _rel_to(prev, run_dir), _rel_to(prmtop, run_dir), self.natom)
            if not cntrl.get("imin"):
                t = t + int(cntrl["nstlim"]) * float(cntrl["dt"])
            prev = f"{run_dir}/{stem}.rst7" if run_dir else f"{stem}.rst7"
        return t


def _rel_to(target: str, run_dir: str) -> str:
    """How a run in ``run_dir`` would name ``target`` (both relative to the root)."""
    if not run_dir:
        return target
    depth = len(Path(run_dir).parts)
    t = Path(target)
    if t.parts[:depth] == Path(run_dir).parts:
        return str(Path(*t.parts[depth:]))
    return str(Path(*([".."] * depth), target))


def prod_chain(
    n: int,
    prefix: str = "prod_",
    width: int = 4,
    skip: set[int] | None = None,
    cntrl: dict[str, Any] | None = None,
) -> list[tuple[str, dict[str, Any]]]:
    skip = skip or set()
    return [
        (f"{prefix}{i:0{width}d}", dict(cntrl or PROD)) for i in range(1, n + 1) if i not in skip
    ]


def full_protocol(n_prod: int = 3) -> list[tuple[str, dict[str, Any]]]:
    return [("min", MIN), ("heat", HEAT), ("equil", EQUIL)] + prod_chain(n_prod)


__all__ = [
    "TreeBuilder",
    "RunSpec",
    "write_prmtop",
    "write_rst7",
    "write_mdcrd",
    "write_nc_traj",
    "write_nc_restart",
    "write_run",
    "mdin_text",
    "mdout_text",
    "prod_chain",
    "full_protocol",
    "PROD",
    "MIN",
    "HEAT",
    "EQUIL",
    "WATER_TIP3P",
    "WATER_OPC",
    "protein_residues",
    "natom_of",
    "math",
]
