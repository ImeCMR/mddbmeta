"""NAMD run settings, from the log's resolved Info lines (engine_report) or the config (stated)."""

from __future__ import annotations

from mddbmeta.engines.namd.config import NamdConfig
from mddbmeta.engines.namd.log import NamdLog
from mddbmeta.model import StepSettings
from mddbmeta.provenance import Confidence, Mined, Source

TEMP_PARAMS = {"langevintemp", "reassigntemp", "rescaletemp", "tcoupletemp",
               "stochrescaletemp", "loweandersentemp"}  # fmt: skip
ON = {"on", "yes", "true", "1"}


def _num(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        return float(text.split()[0])
    except (ValueError, IndexError):
        return None


def settings_from_log(log: NamdLog, rel: str) -> StepSettings:
    s = StepSettings()

    def rep(value, locator: str, unit: str | None = None) -> Mined:
        return Mined(value=value, unit=unit, sources=(Source(rel, locator),), method="engine_report",
                     confidence=Confidence.EXACT)  # fmt: skip

    dyn, mini = log.dynamics_steps, log.minimization_steps
    s.minimization = rep(bool(mini and not dyn), "TCL: Running/Minimizing")
    s.nsteps = rep(dyn or mini, "TCL: Running/Minimizing for N steps")
    if mini and dyn:
        s.minimization = s.minimization.with_value(
            False, note=f"also minimizes for {mini} steps first"
        )
    dt = _num(log.one("TIMESTEP"))
    if dt is not None and dyn:
        s.dt_ps = rep(round(dt / 1000.0, 9), "Info: TIMESTEP (fs)", "ps")
    if not dyn:
        s.ensemble = Mined(value=None, method="not_applicable", note="energy minimization")
        return s

    flags = " | ".join(sorted(log.flags))
    thermo, temp_key = None, None
    for flag, name, key in (
        ("LANGEVIN DYNAMICS ACTIVE", "langevin", "LANGEVIN TEMPERATURE"),
        (
            "STOCHASTIC RESCALING ACTIVE",
            "stochastic-velocity-rescaling",
            "STOCHASTIC RESCALING TEMPERATURE",
        ),
        ("LOWE-ANDERSEN DYNAMICS ACTIVE", "lowe-andersen", "LOWE-ANDERSEN TEMPERATURE"),
        ("TEMPERATURE COUPLING ACTIVE", "berendsen", "TEMPERATURE COUPLING TEMP"),
    ):
        if flag in log.flags:
            thermo, temp_key = name, key
            break
    if thermo is None and log.one("VELOCITY REASSIGNMENT FREQ"):
        thermo, temp_key = "velocity-reassignment", "VELOCITY REASSIGNMENT TEMP"
    s.thermostat = rep(thermo or "none", f"Info flags: {flags[:80]}")
    changes = [(k.lower(), v) for k, v in log.parameter_changes if k.lower() in TEMP_PARAMS]
    if thermo is not None:
        final = _num(changes[-1][1]) if changes else _num(log.one(temp_key)) if temp_key else None
        if final is not None:
            s.temp_K = rep(
                final, "TCL: Setting parameter ... (last)" if changes else f"Info: {temp_key}", "K"
            )
    ramp = len({v for _, v in changes}) > 1 or bool(log.one("VELOCITY REASSIGNMENT INCR"))
    s.temp_ramp = rep(ramp, "TCL: Setting parameter <temperature>")
    initial = _num(log.one("INITIAL TEMPERATURE"))
    if initial is not None:
        s.temp_initial_K = rep(initial, "Info: INITIAL TEMPERATURE", "K")
    s.restart = rep(bool(log.one("VELOCITY FILE")), "Info: VELOCITY FILE (binvelocities)")

    baro = None
    for flag, name in (
        ("LANGEVIN PISTON PRESSURE CONTROL ACTIVE", "langevin-piston"),
        ("BERENDSEN PRESSURE COUPLING ACTIVE", "berendsen"),
        ("MONTE CARLO PRESSURE CONTROL ACTIVE", "monte-carlo"),
    ):
        if flag in log.flags:
            baro = name
            break
    s.barostat = rep(baro or "none", "Info flags")
    s.constant_pressure = rep(baro is not None, "Info flags")
    s.periodic = rep(bool(log.one("PERIODIC CELL BASIS 1")), "Info: PERIODIC CELL BASIS")
    s.restraints = rep(
        "HARMONIC CONSTRAINTS ACTIVE" in log.flags, "Info: HARMONIC CONSTRAINTS ACTIVE"
    )
    freq = _num(log.one("DCD FREQUENCY"))
    if freq is not None:
        s.traj_interval_steps = rep(int(freq), "Info: DCD FREQUENCY")
    cut = _num(
        log.one("CUTOFF") or log.one("SWITCHING OFF")
    )  # with switching, NAMD prints the cutoff this way
    if cut is not None:
        s.cutoff_A = rep(cut, "Info: CUTOFF", "A")
    if thermo is not None:
        ens = "NPT" if baro else "NVT"
    else:
        ens = "NPH" if baro else "NVE"
    s.ensemble = Mined(value=ens, sources=(Source(rel, "Info flags"),), method="derived: thermostat/barostat flags",
                       confidence=Confidence.DERIVED)  # fmt: skip
    return s


def settings_from_config(cfg: NamdConfig, rel: str) -> StepSettings:
    """For a config that has not run: what it states (no Tcl evaluation)."""
    s = StepSettings()

    def st(value, key: str, unit: str | None = None) -> Mined:
        return Mined(value=value, unit=unit, sources=(Source(rel, key),))

    dyn = sum(n for k, n in cfg.runs if k == "run")
    mini = sum(n for k, n in cfg.runs if k == "minimize")
    s.minimization = st(bool(mini and not dyn), "run/minimize")
    if dyn or mini:
        s.nsteps = st(dyn or mini, "run/minimize")
    ts = _num(cfg.get("timestep"))
    if ts is not None and dyn:
        s.dt_ps = st(ts / 1000.0, "timestep", "ps")
    langevin = (cfg.get("langevin") or "").lower() in ON
    if langevin:
        t = _num(cfg.get("langevintemp"))
        if t is not None:
            s.temp_K = st(t, "langevinTemp", "K")
    piston = (cfg.get("langevinpiston") or "").lower() in ON
    if dyn:
        s.ensemble = Mined(value=("NPT" if piston else "NVT") if langevin else ("NPH" if piston else "NVE"),
                           sources=(Source(rel, "langevin/langevinPiston"),), method="derived: config",
                           confidence=Confidence.DERIVED)  # fmt: skip
    s.restraints = st((cfg.get("constraints") or "").lower() in ON, "constraints")
    freq = _num(cfg.get("dcdfreq"))
    if freq is not None:
        s.traj_interval_steps = st(int(freq), "dcdfreq")
    return s
