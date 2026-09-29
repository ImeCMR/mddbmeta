"""Turn AMBER &cntrl parameters into engine-neutral StepSettings.

Resolution order for every parameter:

1. stated in the mdin (``method="stated"``, EXACT);
2. reported by the engine in the mdout "CONTROL DATA" section
   (``method="engine_report"``, EXACT -- this is what actually ran);
3. the documented AMBER default (``method="engine_default"``, DERIVED).

Defaults are spelled out here, in one place, rather than hidden in parsers.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from mddbmeta.model import StepSettings
from mddbmeta.provenance import Confidence, Mined, Source, derived, engine_default

THERMOSTATS = {
    0: "none",
    1: "berendsen",
    2: "andersen",
    3: "langevin",
    9: "optimized-isokinetic-nose-hoover",
    10: "stochastic-isokinetic-nose-hoover-respa",
    11: "bussi",
}
BAROSTATS = {1: "berendsen", 2: "monte-carlo"}


def _default_ntb(ctx: Callable[[str], Any]) -> int:
    # AMBER: ntb defaults to 0 with implicit solvent, 2 with constant pressure, else 1.
    if (ctx("igb") or 0) > 0:
        return 0
    if (ctx("ntp") or 0) > 0:
        return 2
    return 1


def _default_cut(ctx: Callable[[str], Any]) -> float:
    return 9999.0 if (ctx("igb") or 0) > 0 else 8.0


def _default_ntwr(ctx: Callable[[str], Any]) -> int | None:
    return None  # differs between engines/versions; leave unstated


DEFAULTS: dict[str, Callable[[Callable[[str], Any]], Any]] = {
    "imin": lambda c: 0,
    "nstlim": lambda c: 1,
    "maxcyc": lambda c: 1,
    "dt": lambda c: 0.001,
    "ntt": lambda c: 0,
    "temp0": lambda c: 300.0,
    "tempi": lambda c: 0.0,
    "ntp": lambda c: 0,
    "barostat": lambda c: 1,
    "igb": lambda c: 0,
    "ntb": _default_ntb,
    "ntr": lambda c: 0,
    "irest": lambda c: 0,
    "ntx": lambda c: 1,
    "ntwx": lambda c: 0,
    "ntwr": _default_ntwr,
    "cut": _default_cut,
    "nmropt": lambda c: 0,
    "t": lambda c: 0.0,
}


class CntrlResolver:
    def __init__(
        self,
        stated_values: dict[str, Any],
        stated_path: str | None,
        stated_locator_prefix: str = "&cntrl",
        reported: dict[str, Any] | None = None,
        reported_path: str | None = None,
    ):
        self.stated = {k.lower(): v for k, v in (stated_values or {}).items()}
        self.stated_path = stated_path
        self.prefix = stated_locator_prefix
        self.reported = {k.lower(): v for k, v in (reported or {}).items()}
        self.reported_path = reported_path
        self._cache: dict[str, Mined[Any]] = {}

    def value(self, key: str) -> Any:
        return self.get(key).value

    def get(self, key: str, unit: str | None = None) -> Mined[Any]:
        key = key.lower()
        if key in self._cache:
            return self._cache[key]
        if key in self.stated and self.stated[key] is not None and self.stated_path:
            m: Mined[Any] = Mined(
                value=self.stated[key],
                unit=unit,
                sources=(Source(self.stated_path, f"{self.prefix}:{key}"),),
            )
        elif key in self.reported and self.reported_path:
            m = Mined(
                value=self.reported[key],
                unit=unit,
                sources=(Source(self.reported_path, f"CONTROL DATA:{key}"),),
                method="engine_report",
                confidence=Confidence.EXACT,
            )
        elif key in DEFAULTS:
            value = DEFAULTS[key](self.value)
            src = self.stated_path or self.reported_path or ""
            m = (
                engine_default(value, src, f"{self.prefix}:{key}", unit=unit)
                if value is not None
                else Mined.missing()
            )
        else:
            m = Mined.missing()
        self._cache[key] = m
        return m


def settings_from_cntrl(res: CntrlResolver, wt: list[dict[str, Any]] | None = None) -> StepSettings:
    s = StepSettings()
    imin = res.get("imin")
    minimization = imin.with_value(bool(imin.value == 1))
    s.minimization = minimization
    if minimization.value:
        s.nsteps = res.get("maxcyc")
        s.ensemble = Mined(value=None, method="not_applicable", note="energy minimization")
    else:
        s.nsteps = res.get("nstlim")
        s.dt_ps = res.get("dt", unit="ps")

    ntt = res.get("ntt")
    thermostat_name = THERMOSTATS.get(int(ntt.value), f"ntt={ntt.value}") if ntt.is_known else None
    s.thermostat = ntt.with_value(thermostat_name, method=f"map:ntt ({ntt.method})")
    thermo_on = bool(ntt.value) if ntt.is_known else False

    if thermo_on and not minimization.value:
        s.temp_K = res.get("temp0", unit="K")
    tempi = res.get("tempi", unit="K")
    if not minimization.value:
        s.temp_initial_K = tempi

    ntp = res.get("ntp")
    ntb = res.get("ntb")
    igb = res.get("igb")
    s.periodic = ntb.with_value(bool(ntb.value and ntb.value > 0), method=f"map:ntb ({ntb.method})")
    cp = bool((ntp.value or 0) > 0 or ntb.value == 2)
    s.constant_pressure = derived(cp, "ntp>0 or ntb==2", (ntp, ntb))
    if cp:
        baro = res.get("barostat")
        name = BAROSTATS.get(int(baro.value), f"barostat={baro.value}") if baro.is_known else None
        s.barostat = baro.with_value(name, method=f"map:barostat ({baro.method})")
    else:
        s.barostat = derived("none", "no pressure coupling", (ntp, ntb))
    s.implicit_solvent = igb.with_value(
        bool((igb.value or 0) > 0), method=f"map:igb ({igb.method})"
    )

    ntr = res.get("ntr")
    s.restraints = ntr.with_value(bool(ntr.value == 1), method=f"map:ntr ({ntr.method})")
    irest = res.get("irest")
    s.restart = irest.with_value(bool(irest.value == 1), method=f"map:irest ({irest.method})")
    if not minimization.value:
        s.traj_interval_steps = res.get("ntwx")
    ntwr = res.get("ntwr")
    if ntwr.is_known:
        s.restart_interval_steps = ntwr
    s.cutoff_A = res.get("cut", unit="A")

    if not minimization.value:
        if cp:
            ens = "NPT" if thermo_on else "NPH"
        else:
            ens = "NVT" if thermo_on else "NVE"
        s.ensemble = derived(ens, "ensemble from ntt/ntb/ntp", (ntt, ntb, ntp))

    # Temperature ramp through the &wt weight-change namelists (needs nmropt > 0).
    nmropt = res.get("nmropt")
    ramp = False
    for w in wt or []:
        wtype = str(w.get("type", "")).strip().upper()
        if wtype == "TEMP0":
            v1, v2 = w.get("value1"), w.get("value2")
            if isinstance(v1, (int, float)) and isinstance(v2, (int, float)) and v1 != v2:
                ramp = True
    if wt is not None and not minimization.value:
        s.temp_ramp = derived(
            bool(ramp and (nmropt.value or 0) > 0), "&wt TEMP0 value1 != value2", (nmropt,)
        )
    return s
