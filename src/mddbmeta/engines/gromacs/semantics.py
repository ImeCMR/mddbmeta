"""Turn GROMACS run parameters into engine-neutral StepSettings.

Resolution order (GROMACS differs from AMBER here on purpose):

1. the "Input Parameters" mdrun printed in its log (``engine_report``, EXACT):
   these are the values that ran, after grompp and any ``convert-tpr``;
2. the matching .mdp (``stated``): the only source for grompp-level options
   such as ``define`` (position restraints) and ``gen-vel``;
3. the documented GROMACS default (``engine_default``, DERIVED).
"""

from __future__ import annotations

from typing import Any

from mddbmeta.model import StepSettings
from mddbmeta.provenance import Confidence, Mined, Source, derived, engine_default

MINIMIZERS = {"steep", "cg", "l-bfgs", "nm"}
STOCHASTIC = {"sd", "bd"}

DEFAULTS: dict[str, Any] = {
    "integrator": "md",
    "dt": 0.001,
    "nsteps": 0,
    "tinit": 0.0,
    "init-step": 0,
    "tcoupl": "no",
    "pcoupl": "no",
    "nstxout": 0,
    "nstxout-compressed": 0,
    "continuation": "no",
    "gen-vel": "no",
    "gen-temp": 300.0,
    "annealing": "no",
    "pbc": "xyz",
    "define": "",
}
TRUE = {"yes", "true", "on"}


def truthy(value: Any) -> bool:
    if isinstance(value, list):
        return any(truthy(v) for v in value)
    return str(value).strip().lower() in TRUE


class GmxResolver:
    def __init__(
        self,
        reported: dict[str, Any] | None,
        reported_path: str | None,
        stated: dict[str, Any] | None,
        stated_path: str | None,
    ):
        self.reported = reported or {}
        self.reported_path = reported_path
        self.stated = stated or {}
        self.stated_path = stated_path

    def get(self, key: str, unit: str | None = None) -> Mined[Any]:
        if key in self.reported and self.reported_path:
            return Mined(
                value=self.reported[key],
                unit=unit,
                sources=(Source(self.reported_path, f"Input Parameters:{key}"),),
                method="engine_report",
                confidence=Confidence.EXACT,
            )
        if key in self.stated and self.stated_path:
            return Mined(
                value=self.stated[key], unit=unit, sources=(Source(self.stated_path, key),)
            )
        if key in DEFAULTS:
            src = self.reported_path or self.stated_path or ""
            return engine_default(DEFAULTS[key], src, key, unit=unit)
        return Mined.missing()


def _single(values: Any) -> Any:
    """ref-t is per temperature-coupling group; one value when they all agree."""
    if isinstance(values, list):
        distinct = {v for v in values}
        return values[0] if len(distinct) == 1 else None
    return values


def settings_from_params(res: GmxResolver) -> StepSettings:
    s = StepSettings()
    integ = res.get("integrator")
    name = str(integ.value).lower()
    minimization = name in MINIMIZERS
    s.minimization = integ.with_value(minimization, method=f"map:integrator ({integ.method})")
    s.nsteps = res.get("nsteps")
    if minimization:
        s.ensemble = Mined(value=None, method="not_applicable", note="energy minimization")
        return s
    s.dt_ps = res.get("dt", unit="ps")

    tcoupl = res.get("tcoupl")
    if name in STOCHASTIC:
        s.thermostat = integ.with_value(
            f"{name} integrator", method=f"map:integrator ({integ.method})"
        )
        thermo_on = True
    else:
        t = str(tcoupl.value).lower()
        thermo_on = t not in ("no", "false", "")
        s.thermostat = tcoupl.with_value(
            t if thermo_on else "none", method=f"map:tcoupl ({tcoupl.method})"
        )
    if thermo_on:
        ref = res.get("ref-t", unit="K")
        if ref.is_known:
            one = _single(ref.value)
            s.temp_K = ref.with_value(
                float(one) if one is not None else None,
                note=None
                if one is not None
                else f"groups coupled to different temperatures: {ref.value}",
            )

    pcoupl = res.get("pcoupl")
    p = str(pcoupl.value).lower()
    cp = p not in ("no", "false", "")
    s.constant_pressure = pcoupl.with_value(cp, method=f"map:pcoupl ({pcoupl.method})")
    s.barostat = pcoupl.with_value(p if cp else "none", method=f"map:pcoupl ({pcoupl.method})")

    pbc = res.get("pbc")
    s.periodic = pbc.with_value(str(pbc.value).lower() != "no", method=f"map:pbc ({pbc.method})")
    cont = res.get("continuation")
    s.restart = cont.with_value(truthy(cont.value), method=f"map:continuation ({cont.method})")
    gen_vel = res.get("gen-vel")
    if truthy(gen_vel.value):
        s.temp_initial_K = res.get("gen-temp", unit="K")

    define = res.get("define")
    if define.method == "stated":
        s.restraints = define.with_value(
            "POSRES" in str(define.value).upper(), method="map:define (-DPOSRES)"
        )
    ann = res.get("annealing")
    s.temp_ramp = ann.with_value(
        any(
            str(v).lower() in ("single", "periodic")
            for v in (ann.value if isinstance(ann.value, list) else [ann.value])
        ),
        method=f"map:annealing ({ann.method})",
    )

    xtc = res.get("nstxout-compressed")
    trr = res.get("nstxout")
    s.traj_interval_steps = xtc if xtc.value else trr
    rvdw = res.get("rvdw")
    if rvdw.is_known and isinstance(rvdw.value, (int, float)):
        s.cutoff_A = rvdw.with_value(
            float(rvdw.value) * 10.0, unit="A", method=f"{rvdw.method}; nm->A"
        )

    if thermo_on:
        ens = "NPT" if cp else "NVT"
    else:
        ens = "NPH" if cp else "NVE"
    s.ensemble = derived(ens, "ensemble from integrator/tcoupl/pcoupl", (integ, tcoupl, pcoupl))
    return s
