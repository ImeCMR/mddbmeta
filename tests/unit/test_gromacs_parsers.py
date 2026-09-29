from __future__ import annotations

import struct
from pathlib import Path

import pytest

from mddbmeta.engines.gromacs import GromacsAdapter, _same_printed, clean_version
from mddbmeta.engines.gromacs.binary import read_tpr_header, read_trr, read_xtc
from mddbmeta.engines.gromacs.log import read_log
from mddbmeta.engines.gromacs.mdp import parse_mdp_text
from mddbmeta.engines.gromacs.semantics import GmxResolver, settings_from_params
from mddbmeta.engines.gromacs.structure import (
    evidence_from_includes,
    force_field_label,
    lengths_angles,
    read_gro,
    read_top,
)
from mddbmeta.export.mddb import duplicated_join_frames
from mddbmeta.mining.boxtype import classify_box
from mddbmeta.model import FileKind
from mddbmeta.provenance import Confidence

REAL = Path(__file__).resolve().parents[1] / "data" / "gromacs_real"
MDDB = Path(__file__).resolve().parents[3] / "MDDB-workflow" / "test" / "data" / "input"

# ------------------------------------------------------------------ mdp


def test_mdp_normalization_and_values() -> None:
    m = parse_mdp_text(
        "; comment\ninteGrator = md\nref_t = 300 310 ; per group\nnstxtcout = 500\n"
        "define = -DPOSRES -DFLEXIBLE\ntc_grps = Protein SOL\ndt=0.002\n"
    )
    assert m["integrator"] == "md" and m["ref-t"] == [300, 310]
    assert m["nstxout-compressed"] == 500  # pre-2016 alias
    assert m["define"] == "-DPOSRES -DFLEXIBLE" and m["tc-grps"] == ["Protein", "SOL"]
    assert m["dt"] == 0.002


# ------------------------------------------------------------------ semantics


def _s(reported=None, stated=None):
    return settings_from_params(
        GmxResolver(reported, "md.log" if reported else None, stated, "md.mdp" if stated else None)
    )


def test_semantics_log_beats_mdp_and_defaults() -> None:
    s = _s(
        reported={"integrator": "md", "dt": 0.004, "nsteps": 10},
        stated={"dt": 0.002, "define": "-DPOSRES"},
    )
    assert s.dt_ps.value == 0.004 and s.dt_ps.method == "engine_report"
    assert s.restraints.value is True  # define only exists in the mdp
    assert s.ensemble.value == "NVE" and s.thermostat.value == "none"
    assert _s(stated={"nsteps": 5}).dt_ps.method == "engine_default"


@pytest.mark.parametrize(
    "params,ensemble,thermostat",
    [
        (
            {"tcoupl": "V-rescale", "pcoupl": "Parrinello-Rahman", "ref-t": [300, 300]},
            "NPT",
            "v-rescale",
        ),
        ({"tcoupl": "nose-hoover", "ref-t": 310}, "NVT", "nose-hoover"),
        ({"integrator": "sd", "ref-t": [298.15]}, "NVT", "sd integrator"),
        ({"pcoupl": "C-rescale"}, "NPH", "none"),
    ],
)
def test_semantics_ensembles(params, ensemble, thermostat) -> None:
    s = _s(reported={"integrator": "md", **params})
    assert s.ensemble.value == ensemble and s.thermostat.value == thermostat


def test_semantics_temperature_groups_and_minimizer() -> None:
    assert _s(reported={"tcoupl": "v-rescale", "ref-t": [300, 300]}).temp_K.value == 300.0
    split = _s(reported={"tcoupl": "v-rescale", "ref-t": [300, 310]}).temp_K
    assert split.value is None and "different temperatures" in split.note
    em = _s(reported={"integrator": "steep", "nsteps": 500})
    assert em.minimization.value is True and em.ensemble.value is None


def test_semantics_gen_vel_and_annealing() -> None:
    s = _s(reported={"tcoupl": "berendsen", "ref-t": 300, "annealing": ["Single", "No"]},
           stated={"gen-vel": "yes", "gen-temp": 50})  # fmt: skip
    assert s.temp_initial_K.value == 50 and s.temp_ramp.value is True


# ------------------------------------------------------------------ structure files


def test_gro_dodecahedron_and_residues() -> None:
    g = read_gro(REAL / "prep" / "npt.gro")
    assert g.natom == 2139 and g.residue_counts == {"ALA": 3, "SOL": 702}
    assert g.residue_atom_names["SOL"] == ["OW", "HW1", "HW2"]
    lengths, angles = g.box_lengths_angles()
    assert angles == pytest.approx((60.0, 60.0, 90.0), abs=0.01)
    assert classify_box(lengths, angles, 1) == "Dodecahedron"


def test_gro_time_in_title_and_truncation(tmp_path: Path) -> None:
    p = tmp_path / "t.gro"
    p.write_text(
        "frame t=  12.50000 step= 6250\n    2\n    1SOL     OW    1   0.100   0.200   0.300\n"
    )
    with pytest.raises(ValueError, match="1 of 2"):
        read_gro(p)
    p.write_text(
        "frame t=  12.50000\n    1\n    1SOL     OW    1   0.100   0.200   0.300\n   1.0 1.0 1.0\n"
    )
    g = read_gro(p)
    assert g.time_ps == 12.5 and g.box_lengths_angles()[1] == (90.0, 90.0, 90.0)


def test_truncated_octahedron_from_gromacs_vectors() -> None:
    d = 1.0
    vec = ((d, 0, 0), (d / 3, 2 * 2**0.5 * d / 3, 0), (-d / 3, 2**0.5 * d / 3, 6**0.5 * d / 3))
    lengths, angles = lengths_angles(vec)
    assert classify_box(lengths, angles, 1) == "Truncated Octahedron"


def test_top_includes_and_force_field_labels() -> None:
    top = read_top(REAL / "prep" / "topol.top")
    assert ("SOL", 702) in top.molecules and top.system == "Protein in water"
    assert evidence_from_includes(top.includes) == (["Amber ff99SB-ILDN"], ["TIP3P"])
    assert force_field_label("charmm36-jul2022.ff") == "CHARMM36 (jul2022)"
    assert force_field_label("oplsaa.ff") == "OPLS-AA/L"
    assert force_field_label("my_custom.ff") == "my_custom"


# ------------------------------------------------------------------ binary files


def test_xtc_frames_and_truncation(tmp_path: Path) -> None:
    x = read_xtc(REAL / "rep2" / "prod.part0003.xtc")
    assert (x.natom, x.n_frames, x.first_step, x.last_step) == (2139, 6, 5000, 7500)
    assert (x.first_time_ps, x.last_time_ps, x.median_interval_ps) == (10.0, 15.0, 1.0)
    data = (REAL / "rep2" / "prod.part0003.xtc").read_bytes()
    cut = tmp_path / "cut.xtc"
    cut.write_bytes(data[: len(data) - 100])
    c = read_xtc(cut)
    assert c.truncated and c.n_frames == 5


def _trr_frame(natoms: int, step: int, t: float) -> bytes:
    box, x = 36, natoms * 12
    head = struct.pack(">iii", 1993, 13, 12) + b"GMX_trn_file"
    sizes = struct.pack(">13i", 0, 0, box, 0, 0, 0, 0, x, 0, 0, natoms, step, 0)
    body = struct.pack(">ff", t, 0.0) + struct.pack(">9f", 3, 0, 0, 0, 3, 0, 0, 0, 3) + b"\x00" * x
    return head + sizes + body


def test_trr_frames(tmp_path: Path) -> None:
    p = tmp_path / "t.trr"
    p.write_bytes(b"".join(_trr_frame(5, s, s * 0.002) for s in (0, 500, 1000)))
    r = read_trr(p)
    assert r.n_frames == 3 and r.natom == 5 and r.times == [0.0, 1.0, 2.0]
    assert r.box_vectors_A[0] == (30.0, 0.0, 0.0)


@pytest.mark.parametrize(
    "path,version,natom",
    [
        (REAL / "rep2" / "prod.tpr", "2022.3-dev-20220902-47c9856-dirty-unknown", 2139),
        (MDDB / "dummy" / "gromacs" / "ala_ala.tpr", "2025.3-conda_forge", 33),
        (MDDB / "raw_project" / "raw_topology.tpr", "2024.5-conda_forge", 8453),
        (MDDB / "raw_structures" / "A0224.tpr", "2022-conda_forge", 10395),
    ],
)
def test_tpr_header_across_versions(path: Path, version: str, natom: int) -> None:
    if not path.exists():
        pytest.skip("MDDB-workflow test data not available")
    h = read_tpr_header(path)
    assert h.version == version and h.natom == natom


# ------------------------------------------------------------------ logs


def test_log_sessions_checkpoints_and_params() -> None:
    log = read_log(REAL / "rep2" / "prod.log")
    first, second = log.sessions
    assert first.checkpoint_file is None and not first.finished  # truncated when appended to
    assert (
        second.appending
        and second.checkpoint_file == "prod.cpt"
        and second.checkpoint_time_ps == 5.0
    )
    assert second.option("-cpi") == "prod.cpt" and second.finished
    p = log.params
    assert p["nsteps"] == 2500  # appended sessions don't reprint Input Parameters
    assert p["ref-t"] == [300, 300] and p["continuation"] == "true" and p["pcoupl"] == "C-rescale"
    assert log.records()[0].step == 0 and log.records()[1].step == 5000


def test_log_noappend_part_and_minimization() -> None:
    part = read_log(REAL / "rep2" / "prod.part0003.log").sessions[0]
    assert (
        part.has_flag("-noappend") and part.simulation_part == 2 and part.checkpoint_time_ps == 10.0
    )
    em = read_log(REAL / "prep" / "em.log")
    assert em.minimization and em.first.converged and em.first.finished


def test_energy_comparison_tolerates_last_digit() -> None:
    assert _same_printed("4.57306e+03", "4.57305e+03")
    assert not _same_printed("4.57306e+03", "4.57304e+03")
    assert clean_version("2022.3-dev-20220902-47c9856") == "2022.3"
    assert clean_version("2019.6") == "2019.6"


# ------------------------------------------------------------------ merge simulation


def test_duplicated_join_frames_match_observed_mdconvert() -> None:
    # observed with MDTraj mdconvert (MDDB-workflow's merge): every frame is kept in list order,
    # so a -noappend continuation's re-written first frame appears twice (rep2: 11 + 6 -> 17)
    assert duplicated_join_frames([(0.0, 10.0), (10.0, 15.0)]) == [1]
    # clock resets lose nothing and duplicate nothing (rep1: 3 x 6 -> 18 frames)
    assert duplicated_join_frames([(0.0, 5.0), (0.0, 5.0), (0.0, 5.0)]) == []
    # AMBER-style parts start one interval after the previous end
    assert duplicated_join_frames([(1.0, 5.0), (6.0, 10.0)]) == []


# ------------------------------------------------------------------ sniffing


def test_sniff_real_files() -> None:
    ad = GromacsAdapter()
    expected = {
        "prep/npt.log": FileKind.LOG, "prep/npt.tpr": FileKind.TOPOLOGY, "prep/topol.top": FileKind.TOPOLOGY,
        "prep/npt.gro": FileKind.COORDS, "prep/npt.mdp": FileKind.CONTROL, "rep1/md_1.xtc": FileKind.TRAJECTORY,
    }  # fmt: skip
    for rel, kind in expected.items():
        path = REAL / rel
        res = ad.sniff(path, path.read_bytes()[:4096])
        assert res is not None and res.kind is kind, rel
    assert _s(reported={"dt": 0.002}).dt_ps.confidence is Confidence.EXACT
