from __future__ import annotations

import os
from pathlib import Path

import pytest

from amberfiles import (
    PROD,
    WATER_OPC,
    WATER_TIP3P,
    RunSpec,
    mdin_text,
    mdout_text,
    protein_residues,
    write_mdcrd,
    write_nc_restart,
    write_nc_traj,
    write_prmtop,
    write_rst7,
)
from mddbmeta.engines.amber import AmberAdapter
from mddbmeta.engines.amber.leaplog import classify_leaprc, parse_leap_text
from mddbmeta.engines.amber.mdin import MdinFormatError, parse_mdin_text
from mddbmeta.engines.amber.mdout import read_mdout
from mddbmeta.engines.amber.netcdf import read_netcdf
from mddbmeta.engines.amber.prmtop import PrmtopFormatError, read_prmtop
from mddbmeta.engines.amber.rst7 import (
    ascii_traj_frames,
    read_ascii_traj,
    read_rst7,
    sniff_ascii_coords,
)
from mddbmeta.engines.amber.semantics import CntrlResolver, settings_from_cntrl
from mddbmeta.model import FileKind
from mddbmeta.provenance import Confidence

# ------------------------------------------------------------------ prmtop


def test_prmtop_pointers_residues_box(tmp_path: Path) -> None:
    residues = protein_residues(3) + [WATER_TIP3P] * 4
    p = write_prmtop(tmp_path / "s.prmtop", residues, ifbox=2, box=(109.4712206, 50.0, 50.0, 50.0))
    top = read_prmtop(p)
    assert top.natom == 3 * 6 + 4 * 3
    assert top.nres == 7
    assert top.ifbox == 2
    assert top.box == pytest.approx((109.4712206, 50.0, 50.0, 50.0))
    assert top.residue_counts() == {"ALA": 1, "GLY": 1, "SER": 1, "WAT": 4}
    atoms, names = top.residue_templates()
    assert atoms["WAT"] == 3 and names["WAT"] == ["O", "H1", "H2"]
    assert top.solvent_pointers == (3, 7, 4)
    assert "MASS" in top.flags_seen  # skipped but seen


def test_prmtop_extra_points_and_chamber(tmp_path: Path) -> None:
    p = write_prmtop(
        tmp_path / "c.prmtop", [WATER_OPC] * 2, numextra=2, chamber=True, title="charmm sys"
    )
    top = read_prmtop(p)
    assert top.numextra == 2
    assert top.is_chamber and top.title == "charmm sys"
    assert top.force_field_type and "CHARMM" in top.force_field_type


def test_prmtop_rejects_old_format(tmp_path: Path) -> None:
    p = tmp_path / "old.prmtop"
    p.write_text("old style title\n   10    2    0\n")
    with pytest.raises(PrmtopFormatError, match="old-format"):
        read_prmtop(p)


def test_prmtop_bad_format_line(tmp_path: Path) -> None:
    p = tmp_path / "bad.prmtop"
    p.write_text("%VERSION x\n%FLAG POINTERS\n%FORMAT(garbage)\n       1\n")
    with pytest.raises(PrmtopFormatError, match="format"):
        read_prmtop(p)


# ------------------------------------------------------------------ mdin


def test_mdin_namelist_parsing_edge_cases() -> None:
    text = """Production run ! with a comment
 &cntrl
   imin = 0, irest=1, ntx=5,   ! restart
   nstlim=500000 dt=4.0d-3,
   restraintmask=':1-10@CA,C,N', restraint_wt = 2.5E+00,
   ntt=3 gamma_ln=1.0, temp0=310.0,
 /
 &ewald
   dsum_tol=1.0e-6
 &end
 &wt type='TEMP0', istep1=0, istep2=1000, value1=0.0, value2=300.0 /
 &wt type='END' /
DISANG=dist.RST
"""
    m = parse_mdin_text(text)
    c = m.cntrl
    assert m.title == "Production run"
    assert c["dt"] == pytest.approx(0.004)
    assert c["nstlim"] == 500000 and c["temp0"] == 310.0
    assert c["restraintmask"] == ":1-10@CA,C,N"
    assert m.first("ewald") == {"dsum_tol": 1e-6}
    wt = m.all("wt")
    assert wt[0]["type"] == "TEMP0" and wt[0]["value2"] == 300.0
    assert wt[1] == {"type": "END"}
    assert "DISANG=dist.RST" in m.trailer


def test_mdin_multiline_and_case() -> None:
    m = parse_mdin_text("t\n &CNTRL\n  IMIN=1,\n  MAXCYC=\n  1000,\n /\n")
    assert m.cntrl == {"imin": 1, "maxcyc": 1000}


def test_mdin_without_namelist() -> None:
    with pytest.raises(MdinFormatError):
        parse_mdin_text("just a title\nnothing else\n")


# ------------------------------------------------------------------ semantics


def _settings(cntrl: dict, reported: dict | None = None):
    res = CntrlResolver(
        cntrl, "run.mdin", reported=reported, reported_path="run.mdout" if reported else None
    )
    return settings_from_cntrl(res, [])


def test_semantics_stated_vs_default() -> None:
    s = _settings({"imin": 0, "nstlim": 100, "ntt": 3, "ntb": 1})
    assert s.dt_ps.value == 0.001 and s.dt_ps.method == "engine_default"
    assert s.dt_ps.confidence is Confidence.DERIVED
    assert s.temp_K.value == 300.0 and s.temp_K.method == "engine_default"
    assert s.nsteps.value == 100 and s.nsteps.method == "stated"


def test_semantics_engine_report_beats_default() -> None:
    s = _settings({"imin": 0, "nstlim": 100, "ntt": 3}, reported={"dt": 0.002, "temp0": 310.0})
    assert s.dt_ps.value == 0.002 and s.dt_ps.method == "engine_report"
    assert s.dt_ps.confidence is Confidence.EXACT


@pytest.mark.parametrize(
    "cntrl,ensemble",
    [
        ({"ntt": 3, "ntb": 2, "ntp": 1}, "NPT"),
        ({"ntt": 3, "ntp": 1}, "NPT"),  # ntb defaults to 2 when ntp > 0
        ({"ntt": 3}, "NVT"),  # ntb defaults to 1
        ({"ntt": 0}, "NVE"),
        ({"ntp": 1}, "NPH"),
        ({"ntt": 3, "igb": 5}, "NVT"),  # implicit solvent, ntb defaults to 0
    ],
)
def test_semantics_ensemble_matrix(cntrl: dict, ensemble: str) -> None:
    s = _settings({"imin": 0, "nstlim": 10, **cntrl})
    assert s.ensemble.value == ensemble


def test_semantics_minimization() -> None:
    s = _settings({"imin": 1, "maxcyc": 500, "ntr": 1})
    assert s.minimization.value is True
    assert s.nsteps.value == 500
    assert s.ensemble.value is None
    assert not s.dt_ps.is_known and not s.temp_K.is_known
    assert s.restraints.value is True


def test_semantics_nve_has_no_temperature() -> None:
    s = _settings({"imin": 0, "nstlim": 10, "ntt": 0, "temp0": 310.0})
    assert not s.temp_K.is_known


def test_semantics_barostat_and_ramp() -> None:
    res = CntrlResolver(
        {"imin": 0, "nstlim": 10, "ntt": 3, "ntp": 1, "barostat": 2, "nmropt": 1}, "r.mdin"
    )
    s = settings_from_cntrl(
        res, [{"type": "TEMP0", "value1": 0.0, "value2": 300.0}, {"type": "END"}]
    )
    assert s.barostat.value == "monte-carlo"
    assert s.temp_ramp.value is True


# ------------------------------------------------------------------ mdout


def test_mdout_complete_run(tmp_path: Path) -> None:
    spec = RunSpec(cntrl=dict(PROD), start_time_ps=100.0)
    names = {
        "MDIN": "prod.mdin",
        "INPCRD": "equil.rst7",
        "PARM": "sys.prmtop",
        "RESTRT": "prod.rst7",
        "MDCRD": "prod.nc",
    }
    p = tmp_path / "prod.mdout"
    p.write_text(mdout_text(spec, names, natom=42))
    out = read_mdout(p)
    assert out.program == "PMEMD" and out.version == "22"
    assert out.variant == "pmemd.cuda"
    assert out.assignments["INPCRD"] == "equil.rst7" and out.assignments["MDCRD"] == "prod.nc"
    assert out.natom == 42
    assert out.begin_time_ps == 100.0
    assert out.control["dt"] == 0.002 and out.control["ntwx"] == 500
    assert out.completed and out.wall_time_s == 42.0
    # the averages/RMS blocks must not count as the last record
    assert out.last_record.nstep == 5000 and out.last_record.time_ps == pytest.approx(110.0)
    assert out.input_echo and "&cntrl" in out.input_echo


def test_mdout_truncated_run(tmp_path: Path) -> None:
    spec = RunSpec(cntrl=dict(PROD), completed=False)
    p = tmp_path / "prod.mdout"
    p.write_text(mdout_text(spec, {}, natom=42))
    out = read_mdout(p)
    assert not out.completed
    assert out.last_record.nstep == 2000


def test_mdout_minimization(tmp_path: Path) -> None:
    spec = RunSpec(cntrl={"imin": 1, "maxcyc": 500})
    p = tmp_path / "min.mdout"
    p.write_text(mdout_text(spec, {}, natom=42))
    out = read_mdout(p)
    assert out.minimization and out.last_min_step == 500


def test_mdout_tolerates_non_utf8(tmp_path: Path) -> None:
    spec = RunSpec(cntrl=dict(PROD))
    p = tmp_path / "x.mdout"
    p.write_bytes(mdout_text(spec, {}, natom=5).encode() + b"\xff\xfe garbage\n")
    assert read_mdout(p).completed


# ------------------------------------------------------------------ coordinates


def test_rst7_time_box_velocities(tmp_path: Path) -> None:
    p = write_rst7(tmp_path / "a.rst7", natom=7, time_ps=1234.5)
    r = read_rst7(p)
    assert r.natom == 7 and r.time_ps == pytest.approx(1234.5)
    assert r.has_velocities and r.box == pytest.approx((40.0, 40.0, 40.0, 90.0, 90.0, 90.0))


def test_inpcrd_without_time(tmp_path: Path) -> None:
    p = write_rst7(tmp_path / "a.inpcrd", natom=5, time_ps=None, velocities=False, box=None)
    r = read_rst7(p)
    assert r.time_ps is None and not r.has_velocities and r.box is None


def test_ascii_sniff_distinguishes_rst7_and_mdcrd(tmp_path: Path) -> None:
    a = write_rst7(tmp_path / "a.rst7", natom=5, time_ps=1.0)
    b = write_mdcrd(tmp_path / "b.mdcrd", natom=5, n_frames=3)
    assert sniff_ascii_coords(a.read_bytes()[:4096]) == "rst7"
    assert sniff_ascii_coords(b.read_bytes()[:4096]) == "mdcrd"


def test_ascii_traj_frame_count(tmp_path: Path) -> None:
    for box in (True, False):
        p = write_mdcrd(tmp_path / f"t{box}.mdcrd", natom=11, n_frames=4, box=box)
        frames, has_box = ascii_traj_frames(read_ascii_traj(p), 11)
        assert frames == 4 and has_box is box


@pytest.mark.netcdf
def test_netcdf_trajectory_and_restart(tmp_path: Path) -> None:
    t = write_nc_traj(tmp_path / "p.nc", natom=6, times=[1.0, 2.0, 3.0, 4.0])
    info = read_netcdf(t)
    assert not info.is_restart and info.n_frames == 4 and info.natom == 6
    assert info.first_time_ps == 1.0 and info.last_time_ps == 4.0 and info.median_dt_ps == 1.0
    assert info.program == "pmemd.cuda" and info.cell_lengths == (40.0, 40.0, 40.0)
    r = write_nc_restart(tmp_path / "p.ncrst", natom=6, time_ps=55.5)
    rinfo = read_netcdf(r)
    assert rinfo.is_restart and rinfo.first_time_ps == 55.5 and rinfo.natom == 6


# ------------------------------------------------------------------ leap


def test_leap_evidence() -> None:
    text = """source leaprc.protein.ff14SB
source leaprc.water.opc
source leaprc.gaff2
m = loadpdb sys.pdb
saveamberparm m system.prmtop system.inpcrd
"""
    ev = parse_leap_text(text, "tleap.in")
    assert ev.force_fields == ["Amber ff14SB", "GAFF2"]
    assert ev.water_models == ["OPC"]
    assert ev.saved_topologies == ["system.prmtop"]


def test_leap_log_source_lines() -> None:
    text = "Welcome to LEaP!\n----- Source: /opt/amber/dat/leap/cmd/leaprc.DNA.OL15\n----- Source: /x/leaprc.water.tip3p\n"
    ev = parse_leap_text(text, "leap.log")
    assert ev.force_fields == ["Amber OL15"] and ev.water_models == ["TIP3P"]


def test_classify_leaprc_unknown_family() -> None:
    assert classify_leaprc("leaprc.protein.fancyFF") == ("ff", "protein.fancyFF")
    assert classify_leaprc("frcmod.ions") is None


# ------------------------------------------------------------------ adapter sniff/parse


def test_adapter_sniff_by_content_not_extension(tmp_path: Path) -> None:
    ad = AmberAdapter()
    top = write_prmtop(tmp_path / "weird.txt", protein_residues(1))
    mdin = tmp_path / "prod.out"  # misleading extension
    mdin.write_text(mdin_text(PROD))
    for path, kind in ((top, FileKind.TOPOLOGY), (mdin, FileKind.CONTROL)):
        res = ad.sniff(path, path.read_bytes()[:4096])
        assert res is not None and res.kind is kind


def test_adapter_parse_never_raises(tmp_path: Path) -> None:
    ad = AmberAdapter()
    p = tmp_path / "broken.prmtop"
    p.write_text("%VERSION x\n%FLAG POINTERS\n%FORMAT(10I8)\n   abc\n")
    sniffed = ad.sniff(p, p.read_bytes())
    pf = ad.parse(p, "broken.prmtop", sniffed)
    assert not pf.ok and pf.errors[0].error_type == "malformed"


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores file permissions")
def test_adapter_permission_error(tmp_path: Path) -> None:
    ad = AmberAdapter()
    p = write_prmtop(tmp_path / "locked.prmtop", protein_residues(1))
    sniffed = ad.sniff(p, p.read_bytes())
    p.chmod(0)
    try:
        pf = ad.parse(p, "locked.prmtop", sniffed)
    finally:
        p.chmod(0o644)
    assert pf.errors[0].error_type == "permission"


def test_real_mddb_dummy_amber_files(mddb_dummy_amber: Path) -> None:
    ad = AmberAdapter()
    for name, kind in (("ala_ala.prmtop", FileKind.TOPOLOGY), ("ala_ala.inpcrd", FileKind.COORDS)):
        path = mddb_dummy_amber / name
        sniffed = ad.sniff(path, path.read_bytes()[:4096])
        assert sniffed is not None and sniffed.kind is kind, name
        pf = ad.parse(path, name, sniffed)
        assert pf.ok, pf.errors
    top = read_prmtop(mddb_dummy_amber / "ala_ala.prmtop")
    crd = read_rst7(mddb_dummy_amber / "ala_ala.inpcrd")
    assert top.natom == crd.natom and top.natom > 0
