from __future__ import annotations

from pathlib import Path

import pytest

from mddbmeta.engines.openmm import OpenmmAdapter, classify_forcefield_files
from mddbmeta.engines.openmm.files import read_state_log, read_xml_state
from mddbmeta.engines.openmm.script import Template, read_script
from mddbmeta.model import FileKind

REAL = Path(__file__).resolve().parents[1] / "data" / "openmm_real"


def _script(tmp_path: Path, text: str):
    p = tmp_path / "run.py"
    p.write_text(
        "from openmm import *\nfrom openmm.app import *\nfrom openmm.unit import *\n" + text
    )
    return read_script(p)


def test_units_integrators_and_barostats(tmp_path: Path) -> None:
    info = _script(
        tmp_path,
        """
dt = 4*femtoseconds
T = 310
integrator = LangevinMiddleIntegrator(T*kelvin, 1/picosecond, dt)
system.addForce(MonteCarloBarostat(1.01325*bar, T*kelvin, 25))
simulation.reporters.append(StateDataReporter('md.csv', 1000, step=True, time=True))
simulation.reporters.append(DCDReporter('md.dcd', 5000))
nsteps = 250000
simulation.step(nsteps)
""",
    )
    (seg,) = info.segments
    assert seg.integrator["dt"].to("time") == pytest.approx(0.004)
    assert seg.integrator["temperature"].to("temperature") == 310.0
    assert seg.barostat["pressure"].to("pressure") == pytest.approx(1.01325)
    assert seg.steps == 250000 and seg.report_interval == 1000
    assert [(f, t.text, n) for f, t, n in seg.trajectories] == [("dcd", "md.dcd", 5000)]


def test_loops_fstrings_and_loads(tmp_path: Path) -> None:
    info = _script(
        tmp_path,
        """
for i in range(1, 4):
    if i > 1:
        simulation.loadCheckpoint(f'seg{i-1}.chk')
    simulation.reporters = [StateDataReporter(f'seg{i}.log', 100), XTCReporter(f'seg{i}.xtc', 100)]
    simulation.step(1000)
    simulation.saveCheckpoint(f'seg{i}.chk')
""",
    )
    (seg,) = info.segments
    assert seg.in_loop and seg.steps == 1000
    values = seg.log.match("seg3.log")
    assert values == {"i": "3"}
    assert seg.loads[0][0].fill(values) == "seg2.chk" and seg.loads[0][1] is True
    assert seg.trajectories[0][1].fill(values) == "seg3.xtc"


def test_steps_outside_loop_multiply(tmp_path: Path) -> None:
    info = _script(
        tmp_path,
        """
simulation.reporters.append(StateDataReporter('heat.csv', 100))
for t in range(0, 301, 30):
    integrator.setTemperature(t*kelvin)
    simulation.step(500)
""",
    )
    (seg,) = info.segments
    assert seg.steps == 11 * 500


def test_heating_ramp_and_new_velocities(tmp_path: Path) -> None:
    info = _script(
        tmp_path,
        """
integrator = LangevinIntegrator(50*kelvin, 1/picosecond, 0.002*picoseconds)
simulation.context.setVelocitiesToTemperature(50*kelvin)
simulation.reporters.append(StateDataReporter('heat.csv', 100))
integrator.setTemperature(100*kelvin)
simulation.step(500)
integrator.setTemperature(300*kelvin)
simulation.step(500)
""",
    )
    (seg,) = info.segments
    assert seg.temperature_changes == [100.0, 300.0] and seg.new_velocities and seg.steps == 1000


def test_with_open_deserialize_and_modeller(tmp_path: Path) -> None:
    info = _script(
        tmp_path,
        """
with open('system.xml') as f:
    system = XmlSerializer.deserialize(f.read())
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=1*nanometer)
with open('solvated.pdb', 'w') as out:
    PDBFile.writeFile(modeller.topology, modeller.positions, out)
""",
    )
    assert [t.text for t in info.deserialized] == ["system.xml"]
    assert info.modeller and [t.text for t in info.written_structures] == ["solvated.pdb"]


def test_script_is_never_executed(tmp_path: Path) -> None:
    marker = tmp_path / "executed"
    _script(tmp_path, f"open({str(marker)!r}, 'w').write('x')\nsimulation.step(10)\n")
    assert not marker.exists()


def test_template_arithmetic() -> None:
    t = Template("prod_{0}.chk", ["i - 1"])
    assert t.fill({"i": "5"}) == "prod_4.chk" and t.fill({}) is None
    assert Template("r{0}/md_{1}.log", ["rep", "n"]).match("r2/md_10.log") == {
        "rep": "2",
        "n": "10",
    }


def test_state_log_and_xml() -> None:
    log = read_state_log(REAL / "rep2" / "prod.log")
    assert log.first["Step"] == 5500 and log.last["Step"] == 10000 and log.report_interval == 500
    assert log.n_rows == 10 and log.volume_min < log.volume_max
    state = read_xml_state(REAL / "prep" / "eq_state.xml")
    assert (state.kind, state.version, state.step) == ("State", "8.6.1", 5000)
    assert (
        state.time_ps == pytest.approx(10.0)
        and state.natom == 884
        and state.box_vectors_A[0][0] > 20
    )
    system = read_xml_state(REAL / "prep" / "system.xml")
    assert any(f["type"] == "MonteCarloBarostat" for f in system.forces)
    integ = read_xml_state(REAL / "prep" / "integrator.xml")
    assert (
        integ.integrator["type"] == "LangevinMiddleIntegrator"
        and integ.integrator["stepSize"] == ".002"
    )


def test_stdout_log_with_other_output(tmp_path: Path) -> None:
    p = tmp_path / "run.log"
    p.write_text(
        'Minimizing...\n#"Step","Time (ps)","Temperature (K)"\n100,0.2,290\nsome warning\n200,0.4,300\n'
    )
    log = read_state_log(p)
    assert log.n_rows == 2 and log.last["Step"] == 200 and log.mean_temperature == 295


@pytest.mark.parametrize(
    "files,ffs,waters",
    [
        (["amber14-all.xml", "amber14/tip3pfb.xml"], ["Amber ff14SB"], ["TIP3P-FB"]),
        (["amber19-all.xml", "amber19/opc.xml"], ["Amber ff19SB"], ["OPC"]),
        (["charmm36.xml", "charmm36/water.xml"], ["CHARMM36"], ["TIP3P"]),
        (["amber99sbildn.xml", "tip3p.xml"], ["Amber ff99SB-ILDN"], ["TIP3P"]),
        (["amoeba2018.xml"], ["AMOEBA 2018"], []),
        (["amber14-all.xml", "implicit/obc2.xml"], ["Amber ff14SB"], []),
    ],
)
def test_forcefield_file_labels(files, ffs, waters) -> None:
    assert classify_forcefield_files(files) == (ffs, waters)


def test_sniff_real_files() -> None:
    ad = OpenmmAdapter()
    expected = {
        "prep/setup.py": FileKind.CONTROL, "prep/npt.csv": FileKind.LOG, "rep2/prod.log": FileKind.LOG,
        "prep/eq_state.xml": FileKind.COORDS, "prep/solvated.pdb": FileKind.COORDS,
    }  # fmt: skip
    for rel, kind in expected.items():
        path = REAL / rel
        res = ad.sniff(path, path.read_bytes()[:4096])
        assert res is not None and res.kind is kind, rel
