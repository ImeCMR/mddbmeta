"""Production run 1 for replica 2 (fresh velocities from the equilibrated state)."""
from openmm import *
from openmm.app import *
from openmm.unit import *

pdb = PDBFile('../prep/solvated.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                                 constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin, 25))
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
integrator.setRandomNumberSeed(2024)
simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('CPU'))
simulation.loadState('../prep/eq_state.xml')
simulation.context.setVelocitiesToTemperature(300*kelvin, 2024)
simulation.reporters.append(XTCReporter('prod.xtc', 500))
simulation.reporters.append(StateDataReporter('prod.log', 500, step=True, time=True, potentialEnergy=True,
                                              temperature=True, volume=True, speed=True))
simulation.reporters.append(CheckpointReporter('prod.chk', 1000))
simulation.step(5000)
simulation.saveState('prod_final.xml')
