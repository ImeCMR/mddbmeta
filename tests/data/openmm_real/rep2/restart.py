"""Continue replica 2 from its checkpoint."""
from openmm import *
from openmm.app import *
from openmm.unit import *

pdb = PDBFile('../prep/solvated.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                                 constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin, 25))
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('CPU'))
simulation.loadCheckpoint('prod.chk')
simulation.reporters.append(XTCReporter('prod2.xtc', 500))
simulation.reporters.append(StateDataReporter('prod2.log', 500, step=True, time=True, potentialEnergy=True,
                                              temperature=True, volume=True, speed=True))
simulation.step(2500)
simulation.saveState('prod2_final.xml')
