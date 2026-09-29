"""Build, minimize and equilibrate (NVT then NPT) an Ala-Ala-Ala peptide in water."""
from openmm import *
from openmm.app import *
from openmm.unit import *

pdb = PDBFile('ala_ala.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield, padding=1.0*nanometers, model='tip3p', ionicStrength=0.15*molar)
with open('solvated.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
                                 constraints=HBonds, rigidWater=True)
temperature = 300*kelvin
integrator = LangevinMiddleIntegrator(temperature, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, Platform.getPlatformByName('CPU'))
simulation.context.setPositions(modeller.positions)

print('Minimizing...')
simulation.minimizeEnergy(maxIterations=500)
simulation.context.setVelocitiesToTemperature(temperature)

# NVT
simulation.reporters.append(StateDataReporter('nvt.csv', 250, step=True, time=True, potentialEnergy=True,
                                              temperature=True, speed=True))
simulation.reporters.append(DCDReporter('nvt.dcd', 250))
simulation.step(2500)
simulation.saveState('nvt_state.xml')

# NPT
system.addForce(MonteCarloBarostat(1*bar, temperature))
simulation.context.reinitialize(preserveState=True)
simulation.reporters = []
simulation.reporters.append(StateDataReporter('npt.csv', 250, step=True, time=True, potentialEnergy=True,
                                              temperature=True, volume=True, density=True, speed=True))
simulation.reporters.append(DCDReporter('npt.dcd', 250))
simulation.step(2500)
simulation.saveState('eq_state.xml')
with open('system.xml', 'w') as f:
    f.write(XmlSerializer.serialize(system))
with open('integrator.xml', 'w') as f:
    f.write(XmlSerializer.serialize(integrator))
