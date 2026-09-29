"""Production in three segments, each continuing from the previous checkpoint."""
import os
from openmm import *
from openmm.app import *
from openmm.unit import *

pdb = PDBFile('../prep/solvated.pdb')
with open('../prep/system.xml') as f:
    system = XmlSerializer.deserialize(f.read())
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('CPU'))
simulation.loadState('../prep/eq_state.xml')

nsegments = 3
steps_per_segment = 2500
for i in range(1, nsegments + 1):
    if i > 1:
        simulation.loadCheckpoint(f'prod_{i-1}.chk')
    simulation.reporters = [
        DCDReporter(f'prod_{i}.dcd', 500),
        StateDataReporter(f'prod_{i}.csv', 500, step=True, time=True, potentialEnergy=True,
                          temperature=True, volume=True, speed=True),
    ]
    simulation.step(steps_per_segment)
    simulation.saveCheckpoint(f'prod_{i}.chk')
    simulation.saveState(f'prod_{i}.xml')
