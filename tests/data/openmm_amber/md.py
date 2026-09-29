import sys
import openmm as mm
import openmm.app as app
import openmm.unit as u

prmtop = app.AmberPrmtopFile('ala_ala.prmtop')
inpcrd = app.AmberInpcrdFile('ala_ala.inpcrd')
system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
system.addForce(mm.AndersenThermostat(310*u.kelvin, 1/u.picosecond))
integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
simulation = app.Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('Reference'))
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.reporters.append(app.DCDReporter('md.dcd', 100))
simulation.reporters.append(app.StateDataReporter(sys.stdout, 100, step=True, time=True, temperature=True))
nsteps = 2000
simulation.step(nsteps)
