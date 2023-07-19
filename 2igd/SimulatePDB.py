from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import parmed as pmd
u = pmd.unit

pdb = PDBFile('2igd_solvated.pdb')
forcefield = ForceField('../amoeba2018.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

e=simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(' Initial energy = %10.4f kcal/mol' % e.value_in_unit(u.kilocalories_per_mole))

# Now write a serialized state that has coordinates
print('Serializing the System...')
with open('system2.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(system))
print('Writing serialized XML restart file...')
with open('restart.xml', 'w') as f:
    f.write(
            openmm.XmlSerializer.serialize(
                simulation.context.getState(
                    getPositions=True, getVelocities=True, getForces=True,
                    enforcePeriodicBox=system.usesPeriodicBoundaryConditions(),
                    getEnergy=True
                )
            )
    )

print('Done!')
quit()

simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
