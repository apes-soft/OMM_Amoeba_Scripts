#!/usr/bin/env python
from __future__ import print_function, division

from argparse import ArgumentParser
import parmed as pmd
u = pmd.unit

import openmm as mm
import openmm.app as app

parser = ArgumentParser()
parser.add_argument('-p', '--pdb', dest='pdb', metavar='FILE', required=True,
                    help='PDB file with the topology of the system')
parser.add_argument('-s', '--system-xml', dest='xml', metavar='FILE',
                    default='system.xml', help='''XML file with the OpenMM
                    System object serialized. Default is %(default)s''')
parser.add_argument('-x', '--xml', dest='output', metavar='FILE',
                    default='minimized.xml', help='''Output XML with minimized
                    coordinates. Default is %(default)s''')

opt = parser.parse_args()

print('Loading the system [%s] and PDB [%s] files...' % (opt.xml, opt.pdb))
system = pmd.load_file(opt.xml)
parm = pmd.load_file(opt.pdb)

print('Creating the Simulation...')
sim = app.Simulation(parm.topology, system, mm.VerletIntegrator(0.001),
                     platform=mm.Platform.getPlatformByName('CUDA'),
                     platformProperties=dict(CudaPrecision='mixed'))

print('Minimizing...')
sim.context.setPositions(parm.positions)
e=sim.context.getState(getEnergy=True).getPotentialEnergy()
print(' Initial energy = %10.4f kcal/mol' % e.value_in_unit(u.kilocalories_per_mole))
sim.minimizeEnergy()
e=sim.context.getState(getEnergy=True).getPotentialEnergy()
print('   Final energy = %10.4f kcal/mol' % e.value_in_unit(u.kilocalories_per_mole))

# Now write a serialized state that has coordinates
print('Finished. Writing serialized XML restart file...')
with open(opt.output, 'w') as f:
    f.write(
            mm.XmlSerializer.serialize(
                sim.context.getState(
                    getPositions=True, getVelocities=True, getForces=True,
                    enforcePeriodicBox=system.usesPeriodicBoundaryConditions(),
                    getEnergy=True
                )
            )
    )

print('Done!')
