#!/usr/bin/env python
from __future__ import division, print_function

from argparse import ArgumentParser

import simtk.openmm as mm
import simtk.unit as u
import simtk.openmm.app as app

parser = ArgumentParser()
group = parser.add_argument_group('Input/Output file options')
group.add_argument('-p', '--pdb', metavar='<PDB FILE>', required=True,
                   dest='pdb', help='PDB file with target system')
group.add_argument('-f', '--forcefield', dest='ff', default='amoeba2013',
                   metavar='<ForceField>', help='''Target force field to use.
                   Should be one of the XML force field files bundled with
                   OpenMM. Default is %(default)s''')
group.add_argument('-s', '--system-xml', dest='system', metavar='FILE',
                   default='system.xml', help='''Name of the output OpenMM
                   System XML file that will be written by this script. Default
                   is %(default)s''')
group = parser.add_argument_group('Potential energy function parameters')
group.add_argument('-c', '--cutoff', metavar='FLOAT', default=7, type=float,
                   help='''Cutoff in Angstroms to use for nonbonded
                   interactions. Default is %(default)s Angstroms (tailored
                   toward AMOEBA simulations). Where vdW and electrostatic
                   cutoffs can differ (only AMOEBA), this is *only* the
                   electrostatic cutoff. -v/--vdw-cutoff for the van der Waals
                   cutoff option.''', dest='cut')
group.add_argument('-v', '--vdw-cutoff', metavar='FLOAT', default=9,
                   type=float, help='''Cutoff in Angstroms to use for van der
                   Waals interactions. This is only applicable to AMOEBA force
                   fields and will be ignored otherwise. Default is
                   %(default)s''', dest='vdwcut')
group.add_argument('-e', '--epsilon', metavar='FLOAT', default=1e-5, dest='eps',
                   help='''Convergence criteria for mutually induced polarizable
                   dipoles in the AMOEBA force field. Default is %(default)g.
                   This option is ignored for fixed-charge FFs''', type=float)
group = parser.add_argument_group('Integration-related options')
group.add_argument('-m', '--hydrogen-mass', metavar='FLOAT', default=None,
                   type=float, help='''Mass of hydrogen atoms (in daltons) to
                   use in the simulation. The default is to use unmodified
                   H-atom masses.  This implements hydrogen mass repartitioning
                   (so the total mass stays the same, but the hydrogens become
                   heavier), which allows a longer timestep to be used.
                   Suggested value for H-mass repartitioning is 3 daltons''',
                   dest='hmass')
group.add_argument('--shake', dest='shake', action='store_true', default=False,
                   help='''Constrain bonds containing hydrogen (using SETTLE on
                   waters). Default is NOT to apply SHAKE''')
# TODO -- add options for implicit solvent... maybe?

opt = parser.parse_args()

# Parse the PDB file
print('Parsing the PDB file [%s]...' % opt.pdb)
pdb = app.PDBFile(opt.pdb)

# Declare the FF we want to use... in this case amoeba13.xml
# Strip off an xml if users added it
ff = opt.ff[:-4] if opt.ff.endswith('.xml') else opt.ff
print('Loading the force field [%s.xml]...' % ff)
ff = app.ForceField(ff + '.xml')

# Create the System object with it
print('Creating the System...')
if opt.shake:
    print('Constraining bonds with hydrogens')
    constraints = app.HBonds
else:
    print('No constraints applied')
    constraints = None
if opt.hmass:
    print('Repartitioning hydrogen masses to %g daltons' % opt.hmass)
    hmass = opt.hmass * u.dalton
else:
    hmass = None
system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=opt.cut*u.angstroms,
                         rigidWater=opt.shake, constraints=constraints,
                         hydrogenMass=hmass)

# Now scan through our forces and change the cutoff for the van der Waals force
# to 9 angstroms, and change our dipole convergence to 1e-6
for force in system.getForces():
    if isinstance(force, mm.AmoebaVdwForce):
        print('Adjusting the vdW cutoff to %g Angstroms...' % opt.vdwcut)
        force.setCutoff(opt.vdwcut*u.angstroms)
    elif isinstance(force, mm.AmoebaMultipoleForce):
        print('Setting the induced dipole convergence criteria to %g' % opt.eps)
        force.setMutualInducedTargetEpsilon(opt.eps)

# Now we are done creating our system. Let's serialize it by writing an XML file
print('Serializing the System...')
with open(opt.system, 'w') as f:
    f.write(mm.XmlSerializer.serialize(system))

# Finish progress report
print('Done.')
