#!/usr/bin/env python

import parmed as pmd
from parmed import unit as u

from simtk import openmm as mm
from simtk.openmm import app

from argparse import ArgumentParser
import sys

parser = ArgumentParser()
group = parser.add_argument_group('Input File Options')
group.add_argument('--xml', dest='xml', default='system.xml', metavar='FILE',
                   help='''OpenMM System XML file. Default is %(default)s''')
group.add_argument('-p', '--pdb', dest='pdb', metavar='<PDB_FILE>', required=True,
                   help='''PDB file with coordinates for all atoms. Is also the
                   reference coordinates''')
group.add_argument('-s', '--state', dest='state', metavar='FILE', default=None,
                   help='''Restart file (any format)''')
group = parser.add_argument_group('Positional Restraint Options')
group.add_argument('--restrain', dest='restraints', metavar='MASK',
                   help='restraint mask (default None)', default=None)
group.add_argument('-k', '--force-constant', dest='force_constant', type=float,
                   metavar='FLOAT', help='''Force constant for cartesian
                   constraints. Default 10 kcal/mol/A^2''', default=10)
group = parser.add_argument_group('Output File Options')
group.add_argument('-r', '--restart', dest='restart', default='restart.nc',
                   metavar='FILE', help='''NetCDF file with information to
                   restart the simulation with another run''')
group.add_argument('-o' , '--output', dest='output', default=sys.stdout,
                   metavar='FILE', help='''Output file for energies''')
group.add_argument('-x', '--trajectory', dest='trajectory', default='md.nc',
                   metavar='FILE', help='''NetCDF trajectory to generate.
                   Snapshots written every --interval steps.''')
group.add_argument('--checkpoint', dest='checkpoint', metavar='FILE',
                   default=None, help='''Name of a checkpoint file to write
                   periodically throughout the simulation. Primarily useful for
                   debugging intermittent and rare errors.''')
group.add_argument('--interval', dest='interval', default=500, metavar='INT',
                   help='Interval between printing state data. Default 500',
                   type=int)
group = parser.add_argument_group('Simulation Options')
group.add_argument('-n', '--num-steps', dest='num_steps', required=True, type=int,
                   help='Number of MD steps to run. Required', metavar='INT')
group.add_argument('--ntp', dest='ntp', default=False, action='store_true',
                   help='Do constant pressure simulation')
group.add_argument('--aniso', dest='aniso', default=False, action='store_true',
                   help='Do anisotropic pressure scaling')
group.add_argument('--dt', dest='timestep', type=float,
                   metavar='FLOAT', help='''time step for Langevin
                   integrator. Default 1 fs''', default=1.0)
group.add_argument('--gamma_ln', dest='gamma_ln', type=float,
                   metavar='FLOAT', help='''collision frequency for Langevin
                   integrator. Default 1 ps-1''', default=1.0)
group.add_argument('--temp', dest='temp', type=float,
                   metavar='FLOAT', help='''target temperature for NVT
                   simulation. Default %(default)s K''', default=300.0)
group.add_argument('--nve', dest='nve', default=False, action='store_true',
                    help='Do constant energy simulation')

opt = parser.parse_args()

print('Command line:\n\t%s' % ' '.join(sys.argv)); sys.stdout.flush()

print('Parsing XML file %s' % opt.xml); sys.stdout.flush()
with open(opt.xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())

# Remove the andersen thermostat that might exist (e.g. from Tinker)
for i in range(system.getNumForces()):
    if isinstance(system.getForce(i), mm.AndersenThermostat):
        print('Deleting the Andersen thermostat'); sys.stdout.flush()
        system.removeForce(i)
        break
else:
    print('No Andersen thermostat found in the system.xml file'); sys.stdout.flush()

print('Parsing PDB file %s' % opt.pdb); sys.stdout.flush()
pdb = pmd.load_file(opt.pdb)

# Add cartesian restraints if desired
if opt.restraints:
    print('Adding restraints (k=%s kcal/mol/A^2) from %s' %
            (opt.force_constant, opt.restraints)); sys.stdout.flush()
    sel = pmd.amber.AmberMask(pdb, opt.restraints).Selection()
    const = opt.force_constant * u.kilocalories_per_mole/u.angstroms**2
    const = const.value_in_unit_system(u.md_unit_system)
    force = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    force.addGlobalParameter('k', const)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    for i, atom_crd in enumerate(pdb.positions):
        if sel[i]:
            force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
    system.addForce(force)

if opt.ntp:
    if opt.aniso:
        print('Using anisotropic barostat'); sys.stdout.flush()
        baro = mm.MonteCarloAnisotropicBarostat(1*u.bar, opt.temp*u.kelvin)
    else:
        print('Using isotropic barostat'); sys.stdout.flush()
        baro = mm.MonteCarloBarostat(1*u.bar, opt.temp*u.kelvin)
    system.addForce(baro)

if opt.gamma_ln == 0.0 and not opt.nve:
    print('Adding Anderson thermostat at %sK, 0.1 psec^-1' % opt.temp); sys.stdout.flush()
    thermo = mm.AndersenThermostat(opt.temp*u.kelvin, 0.1/u.picosecond)
    system.addForce(thermo)

# Create the simulation
if opt.gamma_ln > 0.0:
    integrator = mm.LangevinIntegrator(opt.temp*u.kelvin,
                 opt.gamma_ln/u.picosecond, opt.timestep*u.femtoseconds)
    print('Langevin: %8.2fK, %8.2f ps-1, %8.2f fs' %
       (opt.temp, opt.gamma_ln, opt.timestep) ); sys.stdout.flush()
else:
    integrator = mm.VerletIntegrator(opt.timestep*u.femtoseconds)
    print('Verlet: %8.2f fs' % opt.timestep )

sim = app.Simulation(pdb.topology, system, integrator,
                     mm.Platform.getPlatformByName('CUDA'),
                     dict(CudaPrecision='mixed') )
sim.reporters.append(
        pmd.openmm.StateDataReporter(opt.output, reportInterval=opt.interval,
                        volume=True,density=True,separator='\t')
)
sim.reporters.append(
        pmd.openmm.ProgressReporter(opt.output + '.info', opt.interval, opt.num_steps)
)
sim.reporters.append(
        pmd.openmm.NetCDFReporter(opt.trajectory, opt.interval*10)
)
sim.reporters.append(
        pmd.openmm.RestartReporter(opt.restart, opt.interval*100, netcdf=True)
)
if opt.checkpoint is not None:
    sim.reporters.append(
            app.CheckpointReporter(opt.checkpoint, opt.interval*100)
    )

if opt.state is not None:
    print('Setting coordinates and velocities from restart file %s' %
        opt.state); sys.stdout.flush()

    if opt.state[-3:] == 'xml':
        with open(opt.state, 'r') as f:
            sim.context.setState(mm.XmlSerializer.deserialize(f.read()))
    elif opt.state[-3:] == 'chk':
        sim.loadCheckpoint(opt.state)
    else:
#       jason's code that is supposed to work for any restart file type:
        rst = pmd.load_file(opt.state)
        sim.context.setPositions(rst.coordinates[-1]*u.angstroms)
        sim.context.setVelocities(rst.velocities[-1]*u.angstroms/u.picoseconds)
        sim.context.setPeriodicBoxVectors(*pmd.geometry.box_lengths_and_angles_to_vectors(*rst.box))
        if hasattr(rst, 'time'):
            try:
                sim.context.setTime(rst.time[-1])
            except TypeError:
                sim.context.setTime(rst.time)

else:
    print('Setting coordinates from PDB file %s' % opt.pdb); sys.stdout.flush()
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(opt.temp)

print('Running the simulation for %d steps!' % opt.num_steps); sys.stdout.flush()
sim.step(opt.num_steps)

# The last step may not have resulted in a restart file being written. Force it
# here
state = sim.context.getState(getPositions=True, getVelocities=True,
        getEnergy=True, getForces=True,
        enforcePeriodicBox=system.usesPeriodicBoundaryConditions())
for rep in sim.reporters:
    if isinstance(rep, pmd.openmm.RestartReporter):
        rep.report(sim, state)
