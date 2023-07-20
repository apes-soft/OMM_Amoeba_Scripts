#!/usr/bin/env python
from __future__ import print_function
import math

import numpy as np
import parmed as pmd
u = pmd.unit

import openmm as mm
import openmm.app as app
from openmm.unit import *
from openmm import *
from openmm.app import *

from argparse import ArgumentParser
import sys

def MTSVVVRIntegrator(temperature, collision_rate, timestep, system, ninnersteps=4):
    """
    Create a multiple timestep velocity verlet with velocity randomization (VVVR) integrator.
    Taken from https://github.com/leeping/forcebalance, by Lee-Ping Wang
    
    ARGUMENTS

    temperature (Quantity compatible with kelvin) - the temperature
    collision_rate (Quantity compatible with 1/picoseconds) - the collision rate
    timestep (Quantity compatible with femtoseconds) - the integration timestep
    system (simtk.openmm.System) - system whose forces will be partitioned
    ninnersteps (int) - number of inner timesteps (default: 4)

    RETURNS

    integrator (openmm.CustomIntegrator) - a VVVR integrator

    NOTES
    
    This integrator is equivalent to a Langevin integrator in the velocity Verlet discretization with a
    timestep correction to ensure that the field-free diffusion constant is timestep invariant.  The inner
    velocity Verlet discretization is transformed into a multiple timestep algorithm.

    REFERENCES

    VVVR Langevin integrator: 
    * http://arxiv.org/abs/1301.3800
    * http://arxiv.org/abs/1107.2967 (to appear in PRX 2013)    
    
    TODO

    Move initialization of 'sigma' to setting the per-particle variables.
    
    """
    # Multiple timestep Langevin integrator.
    for i in system.getForces():
        if i.__class__.__name__ in ["NonbondedForce", "CustomNonbondedForce", "AmoebaVdwForce", "AmoebaMultipoleForce", "MonteCarloBarostat"]:
            # Slow force.
            print('   %s is a Slow Force' % i.__class__.__name__);
            # logger.info(i.__class__.__name__ + " is a Slow Force\n")
            i.setForceGroup(1)
        else:
            # Fast force.
            print('   %s is a Fast Force' % i.__class__.__name__);
            # logger.info(i.__class__.__name__ + " is a Fast Force\n")
            i.setForceGroup(0)

    kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
    kT = kB * temperature
    
    integrator = CustomIntegrator(timestep)
    
    integrator.addGlobalVariable("dt_fast", timestep/float(ninnersteps)) # fast inner timestep
    integrator.addGlobalVariable("kT", kT) # thermal energy
    integrator.addGlobalVariable("a", np.exp(-collision_rate*timestep)) # velocity mixing parameter
    integrator.addGlobalVariable("b", np.sqrt((2/(collision_rate*timestep)) * np.tanh(collision_rate*timestep/2))) # timestep correction parameter
    integrator.addPerDofVariable("sigma", 0) 
    integrator.addPerDofVariable("x1", 0) # position before application of constraints

    #
    # Pre-computation.
    # This only needs to be done once, but it needs to be done for each degree of freedom.
    # Could move this to initialization?
    #
    integrator.addComputePerDof("sigma", "sqrt(kT/m)")

    # 
    # Velocity perturbation.
    #
    integrator.addComputePerDof("v", "sqrt(a)*v + sqrt(1-a)*sigma*gaussian")
    integrator.addConstrainVelocities();
    
    #
    # Symplectic inner multiple timestep.
    #
    integrator.addUpdateContextState(); 
    integrator.addComputePerDof("v", "v + 0.5*b*dt*f1/m")
    for innerstep in range(ninnersteps):
        # Fast inner symplectic timestep.
        integrator.addComputePerDof("v", "v + 0.5*b*dt_fast*f0/m")
        integrator.addComputePerDof("x", "x + v*b*dt_fast")
        integrator.addComputePerDof("x1", "x")
        integrator.addConstrainPositions();        
        integrator.addComputePerDof("v", "v + 0.5*b*dt_fast*f0/m + (x-x1)/dt_fast")
    integrator.addComputePerDof("v", "v + 0.5*b*dt*f1/m") # TODO: Additional velocity constraint correction?
    integrator.addConstrainVelocities();

    #
    # Velocity randomization
    #
    integrator.addComputePerDof("v", "sqrt(a)*v + sqrt(1-a)*sigma*gaussian")
    integrator.addConstrainVelocities();

    return integrator

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
                   Snapshots written every 10 * --interval steps.''')
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
                   metavar='FLOAT', help='''time step for integrator (outer
                   time-step for RESPA integrator) Default 1 fs''', default=1.0)
group.add_argument('--tfreq', dest='tfreq', type=float,
                   metavar='FLOAT', help='''frequency for Andersen thermostat.
                   Default 0.1 ps^-1''', default=0.1)
group.add_argument('--nrespa', dest='nrespa', type=int, metavar='INT',
                   default=1, help='''Number of inner time steps to run
                   (evaluating fast forces) for every outer timestep. Default is
                   1 (no RESPA). Best value to use for AMOEBA is at most 2. Only
                   AMOEBA is supported.''')
group.add_argument('--gamma_ln', dest='gamma_ln', type=float,
                   metavar='FLOAT', help='''collision frequency for Langevin
                   integrator. Default %(default)s ps-1''', default=0.0)
group.add_argument('--temp', dest='temp', type=float,
                   metavar='FLOAT', help='''target temperature for NVT
                   simulation. Default %(default)s K''', default=300.0)
group.add_argument('--nve', dest='nve', default=False, action='store_true',
                    help='Do constant energy simulation')
group.add_argument('--watch-for-errors', dest='hawkeye', default=False,
                   action='store_true',
                   help='''Watch energy every step and if energy becomes large
                   or NaN, print out each component to find where the energy is
                   blowing up. This may slow the simulation down considerably,
                   so it is primarily of use when debugging.''')

opt = parser.parse_args()

print('Command line:\n\t%s' % ' '.join(sys.argv)); sys.stdout.flush()

print('Parsing XML file %s' % opt.xml); sys.stdout.flush()
with open(opt.xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())

if opt.hawkeye:
    if opt.nrespa > 1:
        raise ValueError('Cannot use MTS integrator and watch for errors')
    groups_and_names = []
    for i, force in enumerate(system.getForces()):
#       if isinstance(force, mm.AmoebaMultipoleForce):
#           # Skip the multipole force, since it's so expensive
#           force.setForceGroup(20)
#           continue
        groups_and_names.append((type(force).__name__, i))
        force.setForceGroup(i)
    
    class ErrorDetectionReporter(app.StateDataReporter):
        def __init__(self, groups_and_names):
            self._reportInterval = 1
            self.groups_and_names = groups_and_names

        def describeNextReport(self, simulation):
            """Get information about the next report this object will generate.
    
            Parameters
            ----------
            simulation : Simulation
                The Simulation to generate a report for
    
            Returns
            -------
            tuple
                A five element tuple. The first element is the number of steps
                until the next report. The remaining elements specify whether
                that report will require positions, velocities, forces, and
                energies respectively.
            """
            return 1, False, False, False, True # only need energies

        def report(self, simulation, state):
            ene = state.getPotentialEnergy().value_in_unit(u.kilojoules_per_mole)
            if not math.isnan(ene) and ene < 1e5:
                return
            print('%30s %.6f kcal/mol' % ('Total Energy', ene))
            for name, i in self.groups_and_names:
                ene = simulation.context.getState(getEnergy=True, groups=1<<i).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
                print('%30s %.6f kcal/mol' % (name, ene), file=sys.stderr)
            sys.exit('Bad energy! Look at components')
    
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
    if opt.nrespa > 1:
        integrator = MTSVVVRIntegrator(opt.temp*u.kelvin,
               opt.gamma_ln/u.picosecond, opt.timestep*u.femtoseconds,
               system, opt.nrespa)
        print('MTSVVVR: %8.2fK, %8.2f ps-1, %8.2f fs, %3d inner steps' %
               (opt.temp, opt.gamma_ln, opt.timestep, opt.nrespa) ); sys.stdout.flush()
    else:
        integrator = mm.LangevinIntegrator(opt.temp*u.kelvin,
               opt.gamma_ln/u.picosecond, opt.timestep*u.femtoseconds)
        print('Langevin: %8.2fK, %8.2f ps-1, %8.2f fs' %
               (opt.temp, opt.gamma_ln, opt.timestep) ); sys.stdout.flush()
elif opt.nrespa > 1:
    slow = (mm.AmoebaMultipoleForce, mm.AmoebaVdwForce,
            mm.AmoebaGeneralizedKirkwoodForce, mm.AmoebaWcaDispersionForce)
    found_slow = False
    for force in system.getForces():
        if isinstance(force, slow):
            found_slow = True
            force.setForceGroup(1)
        else:
            force.setForceGroup(0)
    if not found_slow:
        raise ValueError('No slow AMOEBA forces found for MTS integrator!')
    # The list given to MTSIntegrator defining the time steps and force
    # decompositions is a list of 2-element tuples where the first element is
    # a force group and the second element is how many times to evaluate that
    # force group each "outer" time-step. So [(0, opt.nrespa), (1, 1)] means
    # force group 0 is executed nrespa times each time step and force group 1 is
    # executed only once. The slow forces are defined above as force group 1 and
    # all others as force group 0.
    integrator = mm.MTSIntegrator(opt.timestep*u.femtoseconds,
                                  [(0, opt.nrespa), (1, 1)])
    print('RESPA MTS Integrator: %8.2f fs outer time-step with %d inner steps' %
          (opt.timestep, opt.nrespa))
else:
    integrator = mm.VerletIntegrator(opt.timestep*u.femtoseconds)
    print('Verlet: %8.2f fs' % opt.timestep )

sim = app.Simulation(pdb.topology, system, integrator,
                     platform=mm.Platform.getPlatformByName('CUDA'),
                     platformProperties=dict(CudaPrecision='mixed') )

if opt.hawkeye:
    # Watch every step... slow!
    sim.reporters.append(ErrorDetectionReporter(groups_and_names))

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
        pmd.openmm.RestartReporter(opt.restart, 99999999, netcdf=True)
)
if opt.checkpoint is not None:
    sim.reporters.append(
            app.CheckpointReporter(opt.checkpoint, opt.interval*10)
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
