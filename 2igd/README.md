Examples
========

The first example comes from the crystal structure solved in PDB code 2IGD. The
crystal was stripped of water and hydrogen atoms and solvent was added by
``tleap`` in AmberTools 16 to create the file ``2igd_solvated.pdb``.

This example walks through preparing this system (which has only standard
residues) for simulation with OpenMM using the AMOEBA force field.

Preparation
-----------

The first step is system preparation. We run the script ``prep_openmm.py`` for
this. Feel free to use the ``--help`` flag to look at the available options. In
this example, I will specify hydrogen masses of 3 Daltons (for hydrogen-mass
repartitioning):

```
$ ../prep_openmm.py -p 2igd_solvated.pdb -f amoeba2013 -m 3.0
Parsing the PDB file [2igd_solvated.pdb]...
Loading the force field [amoeba2013.xml]...
Creating the System...
No constraints applied
Repartitioning hydrogen masses to 3 daltons
Adjusting the vdW cutoff to 9 Angstroms...
Serializing the System...
Done.
```

Of course, if you did not want to repartition masses, leave out ``-m 3.0``.  If
you look at the top of the ``system.xml`` file, you will see

```
$ head -n 15 system.xml
<?xml version="1.0" ?>
<System openmmVersion="7.0" type="System" version="1">
	<PeriodicBoxVectors>
		<A x="6.7929" y="0" z="0"/>
		<B x="-2.264163559406279" y="6.404455775962287" z="0"/>
		<C x="-2.264163559406279" y="-3.2019384603140684" z="5.54658849047036"/>
	</PeriodicBoxVectors>
	<Particles>
		<Particle mass="8.030999999999999"/>
		<Particle mass="3"/>
		<Particle mass="3"/>
		<Particle mass="3"/>
		<Particle mass="10.018999999999998"/>
		<Particle mass="3"/>
		<Particle mass="8.026999999999997"/>
```

We can see the truncated octahedron periodic box vectors and the first several
atoms. The first atom is a nitrogen with 3 attached hydrogens. The masses have
been repartitioned such that the NH3 group is still 17.031 daltons, but each
H-atom is now 3 daltons. You would see the standard atomic masses if you had
left the ``-m 3.0`` argument out from the ``prep_openmm.py`` command.

Minimization
------------

Now that we have the ``system.xml`` file set up, we can minimize it. In this
example, we will use the command:

```
$ ../minimize.py -p 2igd_solvated.pdb -s system.xml -x 2igd_solvated.min.xml
Loading the system [system.xml] and PDB [2igd_solvated.pdb] files...
Creating the Simulation...
Minimizing...
 Initial energy = -57365.2367 kcal/mol
   Final energy = -86110.7256 kcal/mol
Finished. Writing serialized XML restart file...
Done!
```

After this is finished, you should see a ``2igd_solvated.min.xml`` file. Now
it's time to move on to the next phase... dynamics.

Dynamics
--------

The script we use in this step will work for heating, density equilibration (via
the Monte Carlo barostat), and regular NVT or NVE simulations. Consult
``runmd.py --help`` for more information.

To heat this system, we will use Langevin dynamics (which is obtained by setting
the friction coefficient to a positive number):

```
$ ./runmd.py -p 2igd_solvated.pdb -s 2igd_solvated.min.xml --restrain @CA,C,N \
             --gamma_ln 1.0 -x 2igd_solvated.heat.nc -o 2igd_solvated.heat.mdout \
             -r 2igd_solvated.heat.ncrst -k 10.0 -n 10000
Command line:
	../runmd.py -p 2igd_solvated.pdb -s 2igd_solvated.min.xml --restrain @CA,C,N --gamma_ln 1.0 -x 2igd_solvated.heat.nc -o 2igd_solvated.heat.mdout -r 2igd_solvated.heat.ncrst -k 10.0 -n 10000
Parsing XML file system.xml
No Andersen thermostat found in the system.xml file
Parsing PDB file 2igd_solvated.pdb
Adding restraints (k=10.0 kcal/mol/A^2) from @CA,C,N
Langevin:   300.00K,     1.00 ps-1,     1.00 fs
Setting coordinates and velocities from restart file 2igd_solvated.min.xml
Running the simulation for 10000 steps!
```

You can see that the script prints details about the calculation, including the
command-line used, for logging purposes. You are encouraged to redirect this
output to a file for future reference. While the simulation is running, you can
monitor the ``2igd_solvated.heat.mdout.info`` file to track the speed and
progress of the calculation. For example, it may look like this:

```
-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

  On step 1500 of 10000

  Total Energy     =  -56256.4151
  Potential Energy =  -70983.7228
  Kinetic Energy   =   14727.3077

 Time for last      500 steps:   104.0117 s. (0.415 ns/day)
 Time for all      1500 steps:   313.5034 s. (0.413 ns/day)

 Estimated time to completion: 29.470 min.

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
```

Continuing Dynamics
-------------------

You can always continue dynamics using the NetCDF restart file you wrote in the
previous step.
