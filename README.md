What is this repo?
==================

This is a repository with some preparatory scripts for getting started running
the AMOEBA force field using OpenMM.  Jason Swails originated this project,
and Dave Case added some tweaks and documentation.

A warning: this code is still in a preliminary stage, and is not for the
faint of heart.  It was written to allow people familiar with Amber to run
GPU-accelerated codes with the Amoeba force field.  The expected use-case
is that system preparation and analysis would use AmberTools, and that
the minimization and MD steps would use OpenMM.  (There is little that is
specific to the Amoeba force field, and that package could in principle be
modified to use any other force field supported by OpenMM.)

This workflow utilizes *only* OpenMM, even using its modelling utilities
to parametrize the system. You do *not* need to install the Tinker-OpenMM
interface to use this workflow.

However, this is incapable of handling non-standard residues. You would
need to use procedures available in Tinker to create the force fields
for these, then hand-edit the amoeba2013_dac.xml file to add them.  I
don't have any instructions here, but if you compare amoeba2013_dac.xml
to the amoeba2013.xml that comes with OpenMM, you can see how I added a
non-covalent ligand.  Follow the same procedure for your ligand.

Getting started
===============

There are scripts in here to prepare the system for simulation, minimize the
structure, and the perform standard dynamics (NVE, NVT, and NpT).

Dependencies
------------

This whole workflow is implemented in Python. You will need to have the
following tools installed:

- [Python](https://www.python.org)
- [numpy](http://www.numpy.org/)
- [OpenMM](http://openmm.org/)
- [ParmEd](http://parmed.github.io/ParmEd/html/index.html)

The easiest way to get all of this up and working is to install
[Anaconda](https://www.continuum.io/downloads) or
[Miniconda](http://conda.pydata.org/miniconda.html). You can then use ``conda``
to install OpenMM quickly and easily! Use the following commands:

```
conda create -n openmm  omnia openmm parmed netCDF4 numpy  # only done once!
conda activate openmm
   #  do work here as described below
conda deactivate
```

Preparing your system
---------------------

The ``prep_openmm.py`` script will take a PDB file (that you should have
prepared and solvated previously using something like ``tleap`` from the
AmberTools suite of programs) that has unit cell information (stored as a
``CRYST1`` record) and parametrize it with a target force field.

To use ``prep_openmm.py``, you can ask for help via the ``--help`` flag. A
sample command-line that will create a ``system.xml`` file from
``structure.pdb`` using the AMOEBA 2013 protein force field is shown below:


```
prep_openmm.py -p structure.pdb -f amoeba2013 -s system.xml
```

For standard AMOEBA simulations, most of the defaults are appropriate. This is
also the stage where you would implement H-mass repartitioning or SHAKE to try
and limit high-frequency motions.

Minimizing your system
----------------------

The ``minimize.py`` script is responsible for doing a local geometry
optimization (using OpenMM's LBFGS minimizer). It takes both an XML file of a
serialized OpenMM System and a PDB file as the source of coordinates. It
produces an XML file of a serialized OpenMM State (which contains the minimized
coordinates) for use in the next step.

Running dynamics
----------------

The main engine here for running dynamics is ``runmd.py``. There are many
options, and this can be used to run restrained dynamics, NVT, NpT (with either
anisotropic or isotropic pressure scaling using the Monte-Carlo Barostat), or
NVE. Note that for NVE, you should set an appropriately small default induced
dipole tolerance in the preparatory step.

Like with the other Python scripts here, you can use the ``-h/--help`` flags to
get a full listing of the available options.

This script will generate a tab-delimited output file with energies and other
state variables, a NetCDF trajectory file (that can be analyzed alongside the
original PDB file with any trajectory analysis program, like cpptraj and
pytraj), and a NetCDF restart file for use in continuing the simulation.

BUGS
====

Only amoebapro2013 (and earlier) force fields are currently offered in
OpenMM.  Comments in the OpenMM source indicate that support for more recent
versions of Amoeba, and for things like nucleic acids, are planned.

These scripts appear to work with OpenMM as of January, 2021.  But they are
not associated with that project, and may not track changes with OpenMM.
Users should be sure to compare short simulations with those you get from
Tinker, gem.pmemd (part of AmberTools) or other programs to verify accuracy.
That said, these scripts are short and simple, and just call OpenMM
functionality, so they should benefit from all the testing done by the
OpenMM community.
