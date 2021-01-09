What is this repo?
==================

This is a repository with some preparatory scripts for getting started running
the AMOEBA force field using OpenMM.

This workflow utilizes *only* OpenMM, even using the modelling utilities inside
to parametrize the system. You do *not* need to install the Tinker-OpenMM
interface to use this workflow.

However, this is incapable of handling non-standard residues. You will need to
use the standard Tinker workflow in that case.

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
