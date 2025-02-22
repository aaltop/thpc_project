# Travelling salesman problem, genetic algorithm in serial and parallel

This codebase implements a genetic algorithm solution to the travelling
salesman problem, both in serial and parallel form.

# Instructions

Due to the use of parallel algorithms, OpenMPI will be needed on the
running system such that the `mpif90`compiler wrapper can be used 
for compilation and the `mpirun` command used for execution of parallel 
programs. The base Fortran compiler used in development 
(and which the mpif90 wraps) has been gfortran from the 9.4.0 GCC. The programs
have been developed on a Windows Subsystem for Linux 2 with Ubuntu 20.04.

## Compilation

Refer further to the README in the src/ directory for compilation instructions.
Should `make` be available, compilation can be easily performed from 
the top-level directory by running the command

> make

which simply delegates
to the Makefile in the src/ directory.

## Running

Refer further to the README in the run/ directory for running instructions.
The programs can be both compiled and run from the top-level directory by the command

> make run

 which delegates compilation and running to the respective directories.

Files and directories created during compilation and running can be cleaned
using

> make clean

 in the top-level directory.