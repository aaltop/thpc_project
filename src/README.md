# Compilation instructions

Should all the necessary utilities be available, compilation can
most easily be perfomed by the command `make`. For a clean install, 
this command compiles like so:
> mpif90 -c globals.f90

> mpif90 -c utils.f90

> mpif90 -c tsp.f90

> mpif90 -c tsp_ga.f90

> mpif90 -O2 globals.o utils.o tsp.o tsp_ga.o -o ../run/tsp_ga

> mpif90 -c tsp_ga_parallel.f90

> mpif90 -O2 globals.o utils.o tsp.o tsp_ga_parallel.o -o ../run/tsp_ga_parallel

These steps can be performed individually for manual compilation. 
In either case, the compiler wrapper `mpif90` will be needed, and 
therefore the OpenMPI library available on the system.

The above compiled files can be easily cleaned by the command `make clean`.
This runs the following command:

> rm ../run/tsp_ga ../run/tsp_ga_parallel globals.o utils.o tsp.o tsp_ga.o tsp_ga_parallel.o globals.mod utils.mod tsp.mod