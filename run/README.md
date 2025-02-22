# Running instructions

Should all the necessary utilities be available and compilation already performed, running can
most easily be perfomed by the command 

> make run

An example of running the command is included in the file `example_output.txt`. 
Running the programs will also produce files in the subdirectory `generated_data`. 
These files will contain on each line the best route for each generation of running 
the algorithm. The first value on a line will be the distance of the route, and the 
rest of the line is the indices of the points of the route. `breed.txt` will contain 
the serial algorithm run, files `parallel1_<i>.txt` will contain the parallel 
algorithm run for process i, files `parallel2_<i>.txt` will contain the parallel 
run of the serial algorithm for process i (the serial algorithm is also run in 
each process in the parallel program). Files 
`random_search.txt`, `random_breed.txt`, and `random_breed_better.txt` 
will contain the runs of more naive algorithms (these are run in the serial program). 
Depending on the arguments to the programs, `generated_data` can be quite large, 
but for the default arguments, should be less than 100MB.


Alternatively, the serial program can be run by the shell script 

> ./serial.sh

and the parallel program by running the shell script 

> ./parallel.sh i

where `i` is the number of processes to run with. Both of these are run by 
the `make run`, the parallel one with two processes. When running by either 
of these or `make run`, the arguments to the program can supplied using the 
text files `serial_args.txt` and `parallel_args.txt`. These should already 
contain usable example arguments. Multiple lines of arguments in the files 
will all be run consecutively.

Finally, the serial program can also be run
by

> ./tsp_ga <serial_args>

and the parallel one by

> mpirun --oversubscribe -n i  ./tsp_ga_parallel <parallel_args>

Here, `i` is the number of processes to run the parallel program with.
Note that the flag `--oversubscribe` will create more processes than
the number of cores available on the CPU, for example, if `i` is specified
to be more than the number of cores.

Running either program without any arguments should print out the list
of arguments and their meaning. Example arguments are in the 
files `serial_args.txt` and `parallel_args.txt`. For example, the parallel 
program can be run by supplying the arguments 
(these are the default arguments in the file `parallel_args.txt`).

> circle_dist.txt 10 15 0.95 10000 2 5

- `circle_dist.txt` is provided as the example distance text file. 
This file and any similar "distance matrix" file should contain the row 
or column length of the (square) distance matrix on the first line of the 
file. The rest of the file is the distance matrix itself. Row i column j 
of the distance matrix contains the "distance" between some point i and point j.*

- 10 (P) is the number of parents that breed during a generation of the algorithm.*

- 15 (C) is the number of children bred by the parents. C should be specified such that P <= C <= P*(P-1).*

- 0.95 is the chance that a child mutates.*

- 10000 is the number of generations to run the algorithm for.*

- 2 is the number of migrators.

- 5 is how often there is a migration.

Of these, the ones marked with * are common to 
(that is, supplied to both) the serial and parallel algorithm.

All files created by running the programs can be cleaned easily by 
running the command `make clean`. This runs the following command:

> rm -r generated_data

As such, cleaning involves deleting the `generated_data` folder.