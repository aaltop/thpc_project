compiler = mpif90

modules = globals.f90 utils.f90 tsp.f90
mod_files = $(modules:.f90=.mod)
mod_objects = $(modules:.f90=.o)

programs = tsp_ga.f90 tsp_ga_parallel.f90
program_objects = $(programs:.f90=.o)

objects = $(mod_objects) $(program_objects)

prog = $(programs:.f90= )
to_clean := $(prog) $(objects) $(mod_files)



.PHONY: all
all: $(prog)

.PHONY: run
run: run_serial run_parallel

.PHONY: run_parallel
run_parallel: tsp_ga_parallel
	mpirun --oversubscribe -n 4 ./tsp_ga_parallel wg59_dist_array.txt 10 15 0.95 10000 2 5

tsp_ga_parallel: $(mod_objects) tsp_ga_parallel.o
	$(compiler) -O2 $^ -o $@

tsp_ga_parallel.o: tsp_ga_parallel.f90
	$(compiler) -c $<

.PHONY: run_serial
run_serial: tsp_ga
	./tsp_ga wg59_dist_array.txt 10 15 0.95 10000

tsp_ga: $(mod_objects) tsp_ga.o
	-$(compiler) -O2 $^ -o $@

tsp_ga.o: tsp_ga.f90
	-$(compiler) -c $<

$(mod_objects): %.o: %.f90
	-$(compiler) -c $<

.PHONY: clean
clean:
	-rm $(to_clean)
	-rm -r generated_data