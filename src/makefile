compiler = mpif90

modules = globals.f90 utils.f90 tsp.f90
mod_files = $(modules:.f90=.mod)
mod_objects = $(modules:.f90=.o)

programs = tsp_ga.f90 tsp_ga_parallel.f90
program_objects = $(programs:.f90=.o)

objects = $(mod_objects) $(program_objects)

binaries = $(programs:.f90= )
binaries := $(addprefix ../run/, $(binaries))
to_clean := $(binaries) $(objects) $(mod_files)


.PHONY: all
all: $(binaries)

../run/tsp_ga_parallel: $(mod_objects) tsp_ga_parallel.o
	$(compiler) -O2 $^ -o $@

tsp_ga_parallel.o: tsp_ga_parallel.f90
	$(compiler) -c $<

../run/tsp_ga: $(mod_objects) tsp_ga.o
	-$(compiler) -O2 $^ -o $@

tsp_ga.o: tsp_ga.f90
	-$(compiler) -c $<

$(mod_objects): %.o: %.f90
	-$(compiler) -c $<

.PHONY: clean
clean:
	-rm $(to_clean)