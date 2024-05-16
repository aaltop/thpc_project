compiler = gfortran

modules = globals.f90 utils.f90 tsp.f90
mod_files = $(modules:.f90=.mod)
mod_objects = $(modules:.f90=.o)

programs = tsp_ga.f90
program_objects = $(programs:.f90=.o)

objects = $(mod_objects) $(program_objects)

prog = $(programs:.f90= )
to_clean := $(prog) $(objects) $(mod_files)



.PHONY: all
all: $(prog)

.PHONY: run
run_serial: tsp_ga
	./tsp_ga

tsp_ga: $(mod_objects) tsp_ga.o
	-$(compiler) -O2 $? -o $@

$(program_objects): %.o: %.f90
	-$(compiler) -c $<

$(mod_objects): %.o: %.f90
	-$(compiler) -c $<

.PHONY: clean
clean:
	-rm $(to_clean)