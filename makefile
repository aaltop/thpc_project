compiler = gfortran

modules = globals.f90 tsp.f90
mod_files = $(modules:.f90=.mod)
mod_objects = $(modules:.f90=.o)

programs = tsp_ga.f90
program_objects = $(programs:.f90=.o)

objects = $(mod_objects) $(program_objects)

prog = $(programs:.f90= )
to_clean := $(prog)



.PHONY: all
all: $(prog)

.PHONY: run
run_serial: tsp_ga
	./tsp_ga

tsp_ga: tsp_ga.o $(mod_objects)
	-$(compiler) -O2 $? -o $@

$(program_objects): %.o: %.f90
	-$(compiler) -c $<

$(mod_objects): %.o: %.f90
	-$(compiler) -c $<

.PHONY: clean
clean:
	-rm $(to_clean)