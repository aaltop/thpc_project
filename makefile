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
all:
	$(MAKE) -C src all


run_commands = run_serial run_parallel
.PHONY: run
run: $(run_commands)

.PHONY: $(run_commands)
$(run_commands): all
	$(MAKE) -C run $@


.PHONY: clean
clean:
	$(MAKE) -C run clean
	$(MAKE) -C src clean