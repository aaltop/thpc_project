compiler = gfortran
deps = tsp_ga.f90
prog = $(deps:.f90= )
to_clean := $(prog)



.PHONY: all
all: $(prog)

.PHONY: run
run: $(prog)
	./run_script.sh

$(prog): %: %.f90
	-$(compiler) -O2 $< -o $@


.PHONY: clean
clean:
	-rm $(to_clean)