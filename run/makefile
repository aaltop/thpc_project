
.PHONY: run
run: run_serial run_parallel

.PHONY: run_parallel
run_parallel:
	./parallel.sh 2

.PHONY: run_serial
run_serial:
	./serial.sh

.PHONY: clean
clean:
	-rm -r generated_data