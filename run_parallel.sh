#!/bin/bash

input="parallel_args.txt"
while IFS= read -r line; do
    echo "args: $line"
    mpirun --oversubscribe -n 4 ./tsp_ga_parallel $line
done < $input