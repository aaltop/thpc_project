#!/bin/bash

input="parallel_args.txt"
while IFS= read -r line; do
        echo "args: $line"
        # needs </dev/null because mpirun stops the while from running, apparently
        mpirun --oversubscribe -n $1 ./tsp_ga_parallel $line </dev/null
done < $input