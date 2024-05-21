#!/bin/bash

input="parallel_args.txt"
while IFS= read -r line; do
        for i in {4,}; do
            echo "args: $line"
            # needs </dev/null because mpirun stops the while from running, apparently
            mpirun --oversubscribe -n $i ./tsp_ga_parallel $line </dev/null
        done
done < $input