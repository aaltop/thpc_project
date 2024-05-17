#!/bin/bash

input="serial_args.txt"
while IFS= read -r line; do
    echo "args: $line"
    ./tsp_ga $line
done < $input