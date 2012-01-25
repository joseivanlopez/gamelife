#!/bin/sh

# Loop varying dimensions world dimensions
for dimension in 1000 $(seq 5000 5000 50000) 
do
  # Sequential method
  ./gamelife -f examples/example1.txt -r $dimension -c 1000 -m 0 -t >> out_sequential.csv 
done
