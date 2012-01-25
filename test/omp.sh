#!/bin/sh

# Loop varying dimensions world dimensions
for dimension in 1000 $(seq 5000 5000 50000) 
do
  # OMP method
    for thr in 1 2 4 8 16 32
    do
        OMP_NUM_THREADS=$thr ./gamelife -f examples/example1.txt -r $dimension -c 1000 -m 1 -t >> out_omp$thr.csv
    done
done
