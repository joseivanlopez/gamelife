#!/bin/sh

# Loop varying dimensions world dimensions
for dimension in 1000 $(seq 5000 5000 50000) 
do
    # MPI method
    mpirun -np 4 --hostfile hostfile88.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 2 -t >> out_mpi_np4_hf88.csv
    mpirun -np 4 --hostfile hostfile22.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 2 -t >> out_mpi_np4_hf22.csv
    mpirun -np 8 --hostfile hostfile88.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 2 -t >> out_mpi_np8_hf88.csv
    mpirun -np 8 --hostfile hostfile22.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 2 -t >> out_mpi_np8_hf22.csv
    mpirun -np 16 --hostfile hostfile88.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 2 -t >> out_mpi_np16_hf88.csv
    mpirun -np 16 --hostfile hostfile22.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 2 -t >> out_mpi_np16_hf22.csv

    # MPI optimized method
    mpirun -np 4 --hostfile hostfile88.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 3 -t >> out_mpiopt_np4_hf88.csv
    mpirun -np 4 --hostfile hostfile22.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 3 -t >> out_mpiopt_np4_hf22.csv
    mpirun -np 8 --hostfile hostfile88.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 3 -t >> out_mpiopt_np8_hf88.csv
    mpirun -np 8 --hostfile hostfile22.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 3 -t >> out_mpiopt_np8_hf22.csv
    mpirun -np 16 --hostfile hostfile88.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 3 -t >> out_mpiopt_np16_hf88.csv
    mpirun -np 16 --hostfile hostfile22.txt gamelife -f examples/example1.txt -r $dimension -c 1000 -m 3 -t >> out_mpiopt_np16_hf22.csv
done
