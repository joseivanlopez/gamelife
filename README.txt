To compile:
$ mpicc gamelife.c -o gamelife -fopenmp

To execute:
  With Sequential or openMP mode:
  $ ./gamelife [options] -m [0|1]
  $ OMP_NUM_THREADS=N ./gamelife [options] -m 1
  
  With MPI mode:
  $ mpirun -np N gamelife [options] -m 2
