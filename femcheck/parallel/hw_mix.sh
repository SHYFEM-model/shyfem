#!/bin/sh
#
# simple MPI/OMP fortran test
#
#--------------------------------------------

echo "simple MPI/OMP test..."
echo "compiling..."
mpif90 -fopenmp hw_mix.f
[ $? -ne 0 ] && exit 1

echo "running..."
mpirun -np 3 a.out

echo "finished..."

