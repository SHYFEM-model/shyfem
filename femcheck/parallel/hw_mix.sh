#!/bin/sh
#
# simple MPI/OMP fortran test
#
#--------------------------------------------

dir=/usr/bin

echo "simple MPI/OMP test..."
echo "compiling..."
$dir/mpif90 -fopenmp hw_mix.f
[ $? -ne 0 ] && exit 1

echo "running..."
$dir/mpirun -np 3 a.out

echo "finished..."

