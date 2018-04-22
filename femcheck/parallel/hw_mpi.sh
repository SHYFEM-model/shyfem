#!/bin/sh
#
# simple MPI fortran test
#
#--------------------------------------------

echo "simple MPI test..."
echo "compiling..."
mpif90 hw_mpi.f
[ $? -ne 0 ] && exit 1

echo "running..."
mpirun -np 3 a.out

echo "finished..."

