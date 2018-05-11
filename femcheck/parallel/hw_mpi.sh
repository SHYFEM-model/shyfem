#!/bin/sh
#
# simple MPI fortran test
#
#--------------------------------------------

dir=/usr/bin
comp=$dir/mpif90
run=$dir/mpirun

comp=mpiifort
run=mpirun

echo "simple MPI test..."
echo "compiling..."
$comp hw_mpi.f
[ $? -ne 0 ] && exit 1

echo "running..."
$run -np 3 a.out

echo "finished..."

