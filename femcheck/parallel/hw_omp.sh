#!/bin/sh
#
# simple OMP fortran test
#
#--------------------------------------------

echo "simple OMP test..."
echo "compiling..."
gfortran -fopenmp hw_omp.f
[ $? -ne 0 ] && exit 1

echo "running..."
./a.out

echo "finished..."

