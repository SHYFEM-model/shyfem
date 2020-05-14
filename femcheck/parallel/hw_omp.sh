#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
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

