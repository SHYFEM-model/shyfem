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

routine=timing_omp.f

echo "simple OMP test..."
echo "compiling..."

#gfortran -fopenmp $routine
#ifort -qopenmp $routine
nvfortran -mp $routine

[ $? -ne 0 ] && exit 1

echo "running..."
./a.out

echo "finished..."

