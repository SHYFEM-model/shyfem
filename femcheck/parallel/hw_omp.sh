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

RunOmp()
{
  echo "============================================="
  echo "simple OMP test..."
  echo "compiling..."
  echo "============================================="

  $comp --version > /dev/null 2>&1
  [ $? -ne 0 ] && echo "no such compiler $comp" && exit 0

  $comp hw_omp.f
  [ $? -ne 0 ] && exit 1

  echo "running..."
  ./a.out

  echo "finished..."
}

#--------------------------------------------

comp="gfortran -fopenmp"

RunOmp

#--------------------------------------------

comp="fintel -fopenmp"

RunOmp

#--------------------------------------------

