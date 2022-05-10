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
# simple MPI/OMP fortran test
#
#--------------------------------------------

RunOmpMpi()
{
  echo "============================================="
  echo "simple MPI/OMP test..."
  echo "compiling..."
  echo "============================================="

  $comp --version > /dev/null 2>&1
  [ $? -ne 0 ] && echo "no such compiler $comp" && exit 0

  $comp hw_mix.f
  [ $? -ne 0 ] && exit 1

  echo "running..."
  $run -np 3 a.out

  echo "finished..."
}

#--------------------------------------------

dir=/usr/bin
comp="$dir/mpif90 -fopenmp"
run=$dir/mpirun

RunOmpMpi

#--------------------------------------------

comp="mpiifort -fopenmp"

RunOmpMpi

#--------------------------------------------

