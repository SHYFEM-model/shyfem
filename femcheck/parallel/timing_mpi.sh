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
# simple MPI fortran test
#
#--------------------------------------------

RunMpi()
{
  echo "============================================="
  echo "simple MPI test..."
  echo "using compiler $comp"
  echo "============================================="

  $comp --version > /dev/null 2>&1
  [ $? -ne 0 ] && echo "no such compiler $comp" && exit 1
  
  echo "compiling..."
  $comp timing_mpi.f
  [ $? -ne 0 ] && exit 1

  echo "running..."
  $run -np 1 a.out
  $run -np 3 a.out
  $run -np 5 a.out

  echo "finished..."
}

#--------------------------------------------

dir=/usr/bin
comp=$dir/mpif90
run=$dir/mpirun
text="mpi gfortran"

RunMpi

#--------------------------------------------

comp=mpifort
run=mpirun

RunMpi

#--------------------------------------------

