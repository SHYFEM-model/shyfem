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
  echo "trying compiler $comp"
  $comp --version > /dev/null 2>&1
  if [ $? -ne 0 ]; then
    echo "no such compiler $comp" 
    echo "============================================="
    return
  fi
  echo "simple MPI test... $text"
  echo "compiling as: $comp hw_mpi.f"
  echo "running as: $command"
  echo "============================================="

  echo "compiling..."
  $comp hw_mpi.f
  [ $? -ne 0 ] && exit 1

  echo "running..."
  $command

  echo "finished..."
}

#--------------------------------------------

dir=/usr/bin
comp=$dir/mpif90
run=$dir/mpirun
command="$run -np 3 ./a.out"
text="mpi gfortran"

RunMpi

#--------------------------------------------

comp=mpiifort
run=mpirun
command="$run -np 3 ./a.out"
text="mpi intel"

RunMpi

#--------------------------------------------
exit 0
#--------------------------------------------

