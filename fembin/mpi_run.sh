#!/bin/sh
#
# runs shyfem model in mpi mode
#
#-----------------------------------------

fem3d=~/shyfem/fem3d
shyfem=$fem3d/shyfem

str=$1
proc=$2
if [ $# -lt 2 ]; then
  echo "Usage: mpi_run.sh str-file nproc"
  exit 1
fi

#-----------------------------------------

if [ $proc = 0 ]; then
  $shyfem $str
  status=$?
elif [ $proc = 1 ]; then
  $shyfem -mpi $str
  status=$?
else
  /usr/bin/mpirun -np $proc $shyfem -mpi $str
  status=$?
fi

#-----------------------------------------

if [ $status -eq 99 ]; then
  exit 0
else
  echo "exit status is $status ... error"
  exit 3
fi

#-----------------------------------------

echo "command run with proc=$proc and status=$status"

#-----------------------------------------

