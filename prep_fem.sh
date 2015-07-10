#!/bin/sh
#
#---------------------------------------------------

dirs="fem3d femplot"

if [ $# -eq 0 ]; then
  echo "Usage: prep_fem.sh number"
  exit 1
fi

number=$1

for dir in $dirs
do
  d=${dir}_$number
  if [ ! -d $d ]; then
    echo "directory does not exist... cannot prepare: $d"
    exit 1
  fi
done

for dir in $dirs
do
  d=${dir}_$number
  prog=./prep_$dir.sh
  echo "preparing $d with $prog"
  cd $d
  $prog
  cd ..
done

