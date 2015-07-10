#!/bin/sh
#
#---------------------------------------------------

dirs="fem3d femplot"

if [ $# -eq 0 ]; then
  echo "Usage: cp_fem.sh number"
  exit 1
fi

number=$1

for dir in $dirs
do
  d=${dir}_$number
  if [ -d $d ]; then
    echo "directory exists... cannot copy: $d"
    exit 1
  fi
done

for dir in $dirs
do
  d=${dir}_$number
  echo "copying from $dir to $d"
  mkdir $d
  cp $dir/* $d
done

dir=fem3d/bin
d=fem3d_$number/bin
echo "copying from $dir to $d"
mkdir $d
cp $dir/* $d

