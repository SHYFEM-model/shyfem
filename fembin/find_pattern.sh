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
# finds a pattern in all source code of SHYFEM
#
#------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

if [ $# -eq 0 ]; then
  echo "Usage: find_pattern.sh pattern"
  exit 0
fi

pattern=$1

base=$FEMDIR

dirs="fem3d femplot femadj femspline"
fembfm=femersem
femgotm=femgotm

for dir in $dirs
do
  echo " ----------------- $dir -------------------"
  cd $base/$dir
  grep -i "$pattern" *.f *.F90 *.f90 2> /dev/null | grep -v __genmod
done

dir=$fembfm
cd $base/$dir
echo " ----------------- $dir -------------------"
find . -name "*.F90" -print0 2> /dev/null | xargs --null grep -i $pattern
find . -name "*.f90" -print0 2> /dev/null | xargs --null grep -i $pattern
find . -name "*.f"   -print0 2> /dev/null | xargs --null grep -i $pattern

dir=$femgotm
cd $base/$dir
echo " ----------------- $dir -------------------"
find . -name "*.F90" -print0 2> /dev/null | xargs --null grep -i $pattern
find . -name "*.f90" -print0 2> /dev/null | xargs --null grep -i $pattern
find . -name "*.f"   -print0 2> /dev/null | xargs --null grep -i $pattern

