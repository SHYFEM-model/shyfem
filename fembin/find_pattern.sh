#!/bin/sh
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
  grep -i "$pattern" *.f *.F90 *.f90
done

dir=$fembfm
cd $base/$dir
echo " ----------------- $dir -------------------"
find . -name "*.F90" -print0 | xargs --null grep -i $pattern
find . -name "*.f90" -print0 | xargs --null grep -i $pattern
find . -name "*.f" -print0 | xargs --null grep -i $pattern

dir=$femgotm
cd $base/$dir
echo " ----------------- $dir -------------------"
find . -name "*.F90" -print0 | xargs --null grep -i $pattern
find . -name "*.f90" -print0 | xargs --null grep -i $pattern
find . -name "*.f" -print0 | xargs --null grep -i $pattern

