#!/bin/sh

femdir=${SHYFEMDIR:=$HOME/shyfem}
bindir=$femdir/fem3d
binbindir=$femdir/fem3d/bin

compiler=f77
options=

compiler=ifort
options="-warn nouncalled"

compiler=gfortran
options=

for file
do
  echo $file
  $compiler $options $file 2> ggg
  $binbindir/dependency.pl -main -$compiler ggg
  [ $? -ne 0 ] && exit 1
done

