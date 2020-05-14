#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

femdir=${SHYFEMDIR:=$HOME/shyfem}
bindir=$femdir/fem3d
binbindir=$femdir/fem3d/bin

compiler=f77
options=

compiler=ifort
options="-warn nouncalled"

compiler=gfortran
options=

rm -f out_all.txt

for file
do
  echo $file | tee -a out_all.txt
  $compiler $options $file 2> out.txt
  $binbindir/dependency.pl -main -$compiler out.txt | tee -a out_all.txt
  [ $? -ne 0 ] && exit 1
done

elab_freq.pl out_all.txt

