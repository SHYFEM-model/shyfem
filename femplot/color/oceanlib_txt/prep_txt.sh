#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

gfortran prep_txt.f

list="FreeSurface-rgb.txt Velocity-rgb.txt Vorticity-rgb.txt OxygenVar-rgb.txt"
datfiles=""

for file in $list
do
  name=$( echo $file | sed -e 's/-rgb.txt//' )
  echo "$file -> $name.dat"
  a.out < $file > $name.tmp
  ./prep_txt.pl $name.tmp > $name.dat
  datfiles="$datfiles $name.dat"
done

cat $datfiles > oceanlib_txt.dat

