#!/bin/sh

hfiles=$( ls *.h )
files=$( ls *.f *.h *.f90 )

show=NO
[ "$1" = "-show" ] && show=YES

for h in $hfiles
do
  n=$( grep "$h" $files | grep -i include | wc -l )
  echo "$n   $h"
  if [ $show = "YES" ]; then
    grep "$h" $files | grep -i include
  fi
done

