#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

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

