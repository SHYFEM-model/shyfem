#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

hfiles=$( ls *.h )
files=$( ls *.f *.h *.f90 )

movedir=obsolete

show=NO
move=NO
file=NO
[ "$1" = "-show" ] && show=YES && shift
[ "$1" = "-move" ] && move=YES && shift
[ "$1" = "-file" ] && file=YES && shift

for h in $hfiles
do
  n=$( grep "$h" $files | grep -i include | wc -l )
  echo "$n   $h"
  if [ $show = "YES" ]; then
    grep "$h" $files | grep -i include
  fi
  if [ $file = "YES" ]; then
    grep "$h" $files | grep -i include | sed -e 's/:.*//' | uniq
  fi
  if [ $n -eq 0 ]; then
    if [ $move = "YES" ]; then
      echo "moving $h to $movedir"
      mkdir -p $movedir
      mv $h $movedir
    fi
  fi
done

#------------------------------------------------------------------------

