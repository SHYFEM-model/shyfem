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
# resets bounding box from command line

tmp=tmp.tmp

if [ $# -le 1 ]; then
  echo "Usage: psmakebb.sh bound-box file[s]"
  echo '   example: psmakebb.sh "79 283 516 607" plot.*.eps'
  exit 0
fi

bb=$1
shift

##################################################################

for file
do
  echo "$file"

  sed -e "s/^%%BoundingBox:.*/%%BoundingBox: $bb/" $file > $tmp
  mv $tmp $file
done

##################################################################

