#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# extracts lines from grd file
#
#---------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: extract_lines.sh [-h|-help] grd-file"
  exit 0
elif [ $1 = '-h' -o $1 = '-help' ]; then
  echo "Usage: extract_lines.sh grd-file"
  echo "  output is written to file lines.grd"
  exit 0
fi

exgrd -lES $1

mv new.grd lines.grd

