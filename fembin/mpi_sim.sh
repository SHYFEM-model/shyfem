#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# compines mpi_basin and mpi_run
#
#--------------------------------------------------

if [ $# -lt 2 ]; then
  echo "Usage: mpi_sim.sh nproc [shyfem-options] str-file"
  echo "  routine prepares basin and runs simulation"
  exit 1
fi

strfile=${@: -1}		#last argument (str-file)
np=$1
shift

grdfile=$( cat $strfile | sed -n '/\$title/,$p' | head -4 | tail -1 )
grdfile=$grdfile.grd

echo "nproc: $np"
echo "grd-file: $grdfile"
echo "str-file: $strfile"
echo "rest of line: $@"

#exit 0

#--------------------------------------------------

mpi_basin.sh $np $grdfile

mpi_run.sh $np $@

#--------------------------------------------------

