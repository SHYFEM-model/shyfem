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
# creates bsin for mpi runs
#
#------------------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

fem3d=$FEMDIR/fem3d

if [ $# -ne 2 ]; then
  echo "Usage: mpi_basin.sh npart grd-file"
  exit 1
fi

npart=$1
file=$2

basin=$( basename $file .grd )

newnode=$basin.$npart.node
visual=$basin.$npart.visual

#---------------------------------------------------
# cleaning files
#---------------------------------------------------

[ -f $newnode ] && rm -f $newnode
[ -f npart.grd ] && rm -f npart.grd
[ -f $basin.bas ] && rm -f $basin.bas

#---------------------------------------------------
# no MPI partitioning
#---------------------------------------------------

if [ $npart -le 1 ]; then
  $fem3d/shypre -silent -noopti $basin.grd
  if [ ! -f $basin.bas ]; then
    echo "*** error creating file $basin.bas"
    exit 1
  fi
  echo "created $basin.bas with no mpi partitioning"
  exit 0
fi

#---------------------------------------------------
# create partitioning
#---------------------------------------------------

echo "running shyparts..."
[ -f $newnode.grd ] && rm -f $newnode.grd
$fem3d/shyparts  -nparts $npart $basin.grd

if [ $? -ne 0 -o ! -f $newnode.grd ]; then
  echo "*** error creating file partition..."
  exit 1
fi

#---------------------------------------------------
# create grd-file for visualization - use "grid -ofT grd-file"
#---------------------------------------------------

echo "running shybas..."
$fem3d/shybas -silent -npart $newnode.grd

if [ ! -f npart.grd ]; then
  echo "*** error creating file npart.grd"
  exit 1
fi
mv npart.grd $visual.grd

#---------------------------------------------------
# create bas file
#---------------------------------------------------

echo "running shypre..."
$fem3d/shypre -silent -noopti -partition $newnode.grd $basin.grd

if [ ! -f $basin.bas ]; then
  echo "*** error creating file $basin.bas"
  exit 1
fi

#---------------------------------------------------
# final info
#---------------------------------------------------

echo ""
echo "bas file created : $basin.bas"
echo "file to visualize: $visual"
echo "  (use: grid -fT $visual)"
echo ""

#---------------------------------------------------
# end of routine
#---------------------------------------------------

