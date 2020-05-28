#!/bin/sh
#
# creates bsin for mpi runs
#
#---------------------------------------------------

if [ $# -ne 2 ]; then
  echo "Usage: mpi_basin.sh grd-file npart"
  exit 0
fi

file=$1
npart=$2

basin=$( basename $file .grd )

newnode=$basin.$npart.node
visual=$basin.$npart.visual

#---------------------------------------------------
# create partitioning
#---------------------------------------------------

shyparts  -nparts $npart $basin.grd

if [ $? -ne 0 -o ! -f $newnode.grd ]; then
  echo "*** error creating file partition..."
  exit 1
fi

#---------------------------------------------------
# create grd-file for visualization - use "grid -ofT grd-file"
#---------------------------------------------------

shybas -silent -npart $newnode.grd
mv npart.grd $visual.grd

#---------------------------------------------------
# create bas file
#---------------------------------------------------

shypre -silent -partition $newnode.grd $basin.grd

#---------------------------------------------------
# final info
#---------------------------------------------------

echo ""
echo "bas file created : $basin.bas"
echo "file to visualize: $visual"
echo ""

#---------------------------------------------------
# end of routine
#---------------------------------------------------

