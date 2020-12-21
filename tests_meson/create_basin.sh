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

if [ $# -lt 2 ] || [ $# -gt 3 ] ; then
  echo "Usage: mpi_basin.sh npart grd-file grd-file-path"
  exit 1
fi

npart=$(($1))
file=$2
if [ $# -eq 3 ] ; then
  tests_dir=$3
  rm $file > /dev/null 2>&1
  ln -s $tests_dir$file $file
else
  tests_dir=''
fi

basin=$( basename $file .grd )

if [ $npart -eq 1 ] ; then
  echo "run ../shypre $file"
  ../shypre $file
  echo "done running ../shypre $file"
else

  newnode=$basin.$npart.node
  visual=$basin.$npart.visual

  #---------------------------------------------------
  # cleaning files
  #---------------------------------------------------

  [ -f $newnode ] && rm -f $newnode
  [ -f npart.grd ] && rm -f npart.grd
  [ -f $basin.bas ] && rm -f $basin.bas

  #---------------------------------------------------
  # create partitioning
  #---------------------------------------------------

  echo "running shyparts..."
  ../shyparts  -nparts $npart $basin.grd

  if [ ! -f $newnode.grd ]; then  # -o $? -ne 0
       
    echo "*** error creating file partition $newnode.grd "
    exit 1
  fi

  #---------------------------------------------------
  # create grd-file for visualization - use "grid -ofT grd-file"
  #---------------------------------------------------

  echo "running shybas..."
  ../shybas -silent -npart $newnode.grd

  if [ ! -f npart.grd ]; then
    echo "*** error creating file npart.grd"
    exit 1
  fi
  mv npart.grd $visual.grd

  #---------------------------------------------------
  # create bas file
  #---------------------------------------------------

  echo "running shypre..."
  ../shypre -silent -noopti -partition $newnode.grd $basin.grd

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

fi

#---------------------------------------------------
# end of routine
#---------------------------------------------------
rm $basin.$npart.bas > /dev/null 2>&1
mv $basin.bas $basin.$npart.bas
if [ ! -f $basin.$npart.bas ]; then
  echo "ERROR : bas file ${_shyfem_input} was not created"
  exit 1
fi
exit 0
