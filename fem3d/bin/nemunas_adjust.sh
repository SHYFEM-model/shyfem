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
# adjusts some code for old gfortran compiler
#
#--------------------------------------------------

dir=fem3d/bin

if [ "$#" -eq 0 ]; then
  echo "Usage: nemunas_adjust.sh {-nemunas|-original}"
  exit 0
fi

if [ "$1" = "-original" ]; then
  option="-undo"
elif [ "$1" = "-nemunas" ]; then
  option=""
else
  echo "unknown option: $1"
  exit 1
fi

if [ ! -f VERSION ]; then
  echo "not in base directory... exiting"
  exit 1
fi

$dir/nemunas_adjust.pl $option Rules.make
diff Rules.make.bak Rules.make
$dir/nemunas_adjust.pl $option fem3d/netcdf.f
diff fem3d/netcdf.f.bak fem3d/netcdf.f

