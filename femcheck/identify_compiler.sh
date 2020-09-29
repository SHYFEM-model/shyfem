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
# try to identify compiler version for given executable
#
#------------------------------------------------

actdir=$( dirname $0 )
base=$actdir/..

program=$1
if [ $# -eq 0 ]; then
  program=$base/fem3d/shyfem
fi
#echo "using program $program"

# check intel
intel=$( strings $program | grep -i "intel fortran" | head -1 )

# check gfortran
gfortran=$( strings $program | grep "GFORTRAN" | head -1 )

# check pgi
pgi=$( strings $program | grep "NVFORTRAN" | head -1 )

#echo "intel: $intel"
#echo "gfortran: $gfortran"

if [ -n "$intel" ]; then
  out=$intel
elif [ -n "$gfortran" ]; then
  gcc=$( strings $program | grep "GCC" | head -1 )
  glibc=$( strings $program | grep "GLIBC" | head -1 )
  out="$gfortran  $gcc  $glibc"
elif [ -n "$pgi" ]; then
  out=$pgi
else
  echo "cannot identify compiler..."
fi

echo $out

