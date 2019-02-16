#!/bin/bash

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

femdir=~/shyfem

bindir=$femdir/fem3d

pushd $bindir
make check_debug
popd

if [ $# -ne 2 ]; then
  echo "Usage: check_debug.sh fort.66.1 fort.66.2"
  exit 1
fi

[ -f debug_one.dat ] && rm -f debug_one.dat
[ -f debug_two.dat ] && rm -f debug_two.dat

ln -fs $1 debug_one.dat
ln -fs $2 debug_two.dat

echo "using files: $1 $2"

$bindir/check_debug 

