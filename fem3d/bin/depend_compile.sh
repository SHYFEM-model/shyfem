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
# compiles fortran program to see if standalone
#
#-----------------------------------------------

[ $# -eq 0 ] && exit 0

echo "\tend" > ggg.f
cat $* >> ggg.f

gfortran ggg.f 2>&1 >/dev/null | grep "undefined reference"

