#!/bin/sh

#--------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#--------------------------------------------------------------------------
#
# plots grid in ncfile
#
#--------------------------------------------

xdim=x
ydim=y
xvar=lon
yvar=lat

file=file.nc

#--------------------------------------------

ncdump -v $xvar $file > xheader.txt
ncdump -v $yvar $file > yheader.txt

./ncgrid.pl $xvar xheader.txt > xdata.txt
./ncgrid.pl $yvar yheader.txt > ydata.txt

./ncgrid

