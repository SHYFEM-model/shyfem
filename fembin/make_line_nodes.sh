#!/bin/sh
#
##############################################################

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

BINDIR=$FEMDIR/fembin
PLOTDIR=$FEMDIR/femplot

prog=$PLOTDIR/line_nodes
CR=$BINDIR/CR
memory=$BINDIR/memory

##############################################################

Usage()
{
  exefile=`basename $0`
  echo "Usage: $exefile grd-line-file [grd-basin-file]"
  echo "   grd-line-file    grd file containing line"
  echo "   grd-basin-file   grd file containing basin (default from .memory)"
  exit 1
}

#--------------------------------------------------------------

if [ $# -le 0 ]; then
  Usage
elif [ $# -eq 1 ]; then
  linefile=`basename $1 .grd`
  basfile=`$memory -b`
elif [ $# -eq 2 ]; then
  linefile=`basename $1 .grd`
  basfile=`basename $2 .grd`
else
  Usage
fi

echo "file containing line: $linefile.grd"
echo "file containing basin: $basfile.grd"

if [ -f $basfile.bas ]; then
  echo "using available $basfile.bas"
elif [ -f $basfile.grd ]; then
  echo "creating $basfile.bas"
  vpgrd $basfile
else
  echo "no basin file available: $basfile ... aborting"
  exit 3
fi

$memory -b $basfile

#--------------------------------------------------------------

exgrd -lES $linefile
mv new.grd line_1.grd

grd2bnd.pl line_1.grd > line_1.bnd

$prog < line_1.bnd
mv line_nodes.grd line_2.grd
mv line_nodes.txt line_2.nod

echo "original line is in line_1.grd"
echo "transformed line is in line_2.grd"
echo "node list is in line_2.nod"

#--------------------------------------------------------------

