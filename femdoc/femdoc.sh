#!/bin/sh
#
# creates *.tex documentation files from *.f files

FEMDIR=..
BINDIR=$FEMDIR/fem3d		#directory where fortran files reside

Proc()
{
  ./femdoc.pl $1
  if [ $? -ne 0 ]; then
    echo "*** Error in processing file $1"
    exit 1
  fi
}

#---------------------------------------------------------

Proc $BINDIR/subsys.f
Proc $BINDIR/subbnd.f
#Proc $BINDIR/subwin.f
Proc $BINDIR/submeteo2.f
Proc $BINDIR/subn35.f
Proc $BINDIR/subver.f
Proc $BINDIR/bio3d.f
Proc $BINDIR/sedi3d.f
Proc $BINDIR/subwaves.f

#---------------------------------------------------------

