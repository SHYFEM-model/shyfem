#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# shows regular fem grid
#
#------------------------------------------------------------------------

Usage()
{
  echo "Usage: grd_show_regular.sh fem-file"
}

MakeRegular()
{
  #femelab -info $femfile
  reginfo=$( femelab -info $femfile | sed '1,/regpar:/d' | head -4 )
  echo $reginfo
  nxy=$( echo $reginfo | sed 's/.*ny://' | sed 's/x0.*//' )
  echo "nxy: $nxy"
  xy0=$( echo $reginfo | sed 's/.*y0://' | sed 's/x1.*//' )
  echo "xy0: $xy0"
  dxy=$( echo $reginfo | sed 's/.*dy://' )
  echo "dxy: $dxy"

  make_regular.pl $nxy $xy0 $dxy > $tmpfile

}

ShowGrd()
{
  grid -fT $grdfile $tmpfile
}

#--------------------------------------------------------------

femfile=$1
tmpfile=regular.tmp.grd

if [ $# -lt 1 ]; then
  Usage; exit 1
elif [ ! -f $femfile ]; then
  echo "*** no such file: $femfile"; exit 3
fi

#--------------------------------------------------------------

MakeRegular
ShowGrd

#--------------------------------------------------------------

