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
# shell to allow for global change of area codes
#
#-----------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
BINDIR=$FEMDIR/fembin

#-----------------------------------------------------------

bfile=area_bound
tfile=area_tmp
afile=area_elem
efile=area_type

#-----------------------------------------------------------

Usage()
{
  echo "Usage: grd_area.sh [-h] {-b|-c|-a} grd-file [elem-file]"
  exit 1
}

FullUsage()
{
  echo ""
  echo "Usage: grd_area.sh {-b|-c|-a} grd-file [elem-file]"
  echo ""
  echo "  Allows for global change of area codes"
  echo ""
  echo "  Create boundary file:"
  echo "     grd_area.sh -b grd-file"
  echo "  Edit macro element file:"
  echo "     grd_area.sh -c grd-file [elem-file]"
  echo "     (use with elem-file if already existing, else create new one)"
  echo "  Change element type:"
  echo "     grd_area.sh -a grd-file elem-file"
  echo ""
  exit 1
}

#-------------------------------------------------------------

if [ $# -eq 0 ]; then
  Usage
elif [ $1 = "-h" ]; then
  FullUsage
fi

if [ "$1" = "-b" ]; then
  shift
  [ $# -eq 0 ] && FullUsage

  $BINDIR/makebline.pl $1 > $bfile.grd
  echo "boundary line written to file $bfile.grd"

  exit 0
fi
  
if [ "$1" = "-c" ]; then
  shift
  [ $# -eq 0 ] && FullUsage

  $BINDIR/grid -t888 -O$tfile $*
  $BINDIR/grd_area.pl -e $tfile
  mv modify.grd $afile.grd
  echo "area elements written to file $afile.grd"
  echo "  (please rename if you want to change it again)"

  exit 0
fi

if [ "$1" = "-a" ]; then
  shift
  [ $# -eq 0 ] && FullUsage

  $BINDIR/grd_area.pl $1 $2
  mv modify.grd $efile.grd
  echo "grd with new area codes written to file $efile.grd"

  exit 0
fi





