#!/bin/sh
#
# extracts background grid from file and meshes lines if necessary
#
#------------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
BINDIR=$FEMDIR/fembin

#exgrd=$BINDIR/exgrd
exgrd=$FEMDIR/mesh/exgrd
backb=$BINDIR/background.pl
modify=$BINDIR/grd_modify.pl

#---------------------------------------------------

backtype=-1
if [ "$1" = "-g" ]; then
  shift
  backtype=$1
  shift
fi

if [ $# -le 0 ]; then
  echo "Usage: background.sh [-g #] file(s)"
  echo "   -g   identifies type of background items"
  exit 0
fi

#--------------- extract background items from file

if [ -n "$backtype" ]; then
  $exgrd -t $backtype -T $backtype -elS $*
else
  $exgrd $*
fi
mv new.grd back_items.grd

#--------------- get background elements

$exgrd -eLS back_items.grd
mv new.grd back_elems.grd

#--------------- mesh lines

$backb back_items.grd

#--------------- combine meshed line background grid

$exgrd -u gtmp*.grd
if [ $backtype -ge 0 ]; then
  $modify -e -type=$backtype new.grd
  mv modify.grd back_lines.grd
else
  mv new.grd back_lines.grd
fi

#--------------- combine background grid from elements and lines (only elements)

$exgrd -eLS back_elems.grd back_lines.grd
mv new.grd background.grd

#--------------- write final message and clean up

rm -f M_*.grd ltmp*.grd gtmp*.grd

echo "background grid has been written to background.grd"

#--------------- end of routine

