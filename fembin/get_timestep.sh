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
# gets time-step from .inf file (deleting synconization times)
#
#----------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
fembin=$FEMDIR/fembin

show="NO"
if [ "$1" = "-show" ]; then
  shift
  show="YES"
fi

if [ $# -lt 1 ]; then
  echo "Usage: get_timestep.sh [-show] simulation"
  exit 1
fi

file=$1

name=`basename $file .inf`
file=$name.inf

$fembin/getkey.pl  timestep  $file  > aaa.tmp
$fembin/clean_sync_dt.pl aaa.tmp > bbb.tmp
mv bbb.tmp timestep.txt

gp -u 1:5 timestep.txt
cp out.ps timestep.ps

echo "data is in timestep.txt and plot in timestep.ps"

[ $show = "NO" ] && exit 0

gv timestep.ps

#----------------------------------------------------------

