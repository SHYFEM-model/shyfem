#!/bin/sh
#
# gets time-step from .inf file (deleting synconization times)
#
#----------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
fembin=$FEMDIR/fembin

if [ $# -lt 1 ]; then
  echo "Usage: get_timestep.sh simulation [idtsyn]"
  exit 1
fi

idtsyn=0
file=$1
[ $# -ge 2 ] && idtsyn=$2

name=`basename $file .inf`
file=$name.inf

if [ $idtsyn -eq 0 ]; then
  getkey.pl  set_timestep  $file
else
  echo "$idtsyn" > aaa.tmp
  getkey.pl  set_timestep  $file  >> aaa.tmp
  $fembin/progs/clean_time < aaa.tmp
fi

