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
# shows given coordinates on top of a grd file
#
#------------------------------------------------------------------------

Usage()
{
  echo "Usage: grd_show_coords.sh grd-file coord-file"
}

MakeGrdCoords()
{

  [ -f $tmpfile ] && rm -f $tmpfile

  echo "using $maxnode for minimum node number"

  j=$maxnode
  while read -r line
  do
    j=$(( j + 1 ))
    SetLatLon $line
    [ -z "$lon" ] && continue
    echo "extracting: $line - $lon - $lat"
    echo "1 $j 3 $lon $lat" >> $tmpfile
  done < "$coordfile"
}

SetLatLon()
{
  lon=$1
  lat=$2
}

ShowGrd()
{
  grid -fT $grdfile $tmpfile
}

MakeMaxNode()
{
  max=$( exgrd -i $grdfile 2>/dev/null \
		| grep "Min/Max Node" | sed -e 's/.* //' )
  l=$( log 10 $max )
  l=$(( l + 1 ))

  x=1
  while((l))
  do 
    let l-=1 x*=10
  done

  maxnode=$(( x - max - 1 ))
}

log()
{ 
  # returns floor[log(x)] for x integer
  # call as "log 777" (base 2) or "log 10 777" (base 10)
  # http://phodd.net/gnu-bc/bcfaq.html#bashlog

  local x=$1 n=2 l=-1;
  if [ "$2" != "" ]; then 
    n=$x;x=$2
  fi

  while((x))
  do 
    let l+=1 x/=n
  done

  echo $l;
} 

#--------------------------------------------------------------

grdfile=$1
coordfile=$2
tmpfile=coords.tmp.grd

if [ $# -lt 2 ]; then
  Usage; exit 1
elif [ ! -f $grdfile ]; then
  echo "*** no such file: $grdfile"; exit 3
elif [ ! -f $coordfile ]; then
  echo "*** no such file: $coordfile"; exit 3
fi

#--------------------------------------------------------------

MakeMaxNode
MakeGrdCoords
ShowGrd

#--------------------------------------------------------------

