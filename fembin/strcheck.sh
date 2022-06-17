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
# checks STR file and its files
#
#-----------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
fembin=$FEMDIR/fembin
fem3d=$FEMDIR/fem3d

if [ $# -ne 1 ]; then
  echo "Usage: strcheck.sh str-file"
  exit 0
fi

str=$1

#-----------------------------------------------------------------

ElabFile()
{
  first=$( grep 'first time record' tmp2.tmp | sed -e 's/.*: //' )
  last=$( grep 'last time record' tmp2.tmp | sed -e 's/.*:  //' )
  dt=$( grep 'regular time step' tmp2.tmp | sed -e 's/.*:  * //' )

  [ -z "$dt" ] && dt="irreg"

  echo "$first - $last - $dt - $file"
}

ElabSTR()
{
  itanf=$( echo $1 | sed -E "s/'//g" )
  itend=$( echo $2 | sed -E "s/'//g" )
  idt=$( echo $3 | sed -E "s/'//g" )
  
  echo "$itanf - $itend - $idt - STR"
}

#-----------------------------------------------------------------

itanf=`$fembin/strparse.pl -quiet -value=itanf $str`
itend=`$fembin/strparse.pl -quiet -value=itend $str`
date=`$fembin/strparse.pl -quiet -value=date $str`
idt=`$fembin/strparse.pl -quiet -value=idt $str`
[ "$date" = "" ] && date=0

ElabSTR $itanf $itend $idt

$fembin/strparse.pl -quiet -files $str > tmp.tmp

while read line
do
  where=$( echo $line | sed -e 's/ :.*//' )
  what=$( echo $line | sed -e 's/^.* : *//' | sed -e 's/ = .*//' )
  file=$( echo $line | sed -e 's/.*= //' | sed -E "s/'//g" )

  #echo "$where - $what - $file"
  [ $what = "grid" ] && continue
  [ $what = "gotmpa" ] && continue

  $fem3d/shyelab -quiet $file > tmp2.tmp

  ElabFile

  #echo "checking file $file"
  #$fembin/strcheck.pl $itanf $itend $date tmp2.tmp

done < tmp.tmp

rm -f tmp*.tmp

