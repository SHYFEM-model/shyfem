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

itanf=`$fembin/strparse.pl -quiet -value=itanf $str`
itend=`$fembin/strparse.pl -quiet -value=itend $str`
date=`$fembin/strparse.pl -quiet -value=date $str`
[ "$date" = "" ] && date=0

echo "STR parameters: $itanf $itend $date"

$fembin/strparse.pl -quiet -files $str > tmp.tmp
tail -n +2 tmp.tmp > tmp1.tmp		# get rid of line with grid (first)

while read line
do
  new1=`echo $line | sed -e 's/.*= //'`
  new2=`echo $new1 | sed -E "s/'//g"`
  file=$new2

  $fem3d/fileinf  $file > tmp2.tmp

  echo "checking file $file"
  $fembin/strcheck.pl $itanf $itend $date tmp2.tmp

  status=$?
  if [ $status -eq 1 ]; then
    status="error"
  elif [ $status -eq 2 ]; then
    status="warning"
  else
    status="ok"
  fi

  line=`head -1 tmp2.tmp`
  echo "$line  ($status)"
done < tmp1.tmp

rm -f tmp*.tmp

