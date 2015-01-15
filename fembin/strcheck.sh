#!/bin/sh
#
# checks STR file and its files
#
#-----------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
fembin=$FEMDIR/fembin

if [ $# -ne 1 ]; then
  echo "Usage: strcheck.sh str-file"
  exit 0
fi

str=$1

itanf=`strparse.pl -quiet -value=itanf $str`
itend=`strparse.pl -quiet -value=itend $str`
date=`strparse.pl -quiet -value=date $str`
[ "$date" = "" ] && date=0

echo "STR parameters: $itanf $itend $date"

strparse.pl -quiet -files $str > tmp.tmp
tail -n +2 tmp.tmp > tmp1.tmp		# get rid of line with grid (first)

while read line
do
  new1=`echo $line | sed -e 's/.*= //'`
  new2=`echo $new1 | sed -E "s/'//g"`
  file=$new2

  $fembin/fileinf  $file > tmp2.tmp

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

