#!/bin/sh

femdir=$FEMDIR

program=$femdir/femersem/fem_extra/bfm_change.pl

option=$1

#echo "$femdir  $program  $option"

files=`ls *.F90 *.f90 *.h    2>/dev/null`

for file in $files
do
  echo "checking file.... $file"
  chmod -x $file
  $program $file > tmp.tmp
  if [ $? -ne 0 ]; then
    if [ "$option" = "-change" ]; then
      echo "*** changing..."
      mv -f tmp.tmp $file
    else
      echo "*** only checking..."
      rm -f tmp.tmp
    fi
  else
    rm -f tmp.tmp
  fi
done

