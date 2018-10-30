#!/bin/sh
#
# checks if copyright has been inserted
#
#------------------------------------------

nfiles=0
files=$( ls *.[fhc] *.tex *.f90 2> /dev/null )

for file in $files
do
  grep "Copyright" $file | grep "Umgiesser" | grep "2018" > /dev/null
  status=$?
  #echo "$status  $file"
  if [ $status -eq 0 ]; then
    echo "ok   $file"
  else
    echo "***  $file"
    nfiles=$(($nfiles+1))
  fi
done

if [ $nfiles -gt 0 ]; then
  echo "$nfiles without copyright found"
fi

