#!/bin/sh
#
# checks if copyright has been inserted
#
#--------------------------------------------------------

progdir=$HOME/shyfem/femcheck/copyright
prog="$progdir/check_copyright.pl -quiet"

#--------------------------------------------------------

Usage()
{
  echo "Usage: check_copyright.sh [-ok|-no|-quiet] dir"
  exit 0
}

#--------------------------------------------------------

ok="YES"
no="YES"

if [ "$1" = '-ok' ]; then
  no="NO"
  shift
elif [ "$1" = '-no' ]; then
  ok="NO"
  shift
elif [ "$1" = '-quiet' ]; then
  ok="NO"
  no="NO"
  shift
fi

if [ $# -eq 0 ]; then
  Usage
fi
cd $1

#--------------------------------------------------------
  
nofiles=0
files=$( ls *.[fhc] *.tex *.f90 *.F90 2> /dev/null )
files=$( ls *.[f] 2> /dev/null )

for file in $files
do
  $prog $file
  status=$?
  #echo "$status  $file"
  if [ $status -eq 0 ]; then
    [ $ok = "YES" ] && echo "ok   $file"
    okfiles=$(($okfiles+1))
  else
    [ $no = "YES" ] && echo "***  $file"
    nofiles=$(($nofiles+1))
  fi
done

echo "$okfiles with copyright found"
echo "$nofiles without copyright found"

