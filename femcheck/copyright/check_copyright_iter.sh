#!/bin/sh
#
#----------------------------------------------------------

progdir=$HOME/shyfem/femcheck/copyright
prog="$progdir/check_copyright.sh"
options="-no -silent"

Usage()
{
  echo "Usage: iterate.sh [-options] dir"
  exit 1
}

#----------------------------------------------------------

if [ $# -ne 1 ]; then
  Usage
fi

mdir=$1

if [ ! -d $mdir ]; then
  echo "no such directory: $mdir"
  Usage
fi

cd $mdir
act=$( pwd )
dirs=$( find -type d )

for dir in $dirs
do
  echo "--------------------------------"
  echo $dir
  echo "--------------------------------"
  cd $dir
  $prog $options .
  cd $act
done

