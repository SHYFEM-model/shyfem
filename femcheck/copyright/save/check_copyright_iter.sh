#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#----------------------------------------------------------

progdir=$HOME/shyfem/femcheck/copyright
prog="$progdir/check_copyright.sh"
options="-no -silent"

Usage()
{
  echo "Usage: check_copyright_iter.sh [-options] dir"
  exit 1
}

FullUsage()
{
  Usage
}
#----------------------------------------------------------

while [ -n "$1" ]
do
   case $1 in
        -h|-help)       FullUsage; exit 0;;
        -*)             options="$options $1";;
        *)              break;;
   esac
   shift
done

if [ $# -le 0 ]; then
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

