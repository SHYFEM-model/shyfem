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
# iterate through dirs and execute command
#
#------------------------------------------------------------------------

if [ $# -lt 1 ]; then
  echo "Usage: iterdir.sh [-dir] command dir(s)"
  exit 1
fi

elabdir="NO"
if [ "$1" = "-dir" ]; then
  elabdir="YES"
  shift
fi

command="$1"
shift

home=`pwd`
if [ $# -eq 0 ]; then
  dirs=.
else
  dirs=$*
fi

#echo "home = $home"
#echo "dirs = $dirs"
#echo "command: $command"

#------------------------------------------------------------------------

Iterate() {
  for d in *; do
    if [ -d "$d" ]; then
      (cd -- "$d" && Iterate)
    fi
    [ $elabdir = "NO" ] && $command
  done
  [ $elabdir = "YES" ] && $command
}

#------------------------------------------------------------------------

for dir in $dirs
do
  echo "starting directory: $dir"
  cd $dir
  Iterate
  cd $home
done

#------------------------------------------------------------------------

