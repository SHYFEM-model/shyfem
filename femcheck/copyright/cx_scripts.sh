#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

Usage()
{
  echo "Usage: cx_scripts.sh [-cx] dir"
}

cx="NO"
if [ "$1" = "-cx" ]; then
  cx="YES"
  shift
fi

if [ $# -eq 0 ]; then
  Usage
  exit 1
fi

if [ ! -d $1 ]; then
  echo "not a directory: $1"
  Usage
  exit 1
fi

cd $1

files_sh=$( findf '*.sh' )
files_pl=$( findf '*.pl' )
files="$files_sh $files_pl"
n=0

for file in $files
do
  if [ ! -x $file ]; then
    if [ $cx = "YES" ]; then
      echo "making file executable: $file"
      cx $file
    else
      echo "file is not executable: $file"
      n=$(( n+1 ))
    fi
  fi
done

if [ $n -gt 0 ]; then
  echo "$n non executable files found"
  echo "use -cx to make files executable"
fi

