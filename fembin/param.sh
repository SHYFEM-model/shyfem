#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

dir=$1;
file=$2;

#echo "param.sh: $dir  $file";

if [ ! -f $dir/$file ]; then
  echo "No file $file in directory $dir"
  exit 0
fi

if [ ! -f $file ]; then
  echo "No file $file in current directory"
  exit 1
fi

cmp $dir/$file $file
status=$?

if [ $status -ne 0 ]; then
  echo "files $file differ ... copying from $dir"
  cp -f $dir/$file .
fi

