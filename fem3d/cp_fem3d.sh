#!/bin/sh

if [ $# -eq 0 ]; then
  echo "Usage: cp_fem3d.sh new-dir"
  exit 1
fi

dir=$1

if [ -d $dir ]; then
  echo "directory exists... cannot copy: $dir"
  exit 1
fi

mkdir $dir
cp -a * $dir

