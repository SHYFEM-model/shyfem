#!/bin/sh
#
#---------------------------------------------------

CheckDir()
{
  if [ ! -d $1 ]; then
    echo "No such dir: $1"
    exit 1
  fi
}

#---------------------------------------------------

undo="NO"
if [ "$1" = "-u" ]; then
  undo="YES"
  shift
fi

if [ $# -eq 0 ]; then
  echo "Usage: ./prep_dirs.sh [-u] number"
  exit 1
fi

number=$1
echo "linking to $number"

CheckDir fem3d_orig
CheckDir femplot_orig
CheckDir fem3d_$number
CheckDir femplot_$number

if [ $undo = "YES" ]; then
  rm fem3d
  rm femplot
  ln -s fem3d_orig fem3d
  ln -s femplot_orig femplot
else
  rm fem3d
  rm femplot
  ln -s fem3d_$number fem3d
  ln -s femplot_$number femplot
fi

