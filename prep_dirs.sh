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
  echo "linking to original dir"
  [ -L fem3d ] && rm fem3d && ln -s fem3d_orig fem3d
  [ -L femplot ] && rm femplot && ln -s femplot_orig femplot
  exit 0
elif [ "$1" = "-o" ]; then	#make original
  echo "making original situation"
  [ -L fem3d ] && rm fem3d 
  [ -L femplot ] && rm femplot 
  mv fem3d_orig fem3d
  mv femplot_orig femplot
  exit 0
elif [ $# -eq 0 ]; then
  echo "Usage: ./prep_dirs.sh [-u|-o] number"
  exit 1
fi

number=$1
echo "linking to $number"

CheckDir fem3d_$number
CheckDir femplot_$number

if [ ! -L fem3d ]; then
  echo "fem3d is not a symbolic link... moving to fem3d_orig"
  mv fem3d fem3d_orig
fi
if [ ! -L femplot ]; then
  echo "femplot is not a symbolic link... moving to femplot_orig"
  mv femplot femplot_orig
fi

CheckDir fem3d_orig
CheckDir femplot_orig

[ -L fem3d ] && rm fem3d
[ -L femplot ] && rm femplot
ln -s fem3d_$number fem3d
ln -s femplot_$number femplot

