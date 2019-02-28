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
#---------------------------------------------------------

shydir="$HOME/shyfem"
copydir="$shydir/femcheck/copyright"
copyfile="$copydir/copy.log"

#---------------------------------------------------------

CheckAll()
{
  CheckDir $shydir/fem3d
  CheckDir $shydir/femplot
  CheckDir $shydir/femadj

  $copydir/count_developers.pl -count $copyfile
}

OldFiles()
{
  files=$*

  for file in $files
  do
    $copydir/count_developers.pl -old $file
  done
}

SubstFiles()
{
  files=$*

  for file in $files
  do
    [ -f tmp.tmp ] && rm -f tmp.tmp
    $copydir/count_developers.pl -copy $file    > tmp.tmp
    if [ -s tmp.tmp ]; then
      echo $file
      cat tmp.tmp
      #$copydir/count_developers.pl -subst $file    > new.tmp
    else
      echo $file
      echo "*** file has no revision log..."
    fi
  done
}

CheckDir()
{
  if [ -d $1 ]; then
    cd $1
    files=$( ls *.f *.f90 2> /dev/null )
  else
    files=$*
  fi

  echo "checking in $( pwd )"

  for file in $files
  do
    $copydir/count_developers.pl $file    >> $copyfile
  done

  cd $actdir
}

#---------------------------------------------------------

if [ "$1" = "-subst" ]; then
  shift
  SubstFiles $*
  exit 0
fi

if [ "$1" = "-old" ]; then
  shift
  OldFiles $*
  exit 0
fi

actdir=$( pwd )

[ -f $copyfile ] && rm -f $copyfile

if [ $# -eq 0 ]; then
  CheckAll
  exit 0
fi

dir=$1

CheckDir $dir
$copydir/count_developers.pl -count $copyfile

#---------------------------------------------------------

