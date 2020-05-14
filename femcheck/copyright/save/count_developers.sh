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
#---------------------------------------------------------

shydir="$HOME/shyfem"
copydir="$shydir/femcheck/copyright"
copyfile="$copydir/copy.log"
tmpfile="copy.tmp"

#---------------------------------------------------------

OldFiles()
{
  files=$*

  for file in $files
  do
    $copydir/count_developers.pl -old $file
  done
}

CheckRevision()
{
  files=$*

  for file in $files
  do
    echo "================================================= $file"
    $copydir/count_developers.pl -copy -revision $file
  done
}

ElabCopy()
{
  files=$*

  for file in $files
  do
    echo "================================================= $file"
    $copydir/count_developers.pl -copy $options $file
  done
}

SubstFiles()
{
  files=$*

  for file in $files
  do
    [ -f $tmpfile ] && rm -f $tmpfile
    echo "================================================= $file"
    $copydir/count_developers.pl -copy $file    > $tmpfile
    $status = $?
    if [ $status -ne 0 ]; then
      echo "*** error handling copyright... not substituting"
    elif [ -s $tmpfile ]; then
      cat $tmpfile
      echo "substituting copyright...";
      cp $file $file.old
      $copydir/count_developers.pl -subst $file.old    > $file
    else
      echo "keeping old copyright...";
    fi
  done

  [ -f $tmpfile ] && rm -f $tmpfile
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
    echo "$file"
    $copydir/count_developers.pl -new $file    >> $copyfile
  done

  cd $actdir
}

CheckSingle()
{
  CheckDir $*

  $copydir/count_developers.pl -count $copyfile
}

CheckAll()
{
  CheckDir $shydir/fem3d
  CheckDir $shydir/femplot
  CheckDir $shydir/femadj

  $copydir/count_developers.pl -count $copyfile
}

#---------------------------------------------------------

Usage()
{
  echo "Usage: count_developers.sh [-h|-help] [options] {dir|files}"
}

FullUsage()
{
  Usage
}

#---------------------------------------------------------


actdir=$( pwd )
[ -f $copyfile ] && rm -f $copyfile
what=""
options=""

while [ -n "$1" ]
do
   case $1 in
        -revision)      options="-revision";;
        -copy)          what="copy";;
        -subst)         what="subst";;
        -old)           what="old";;
        -check)         what="check";;
        -checkall)      what="checkall";;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

if [ $# -eq 0 ]; then
  Usage
elif [ -z "$what" ]; then
  echo "*** no action given... exiting"
  Usage
elif [ "$what"  = "copy" ]; then
  ElabCopy $*
elif [ "$what"  = "subst" ]; then
  SubstFiles $*
elif [ "$what"  = "old" ]; then
  OldFiles $*
elif [ "$what"  = "check" ]; then
  CheckSingle $1
elif [ "$what"  = "checkall" ]; then
  CheckAll $1
fi

#---------------------------------------------------------

