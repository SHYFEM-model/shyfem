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
# skeleton for bash scripts
#
#-------------------------------------------------------------

FEMDIR=~/fem
tests="test1 hakata taranto_h taranto_z"

files="Makefile README"
filesbin="checkts.pl compile.sh ggunzip progrg.sh"

error=0

#-------------------------------------------------------------

Check()
{
    file1=$1
    file2=$2

    cmp --quiet $file1 $file2
    status=$?
    if [ $status -eq 0 ]; then
      status="ok "
      if [ $verbose = "YES" ]; then
        echo "status: $status   file: $file2   test: $test"
      fi
    else
      status="***"
      error=1
      echo "status: $status   file: $file2   test: $test"
    fi
}

#-------------------------------------------------------------

HandleOptions()
{
  dist="NO"
  check="NO"
  verbose="NO"

  while [ -n "$1" ]
  do
    case $1 in
        -dist)          dist="YES";;
        -check)         check="YES";;
        -verbose)       verbose="YES";;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
    esac
    shift
  done

  file=$1
}

Usage()
{
  echo "Usage: diff.sh [-h|-help] [-options] [file]"
}

FullUsage()
{
  Usage

  echo "  options: (one of -dist or -check must be given)"
  echo "     -dist        distribute file (file name must be given)"
  echo "     -check       check differences of files"
  echo "     -verbose     be verbose about checks"
}

ErrorOption()
{
  echo "No such option: $1"
}

#-------------------------------------------------------------

Diff()
{
  subdir=.
  for file in $files
  do
    echo "checking file $subdir/$file"
    for test in $tests
    do
      base="../$test"
      Check $base/$file $file
    done
  done

  subdir=bin
  for file in $filesbin
  do
    echo "checking file $subdir/$file"
    for test in $tests
    do
      base="../$test/$subdir"
      Check $base/$file $subdir/$file
    done
  done

  if [ $error -ne 0 ]; then
    echo "*** some files are not up to date"
  else
    echo "All files are updated"
  fi
}

Dist()
{
  file=$1

  if [ -z "$file" ]; then
    echo "option -dist needs also file name"
    FullUsage
    exit 1
  elif [ ! -f "$file" ]; then
    echo "No such file: $file"
    exit 1
  fi

  echo "distributing file $file"

  for test in $tests
  do
    base="../$test"
    Check $base/$file $file
    if [ "$status" = "***" ]; then
      echo "copying file $file to $base"
      cp -f $file $base/$file
    fi
  done
}

#-------------------------------------------------------------

HandleOptions $*

if [ "$dist" = "YES" ]; then
  Dist $file
elif [ "$check" = "YES" ]; then
  Diff
else
  Usage
fi

#-------------------------------------------------------------

