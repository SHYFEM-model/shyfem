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
# returns major version of compilers (gfortran and intel)
#
#-----------------------------------------------------------

Usage()
{
  echo "Usage: cmv.sh [-info] {gfortran|intel}"
}

info="NO"
if [ "$1" = "-info" ]; then
  info="YES"
  shift
fi

if [ $# -eq 0 ]; then
  Usage
  exit 1
fi

#-----------------------------------------------------------

version=unknown
mversion=0

compiler=$1
compexe=$compiler
[ $compiler = intel ] && compexe=ifort

prog=$( which $compexe )

if [ -n "$prog" ]; then
  if [ $compiler = gfortran ]; then
    version=$( $prog -v 2>&1 | tail -1 | cut -d " " -f 3 )
    if [ -n "$version" ]; then
      mversion=$( echo $version | sed -e 's/\..*//' )
    fi
  elif [ $compiler = intel ]; then
    version=$( $prog -v 2>&1 | sed -e 's/.*version *//' )
    if [ -n "$version" ]; then
      mversion=$( echo $version | sed -e 's/\..*//' )
    fi
  else
    echo "compiler not supported: $compiler"
    Usage
    exit 1
  fi
else
  echo "cannot find compiler: $compiler"
  Usage
  exit 1
fi

if [ $info = "YES" ]; then
    echo "$compiler  version=$version  major_version=$mversion"
else
    echo "$mversion"
fi

#-----------------------------------------------------------

