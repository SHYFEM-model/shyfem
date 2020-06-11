#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
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
  echo "Usage: cmv.sh [-info|-quiet] {gfortran|intel}"
}

info="NO"
if [ "$1" = "-info" ]; then
  info="YES"
  shift
fi

quiet="NO"
if [ "$1" = "-quiet" ]; then
  quiet="YES"
  shift
fi

if [ $# -eq 0 ]; then
  Usage
  exit 1
fi

#-----------------------------------------------------------

version=unknown
mversion=0
exitcode=0

compiler=$1
compexe=$compiler
[ $compiler = intel ] && compexe=ifort

prog=$( which $compexe 2> /dev/null )

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
    [ $quiet = NO ] && echo "compiler not supported: $compiler" >> /dev/stderr
    exitcode=1
  fi
else
  [ $quiet = NO ] && echo "cannot find compiler: $compiler" >> /dev/stderr
  exitcode=1
fi

[[ $mversion != +([0-9]) ]] && mversion=0   #handle error in running compiler

if [ $info = "YES" ]; then
    echo "$compiler  version=$version  major_version=$mversion"
else
    echo "$mversion"
fi

exit $exitcode

#-----------------------------------------------------------

