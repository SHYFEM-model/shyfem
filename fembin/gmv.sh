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
# returns major version of gfortran
#
#-----------------------------------------------------------

prog=$( which gfortran )

if [ -n "$prog" ]; then
  version=$( $prog -v 2>&1 | tail -1 | cut -d " " -f 3 )
  if [ -z "$version" ]; then
    mversion=""
  else
    mversion=$( echo $version | sed -e 's/\..*//' )
  fi
fi

if [ "$1" = "-info" ]; then
    echo "gfortran  version=$version  major_version=$mversion"
else
    echo "$mversion"
fi

#if [ -n "$prog" ]; then
#  $prog --version
#               | head -1               \
#                | sed -e 's/.*) *//'   \
#               | sed -e 's/\..*//'
#fi

