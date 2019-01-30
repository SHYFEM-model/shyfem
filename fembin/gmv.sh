#!/bin/sh
#
# returns major version of gfortran
#
#-----------------------------------------------------------

prog=$( which gfortran )

if [ -n "$prog" ]; then
  version=$( $prog -v 2>&1 | tail -1 | cut -d " " -f 3 )
  if [ -z "$version" ]; then
    #echo "*** cannot determine gfortran version... (gmv.sh)"
    #$prog -v
    #exit 1
    echo ""
  else
    mversion=$( echo $version | sed -e 's/\..*//' )
    #echo "version=$version  major_version=$mversion"
    echo "$mversion"
  fi
fi

#if [ -n "$prog" ]; then
#  $prog --version
#               | head -1               \
#                | sed -e 's/.*) *//'   \
#               | sed -e 's/\..*//'
#fi

