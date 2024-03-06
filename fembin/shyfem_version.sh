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
# shyfem_version.sh
#
# extracts version of SHYFEM
#
# called without argument returns version of actual SHYFEM distribution
# called with directory returns version of this fem directory
# we must be in the base directory of shyfem
# we cannot use exit because this script might be called  through sourcing
#
# Usage: shyfem_dir.sh [femdir]
#
#------------------------------------------------------------------------

dist=shyfem

FEMDIR=${SHYFEMDIR:=$HOME/$dist}
script=$(realpath $0)
FEMBIN=$(dirname $script)

if [ $# -eq 1 ]; then
  dir=$1
else
  dir=$FEMDIR
fi

version_file=$dir/VERSION
version="unknown"

if [ ! -f $version_file ]; then
  echo "Cannot find file: VERSION ...aborting" 1>&2
  exit 1
elif [ ! -x $FEMBIN/shyfem_version.pl ]; then
  echo "Cannot find $FEMBIN/shyfem_version.pl ...aborting" 1>&2
  exit 3
else
  version=`$FEMBIN/shyfem_version.pl -tag_extra $version_file`
  if [ $? -ne 0 ]; then
    echo "file $version_file does not contain SHYFEM version ...aborting" 1>&2
    exit 5
  fi
fi

echo $version
  
#------------------------------------------------------------------------

