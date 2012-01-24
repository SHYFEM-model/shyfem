#!/bin/sh
#
# shyfem_version.sh
#
# extracts version of SHYFEM
#
# called without argument returns version of actual SHYFEM distribution
# called with directory returns version of this fem directory
#
# Usage: shyfem_dir.sh [femdir]
#
#------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

if [ $# -eq 1 ]; then
  dir=$1
else
  dir=$FEMDIR
fi

version_file=$dir/VERSION
if [ ! -f $version_file ]; then
  version_file=$dir/fem3d/VERSION		#try this for old versions
fi

if [ ! -f $version_file ]; then
  echo "Cannot find file: VERSION ...aborting" 1>&2
  echo "unknown"
  exit 1
fi

version=`$FEMDIR/fembin/shyfem_version.pl -tag_extra $version_file`

if [ $? -ne 0 ]; then
  echo "file $version_file does not contain SHYFEM version...aborting" 1>&2
  echo "$line" 1>&2
  echo "unknown"
  exit 1
fi

echo $version
  
