#!/bin/sh
#
# shyfem_dir.sh
#
# resets path and and environment variables
#
# called without argument returns info on actual SHYFEM distribution
# called with argument sets the new shyfem directory
#
# in order to work this script must be called through shyfemdir
# or must be sourced: ". shyfem_dir.sh" or "source shyfem_dir.sh"
#
# Usage: 
#    . shyfem_dir.sh [femdir]
#    source shyfem_dir.sh [femdir]
#    shyfemdir [femdir]
#
#------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
fembin=$FEMDIR/fembin

FEMDIR_INSTALL=${SHYFEM_INSTALL:=$HOME/shyfem}
fembin_install=$FEMDIR_INSTALL/fembin

# command line options ----------------------------

write="normal"
if [ "$1" = "-quiet" ]; then
  write="quiet"
  shift
elif [ "$1" = "-debug" ]; then
  write="debug"
  shift
fi

# change path and environment variables ---------------

if [ -n "$1" ]; then

  dir=`readlink -f $1`	# get full path name

  #echo "debug message: using dir as $dir"

  version=`$fembin_install/shyfem_version.sh $dir`
  if [ -z "$version" -o "$version" = "unknown" ]; then
    echo "cannot get version for $dir ... aborting" 1>&2
    exit 1
  fi

  export SHYFEMDIR=$dir
  FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
  fembin=$FEMDIR/fembin

  path=`$fembin_install/shyfem_path.pl $PATH`
  export PATH=$path:$fembin:$HOME/shyfem/fembin

fi

# show path and environment variables ---------------

version=`$fembin_install/shyfem_version.sh`

if [ $write != "quiet" ]; then

  echo "SHYFEM version: $version"
  echo "SHYFEM directory: $FEMDIR"
  echo "SHYFEM install directory: $FEMDIR_INSTALL"

  if [ $write = "debug" ]; then
    echo "SHYFEM path: $PATH"
  fi
fi

# end of routine ----------------------------------------

