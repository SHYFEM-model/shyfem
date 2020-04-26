#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# for galileo server - gfortran
#
# this file must be sourced (use command "source" or ".")
#
#--------------------------------------------------------------

AddLibrary()
{
  libdir=$1
  pre=$2

  [ -n "$pre" ] && libdir=$libdir/$pre

  [ -z "$LD_LIBRARY_PATH_ORIG" ] && LD_LIBRARY_PATH_ORIG=$LD_LIBRARY_PATH

  echo $LD_LIBRARY_PATH | grep $libdir > /dev/null
  if [ $? -ne 0 ]; then
    echo "adding $libdir to library path"
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$libdir
  fi
}

ShowDirs()
{
  echo "NETCDF_SHYFEM_LIBCDIR = $NETCDF_SHYFEM_LIBCDIR"
  echo "NETCDF_SHYFEM_LIBFDIR = $NETCDF_SHYFEM_LIBFDIR"
  echo "NETCDF_SHYFEM_INCDIR = $NETCDF_SHYFEM_INCDIR"
  echo "NETCDF_SHYFEM_MODDIR = $NETCDF_SHYFEM_MODDIR"
  echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
}

ResetDirs()
{
  export NETCDF_SHYFEM_LIBCDIR=
  export NETCDF_SHYFEM_LIBFDIR=
  export NETCDF_SHYFEM_INCDIR=
  export NETCDF_SHYFEM_MODDIR=
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ORIG
  module purge
}

SetDirs()
{
  basedir=/cineca/prod/opt/libraries

  export NETCDF_SHYFEM_LIBCDIR=$basedir/netcdf/4.6.1/intel--pe-xe-2018--binary
  export NETCDF_SHYFEM_LIBFDIR=$basedir/netcdff/4.4.4/intel--pe-xe-2018--binary
  export NETCDF_SHYFEM_INCDIR=$basedir/netcdff/4.4.4/intel--pe-xe-2018--binary
  export NETCDF_SHYFEM_MODDIR=$basedir/netcdff/4.4.4/intel--pe-xe-2018--binary

  AddLibrary $NETCDF_SHYFEM_LIBCDIR lib
  AddLibrary $NETCDF_SHYFEM_LIBFDIR lib

  module purge
  module load profile/advanced
  module load autoload
  module load intel/pe-xe-2018--binary
  module load intelmpi/2018--binary
  module load netcdf/4.6.1--intel--pe-xe-2018--binary
  module load netcdff/4.4.4--intel--pe-xe-2018--binary
}

#--------------------------------------------------------------

if [ "$1" = "-show" ]; then
  ShowDirs
elif [ "$1" = "-reset" ]; then
  ResetDirs
else
  SetDirs
fi

#--------------------------------------------------------------

