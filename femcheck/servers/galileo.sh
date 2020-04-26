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
# for galileo server
#
# this file must be sourced (use command "source" or ".")
#
# Usage: "source galileo.sh {-show|-reset|-gfortran|-intel}"
#
#--------------------------------------------------------------

server=galileo

(return 0 2>/dev/null) && sourced=1 || sourced=0

#--------------------------------------------------------------

CheckSourced()
{
  if [ $sourced -eq 0 ]; then
    echo "this script must be sourced... please run:"
    echo "  source $0 $*"
    exit 1
  fi
}

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
  echo "SERVER_SHYFEM_HOSTNAME = $SERVER_SHYFEM_HOSTNAME"
  echo "SERVER_SHYFEM_COMPILER = $SERVER_SHYFEM_COMPILER"
  echo "NETCDF_SHYFEM_LIBCDIR = $NETCDF_SHYFEM_LIBCDIR"
  echo "NETCDF_SHYFEM_LIBFDIR = $NETCDF_SHYFEM_LIBFDIR"
  echo "NETCDF_SHYFEM_INCDIR = $NETCDF_SHYFEM_INCDIR"
  echo "NETCDF_SHYFEM_MODDIR = $NETCDF_SHYFEM_MODDIR"
  echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
}

ResetDirs()
{
  echo "resetting server variables on $server"
  export SERVER_SHYFEM_HOSTNAME=
  export SERVER_SHYFEM_COMPILER=
  export NETCDF_SHYFEM_LIBCDIR=
  export NETCDF_SHYFEM_LIBFDIR=
  export NETCDF_SHYFEM_INCDIR=
  export NETCDF_SHYFEM_MODDIR=
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ORIG
  module purge
}

SetDirs_intel()
{
  basedir=/cineca/prod/opt/libraries

  [ "$SERVER_SHYFEM_COMPILER" = "INTEL" ] && return

  ResetDirs
  echo "setting server variables on $server for compiler intel"

  export SERVER_SHYFEM_HOSTNAME=$server
  export SERVER_SHYFEM_COMPILER=INTEL

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

SetDirs_gfortran()
{
  basedir=/cineca/prod/opt/libraries

  [ "$SERVER_SHYFEM_COMPILER" = "GNU_GFORTRAN" ] && return

  ResetDirs
  echo "setting server variables on $server for compiler gfortran"

  export SERVER_SHYFEM_HOSTNAME=$server
  export SERVER_SHYFEM_COMPILER=GNU_GFORTRAN

  export NETCDF_SHYFEM_LIBCDIR=$basedir/netcdf/4.4.1/gnu--6.1.0
  export NETCDF_SHYFEM_LIBFDIR=$basedir/netcdff/4.4.4/gnu--6.1.0
  export NETCDF_SHYFEM_INCDIR=$basedir/netcdff/4.4.4/gnu--6.1.0
  export NETCDF_SHYFEM_MODDIR=$basedir/netcdff/4.4.4/gnu--6.1.0

  AddLibrary $NETCDF_SHYFEM_LIBCDIR lib
  AddLibrary $NETCDF_SHYFEM_LIBFDIR lib
}

#--------------------------------------------------------------

if [ "$1" = "-show" ]; then
  ShowDirs
elif [ "$1" = "-reset" ]; then
  CheckSourced $*
  ResetDirs
elif [ "$1" = "-gfortran" ]; then
  CheckSourced $*
  SetDirs_gfortran
elif [ "$1" = "-intel" ]; then
  CheckSourced $*
  SetDirs_intel
else
  echo "no command or unknown command given: $1"
fi

#--------------------------------------------------------------

