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
# for server fluxus
#
# this file must be run through check_server.sh
#
#--------------------------------------------------------------

this_server=fluxus

ResetDirs()
{
  echo "resetting server variables on $this_server"
  module purge
}

SetDirs_intel()
{
  [ "$SERVER_SHYFEM_COMPILER" = "INTEL" ] && return
  echo "setting server variables on $this_server for compiler intel"

  basedir=/usr/local/intel

  export SERVER_SHYFEM_HOSTNAME=$this_server
  export SERVER_SHYFEM_COMPILER=INTEL

  export NETCDF_SHYFEM_LIBCDIR=$basedir
  export NETCDF_SHYFEM_LIBFDIR=$basedir
  export NETCDF_SHYFEM_INCDIR=$basedir
  export NETCDF_SHYFEM_MODDIR=$basedir

  AddLibrary $NETCDF_SHYFEM_LIBCDIR lib
  AddLibrary $NETCDF_SHYFEM_LIBFDIR lib

  module load intel-openmpi
}

SetDirs_gfortran()
{
  [ "$SERVER_SHYFEM_COMPILER" = "GNU_GFORTRAN" ] && return
  echo "setting server variables on $this_server for compiler gfortran"

  basedir=/usr/lib64

  export SERVER_SHYFEM_HOSTNAME=$this_server
  export SERVER_SHYFEM_COMPILER=GNU_GFORTRAN

  export NETCDF_SHYFEM_LIBCDIR=$basedir
  export NETCDF_SHYFEM_LIBFDIR=$basedir
  export NETCDF_SHYFEM_INCDIR=$basedir
  export NETCDF_SHYFEM_MODDIR=$basedir/gfortran/modules

  AddLibrary $NETCDF_SHYFEM_LIBCDIR lib
  AddLibrary $NETCDF_SHYFEM_LIBFDIR lib

  module load openmpi-x86_64
}

#--------------------------------------------------------------

