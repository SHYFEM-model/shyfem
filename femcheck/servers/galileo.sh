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
# for server galileo
#
# this file must be run through check_server.sh
#
#--------------------------------------------------------------

this_server=galileo

ResetDirs()
{
  echo "resetting server variables on $this_server"
  module purge
}

SetDirs_intel()
{
  [ "$SERVER_SHYFEM_COMPILER" = "INTEL" ] && return
  echo "setting server variables on $this_server for compiler intel"

  basedir=/cineca/prod/opt/libraries

  export SERVER_SHYFEM_HOSTNAME=$this_server
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
  [ "$SERVER_SHYFEM_COMPILER" = "GNU_GFORTRAN" ] && return
  echo "setting server variables on $this_server for compiler gfortran"

  basedir=/cineca/prod/opt/libraries

  export SERVER_SHYFEM_HOSTNAME=$this_server
  export SERVER_SHYFEM_COMPILER=GNU_GFORTRAN

  export NETCDF_SHYFEM_LIBCDIR=$basedir/netcdf/4.4.1/gnu--6.1.0
  export NETCDF_SHYFEM_LIBFDIR=$basedir/netcdff/4.4.4/gnu--6.1.0
  export NETCDF_SHYFEM_INCDIR=$basedir/netcdff/4.4.4/gnu--6.1.0
  export NETCDF_SHYFEM_MODDIR=$basedir/netcdff/4.4.4/gnu--6.1.0

  AddLibrary $NETCDF_SHYFEM_LIBCDIR lib
  AddLibrary $NETCDF_SHYFEM_LIBFDIR lib
}

#--------------------------------------------------------------

