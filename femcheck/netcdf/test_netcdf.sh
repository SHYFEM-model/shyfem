#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#-------------------------------------------------------------

#LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
#LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

Compile_intel()
{
  echo "compiling with ifort"
  f77=ifort

  NETCDF_INCDIR=/usr/local/intel/include/
  NETCDF_MODDIR=/usr/local/intel/include/
  NETCDF_LIBDIR=/usr/local/intel/lib/
  NETCDF_INCDIR=$( ./set_netcdf.sh netcdf.inc INTEL )
  NETCDF_MODDIR=$( ./set_netcdf.sh netcdf.mod INTEL )
  NETCDF_LIBDIR=$( ./set_netcdf.sh lib INTEL )
  NETCDF_LIBS="-lnetcdf -lnetcdff"

  FLAGS="-I$NETCDF_INCDIR -I$NETCDF_MODDIR"
  FLAGS="$FLAGS $NETCDF_LIBDIR $NETCDF_LIBS"

  #echo "  using flags: $FLAGS"

  file=nc_include.f
  echo "  compiling $file"
  echo "    FLAGS: $FLAGS"
  #$f77 $FLAGS $file
  $f77 $file $FLAGS
  [ $? -eq 0 ] || exit 1

  file=nc_mod.f
  echo "  compiling $file"
  echo "    FLAGS: $FLAGS"
  $f77 $FLAGS $file
  [ $? -eq 0 ] || exit 1
}

Compile_gfortran()
{
  echo "compiling with gfortran"
  f77=gfortran

  NETCDF_INCDIR=/usr/include/
  NETCDF_MODDIR=/usr/lib64/gfortran/modules
  NETCDF_LIBDIR=/usr/lib
  NETCDF_INCDIR=$( ./set_netcdf.sh netcdf.inc )
  NETCDF_MODDIR=$( ./set_netcdf.sh netcdf.mod )
  NETCDF_LIBDIR=$( ./set_netcdf.sh lib )
  NETCDF_LIBS="-lnetcdf -lnetcdff"

  echo "NETCDF_INCDIR: $NETCDF_INCDIR"
  echo "NETCDF_MODDIR: $NETCDF_MODDIR"
  echo "NETCDF_LIBDIR: $NETCDF_LIBDIR"

  FLAGS="-I$NETCDF_INCDIR -I$NETCDF_MODDIR"
  FLAGS="$FLAGS $NETCDF_LIBDIR"

  #echo "  using flags: $FLAGS"

  file=nc_include.f
  echo "  compiling $file"
  echo "    FLAGS: $FLAGS"
  $f77 $FLAGS $file
  [ $? -eq 0 ] || exit 1

  file=nc_mod.f
  echo "  compiling $file"
  echo "    FLAGS: $FLAGS"
  $f77 $FLAGS $file
  [ $? -eq 0 ] || exit 1
}

#----------------------------------------------------------

#Compile_gfortran
Compile_intel

#----------------------------------------------------------

