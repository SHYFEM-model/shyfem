#!/bin/sh
#
# test compile for file dependencies
#
#-------------------------------------------------------

F77=gfortran

DIRLIB_NETCDF=/usr/lib/x86_64-linux-gnu
LIBNETCDF=netcdff

LIBG_NETCDF="-L$DIRLIB_NETCDF -l$LIBNETCDF"
LIBF_NETCDF="$DIRLIB_NETCDF/lib$LIBNETCDF.a"

LIBGS=$LIBG_NETCDF
LIBFS=$LIBF_NETCDF

#-------------------------------------------------------

MakeFake()
{
  echo "\tend" > fake.f
}

#-------------------------------------------------------

MakeFake

$F77 fake.f
$F77 fake.f subclo.f
$F77 fake.f subdts.f 
$F77 fake.f subgeo.f
$F77 fake.f subiso8601.f
$F77 fake.f subscn.f
 
$F77 fake.f ncf_util.f $LIBGS
$F77 fake.f ncf_util.f ncf_tutil.f subdts.f subiso8601.f $LIBGS
$F77 fake.f ncf_dim_coords.f ncf_util.f $LIBGS
#$F77 fake.f ncf_util.f ncf_tutil.f subiso8601.f $LIBGS

#-------------------------------------------------------

