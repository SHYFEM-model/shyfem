#!/bin/sh
#
# finds netcdf libraries
#
#-----------------------------------------------------
#   ldconfig -p | grep libnetcdff
#   whereis libnetcdff
#   locate libnetcdff.a
#-----------------------------------------------------

debug="YES"
debug="NO"

what=$1			#can be lib or include
compiler=$2
usrdir=$3

[ -z "$what" ] && what=lib
[ -z "$compiler" ] && compiler=gfortran

#-----------------------------------------------------

SetupDirs()
{
  dirs="/usr /opt/sw/netcdf /usr/local/netcdf /usr/local"

  if [ $compiler = "INTEL" ]; then
    dirs="$usrdir /usr/local/intel $dirs"
  else
    dirs="$usrdir $dirs"
  fi
}

LookUp()
{
  pre=$1
  file=$2

  nc_out=""

  for dir in $dirs
  do
    if [ -f $dir/$pre/$file ]; then
      nc_out=$dir/$pre
      return 0
    fi
  done

  dir=/usr/lib/x86_64-linux-gnu
  if [ -f $dir/$file ]; then
    nc_out=$dir
    return 0
  fi

  return 1
}

SetLibDir()
{
  SetDir lib $1 a so
}

SetIncDir()
{
  SetDir include $1
}

SetDir()
{
  pre=$1
  lib=$2
  shift 2

  file=$lib
  if [ $# -gt 0 ]; then
    file=$lib.$1
    shift
  fi

  LookUp $pre $file
  if [ $? -ne 0 ]; then
    file=$lib.$1
    LookUp $pre $file
    if [ $? -ne 0 ]; then
      echo "*** cannot find directory for $lib... aborting" >> /dev/stderr
      exit 1
    fi
  fi

  if [ $debug = "YES" ]; then
    echo "directory found for $file: $nc_out" >> /dev/stderr
  fi

  if [ $pre = "lib" ]; then
    #name=$( echo $file | sed -e 's/^lib//' )
    name=$( echo $lib | sed -e 's/^lib//' )
    nc_libdir="$nc_libdir -L$nc_out"
    nc_libs="$nc_libs -l$name"
  else
    nc_inc="$nc_inc -I$nc_out"
  fi
}

Unique()
{
  echo "$*" | xargs -n1 | sort -u | xargs
}

#-----------------------------------------------------

SetupDirs

if [ $what = "lib" ]; then
  SetLibDir libnetcdff
  SetLibDir libnetcdf
  nc_libdir=$( Unique $nc_libdir )
  nc_libds=$( Unique $nc_libs )
  echo "$nc_libdir $nc_libs"
else
  SetIncDir netcdf.inc
  SetIncDir netcdf.mod
  nc_inc=$( Unique $nc_inc )
  echo "$nc_inc"
fi

#echo "$nc_libdir" | xargs -n1 | sort -u | xargs

#-----------------------------------------------------

