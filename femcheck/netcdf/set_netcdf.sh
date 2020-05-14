#!/bin/bash

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#-----------------------------------------------------
#
# finds netcdf libraries
#
#-----------------------------------------------------
#   ldconfig -p | grep libnetcdff
#   whereis libnetcdff
#   locate libnetcdff.a
#-----------------------------------------------------

SetupVars()
{
  debug="YES"
  debug="NO"

  what=$1			#can be lib or include or file name
  compiler=$2
  usrdir=$3

  [ -z "$what" ] && what=lib
  [ -z "$compiler" ] && compiler=gfortran

  Info "set_netcdf.sh what=$what compiler=$compiler usrdir=$usrdir"
}

#-----------------------------------------------------

Stderr()
{
  echo "$1: $2" >> /dev/stderr
}

Info()
{
  [ $debug = "NO" ] && return
  Stderr "Info" "$1"
}

Warning()
{
  Stderr "Warning" "$1"
}

Error()
{
  Stderr "*** Error" "$1"
}

#-----------------------------------------------------

SetupDirs()
{
  dirs="/usr /opt/sw/netcdf /usr/local/netcdf /usr/local"
  dirs_intel=" \
		/usr/local/intel \
		/usr/lib/x86_64-linux-gnu \
		"
  dirs_gfortran=" \
		/usr/lib64 /usr/lib64/gfortran/modules \
		/usr/lib/x86_64-linux-gnu \
		"

  [ -n $NETCDF_SHYFEM_LIBCDIR ] && usrdir="$usrdir $NETCDF_SHYFEM_LIBCDIR"
  [ -n $NETCDF_SHYFEM_LIBFDIR ] && usrdir="$usrdir $NETCDF_SHYFEM_LIBFDIR"
  [ -n $NETCDF_SHYFEM_INCDIR ] && usrdir="$usrdir $NETCDF_SHYFEM_INCDIR"
  [ -n $NETCDF_SHYFEM_MODDIR ] && usrdir="$usrdir $NETCDF_SHYFEM_MODDIR"

  #echo "usrdir = $usrdir"

  if [ $compiler = "INTEL" ]; then
    dirs="$usrdir $dirs_intel $dirs"
  else
    dirs="$usrdir $dirs_gfortran $dirs"
  fi
}

TestFile()
{
   nc_out=$1
   Info "checking $nc_out/$2"
   [ -f $nc_out/$2 ] && return 0
   return 1
}

LookUp()
{
  pre=$1
  file=$2

  nc_out=""

  for dir in $dirs
  do
    [ -n "$pre" ] && TestFile $dir/$pre $file && return 0
    TestFile $dir $file && return 0
  done

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

SetFileDir()
{
  pre=""
  [[ $1 == *.mod ]] && pre=include
  [[ $1 == *.inc ]] && pre=include
  [[ $1 == *.a ]] && pre=lib
  [[ $1 == *.so ]] && pre=lib

  if [ -z "$pre" ]; then
    if [[ $1 == lib* ]]; then
      SetLibDir $1
    else
      Error "unknown file type $1 ...aborting"
      exit 7
    fi
  else
    SetDir $pre $1
  fi
}

SetDir()
{
  pre=$1
  search=$2
  shift 2

  files=$search
  for ext
  do
    files="$files $search.$ext"
  done
  Info "looking for $files"

  for file in $files
  do
    LookUp "$pre" $file
    status=$?
    [ $status -eq 0 ] && break
  done
  
  if [ $status -ne 0 ]; then
    Error "cannot find directory for $search ...aborting"
    exit 1
  fi

  Info "directory found for $file: $nc_out"

  if [ "$pre" = "lib" ]; then
    #name=$( echo $file | sed -e 's/^lib//' )
    name=$( echo $search | sed -e 's/^lib//' )
    nc_libdir="$nc_libdir -L$nc_out"
    nc_libs="$nc_libs -l$name"
  elif [ "$pre" = "include" ]; then
    nc_inc="$nc_inc -I$nc_out"
  else
    nc_inc="$nc_out"
  fi
}

Unique()
{
  echo "$*" | xargs -n1 | sort -u | xargs
}

#-----------------------------------------------------

#-----------------------------------------------------

SetupVars $*
SetupDirs

if [ $what = "lib" ]; then
  SetLibDir libnetcdf
  SetLibDir libnetcdff
  nc_libdir=$( Unique $nc_libdir )
  nc_libds=$( Unique $nc_libs )
  echo "$nc_libdir $nc_libs"
elif [ $what = "include" ]; then
  SetIncDir netcdf.inc
  SetIncDir netcdf.mod
  nc_inc=$( Unique $nc_inc )
  echo "$nc_inc"
else
  SetFileDir $what
  echo "$nc_out"
fi

#echo "$nc_libdir" | xargs -n1 | sort -u | xargs

#-----------------------------------------------------

