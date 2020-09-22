#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# this script handles server specific settings
#
#--------------------------------------------------------------

debug="YES"
debug="NO"

#--------------------------------------------------------------

FindServer()
{
  if [ "$host" = galileo -o "$HPC_SYSTEM" = galileo ]; then
    server=galileo
  elif [ "$host" = fluxus -o "$SGE_CLUSTER_NAME" = fluxus ]; then
    server=fluxus
  fi
}

SetCompilerOption()
{
  if [ -z "$compiler" ]; then	#no compiler given
    option="{gfortran|intel}"
  elif [ $compiler = GNU_GFORTRAN ]; then
    option="gfortran"
  elif [ $compiler = gfortran ]; then
    option="gfortran"
  elif [ $compiler = INTEL ]; then
    option="intel"
  elif [ $compiler = intel ]; then
    option="intel"
  elif [ $compiler = PGI ]; then
    option="pgi"
  elif [ $compiler = pgi ]; then
    option="pgi"
  else
    echo "*** (check_server.sh) unknown compiler: $compiler"
  fi
}

CheckServer()
{
  [ $debug = "YES" ] && echo "checking server $server"

  if [ "$debug" = "YES" ]; then
    echo "host = $host   HPC_SYSTEM = $HPC_SYSTEM"
    echo "SERVER_SHYFEM_HOSTNAME = $SERVER_SHYFEM_HOSTNAME"
    echo "FORTRAN_COMPILER = $compiler"
    echo "SERVER_SHYFEM_COMPILER = $SERVER_SHYFEM_COMPILER"
  fi

  SetCompilerOption
  script="$BASH_SOURCE -load $option"

  if [ -n "$server" ]; then
    if [ "$SERVER_SHYFEM_HOSTNAME" != $server ]; then
      echo "initialization routine for $server has not been run..."
      [ -n "$compiler" ] && echo "shyfem compiler is set to $compiler"
      echo "please run the following command:"
      echo "  source $script"
      exit 1
    fi
    if [ -z "$compiler" ]; then
      echo "this script must be run from within Makefile:"
      echo "  make check_server"
      exit 1
    elif [ "$SERVER_SHYFEM_COMPILER" != "$compiler" ]; then
      echo "compiler has changed..."
      echo "server prepared for compiler: $SERVER_SHYFEM_COMPILER"
      echo "shyfem compiler is set to $compiler"
      echo "please run the following command:"
      echo "  source $script"
      exit 1
    fi
  fi
}

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

ResetDirsGeneric()
{
  echo "resetting server variables (generic)"
  export SERVER_SHYFEM_HOSTNAME=
  export SERVER_SHYFEM_COMPILER=
  export NETCDF_SHYFEM_LIBCDIR=
  export NETCDF_SHYFEM_LIBFDIR=
  export NETCDF_SHYFEM_INCDIR=
  export NETCDF_SHYFEM_MODDIR=
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ORIG
}

#--------------------------------------------------------------

Usage()
{
  echo "Usage: check_server.sh [-h|-help] [-action] [compiler]"
}

FullUsage()
{
  Usage
}

ParseOptions()
{
  while [ -n "$1" ]
  do
    case $1 in
        -server)        what="server";;
        -check)         what="check";;
        -show)          what="show";;
        -reset)         what="reset";;
        -load)          what="load"; setting=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
    esac
    shift
  done

  compiler=$1
}

#--------------------------------------------------------------

(return 0 2>/dev/null) && sourced=1 || sourced=0

ParseOptions $*
FindServer

#echo "$server $what"

if [ -z "$what" ]; then
  Usage
elif [ "$what" = "server" ]; then
  echo $server
elif [ "$what" = "check" ]; then
  CheckServer
elif [ "$what" = "show" ]; then
  ShowDirs
elif [ "$what" = "reset" ]; then
  CheckSourced $*
  ResetDirsGeneric
elif [ "$what" = "load" ]; then
  if [ -z "$server" ]; then
    echo "no setting for this server existing... cannot load"
    setting="none"
  elif [ -z "$setting" ]; then
    echo "no setting for this server existing... cannot load"
    setting="none"
  else
    dir=$( dirname $BASH_SOURCE )
    [ -z "$dir" ] && dir="."
    server_script=$dir/$server.sh
    if [ -f $server_script ]; then
      #echo "running... $server_script (from dir $dir)"
      . $server_script
    else
      echo "no such script: $server_script"
      setting=unknown
    fi
  fi
  if [ "$setting" = "gfortran" ]; then
    CheckSourced $*
    ResetDirsGeneric
    ResetDirs
    SetDirs_gfortran
  elif [ "$setting" = "intel" ]; then
    CheckSourced $*
    ResetDirsGeneric
    ResetDirs
    SetDirs_intel
  elif [ "$setting" = "none" ]; then
    :
  else
    echo "cannot load this setting: $setting"
  fi
fi

#--------------------------------------------------------------

