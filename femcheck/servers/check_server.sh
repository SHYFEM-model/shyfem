#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#--------------------------------------------------------------

debug="NO"

#--------------------------------------------------------------

FindServer()
{
  if [ "$host" = galileo -o "$HPC_SYSTEM" = galileo ]; then
    server=galileo
  fi
}

CheckServer()
{
  #echo "checking server $server"

  if [ "$debug" = "YES" ]; then
    echo "host = $host   HPC_SYSTEM = $HPC_SYSTEM"
    echo "SERVER_SHYFEM_HOSTNAME = $SERVER_SHYFEM_HOSTNAME"
    echo "FORTRAN_COMPILER = $compiler"
    echo "SERVER_SHYFEM_COMPILER = $SERVER_SHYFEM_COMPILER"
  fi

  if [ "$host" = $server -o "$HPC_SYSTEM" = $server ]; then
    if [ "$SERVER_SHYFEM_HOSTNAME" != $server ]; then
      echo "initialization routine for $server has not been run..."
      echo "please run one of the following commands (depending on compiler)"
      echo "  source femcheck/servers/galileo.sh {-gfortran|-intel}"
      exit 1
    fi
    if [ -z "$compiler" ]; then
      echo "this script must be run from within Makefile:"
      echo "  make check_server"
      exit 1
    elif [ "$SERVER_SHYFEM_COMPILER" != "$compiler" ]; then
      echo "compiler has changed..."
      echo "server prepared for compiler: $SERVER_SHYFEM_COMPILER"
      echo "fortran compiler chosen: $compiler"
      echo "please run one of the following commands (depending on compiler)"
      echo "  source femcheck/servers/galileo.sh {-gfortran|-intel}"
      exit 1
    fi
  fi
}

#--------------------------------------------------------------

compiler=$1

FindServer
if [ "$compiler" = "-print" ]; then
  echo $server
  exit 0
elif [ -n "$server" ]; then
  CheckServer
fi

#--------------------------------------------------------------

