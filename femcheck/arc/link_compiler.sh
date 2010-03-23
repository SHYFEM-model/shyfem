#!/bin/sh

log=COMPILER_LOG
fembin=$HOME/fem/fembin

CheckCommand()
{
  command=$1

  ($command) >> $log 2>&1

  status=$?

}

#---------------------------------------------------------

GetFortranCompiler()
{
  CheckCommand "g77 -v"
  if [ $status -eq 0 ]; then
    fortran_compiler=`which g77`
    return
  fi

  CheckCommand "f77 -v"
  if [ $status -eq 0 ]; then
    fortran_compiler=`which f77`
    return
  fi

  echo "Cannot find Fortran compiler... exiting"
  exit 1
}

SetFortranCompiler()
{
  CheckCommand "g77 -v"
  if [ $status -ne 0 ]; then
    echo "linking $fortran_compiler $fembin/g77"
    ln -s $fortran_compiler $fembin/g77
  fi

  CheckCommand "f77 -v"
  if [ $status -ne 0 ]; then
    echo "linking $fortran_compiler $fembin/f77"
    ln -s $fortran_compiler $fembin/f77
  fi
}

#---------------------------------------------------------

GetCCompiler()
{
  CheckCommand "gcc -v"
  if [ $status -eq 0 ]; then
    c_compiler=`which gcc`
    return
  fi

  CheckCommand "cc -v"
  if [ $status -eq 0 ]; then
    c_compiler=`which cc`
    return
  fi

  echo "Cannot find C compiler... exiting"
  exit 1
}

SetCCompiler()
{
  CheckCommand "gcc -v"
  if [ $status -ne 0 ]; then
    echo "linking $c_compiler $fembin/gcc"
    ln -s $c_compiler $fembin/gcc
  fi

  CheckCommand "cc -v"
  if [ $status -ne 0 ]; then
    echo "linking $c_compiler $fembin/cc"
    ln -s $c_compiler $fembin/cc
  fi
}

#---------------------------------------------------------

GetFortranCompiler 
SetFortranCompiler

GetCCompiler 
SetCCompiler

#---------------------------------------------------------

