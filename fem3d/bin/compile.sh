#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

comp=gfortran
comp="gfortran -cpp"
dmain=dummy_main.f

actdir=$( pwd )
tmpdir=tmp_compile

#--------------------------------------

MakeDummy()
{
  echo "	end" > $dmain
}

DoOutput()
{
  if [ $status -eq 0 ]; then
    echo "$status  $file"
  else
    [ $only0 = "NO" ] && echo "$status  $file"
    if [ $show = "YES" -o $verbose = "YES" ]; then
      cat logfile.log
    fi
  fi
}

#--------------------------------------

Usage()
{
  echo "Usage: compile.sh [-h|-help] [-options] fortran-file(s)"
  exit 0
}

FullUsage()
{
  echo "Usage: compile.sh [-h|-help] [-options] fortran-file(s)"
  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -show            show error messages"
  echo "  -verbose         be verbose"
  echo "  -quiet           be quiet"
  echo "  -single          compile files one by one"
  echo "  -keep            do not delete directory and objetcs"
  echo "  -only0           only show files with status 0"
  echo ""

  exit 0
}

#--------------------------------------

show=NO
verbose=NO
show=NO
quiet=NO
single=NO
keep=NO
only0=NO

while [ -n "$1" ]
do
   case $1 in
        -quiet)         quiet="YES";;
        -verbose)       verbose="YES";;
        -show)          show="YES";;
        -single)        single="YES";;
        -keep)          keep="YES";;
        -only0)         only0="YES";;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   shift
done

[ $# -eq 0 ] && Usage

#--------------------------------------

mkdir -p $tmpdir
cp $* $tmpdir
cd $tmpdir

MakeDummy

if [ $single = "YES" ]; then
  for file
  do
    $comp $dmain $file > logfile.log 2>&1
    status=$?
    DoOutput
  done
else
  $comp $dmain $* > logfile.log 2>&1
  status=$?
  DoOutput
fi

cd $actdir
if [ $keep = "NO" ]; then
  rm -rf $tmpdir
fi

#------------------------------------------------------------------------

