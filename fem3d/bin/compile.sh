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

if [ $# -eq 0 ]; then
  echo "Usage: compile.sh [-single|-show] fortran-file(s)"
  exit 0
fi

single=NO
show=NO
keep=NO
if [ "$1" = "-single" ]; then
  single="YES"
  shift
fi
if [ "$1" = "-show" ]; then
  show="YES"
  shift
fi
if [ "$1" = "-keep" ]; then
  keep="YES"
  shift
fi

#--------------------------------------

MakeDummy()
{
  echo "	end" > $dmain
}

#--------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: compile.sh fortran-file(s)"
  exit 0
fi

mkdir -p $tmpdir
cp $* $tmpdir
cd $tmpdir

MakeDummy

if [ $single = "YES" ]; then
  for file
  do
    $comp $dmain $file > logfile.log 2>&1
    status=$?
    echo "$status  $file"
    [ $status -ne 0 -a $show = "YES" ] && cat logfile.log
  done
else
  $comp $dmain $* > logfile.log 2>&1
  status=$?
  echo "status = $status"
  [ $status -ne 0 -a $show = "YES" ] && cat logfile.log
fi

cd $actdir
if [ $keep = "NO" ]; then
  rm -rf $tmpdir
fi

#------------------------------------------------------------------------

