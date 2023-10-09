#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

compiler=g77
compiler=gfortran

tmpdir=tmp_compile

mkdir $tmpdir

echo "	end" > main_dummy.f

if [ $# -eq 0 ]; then
  echo "Usage: compile_test.sh [-single] fortran_files.f"
  exit 0
fi

single="NO"
if [ $1 = "-single" ]; then
  single="YES"
  shift
fi

status=0

if [ $single = "YES" ]; then	# compile files one by one
  for file in $*
  do
    echo "--------------------------------------------"
    echo "   $file"
    echo "--------------------------------------------"
    $compiler main_dummy.f $file
    [ $? -ne 0 ] && $(( status = status + 1 ))   
  done
else				# compile all files together
  $compiler main_dummy.f $*
  status=$?
fi

rm -f main_dummy.f

exit $status

