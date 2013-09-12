#!/bin/sh

compiler=g77
compiler=gfortran

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

if [ $single = "YES" ]; then	# compile files one by one
  for file in $*
  do
    echo "--------------------------------------------"
    echo "   $file"
    echo "--------------------------------------------"
    $compiler main_dummy.f $file
  done
else				# compile all files together
  $compiler main_dummy.f $*
fi

rm -f main_dummy.f

