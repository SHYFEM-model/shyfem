#!/bin/bash
#
# check corresponance between objects and fortran files
#
#-------------------------------------------------

debug=NO
debug=YES

nf90=0
nf=0
nno=0
nobj=0
nnoobj=0

#------------------------------------------
# set directory
#------------------------------------------

dir=.
[ $# -gt 0 ] && dir=$1
cd $dir

#------------------------------------------
# check object files
#------------------------------------------

ofiles=$( ls *.o )

for obj in $ofiles
do
  name=$( basename $obj .o )
  #echo "$name"
  if [ -f $name.f90 ]; then
    nf90=$(( nf90 + 1 ))
  elif [ -f $name.f ]; then
    nf=$(( nf + 1 ))
  else
    echo "*** no such file: $name"
    nno=$(( nno + 1 ))
  fi
done

nfall=$(( nf90 + nf ))
[ $debug = YES ] && echo "$nf90  $nf  $nfall  $nno"

if [ $nf -gt 0 ]; then
  echo "files with ending .f present: $nf"
elif [ $nf -gt 0 ]; then
  echo "files not found: $nno"
fi

#------------------------------------------
# check fortran files
#------------------------------------------

if [ $nf -gt 0 ]; then
  ffiles=$( ls *.f90 *.f )
else
  ffiles=$( ls *.f90 )
fi

for ffile in $ffiles
do
  name=${ffile%.*}
  #echo "$name"
  if [ -f $name.o ]; then
    nobj=$(( nobj + 1 ))
  else
    nnoobj=$(( nnoobj + 1 ))
    echo "*** no object to file: $ffile"
  fi
done

[ $debug = YES ] && echo "$nobj  $nnoobj"

#------------------------------------------
# summary
#------------------------------------------

if [ $nnoobj -gt 0 ]; then
  echo "files without object present: $nnoobj"
elif [ $nfall -ne $nobj ]; then
  echo "number of fortran files and objects are different: $nfall $nobj"
else
  echo "successful check"
fi

#-------------------------------------------------
# end of routine
#------------------------------------------

