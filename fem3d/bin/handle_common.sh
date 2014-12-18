#!/bin/sh
#
# handles conversion of common blocks
#
#-------------------------------------------------------

bindir=/home/georg/fem/fem3d/bin

what=$1

files=`ls *.f`
#files=`ls sub*.f`
#files=`ls ht.f 2> /dev/null`

#---------------------------------------------------------------

if [ -z "$what" ]; then
  echo "Usage: handle_common.sh [-h|-help] [-option]"
  exit 1
elif [ $what = "-h" -o $what = "-help" ]; then
  echo "Usage: handle_common.sh [-h|-help] [-option]"
  echo "  options:"
  echo "    -h|-help     this help screen"
  echo "    -subst       substitutes common and creates .new"
  echo "    -clean       cleans comments from substitution and creates .clean"
  echo "    -cleanall    cleans comments from substitution, no revert possible"
  echo "    -new         copies .new files to .f, creates .old"
  echo "    -old         copies .old files to .f"
  exit 1
else
  echo "option: $what"
fi

#---------------------------------------------------------------

if [ $what = "-subst" ]; then
  echo "substituting in files"
  $bindir/newsubst.pl -subst -write $files
elif [ $what = "-clean" ]; then
  echo "cleaning in files"
  $bindir/newsubst.pl -clean $files
elif [ $what = "-cleanall" ]; then
  echo "cleaning all in files"
  $bindir/newsubst.pl -clean $files
  rm -f *.new *.old
  files=`ls *.clean 2> /dev/null`
  for file in $files
  do
    orig=`basename $file .clean`
    echo "moving $file -> $orig"
    mv -f $file $orig
  done
elif [ $what = "-new" ]; then
  files=`ls *.new 2> /dev/null`
  for file in $files
  do
    new=$file
    orig=`basename $new .new`
    old=$orig.old
    echo "copying $new to $orig, creating $old"
    cp -f $orig $old
    cp -f $new $orig 
  done
elif [ $what = "-old" ]; then
  files=`ls *.old 2> /dev/null`
  for file in $files
  do
    old=$file
    orig=`basename $old .old`
    echo "copying $old to $orig"
    cp -f $old $orig 
  done
else
  echo "*** unknown option: $what"
  exit 1
fi

#---------------------------------------------------------------

