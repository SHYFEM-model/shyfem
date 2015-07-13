#!/bin/sh
#
# handles conversion of common blocks
#
#-------------------------------------------------------

bindir=/home/georg/fem/fem3d/bin

what=$1
[ -n "$what" ] && shift
files=$*

if [ -z "$files" ]; then
  files=`ls *.f`
fi

#echo $files
#exit 1

#---------------------------------------------------------------

if [ -z "$what" ]; then
  echo "Usage: handle_common.sh [-h|-help] [-option] [files]"
  exit 1
elif [ $what = "-h" -o $what = "-help" ]; then
  echo "Usage: handle_common.sh [-h|-help] [-option] [files]"
  echo "  options:"
  echo "    -h|-help     this help screen"
  echo "    -subst       substitutes common and creates .new"
  echo "    -include     cleans include files and creates .new"
  echo "    -inc2use     substitutes include with use"
  echo "    -revert      revert to original version"
  echo "    -check       checks remaining common blocks"
  echo "    -clean       deletes .old and .new files"
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
elif [ $what = "-include" ]; then
  echo "cleaning include in files"
  $bindir/newsubst.pl -include -write $files
elif [ $what = "-check" ]; then
  echo "checking remaining common blocks"
  $bindir/newsubst.pl -check $files
elif [ $what = "-inc2use" ]; then
  echo "substituting include with use"
  $bindir/newsubst.pl -inc2use -write $files
elif [ $what = "-revert" ]; then
  echo "reverting to old version"
  $bindir/newsubst.pl -revert $files
elif [ $what = "-clean" ]; then
  echo "cleaning files"
  rm -f *.new *.old
elif [ $what = "-cleanall" ]; then
  echo "cleaning all in files"
  rm -f *.new *.old
  $bindir/newsubst.pl -clean $files
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

