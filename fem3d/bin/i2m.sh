#!/bin/sh
#
# changes includes to modules
#
#----------------------------------------------------

bin=~/shyfem/fem3d/bin

what=$1
if [ $# -eq 0 ]; then
  echo "Usage: i2m.sh what"
  echo "  include what.h will be substituted with module what"
  exit 0
fi

include=$what.h

files=$( grep -l $include *.f90 *.f )

for file in $files
do
  echo "  $file"
  $bin/i2m.pl $what $file > $file.mmm
  #cp $file $file.old
  mv $file.mmm $file
done

#----------------------------------------------------

