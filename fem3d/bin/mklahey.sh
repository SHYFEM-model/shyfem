#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

################################################## make directory

echo "....................................................... making directory"

if [ ! -d lahey ]; then
  mkdir lahey
fi

################################################## .f

list=`ls *.f`
echo "........................................................ copying f files"

for file in $list
do
  bfile=`basename $file .f`
  nfile=lahey/$bfile.for
  echo "$file -> $nfile"
  cp $file $nfile
done

################################################## .h

list=`ls *.h`
echo "........................................................ copying h files"

for file in $list
do
  nfile=lahey/$file
  echo "$file -> $nfile"
  cp $file $nfile
done

################################################## .F

list=`ls *.F`
echo "........................................................ copying F files"

for file in $list
do
  bfile=`basename $file .F`
  nfile=lahey/$bfile.for
  echo "$file -> $nfile"
  cpp -P $file $nfile
done

################################################## end

echo "................................................................... done"

