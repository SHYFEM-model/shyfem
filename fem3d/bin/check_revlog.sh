#!/bin/sh

bindir=$HOME/fem/fem3d/bin

for file
do
  echo "$file"
  $bindir/check_revlog.pl $file
done

