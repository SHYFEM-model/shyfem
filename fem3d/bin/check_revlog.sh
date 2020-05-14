#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

bindir=$HOME/fem/fem3d/bin

for file
do
  echo "$file"
  $bindir/check_revlog.pl $file
done

