#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

files=`ls *.py`

echo " " > cm.dat

for file in $files
do
  echo "$file"
  ./prep.pl $file >> cm.dat
  echo " " >> cm.dat
done

mv cm.dat oceanlib.dat

