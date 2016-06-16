#!/bin/sh

files=`ls *.py`

echo " " > cm.dat

for file in $files
do
  echo "$file"
  ./prep.pl $file >> cm.dat
  echo " " >> cm.dat
done

mv cm.dat oceanlib.dat

