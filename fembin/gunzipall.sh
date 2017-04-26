#!/bin/sh
#
# unzips file(s) until we have original file
# this resolves cases like grid.grd.gz.gz.gz -> grid.grd
#
#-----------------------------------------------------

files=$( ls *.gz 2>/dev/null )

for file in $files
do
  orig=$( echo $file | sed -e 's/.gz.*//' )
  echo "$file -> $orig"
  while [ -n "$file" ]
  do
    gunzip $file
    file=$( ls $orig.gz* 2>/dev/null )
    echo "...$file"
  done
done

