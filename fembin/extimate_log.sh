#!/bin/sh
#
# estimates running time in all subdirectories
#
#------------------------------------------------

files=$( findf '*.log' 2>/dev/null )


for file in $files
do
  #echo $file
  extimate.sh $file
done

#------------------------------------------------

