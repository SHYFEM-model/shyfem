#!/bin/sh
#
# resets bounding box from plot hints in eps file

tmp=tmp.tmp

if [ $# -eq 0 ]; then
  echo "Usage: psresetbb.sh file[s]"
  exit 0
fi

##################################################################

for file
do

echo "$file"

##################### extract bounding box information

bb=`grep InternalBoundingBox $file | head -1 | sed -e 's/^.*Box//'`
bb=`echo $bb | sed -e 's/  */ /g'`
bbb="%%BoundingBox: $bb"

#echo "$bb"
#echo "$bbb"

##################### substitute in ps file and rename files

if [ -n "$bb" ]; then
  sed -e "s/^%%BoundingBox:.*/$bbb/" $file > $tmp

  #mv $file $file.bak
  mv $tmp $file
fi

done

##################################################################

