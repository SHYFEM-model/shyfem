#!/bin/sh
#
# rename files
#
#----------------------------------------------------------------

copy="NO"
if [ "$1" = "-copy" ]; then
  copy="YES"
  shift
fi

if [ $# -lt 3 ]; then
  echo "Usage: rename.sh [-copy] from_pattern to_pattern file(s)"
  exit 0
fi

pattern1=$1
pattern2=$2
shift; shift

for file
do
  new=$( echo $file | sed -e "s/$pattern1/$pattern2/" )
  echo "$file -> $new"
  [ $copy = "YES" ] && mv $file $new
done

if [ $copy = "YES" ]; then
  echo "all files have been renamed."
else
  echo "no files have been renamed. To actually rename use -copy."
fi

#----------------------------------------------------------------

