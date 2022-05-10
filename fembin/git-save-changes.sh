#!/bin/sh
#
# saves all changes of local files to a new directory
#
#-----------------------------------------------------------

savedir=save/save.$$
mkdir -p $savedir

files=$( git s | grep modified: | sed -e 's/.*modified:   //' )

for item in $files
do
  dir=$( dirname $item )
  file=$( basename $item )
  echo "file: $item - $dir - $file"
  mkdir -p $savedir/$dir
  cp $item $savedir/$item
done

echo "files saved to dir: $savedir"

#-----------------------------------------------------------

