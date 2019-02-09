#!/bin/sh

home=$HOME
copydir="$home/shyfem/femcheck/copyright"

Usage()
{
  echo "Usage: include_copyright.sh [-dir subdir] file(s)"
  exit 1
}

if [ "$1" = "-dir" ]; then
  dir=$2
  if [ -z "$dir" -o ! -d "$dir" ]; then
    echo "not a directory: $dir"
    Usage
  fi
  shift; shift
  cd $dir
  files=$( ls $1 )
else
  files=$*
fi

[ -z "$files" ] && Usage

for file in $files
do
  [ ! -f $file ] && continue			#no regular file
  $copydir/check_copyright.pl -quiet $file
  [ $? -eq 0 ] && continue			#file has already copyright
  echo "inserting copyright in $file"
  $copydir/include_copyright.pl $file > $file.copy
  mv -f $file.copy $file
done

