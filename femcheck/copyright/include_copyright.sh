#!/bin/sh

copydir="/home/georg/shyfem/femcheck/copyright"

Usage()
{
  echo "Usage: include_copyright.sh {f|c|t} file(s)"
  exit 1
}

[ $# -eq 0 ] && Usage

lang=$1
shift

for file
do
  echo "$file"
  $copydir/include_copyright.pl $lang $file > $file.copy
  mv -f $file.copy $file
done

