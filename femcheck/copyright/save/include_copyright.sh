#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# inserts copyright into files
#
#------------------------------------------------------------

home=$HOME
copydir="$home/shyfem/femcheck/copyright"

#------------------------------------------------------------

Usage()
{
  echo "Usage: include_copyright.sh [-h|-help] [-options] {dir|file(s)}"
}

FullUsage()
{
  Usage
  echo "  inserts copyright into files"
  echo ""
  echo "  -h|-help      this help screen"
  echo "  -type type    use this type for files"
  echo ""
}

#------------------------------------------------------------

type=""

while [ -n "$1" ]
do
   case $1 in
        -type)          type="-type=$2"; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

if [ $# -eq 0 ]; then
  Usage; exit 1
fi

if [ -d "$1" ]; then
  cd $1
  files=$( ls )
else
  files=$*
fi

if [ -z "$files" ]; then
  Usage; exit 1
fi

#------------------------------------------------------------

for file in $files
do
  [ ! -f $file ] && continue			#no regular file
  $copydir/check_copyright.pl -quiet $file
  [ $? -eq 0 ] && continue			#file has already copyright
  echo "inserting copyright in $file"
  $copydir/include_copyright.pl $type $file > $file.copy
  status=$?
  if [ $status -eq 0 ]; then
    if [ -x $file ]; then
      chmod +x $file.copy
    fi
    mv -f $file.copy $file
  else
    echo "error in command...not copying - status=$status"
    rm -f $file.copy
  fi
done

#------------------------------------------------------------

