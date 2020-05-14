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
# shell for hsb to gray transformation

Usage()
{
  echo "Usage: hsb2gray.sh [ -h | -help ] [ -options ] files[s]"
}

FullUsage()
{
  Usage

  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -c               convert to color output"
  echo "  -i               invert color scale"
}

ErrorOption()
{
  echo "No such option : $1"
}

while [ -n "$1" ]
do
   case $1 in
        -c)             color="-color";;
        -i)             invert="-invert";;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   shift
done

if [ $# -eq 0 ]; then
  Usage
  exit 0
fi

for file
do
  f=`basename $file .ps`
  new=$f-c.ps
  echo "$file -> $new"
  hsb2gray $color $invert $file > $f-c.ps
done
