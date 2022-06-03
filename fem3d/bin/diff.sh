#!/bin/bash
#
#--------------------------------------------------

rdir=~/shyfem/fem3d
rdir=../../shyfem/fem3d

#--------------------------------------------------

Usage()
{
  echo "Usage: diff.sh [-h|-help] [-options]"
}

FullUsage()
{
  echo ""

  Usage

  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -compare         compare files interactively"
  echo "  -compile         compile files locally and remotely"
  echo "  -auto            automatic accept copying files"
}

HandleOptions()
{
  compare="NO"
  compile="NO"
  auto="NO"

  while [ -n "$1" ]
  do
    case $1 in
        -h|-help)       FullUsage; exit 0;;
        -compare)       compare="YES";;
        -compile)       compile="YES";;
        -auto)          auto="YES";;
        -*)             echo "No such option: $1"; exit 1;;
        *)              break;;
    esac
    shift
  done
}

#--------------------------------------------------

YesNo()
{
  [ $auto = "YES" ] && echo "y" && return
  echo -n "$1 ? (y/n) : " | cat >&2
  read yesno
  echo "$yesno"
}

#--------------------------------------------------

HandleOptions $*

files=$( diffs -d $rdir *.f )

nf=$( echo "$files" | wc -w )
echo "$nf files found"
[ -z "$files" ] && exit 0

#--------------------------------

for file in $files
do
  echo $file
done

#--------------------------------

if [ "$compare" = "YES" ]; then
  echo "comparing files with $rdir..."
  for file in $files
  do
    echo "comparing $file with $rdir"
    tkdiff $file $rdir
  done
fi

#--------------------------------

answer=$( YesNo "Copy all files to $rdir" )
[ "$answer" = "y" ] || exit 0
echo "copying all files..."

for file in $files
do
  echo "copying $file to $rdir"
  cp $file $rdir
done

[ $compile = "NO" ] && exit 0

echo "compiling locally and in $rdir..."
make
make
make
cd $rdir
make
make
make

#--------------------------------------------------

