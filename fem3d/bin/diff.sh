#!/bin/bash
#
#--------------------------------------------------

#rdir=~/shyfem/fem3d
#rdir=../../shyfem/fem3d
# rdir is set in FindRemote()

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
  echo "  -auto            automatic accept copying and compiling files"
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

FindRemote()
{
  this_dir=$( pwd )
  this_shyfem=$( echo $this_dir | sed -e 's/.*shyfem/shyfem/' \
		| sed -e 's/\/.*//' )

  if [ $this_shyfem = "shyfem" ]; then
    other_shyfem=shyfem-mpi
  elif [ $this_shyfem = "shyfem-mpi" ]; then
    other_shyfem=shyfem
  else
    echo "*** cannot determine this shyfem dir: $this_shyfem"
    exit 1
  fi
  
  other_dir=$( echo $this_dir | sed -e "s/$this_shyfem/$other_shyfem/" )

  echo "this shyfem: $this_shyfem   other shyfem: $other_shyfem"
  echo "this dir: $this_dir"
  echo "other dir: $other_dir"

  rdir=$other_dir
  #exit 1
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

FindRemote

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

if [ $compile = "NO" ]; then
  answer=$( YesNo "compile files locally and in $rdir?" )
  [ "$answer" = "y" ] || exit 0
fi

echo "compiling locally and in $rdir..."
make
make
make
cd $rdir
make
make
make

#--------------------------------------------------

