#!/bin/sh
#
# finds file or pattern in source directories
#
#-----------------------------------------------

base="$HOME/shyfemcm/shyfemcm"

src="src src/util src/shyutil src/shympi src/shyutilmpi"
tools="tools/shybas tools/shyelab tools/shyplot tools/shypre"

alldir="$src $tools"

actdir=$( pwd )

#-----------------------------------------------

Usage()
{
  echo "Usage: shyfind.sh [-h|-help] [-options] pattern"
}

FullUsage()
{
  Usage
  
  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -verbose         be verbose"
  echo "  -quiet           be quiet"
  echo "  -file            look for pattern in file names"
  echo "  -pattern         look for pattern in files"
  echo ""
}

#-----------------------------------------------

do_file=NO
do_pattern=NO
quiet=NO
verbose=NO
while [ -n "$1" ]
do
   case $1 in
        -quiet)         quiet="YES";;
        -verbose)       verbose="YES";;
        -file)          do_file="YES";;
        -pattern)       do_pattern="YES";;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "No such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

[ $# -eq 0 ] && Usage && exit 0

echo "looking for $*"

#-----------------------------------------------

if [ $do_file = "YES" ]; then
  for dir in $alldir
  do
    echo "---- $dir ----" >&2
    cd $base/$dir
    ls $* 2> /dev/null
    cd $actdir
  done
elif [ $do_pattern = "YES" ]; then
  for dir in $alldir
  do
    echo "---- $dir ----" >&2
    cd $base/$dir
    grep $* *.f90 *.f 2> /dev/null
    cd $actdir
  done
elif [ $do_pattern = "YES" ]; then
  :
else
  echo "do not know what to do..."
fi

#-----------------------------------------------

