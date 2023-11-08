#!/bin/sh
#
# finds file or pattern in source directories
#
#-----------------------------------------------

base="$HOME/shyfemcm/shyfemcm"

src="src/shyfem src/utils/generic src/utils/shyutil \
		src/utils/shympi src/utils/shyutilmpi"
tools="src/tools/shybas \
	src/tools/shypre \
	src/tools/shyadj src/tools/shynetcdf src/tools/shyparts \
	src/tools/shyelab src/tools/shyplot \
	src/tools/shyplot/util \
	"
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

Do_File()
{
  echo "looking for file $*"
  for dir in $alldir
  do
    echo "---- $dir ----" >&2
    cd $base/$dir
    ls $* 2> /dev/null
    cd $actdir
  done
}

Do_Pattern()
{
  echo "looking for pattern $*"
  for dir in $alldir
  do
    echo "---- $dir ----" >&2
    cd $base/$dir
    grep $* *.f90 *.f 2> /dev/null
    cd $actdir
  done
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
  Do_File $*
elif [ $do_pattern = "YES" ]; then
  Do_Pattern $*
else
  Do_File $*
  Do_Pattern $*
fi

#-----------------------------------------------

