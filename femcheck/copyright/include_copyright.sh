#!/bin/sh
#
# inserts copyright into files
#
#------------------------------------------------------------

home=$HOME
copydir="$home/shyfem/femcheck/copyright"

#------------------------------------------------------------

Usage()
{
  echo "Usage: include_copyright.sh [-dir subdir] file(s)"
}

FullUsage()
{
  Usage
  echo "  inserts copyright into files"
  echo ""
  echo "  -h|-help      this help screen"
  echo "  -dir subdir   include copyright in all files in this subdir"
  echo "  -type type    use this type for files"
  echo ""
}

#------------------------------------------------------------

type=""

while [ -n "$1" ]
do
   case $1 in
        -dir)           dir=$2; shift;;
        -type)          type="-type=$2"; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

if [ -n "$dir" ]; then
  if [ ! -d "$dir" ]; then
    echo "not a directory: $dir"
    Usage; exit 1
  fi
  cd $dir
  files=$( ls $1 )
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
    mv -f $file.copy $file
  else
    echo "error in command...not copying - status=$status"
    rm -f $file.copy
  fi
done

#------------------------------------------------------------

