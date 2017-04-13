#!/bin/sh
#
# processes revision log from fortran files
#
#-------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
export FEMDIR

tmpfile=tmp.tmp
after=0
befor=30000000
recurse=NO

#-------------------------------------------------------------
# function definition
#-------------------------------------------------------------

Usage()
{
  echo 'Usage: revisionlog.sh [-h|-help] [options] files'
}

Help()
{
  Usage
  echo '  -help           show this help screen'
  echo '  -after date     show only dates after date'
  echo '  -befor date     show only dates befor date'
  echo '  -noname         do not show file name'
  echo '  -sepname        seperate visually file name'
  echo '  -recurse        recurse into subdirs'
  echo '  -check          check header of fortran files'
  echo 'give date in format YYYYMMDD'
}

#-------------------------------------------------------------
# treat command line
#-------------------------------------------------------------

if [ $# -eq 0 ]; then Usage; exit 0; fi

options=""

while [ $# -gt 0 ]
do
  case $1 in
        -h|-help) help; exit 0;;
        -after) after=$2; shift;;
        -befor) befor=$2; shift;;
        -noname) options="$options -noname";;
        -sepname) options="$options -sepname";;
        -check) options="$options -check";;
        -recurse) recurse="YES";;
        -*) echo "Unknown option: $1"; exit 1;;
        *) break;;
  esac
  shift
done

if [ $recurse = "YES" ]; then
  echo "looking in subdirs..."
  files=`find . -name "$*"`
else
  files=$*
fi

#echo "$options -file $file -after $after -befor $befor"

#-------------------------------------------------------------
# extract information from files
#-------------------------------------------------------------

for file in $files
do
  faux=$(echo $file | sed '/\/tmp\//'d)
  [ -z "$faux" ] && continue

  getheader.pl $file | revisionlog.pl \
		$options -file $file -after $after -befor $befor
done

#-------------------------------------------------------------
# end of routine
#-------------------------------------------------------------

