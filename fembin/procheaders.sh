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
# processes headers from fortran files

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
export FEMDIR

###############################################################

usage()
{
  echo Usage: procheaders [-help] [-what] files
}

help()
{
  usage
  echo "  -help          show this help screen"
  echo "  -keywords      show keywords used in headers"
}

###############################################################

keywords()
{
rm -f $tmpfile

for file in $files
do
  getheader.pl $file | keyword.pl >> $tmpfile
done

sort -u $tmpfile
rm -f $tmpfile
}

revisions()
{
rm -f $tmpfile

for file in $files
do
#  echo "" >> $tmpfile
#  echo "file: $file" >> $tmpfile
#  echo "" >> $tmpfile
  getheader.pl $file | revisionlog.pl -file $file -after 19980000 >> $tmpfile
done

}

###############################################################

tmpfile=tmp.tmp

if [ $# -eq 0 ]; then
  usage
  exit 0
fi

while [ $# -gt 0 ]
do

  opt=$1
  shift

  case $opt in
        -h|-help) help; exit 0;;
        -k|-key|-keywords) what=keywords; break;;
        -r|-rev|-revisions) what=revisions; break;;
        -*) echo "Unknown option: $opt"; exit 1;;
        *) echo "No action specified"; exit 2;;
  esac

done

files=$*

case $what in
      keywords) keywords;;
      revisions) revisions;;
      *) echo "Unknown action"; exit 2;;
esac

