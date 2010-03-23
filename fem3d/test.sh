#!/bin/sh

date2num()
{
#  echo $1 | sed 's/\(.*\)-\(.*\)-\(.*\)/\3\2\1/'
  echo $1 | sed 's/\(.*\)-\(.*\)-\(.*\)/\3\2\1/'
}

tmpfile=tmp0.tmp
tmpfile1=tmp1.tmp
tmpfile2=tmp2.tmp

firstline=`head -1 VERSION | sed 's/ \{1,\}/ /g'`
actdate=`echo $firstline | cut -d" " -f3`
comparedate=`date2num $actdate`

revisionlog -after $comparedate *.[cfFh] > $tmpfile
if [ -f VERSION ]; then
  revisionlog_adjust.pl $tmpfile VERSION > $tmpfile2
fi
echo "============================="
cat $tmpfile
echo "============================="
cat $tmpfile2
echo "============================="

