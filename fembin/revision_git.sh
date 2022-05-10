#!/bin/bash
#
# shows revision log of files changed in git
#
# in nshow is maximum number of lines to check in revision log
#
#------------------------------------------------

nshow=5

IFS=$'\n'

#------------------------------------------------

CheckDate()
{
  for line in $revlog
  do
    origdate=$( echo $line | sed -E 's/^[!c] *//' | sed -E 's/\s+.*// ' )
    [ -z "$origdate" ] && continue
    #date=$( echo $date | sed -e 's/\./-/g' )
    day=$( echo $origdate | cut -f 1 -d '.' )
    month=$( echo $origdate | cut -f 2 -d '.' )
    year=$( echo $origdate | cut -f 3 -d '.' )
    date="$year-$month-$day"
    sdate=$( date --date=$date +%s )
    if [ $scdate -le $sdate ]; then
      [ -n "$before" ] && echo "$before"
      before=""
      echo "$line"
    else
      before=$line
    fi
    #echo "$date - $origdate - $sdate"
  done
}

#------------------------------------------------
# get last commit date
#------------------------------------------------

last_commit=$( git-tags | tail -1 )
cdate=$( echo $last_commit | sed -E 's/^\w+ //' | sed -e 's/ .*//' )
#echo "last commit: $cdate"
scdate=$( date --date=$cdate +%s )
echo "last commit date: $cdate - $scdate"

#------------------------------------------------
# get changed files
#------------------------------------------------

changed=$( git s | grep modified: )

#------------------------------------------------
# check revision log in changed files
#------------------------------------------------

for line in $changed
do
  file=$( echo $line | sed -e 's/.*modified: *//' )
  echo "----------------------"
  echo $file
  echo "----------------------"
  revlog=$( revisionlog.pl $file | tail -$nshow )
  CheckDate
done

#------------------------------------------------
# end of routine
#------------------------------------------------

