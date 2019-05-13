#!/bin/sh
#
# check if files have revision log
#
#----------------------------------------------

for file
do
  grep -i "revision log :" $file > /dev/null
  status=$?
  if [ $status -eq 0 ]; then
    true
    #echo "file has revision log: $file"
  else
    true
    echo "*** file has no revision log: $file"
  fi
done

#----------------------------------------------

