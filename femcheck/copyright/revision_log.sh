#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# check if files have revision log
#
#----------------------------------------------

shydir=$HOME/shyfem
copydir=$shydir/femcheck/copyright
 
write=NO

while [ -n "$1" ]
do
   case $1 in
        -write)   write=YES;;
        -*)       option="$option $1";;
	*)        break
  esac
  shift
done
[ -z "$option" ] && option="-check"
#echo "option: $option"

for file
do
  newfile=$file.new
  $copydir/revision_log.pl $option $file
  status=$?
  cmp $file $newfile
  changed=$?
  if [ $changed -ne 0 ]; then
    echo "*** $file has been changed..."
  else
    if [ $status -eq 0 ]; then
      echo "$file has no revision log"
    fi
    rm -f $newfile
  fi
  if [ $changed -ne 0 -a $write = YES ]; then
    echo "$file has been written"
    mv -f $newfile $file
  fi
done

#----------------------------------------------

