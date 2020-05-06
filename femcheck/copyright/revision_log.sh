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
stats=NO
gui=NO

while [ -n "$1" ]
do
   case $1 in
        -write)   write=YES;;
        -stats)   stats=YES; option="$option $1";;
        -gui)     gui=YES;;
        -*)       option="$option $1";;
	*)        break
  esac
  shift
done
[ -z "$option" ] && option="-check"
#echo "option: $option"

for file
do
  [ -d $file ] && continue
  newfile=$file.new
  $copydir/revision_log.pl $option $file
  status=$?
  changed=0
  if [ -f $newfile ]; then
    cmp $file $newfile > /dev/null
    changed=$?
  fi
  if [ $changed -ne 0 ]; then
    echo "   $file has been changed..."
    if [ $gui = YES ]; then
      tkdiff $file $newfile
    fi
  else
    if [ $stats = NO -a $status -eq 0 ]; then
      echo "   $file has no revision log"
    fi
    rm -f $newfile
  fi
  if [ $changed -ne 0 -a $write = YES ]; then
    echo "   $file has been written"
    mv -f $newfile $file
  fi
done

#----------------------------------------------

