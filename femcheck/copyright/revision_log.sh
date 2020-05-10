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
#------------------------------------------------------------------------

shydir=$HOME/shyfem
copydir=$shydir/femcheck/copyright
 
write=NO
keep=NO
gui=NO

while [ -n "$1" ]
do
   case $1 in
        -write)   write=YES;;
        -keep)    keep=YES;;
        -gui)     gui=YES;;
        -*)       option="$option $1";;
	*)        break
  esac
  shift
done
[ -z "$option" ] && option="-check"

for file
do
  [ -d $file ] && continue
  [ -L $file ] && continue
  newfile=$file.new
  #echo "revision_log.sh: treating file $file ($newfile)"
  [ -f $newfile ] && rm -f $newfile
  $copydir/revision_log.pl $option $file
  changed=0
  if [ -f $newfile ]; then
    cmp $file $newfile > /dev/null 2>&1
    changed=$?
  fi
  if [ $changed -ne 0 ]; then
    echo "    $file has been changed..."
    if [ $gui = YES ]; then
      tkdiff $file $newfile
    fi
  else
    [ -f $newfile ] && rm -f $newfile
  fi
  if [ $changed -ne 0 -a -f $newfile ]; then
    if [ $write = YES ]; then
      echo "    $file has been written"
      mv -f $newfile $file
    elif [ $keep = NO ]; then
      rm -f $newfile
    fi
  fi
done

#------------------------------------------------------------------------

