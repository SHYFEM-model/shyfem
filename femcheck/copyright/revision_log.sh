#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
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

changes=0

for file
do
  [ -d $file ] && continue
  [ -L $file ] && continue
  newfile=$file.revnew
  #echo "revision_log.sh: treating file $file ($newfile)"
  [ -f $newfile ] && rm -f $newfile
  $copydir/revision_log.pl $option $file
  if [ $? != 0 ]; then
    echo "    error running revision_log.pl on file $file"
    [ $keep = NO ] && rm -f $newfile
    continue
  fi
  changed=0
  if [ -f $newfile ]; then
    cmp $file $newfile > /dev/null 2>&1
    changed=$?
  fi
  if [ $changed -ne 0 ]; then
    echo "    $file is changed..."
    changes=$(( changes + 1 ))
    if [ $gui = YES ]; then
      tkdiff $file $newfile
    fi
  else
    [ -f $newfile ] && rm -f $newfile
  fi
  if [ $changed -ne 0 -a -f $newfile ]; then
    if [ $write = YES ]; then
      type=$( $copydir/find_file_type.pl $file )
      if [ -z "$type" ]; then
        echo "    cannot determine type of file $file"
        [ $keep = NO ] && rm -f $newfile
      else
        mv -f $newfile $file
        [ $type = "script" ] && chmod +x $file
        echo "    $file has been written"
      fi
    elif [ $keep = NO ]; then
      rm -f $newfile
    fi
  fi
done

if [ $changes -ne 0 ]; then
  if [ $write = YES ]; then
    echo "$changes file(s) have been written"
  else
    echo "$changes file(s) are changed"
    echo "but changes have not been written to file"
    echo "use --write to really change files"
  fi
fi

#------------------------------------------------------------------------

