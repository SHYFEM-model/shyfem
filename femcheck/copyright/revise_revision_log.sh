#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#---------------------------------------------------------

shydir="$HOME/shyfem"
copydir="$shydir/femcheck/copyright"
copyfile="$copydir/copy.log"

#---------------------------------------------------------

SubstFiles()
{
  files=$*

  for file in $files
  do
    [ -f tmp.tmp ] && rm -f tmp.tmp
    $copydir/revise_revision_log.pl $file  > new.tmp  2> tmp.tmp
    status=$?
    if [ -s tmp.tmp ]; then	#some chages have been found or error
      echo $file
      cat tmp.tmp
      mv $file $file.old
      mv new.tmp $file
      #$copydir/count_developers.pl -subst $file    > new.tmp
    elif [ $status -eq 0 ]; then
      #echo $file
      echo "file $file has no revision log"
      true
    fi
  done
}

#---------------------------------------------------------

  SubstFiles $*

#---------------------------------------------------------

