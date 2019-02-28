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
    if [ -s tmp.tmp ]; then
      echo $file
      cat tmp.tmp
      mv new.tmp $file.new
      #$copydir/count_developers.pl -subst $file    > new.tmp
    else
      #echo $file
      echo "file $file has no old revision..."
    fi
  done
}

#---------------------------------------------------------

  SubstFiles $*

#---------------------------------------------------------

