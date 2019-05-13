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

what=check
if [ $# -eq 0 ]; then
  what=check
elif [ $1 = "-check" ]; then
  what=check
  shift
elif [ $1 = "-change" ]; then
  what=change
  shift
else
  echo "unknown option: $1"
  exit 1
fi

#---------------------------------------------------------

FindRevisionLog()
{
  files=$*

  for file in $files
  do
    [ -f tmp.tmp ] && rm -f tmp.tmp
    $copydir/revise_revision_log.pl $file  > $file.new  2> tmp.tmp
    status=$?
    if [ -s tmp.tmp ]; then	#some changes have been found or error
      echo $file
      cat tmp.tmp
      if [ $what = "change" ]; then
        mv $file $file.old
        mv $file.new $file
      elif [ $what = "check" ]; then
	rm -f $file.new
      fi
      #$copydir/count_developers.pl -subst $file    > new.tmp
    elif [ $status -eq 0 ]; then
      #echo $file
      echo "file $file has no revision log"
      true
    fi
  done
}

SubstRevisionLog()
{
  files=$*

  for file in $files
  do
    [ -f tmp.tmp ] && rm -f tmp.tmp
    $copydir/revise_revision_log.pl -extract $file  > $file.new  2> tmp.tmp
    status=$?
    git-file -revlog $file > aux.tmp
    sed -n '/revision log :/,$p' aux.tmp | tail -n +2 > revlog.git.tmp
  done
}

#---------------------------------------------------------

FindRevisionLog $*
#SubstRevisionLog $*

#---------------------------------------------------------

