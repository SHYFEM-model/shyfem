#!/bin/sh
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
#-----------------------------------------------------------------------
#
# still needed:
#
# -integrate
# -combine
# -update_copyright
#
#-----------------------------------------------------------------------

copydir=$HOME/shyfem/femcheck/copyright

#-----------------------------------------------------------------------

GetGitLog()
{
  local file=$1

  git-file -revlog  $file > aux.tmp
  sed -n '/revision log :/,$p' aux.tmp | tail -n +3 | head --lines=-1 \
      		> revlog_git.tmp
}

#-----------------------------------------------------------------------

Usage()
{
  echo "Usage: check_file.sh [-h|-help] [options] file(s)"
}

FullUsage()
{
  Usage
  echo "  deal with revision logs and copyright"
  echo ""
  echo "  -h|-help      this help screen"
  echo "  -gitrevlog    integrate git revision log if revision log is missing"
  echo "  -update       *** combines two revision logs"
  echo "  -change       really change in the file"
  echo "  -diff         compares .old and .new files"
  echo "  -clean        cleans dir from auxiliary files"
  echo "  -debug        writes more information"
  echo ""
}

#------------------------------------------------------------

options="-warn -obsolete"
options="-obsolete"
options="-warn"

gitrevlog="NO"
update="NO"
change="NO"
diff="NO"
clean="NO"
debug="NO"

while [ -n "$1" ]
do
   case $1 in
        -h|-help)       FullUsage; exit 0;;
        -gitrevlog)     gitrevlog="YES"; options="$options -extract";;
        -update)        update="YES"; options="$options -extract";;
        -change)        change="YES";;
        -diff)          diff="YES";;
        -clean)         clean="YES";;
        -debug)         debug="YES";;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

[ $# -eq 0 ] && Usage

#------------------------------------------------------------

MakeDiff()
{
  for file
  do
    echo "========= $file =============================================="
    [ -f $file.old ] && tkdiff -w $file.old $file.new
    #[ -f $file.old ] && diff -w $file.old $file.new
  done
}

MakeList()
{
  for file
  do
    echo $file
  done
}

MakeClean()
{
  for file
  do
    echo $file
    [ -f $file.old ] && rm $file.old
    [ -f $file.new ] && rm $file.new
  done
  rm -f *.tmp
}

#------------------------------------------------------------

MakeCheck()
{
 for file
 do
  [ -f $file ] || continue
  $copydir/check_file.pl $options $file
  status=$?
  if [ $status -eq 0 ]; then
    true
    echo "*** file has no revision log: $file"
    if [ "$gitrevlog" = "YES" ]; then
      echo "...integrating revision log from git..."
      GetGitLog $file
      $copydir/check_file.pl -integrate=revlog_git.tmp $file > $file.new
      cp $file $file.old
      [ $change = "YES" ] && cp $file.new $file
    fi
  elif [ $status -eq 2 ]; then
    echo "file has manual copyright: $file"
  elif [ $update = "YES" ]; then
      echo "...updating revision log from git: $file"
      GetGitLog $file
      $copydir/check_file.pl -combine revlog.tmp revlog_git.tmp
      $copydir/check_file.pl -substitute=revlog_new.tmp $file > $file.new
      cp $file $file.old
      wold=$( wc -l revlog.tmp | sed -e 's/ .*//' )
      wnew=$( wc -l revlog_new.tmp | sed -e 's/ .*//' )
      if [ $change = "YES" ]; then
        cmp revlog.tmp revlog_new.tmp > /dev/null 2>&1
	if [ $? -eq 0 ]; then
	  echo "  no changes in revlog: $wold $wnew"
	  rm -f $file.old $file.new
	else
          echo "  changes in revlog: $wold $wnew ... copying"
          cp $file.new $file
	fi
      else
        echo "  changes in revlog: $wold $wnew"
      fi
  elif [ $status -eq 1 ]; then
    true
    [ $debug = "YES" ] && echo "file has revision log: $file"
  else
    echo "*** unknown resturn status: $file"
  fi
 done
}

#-----------------------------------------------------------------------

if [ $diff = "YES" ]; then
  MakeDiff $*
elif [ $clean = "YES" ]; then
  MakeClean $*
else
  MakeCheck $*
fi

#-----------------------------------------------------------------------

