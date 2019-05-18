#!/bin/sh
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
  echo "  -change       really change in the file"
  echo ""
}

#------------------------------------------------------------

options="-warn -obsolete"
options="-obsolete"
options="-warn"

gitrevlog="NO"
change="NO"
diff="NO"

while [ -n "$1" ]
do
   case $1 in
        -h|-help)       FullUsage; exit 0;;
        -gitrevlog)     gitrevlog="YES"; options="$options -extract";;
        -change)        change="YES";;
        -diff)          diff="YES";;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

[ $# -eq 0 ] && Usage

#------------------------------------------------------------

if [ $diff = "YES" ]; then
  for file
  do
    [ -f $file.old ] && tkdiff $file.old $file.new
  done
  exit 0
fi

#------------------------------------------------------------

for file
do
  $copydir/check_file.pl $options $file
  status=$?
  if [ $status -eq 0 ]; then
    true
    echo "*** file has no revision log: $file"
    if [ "$gitrevlog" = "YES" ]; then
      echo "  ...integrating revision log from git..."
      git-file -revlog  $file > aux.tmp
      sed -n '/revision log :/,$p' aux.tmp | tail -n +3 | head --lines=-1 \
      		> revlog_git.tmp
      $copydir/check_file.pl -integrate=revlog_git.tmp $file > $file.new
      cp $file $file.old
      [ $change = "YES" ] && cp $file.new $file
    fi
  elif [ $status -eq 2 ]; then
    echo "file has manual copyright: $file"
  else
    true
    #echo "file has revision log: $file"
  fi
done

#-----------------------------------------------------------------------

