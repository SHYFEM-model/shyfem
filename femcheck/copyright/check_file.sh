#!/bin/sh
#
# check if files have revision log
#
#----------------------------------------------
#
# still needed:
#
# -integrate
# -combine
# -update_copyright
#
#----------------------------------------------

copydir=$HOME/shyfem/femcheck/copyright

norevlog=
norevlog=1

gitrevlog=
gitrevlog=1

options="-warn -obsolete"
options="-warn"
options="-obsolete"
options=
options="-extract"

for file
do
  $copydir/check_file.pl $options $file
  status=$?
  if [ $status -eq 0 ]; then
    true
    [ -n "$norevlog" ] && echo "*** file has no revision log: $file"
    if [ -n "$gitrevlog" ]; then
      echo "  ...integrating revision log from git..."
      git-file -revlog  $file > aux.tmp
      sed -n '/revision log :/,$p' aux.tmp | tail -n +3 | head --lines=-1 \
      		> revlog_git.tmp
      $copydir/check_file.pl -integrate=revlog_git.tmp $file > $file.new
      cp $file $file.old
      #cp $file.new $file
    fi
  else
    true
    #echo "file has revision log: $file"
  fi
done

#----------------------------------------------

