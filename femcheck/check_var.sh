#!/bin/bash
#
# checks various things
#
#--------------------------------------------------


CheckExe()
{
  pushd $1 > /dev/null
  echo "checking directory `pwd`"

  files=`ls`

  for file in $files
  do
    if [ -x $file ]; then
      :
    elif [ -d $file ]; then
      :
    elif [ $file = CR -o $file = Makefile ]; then
      :
    elif [ $file = INSTALL-LIST -o $file = README ]; then
      :
    elif [ $file = Rules.dist ]; then
      :
    else
      echo "*** file is not executable: $file"
    fi
  done

  popd > /dev/null
}

#pwd

CheckExe fembin
CheckExe femcheck

