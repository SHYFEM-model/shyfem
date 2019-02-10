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
# checks various things
#
#--------------------------------------------------

except="Makefile INSTALL-LIST README copyright_notice.txt"

CheckExe()
{
  pushd $1 > /dev/null
  echo "checking directory `pwd`"

  files=`ls`

  for file in $files
  do
    [ -x $file ] && continue
    [ -d $file ] && continue
    is_special="NO"
    for special in $except
    do
      [ $file = $special ] && is_special="YES"
    done
    if [ $is_special = "NO" ]; then
      echo "*** file is not executable: $file"
    fi
  done

  popd > /dev/null
}

#pwd

CheckExe fembin
CheckExe femcheck

