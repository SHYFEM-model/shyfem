#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: install-all.sh [-info|-install]"
  exit 0
fi

option=$1
if [ $option = "-info" ]; then
  install="NO"
elif [ $option = "-install" ]; then
  install="YES"
else
  echo "unknown option: $option"
  exit 1
fi

pcks=`cat INSTALL-LIST`

for pck in $pcks
do
  stat=`dpkg -s $pck  2>&1 | grep '^Status'`
  #echo "$pck: $stat"
  if [ -z "$stat" ]; then
    echo "*** package $pck is not known"
    if [ $install = "YES" ]; then
      echo "trying to install anyway ..."
      apt-get install $pck
    fi
  elif [ "$stat" = "Status: install ok installed" ]; then
    echo "package $pck is installed"
  else
    echo -n "*** package $pck is not installed"
    if [ $install = "YES" ]; then
      echo " ... installing"
      apt-get install $pck
    else
      echo ""
    fi
  fi
done

