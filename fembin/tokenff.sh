#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

nocom=NO

while [ -n "$1" ]
do
  case $1 in
	-nocom)		nocom=YES;;
	*)		break;;
  esac
  shift
done


if [ $nocom = YES ]; then
  token $* | grep -vi ': c' | grep -v '\s*!'
else
  token $*
fi


