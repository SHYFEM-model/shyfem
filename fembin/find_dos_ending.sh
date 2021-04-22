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
# finds occurence of dos end of line in files
#
#------------------------------------------------------------------------


if [ $# -eq 0 ]; then
  echo "Usage: find_dos_ending.sh [-detail] file(*)"
  exit 1
fi

option="-r -c"

if [ $1 = "-detail" ]; then
  option="-r"
  shift
fi

grep $option -P '\x0d' $* | grep -v ":0$" | grep -v "^Binary file"

