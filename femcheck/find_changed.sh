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
# finds changed files ignoring compiler generated ones
#
#---------------------------------------------------------

files=$(find . -newer VERSION -type f | 
	grep -v '/arc/' | \
	grep -v '/.git/' | \
	grep -v '/femlib/mod/' | \
	grep -v '__genmod\.f90' | \
	grep -v '\.o' | \
	grep -v '\.a' | \
	grep -v '\.swp' | \
	grep -v 'param\.h' | \
	grep -v 'CHECKLOG\.tmp' | \
	grep -v '\.mod')

for file in $files
do
  type=$(file $file | sed -e 's/^.*: //' | sed -e 's/, x86-64.*$//')
  #echo $type
  if [ "$type" = "ELF 64-bit LSB executable" ]; then
    :
  else
    #echo "$file $type"
    echo "$file"
  fi
done

