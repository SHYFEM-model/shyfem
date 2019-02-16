#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

echo "removing symbolic links..."

files=$( ls )

for file in $files
do
  if [ -h $file ]; then
    echo "removing symbolic link $file"
    #ls -la $file
    rm -f $file
  fi  
done

