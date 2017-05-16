#!/bin/sh

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

