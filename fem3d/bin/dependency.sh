#!/bin/sh

for file
do
  echo $file
  f77 $file 2> ggg
  dependency.pl ggg
  [ $? -ne 0 ] && exit 1
done

