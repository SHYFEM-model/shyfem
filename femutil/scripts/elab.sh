#!/bin/sh

for file
do
  echo $file
  ./elab.pl $file
done

