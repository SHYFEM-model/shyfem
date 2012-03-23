#!/bin/sh

for name in $(git diff --name-only $1)
do 
  echo "============================================================"
  echo "        $name                                               "
  echo "============================================================"
  git --no-pager diff $1 $name
done

