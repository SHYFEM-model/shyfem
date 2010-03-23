#!/bin/sh

for file
do
  echo "checking dependencies in $file"

  procs=`get_procs.pl $file`
  echo "---------------------------------------------------"
  echo $procs
  echo "---------------------------------------------------"

  for proc in $procs
  do
    echo "  -----------------------------------------------"
    echo "  checking proceedure $proc"
    echo "  -----------------------------------------------"
    token -i $proc *.f *.h
  done
done
