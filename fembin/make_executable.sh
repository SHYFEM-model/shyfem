#!/bin/sh

# make scripts executable

dirs="femanim femdoc grid femcheck femplot fem3d fem3d/bin \
	femersem femspline"

CHX()
{
  dir=$1
  echo "making scripts executable in directory $dir"
  shift
  chmod +x $* 2> /dev/null
  #chmod +x $*
}

for dir in $dirs
do
  CHX $dir $dir/*.sh $dir/*.pl
done

dir=femlib/perl
CHX $dir $dir/*.sh $dir/*.pl $dir/*.pm

dir=femlib/python
CHX $dir $dir/*.sh $dir/*.py

exit 0

