#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# changes paralution source code to be compiled with SHYFEM
#
#----------------------------------------------------------------

parascripts=$( pwd )
paradir=$1

#----------------------------------------------------------------

if [ ! -d $paradir ]; then
  echo "*** cannot find directory $paradir"
  exit 1
fi

cd $paradir/src/utils

file=def.hpp
cp $file $file.tmp
cat $file.tmp | sed -e 's/VERBOSE_LEVEL 2/VERBOSE_LEVEL 0/' > $file
echo "file $file adjusted..."

file=allocate_free.cpp
cp $file $file.tmp
cat $file.tmp | sed -e 's/false/NULL/' > $file
echo "file $file adjusted..."

file=log.hpp
cp $file $file.tmp
$parascripts/log_info.pl $paradir/src/utils/$file.tmp > $paradir/src/utils/$file
echo "file $file adjusted..."

cd $paradir/src/solvers

file=iter_ctrl.cpp
cp $file $file.tmp
cat $file.tmp | sed -e 's/this->verb_ = 1/this->verb_ = 0/' > $file
echo "file $file adjusted..."

#----------------------------------------------------------------

cp $paradir/src/Makefile $paradir/src/Makefile.old

dir1=$parascripts/../fembin
file1=$parascripts/../Rules.make

cat $parascripts/Makefile.src | sed -e "s=DIR1=$dir1=" | \
			        sed -e "s=FILE1=$file1=" \
				> $paradir/src/Makefile

#----------------------------------------------------------------

