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
# para_get.sh - gets paralution code
#
#----------------------------------------------------------------

parafile=http://www.paralution.com/downloads/paralution-1.1.0.zip

#----------------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "para_setup.sh needs install dir... no directory given"
  exit 1
fi

date=$( date +%Y-%m-%dT%H_%M_%S )
actdir=$( pwd )
paradir=$1
echo "downloading paralution code to $paradir"

#----------------------------------------------
# create temporary directory
#----------------------------------------------

tmpdir=tmp-$date-$$
mkdir -p $tmpdir

#----------------------------------------------
# get paralution code
#----------------------------------------------

cd $tmpdir
wget -q $parafile
if [ $? -ne 0 ]; then
  echo "error downloading paralution code... exiting"
  cd $actdir
  rm -rf $tmpdir
  exit 3
fi

#----------------------------------------------
# if $paradir exists move it out of the way
#----------------------------------------------

if [ -e $paradir ]; then
  paranew=$paradir-$date-$$
  mv $paradir $paranew
  echo "$paradir exists... moving to $paranew"
fi

#----------------------------------------------
# unzip archive, move to $paradir, and clean up
#----------------------------------------------

unzip -q paralution-1.1.0.zip
mv paralution-1.1.0 $paradir
cd $actdir
rm -rf $tmpdir
echo "paralution code has been downloaded to $paradir"

#----------------------------------------------
# end of routine
#----------------------------------------------

