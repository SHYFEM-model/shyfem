#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# installs shyfem through hard references
#
#-------------------------------------------------

if [ ! -f VERSION ]; then
  echo "this routine must be run in main fem directory... aborting"
  exit 1
fi

if [ "$1" = "-reset" ]; then
  option="-reset"
  extension="ggu_reset"
  echo "resetting: $option"
else
  option="-install"
  extension="ggu_hard"
  echo "installing: $option"
fi

dir=`pwd`

echo "========================================================="
echo "running shyfem_install_hard.sh"
echo "      using directory: $dir"
echo "      using option:    $option"
echo "========================================================="

cd fembin

files=`ls -d *`

for file in $files
do
  ./shyfem_install_hard.pl $option $dir $file > tmp.tmp
  status=$?
  #echo "status: $status  $file"
  if [ $status -eq 1 ]; then
    new=$file.$extension	# for test
    new=$file			# for real
    echo "changed... creating new file version: $file $new"
    touch -r $file tmp.tmp
    chmod +x tmp.tmp
    mv -f tmp.tmp $new
  elif [ $status -eq 0 ]; then	# no change
    rm -f tmp.tmp
  else
    echo "*** unknown error code... aborting...  $status"
    exit 1
  fi
done

