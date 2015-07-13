#!/bin/sh
#
#----------------------------------------

force="NO"
if [ "$1" = "-force" ]; then
  force="YES"
  shift
fi

files=$*

act=`pwd`
echo $act

dist="/home/georg/fem/femplot"
handle=/home/georg/fem/fem3d/bin/handle_common.sh

if [ $force != "YES" ]; then
 if [ $act = $dist ]; then
  echo "we are in distribution ... aborting"
  exit 1
 fi
fi

$handle -inc2use $files
[ $? -ne 0 ] && echo "error running handle_common.sh" && exit 1
$handle -new
$handle -clean

#rm -f supdum.f
#touch supdum.f

#if [ -f subnev_new.f ]; then
#  mv -f subnev.f subnev_old.f
#  mv -f subnev_new.f subnev.f
#fi

make depend

