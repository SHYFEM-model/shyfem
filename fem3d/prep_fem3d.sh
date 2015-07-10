#!/bin/sh
#
#----------------------------------------

files=$*

act=`pwd`
echo $act

dist="/home/georg/fem/fem3d"
handle=/home/georg/fem/fem3d/bin/handle_common.sh

if [ $act = $dist ]; then
  echo "we are in distribution ... aborting"
  exit 1
fi

$handle -inc2use $files
[ $? -ne 0 ] && echo "error running handle_common.sh" && exit 1
$handle -new
$handle -clean

rm -f newdum.f newdum1.f
touch newdum.f newdum1.f

if [ -f subnev_new.f ]; then
  mv -f subnev.f subnev_old.f
  mv -f subnev_new.f subnev.f
fi

make depend

