#!/bin/sh

act=`pwd`
echo $act

dist="/home/georg/fem/fem3d"

if [ $act = $dist ]; then
  echo "we are in distribution ... aborting"
  exit 1
fi

bin/handle_common.sh -inc2use
bin/handle_common.sh -new
bin/handle_common.sh -clean

rm -f newdum.f
touch newdum.f

make depend

