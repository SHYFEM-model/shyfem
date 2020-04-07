#!/bin/sh
#
# finds netcdf libraries
#
#-----------------------------------------------------

debug="YES"
debug="NO"

compiler=$1
usrdir=$2

file=libnetcdff.a
dirs="$usrdir /usr /opt/sw/netcdf /usr/local/netcdf /usr/local"

if [ $compiler = "INTEL" ]; then
  dirs="$dirs /usr/local/intel"
fi

if [ $debug = "YES" ]; then
  echo "vars - $compiler $usrdir - $dirs"
  exit 0
fi

for dir in $dirs
do
  if [ -f $dir/lib/$file ]; then
    echo $dir/lib
    exit 0
  fi
done

dir=/usr/lib/x86_64-linux-gnu
if [ -f $dir/$file ]; then
  echo $dir
  exit 0
fi

echo "unknown"

