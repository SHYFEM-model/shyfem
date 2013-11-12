#!/bin/sh

cd tmp

archive=`ls -1tr ../arc/*_stable.tar.gz | tail -1`
dir=`basename $archive .tar.gz`

actdir=`pwd`
echo "We are in dir: $actdir"
echo "Using archive: $archive"
echo "Using directory: $dir"

tar xvzf $archive > /dev/null

cd $dir
#version=`make version`
version=`head -1 VERSION`
echo "Using version: $version"
actdir=`pwd`
echo "Compiling in: $actdir"

make test_compile

cd ..
actdir=`pwd`
echo "We are in dir: $actdir"

#echo "deleting directory: $dir"
#rm -rf $dir

