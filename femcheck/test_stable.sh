#!/bin/sh


commit=`gittags | tail -1`
commit_number=`echo $commit | cut -d ' ' -f 1`
echo "commit_line  : $commit"
echo "commit number: $commit_number"

name=`gittar -info $commit_number`
if [ $? -ne 0 ]; then
  echo $name
  exit 1
fi

echo "using name: $name"
rm -f $name.tar.gz

gittar $commit_number
mv -f $name.tar.gz tmp
cd tmp

exit

archive=`ls -1tr ../arc/*_stable.tar.gz | tail -1`
dir=`basename $archive .tar.gz`

actdir=`pwd`
echo "We are in dir: $actdir"
echo "Using archive: $archive"
echo "Using directory: $dir"

tar xvzf $archive > /dev/null

cd $dir
version=`head -1 VERSION`
echo "Using version: $version"
actdir=`pwd`
echo "Compiling in: $actdir"

make test_compile

cd ..
actdir=`pwd`
echo "We are in dir: $actdir"

