#!/bin/sh
#
# makes a tarball of last committed version
# substitutes stable versions of files
# unzips and test compiles this version
#
#-----------------------------------------------------

femdir=$SHYFEMDIR

commit_number=`gittags | tail -1 | cut -d " " -f 1`

echo '============================================'
echo ' preparing stable version from last commit'
echo '============================================'

echo "commit_number: $commit_number"

name=`gittar -info $commit_number`
if [ $? -ne 0 ]; then
  echo "error: $name"
  exit 1
fi

echo "using name for tar: $name"
rm -f $name.tar.gz tmp/$name.tar.gz

gittar $commit_number
mv -f $name.tar.gz tmp
cd tmp

$femdir/stable/make_stable.sh $name.tar.gz

echo '============================================'
echo ' unpacking stable distribution '
echo '============================================'

archive=`ls -1tr ./*_stable.tar.gz | tail -1`
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

echo '============================================'
echo ' test compiling the distribution'
echo '============================================'

make test_compile

echo '============================================'
echo ' end of test compiling the distribution'
echo '============================================'

cd ..
actdir=`pwd`
echo "Test compilation was done in: $actdir"

#-----------------------------------------------------

