#!/bin/sh
#
# compares two regression runs

#----------------------------------------------------------

Usage()
{
  echo "Usage: regdiff.sh version1 [version2]"
  echo "   example: regress.sh 4.74 4.75"
  echo "   regression.sh must have already been run on those versions"
  echo "   if version2 is ommitted the actual version is used"
}

FullUsage()
{
  Usage
}

ErrorOption()
{
  echo "no such option: $1"
}

#----------------------------------------------------------

tests="test1 test2"
tests="test3"
tests="test1 test2 test3 test4"

#----------------------------------------------------------

if [ $# -le 0 ]; then
  Usage
  exit 1
elif [ $# -eq 1 ]; then
  version1=$1
  version2=actual
else
  version1=$1
  version2=$2
fi


aux1=`echo $version1 | sed -e 's/\./_/'`
aux2=`echo $version2 | sed -e 's/\./_/'`
test1dir=test_$aux1
test2dir=test_$aux2

echo "using directories $test1dir and $test2dir"

if [ ! -d $test1dir ]; then
  echo "no such directory: $test1dir"
  exit 1
fi

if [ ! -d $test2dir ]; then
  echo "no such directory: $test2dir"
  exit 1
fi

#----------------------------------------------------------

cd $test1dir

for test in $tests
do

if [ ! -d $test ]; then
  echo "** no such directory: $test ... skipping"
  continue
fi

cd $test
compdir=../../$test2dir/$test

if [ ! -d $compdir ]; then
  echo "** no such directory: $compdir ... skipping"
  continue
fi


echo "actual directory: $PWD"
echo "comparing $test between here and $compdir"

echo "difference in z ..."
diffs -d -n $compdir z.*

echo "difference in u ..."
diffs -d -n $compdir u.*

echo "difference in v ..."
diffs -d -n $compdir v.*

echo "difference in m ..."
diffs -d -n $compdir m.*

#echo "difference in .inf ..."
#diff  $compdir/$test.inf $test.inf

#echo "difference in .log ..."
#diff  $compdir/$test.log $test.log

vers1=`head -20 $test.log | grep version | grep FEM`
vers2=`head -20 $compdir/$test.log | grep version | grep FEM`

echo  "version1: $vers1   directory: $test1dir"
echo  "version2: $vers2   directory: $test2dir"

cd ..

done

#----------------------------------------------------------






