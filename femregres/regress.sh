#!/bin/sh
#
# exports version of fem3d and runs regression test

#----------------------------------------------------------

Usage()
{
  echo "Usage: regress.sh [-e|-c|-t|-r|-p|-a] [-h|-help] version"
  echo "   example: regress.sh -a 4.74"
  echo "   the special keyword actual is recognized as the actual version"
}

FullUsage()
{
  Usage
  echo "   -e   export"
  echo "   -c   compile"
  echo "   -t   make test"
  echo "   -r   run test"
  echo "   -p   post processing"
  echo "   -a   all of the above"
}

ErrorOption()
{
  echo "no such option: $1"
}

#----------------------------------------------------------

tests="test1 test2"
tests="test3"
tests="test2"
tests="test1 test2 test3"
tests="test4"
tests="test1"
tests="test1 test2 test3 test4"

#----------------------------------------------------------

export="NO"
compile="NO"
test="NO"
run="NO"
post="NO"
all="NO"

while [ -n "$1" ]
do
   case $1 in
        -e)         export="YES";;
        -c)         compile="YES";;
        -t)         test="YES";;
        -r)         run="YES";;
        -p)         post="YES";;
        -a)         all="YES";;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   shift
done

if [ $all = "YES" ]; then
  export="YES"
  compile="YES"
  test="YES"
  run="YES"
  post="YES"
fi

#----------------------------------------------------------

if [ $# -eq 0 ]; then
  Usage
  exit 1
fi

version=$1
aux=`echo $version | sed -e 's/\./_/'`
tag=VERS_$aux
act=`pwd`
femdir=$act/fem3d_$aux
testdir=$act/test_$aux

basin=venlag62

if [ $version = "actual" ]; then
  aux=actual
  femdir=$HOME/fem/fem3d
  testdir=test_actual
  export="NO"
fi

echo $version $aux $tag $femdir $testdir
echo $basin $simul

#----------------------------------------------------------

if [ $export = "YES" ]; then
  echo "...exporting"
  if [ $femdir = "$HOME/fem/fem3d" ]; then
    echo "Cannot remove $femdir"
    exit 1
  else
    rm -rf $femdir
  fi
  cvs export -r $tag -d $femdir fem3d
  status=$?
  if [ $status -ne 0 ]; then
    echo "Status: $status"
    echo "No such tag $tag ... exiting"
    exit 1
  fi
fi

#----------------------------------------

if [ $compile = "YES" ]; then
  echo "...compiling in $femdir"
  if [ $version != "actual" ]; then
    cmp skel/param.h $femdir/param.h
    status=$?
    echo "comparing param.h ... $status"
    if [ $status -ne 0 ]; then
      echo "copying param.h from skel to source..."
      cp skel/param.h $femdir
    fi
  fi
  cd $femdir
  make fem
  cd $act
fi

#----------------------------------------

if [ $test = "YES" ]; then
  echo "...making test in $testdir"
  #rm -rf $testdir
  mkdir -p $testdir
  cd $testdir
  for test in $tests
  do
    echo "making directory $test in $testdir"
    mkdir -p $test
    cp ../skel/$test/* ./$test
    cd $test
    make cleanall
    make all FEMDIR=$femdir
    cd ..
  done
  cd $act
fi

#----------------------------------------

if [ $run = "YES" ]; then
  echo "...running test in $testdir with ht in $femdir"
  cd $testdir
  for test in $tests
  do
    cd $test
    $femdir/ht < $test.str | tee $test.log
    cd ..
  done
  cd $act
fi

#----------------------------------------

if [ $post = "YES" ]; then
  echo "...postprocessing test in $testdir"
  cd $testdir
  for test in $tests
  do
    cd $test
    make memory
    make post
    cd ..
  done
  cd $act
fi

#----------------------------------------------------------






