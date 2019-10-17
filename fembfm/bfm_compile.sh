#!/bin/bash
#
# this script sets up the BFM model for SHYFEM
#
# this script has to be executed in the main directory of SHYFEM
#
#---------------------------------------------------------------

if [ ! -f VERSION ]; then
  echo "*** you must be in the main directory of SHYFEM to run this script"
  exit 1
fi

BFMDIR=$1
if [ -z "$BFMDIR" ]; then
  echo "*** BFMDIR is not given - please specify in Rules.make"
  echo "    Please see Rules.make for more details"
  exit 3
fi

if [ ! -d $BFMDIR ]; then
  echo "*** BFMDIR is not exisiting - please correct in Rules.make"
  exit 5
fi

shydir=$( pwd )
echo "using BFM directory: $BFMDIR"
echo "this creates the BFM library to be linked with SHYFEM"
echo "this has to be done only once"

#---------------------------------------------------------------
# make new profile to build library
#---------------------------------------------------------------

cd $BFMDIR/build/configurations
[ -d BFMLIB ] && rm -rf BFMLIB
cp -a STANDALONE_PELAGIC BFMLIB
cd BFMLIB

# here we add an extra line to specify BFMEXE
cp configuration configuration.tmp
cat configuration.tmp | sed -e 's/\//        BFMEXE  = 'libbfm.a',\n\//' \
	> configuration

#---------------------------------------------------------------
# add extra files for SHYFEM in directory shyfem
#---------------------------------------------------------------

cd $BFMDIR/src
[ -d shyfem ] && rm -rf shyfem
mkdir -p shyfem
cp $shydir/fembfm/*.F90 shyfem

#---------------------------------------------------------------
# insert new files into configuration
#---------------------------------------------------------------

cd $BFMDIR/build
cp bfm_configure.sh bfm_configure.sh.tmp
$shydir/fembfm/bfm_insert.pl bfm_configure.sh.tmp > bfm_configure.sh
chmod +x bfm_configure.sh

#---------------------------------------------------------------
# start confguration
#---------------------------------------------------------------

cd $BFMDIR/build

./bfm_configure.sh -gcd -p BFMLIB

#---------------------------------------------------------------
# copy library to lib
#---------------------------------------------------------------

cd $BFMDIR
mkdir -p $BFMDIR/lib
cp bin/libbfm.a lib

#---------------------------------------------------------------
# end of routine
#---------------------------------------------------------------

