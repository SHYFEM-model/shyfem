#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# prepares VERSION file for beta version
#
#--------------------------------------------------------------

fembin=./fembin

version_file=VERSION

#-----------------------------------------

line=`$fembin/shyfem_version.pl -noextra $version_file`
extra=beta_`date +"%Y-%m-%d"`

echo "$line   $extra"                    > ver.tmp
echo ""                                 >> ver.tmp
echo "================================" >> ver.tmp
echo ""                                 >> ver.tmp
date                                    >> ver.tmp
echo ""                                 >> ver.tmp
echo "beta version"                     >> ver.tmp
echo ""                                 >> ver.tmp
git diff --cached --stat                >> ver.tmp
echo ""                                 >> ver.tmp
git diff --stat                         >> ver.tmp
echo ""                                 >> ver.tmp

cat $version_file                       >> ver.tmp

mv ver.tmp $version_file

