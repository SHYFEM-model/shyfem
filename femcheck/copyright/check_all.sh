#!/bin/sh
#
# checks all folders for copyright
#
#---------------------------------------------------

basedir=$HOME/shyfem
copydir=$HOME/shyfem/femcheck/copyright
check=$copydir/check.log

#---------------------------------------------------

CheckDir()
{
  echo "========================================"
  echo "./$1"
  echo "========================================"
  $copydir/check_copyright_iter.sh $basedir/$1
}

ElabLog()
{
  echo "+++++++++++++++++++++++++++++++++++++++"
  cat $check | sort | uniq \
	| grep -v -- '-----' \
	| grep -v '=====' \
	| grep -v '^\.' \
	| grep -v '\.grd$' 
  echo "+++++++++++++++++++++++++++++++++++++++"
}

#---------------------------------------------------

[ -f $check ] && rm $check

CheckDir femadj			| tee -a $check
CheckDir femanim		| tee -a $check
CheckDir femcheck		| tee -a $check
CheckDir femdummy		| tee -a $check
CheckDir femlib			| tee -a $check
CheckDir femplot		| tee -a $check

ElabLog
