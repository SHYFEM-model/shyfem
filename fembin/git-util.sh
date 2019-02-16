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
# git utilities... not to be used standalone
#
# use with ". $fembin/git-util.sh"
#
#-----------------------------------------------------------

ParseVersion()
{
  local version=$1
  local first=`echo $version | cut -b 1`
  local pre

  vers=""
  pre=""

  if [ "$first" = "V" ]; then
    vers=`echo $version | sed -e s/^VERS_//`
    pre="VERS_"
  elif [ "$first" = "v" ]; then
    vers=`echo $version | sed -e s/^v//`
    pre="v"
  fi

  tag=$pre$vers		#this will be globaly available
}

ParseCommit()
{
  local version=$1
  git show $version > /dev/null 2>&1
  local status=$?

  if [ $status -eq 0 ]; then
    #echo "commit found ... using this version"
    local date=`git show --pretty=format:"%ci" $version|head -1|cut -d " " -f 1`
    #echo "commit date: $date"

    vers=${version}_$date
    tag=${version}	#this will be globaly available
  else
    vers=""
    tag=""
  fi
}

InBaseDir()
{
  if [ -f VERSION ]; then
    echo "YES"
  else
    echo "NO"
  fi
}

#-----------------------------------------------------------

prog=$( echo $0 | sed -e 's/.*\///' )
#echo $0 $prog

if [ $prog = git-util.sh ]; then
  echo "utility routines for shell scripts"
  echo "no standalone invocation of this program"
  echo "please import in shell scripts"
fi

#-----------------------------------------------------------

