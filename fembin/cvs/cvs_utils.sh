#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# utility for cvs routines
# could also use IFS="\n" in loop
#
#---------------------------------------------------------------

Find_rev_to_version()
{
  local tag

  if [ $# -lt 2 ]; then
    echo "Missing version information... aborting"
    exit 7
  fi

  tag=`echo $2 | sed -e 's/\./_/g'`
  tag="VERS_$tag"

  Find_rev_to_tag $1 $tag
}

#---------------------------------------------------------------

Find_rev_to_tag()
{
  local file tag in_tag found 

  file=$1
  tag=$2

  lines=`rlog -r $file`

  in_tag="NO"
  found="NO"
  revision=""

  for line in $lines
  do
    if [ $found = "YES" ]; then
      revision=$line
      break
    fi
    if [ $in_tag = "YES" ]; then
      [ "$line" = "$tag:" ] && found="YES"
    fi
    [ "$line" = "names:" ] && in_tag="YES"
    #echo "new line: $line"
  done

  return

  if [ -n "$revision" ]; then
    echo "revision found: $revision"
  else
    echo "no revision found to tag: $tag"
    exit 5
  fi
}

#---------------------------------------------------------------

TestRevision()
{
echo "rev: $revision"
Find_rev_to_tag /home/georg/CVS/fem3d/newbcl.f VERS_4_98
echo "rev: $revision"

Find_rev_to_version /home/georg/CVS/fem3d/newbcl.f 4.90
echo "rev: $revision"
}

#---------------------------------------------------------------
#TestRevision
#---------------------------------------------------------------

