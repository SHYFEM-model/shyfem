#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# tries to automerge some files
#
#--------------------------------------------------------------------

git-automerge.pl COMMIT > COMMIT.tmp
git-automerge.pl VERSION > VERSION.tmp
git-automerge.pl fem3d/subver.f > fem3d/subver.f.tmp

echo "following files have been merged:"

echo "    COMMIT         -> COMMIT.tmp"
echo "    VERSION        -> VERSION.tmp"
echo "    fem3d/subver.f -> fem3d/subver.f.tmp"

if [ "$1" = "-write" ]; then
  mv -f COMMIT.tmp COMMIT
  mv -f VERSION.tmp VERSION
  mv -f fem3d/subver.f.tmp fem3d/subver.f
  echo "files have been written"
elif [ "$1" = "-diff" ]; then
  echo "------------- diffing COMMIT -----------------"
  diff COMMIT.tmp COMMIT
  echo "------------- diffing VERSION ----------------"
  diff VERSION.tmp VERSION
  echo "---------- diffing fem3d/subver.f ------------"
  diff fem3d/subver.f.tmp fem3d/subver.f
  echo "----------------------------------------------"
elif [ "$1" = "-tkdiff" ]; then
  tkdiff COMMIT.tmp COMMIT
  tkdiff VERSION.tmp VERSION
  tkdiff fem3d/subver.f.tmp fem3d/subver.f
else
  echo "in order to write changes use option -write"
fi

#--------------------------------------------------------------------

