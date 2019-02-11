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
# finds executable files
#
# -R	recursive
# -l	only file name
#
#---------------------------------------------------------

dirs="fembin fem3d femplot femanim femcheck femdoc femregres"
dirs=femlib

#---------------------------------------------------------

Find_sh()
{
  echo "========= looking for sh =========="
  grep -R '^\#\!.*\/sh' $*
}

Find_csh()
{
  echo "========= looking for csh =========="
  grep -R '^\#\!.*\/csh' $*
}

Find_bash()
{
  echo "========= looking for bash =========="
  grep -R -l '^\#\!.*\/bash' $*
}

Find_perl()
{
  echo "========= looking for perl =========="
  grep -R -l '^\#\!\S*\/perl' $*
}

Find_perl_nousr()
{
  echo "========= looking for perl (no usr)=========="
  grep -R '^\#\!\S*\/perl' $* | grep -v usr
}

#---------------------------------------------------------

#Find_sh $dirs
#Find_bash $dirs
#Find_csh $dirs
Find_perl $dirs
#Find_perl_nousr $dirs

#---------------------------------------------------------
