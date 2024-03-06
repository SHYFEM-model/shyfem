#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# shyfem_profile.sh
#
# sets path and aliases
#
# is called from within .bashrc, .bash_profile, .profile
#
#------------------------------------------------------

dist=shyfem

FEMDIR=${SHYFEMDIR:=$HOME/$dist}
fembin=$FEMDIR/bin

FEMDIR_INSTALL=${SHYFEM_INSTALL:=$HOME/$dist}
fembin_install=$FEMDIR_INSTALL/bin

# set PATH ----------------------------------------

binutil=$fembin_install
[ -d $fembin_install/shyfem_util ] && binutil=$fembin_install/shyfem_util

path=`$binutil/shyfem_path.pl $PATH`
export PATH=$path:$fembin:$HOME/$dist/bin

# set aliases ----------------------------------------

alias shyfemdir=". $binutil/shyfem_dir.sh"
alias shyfeminstall=". $binutil/shyfem_dir.sh $FEMDIR_INSTALL"
alias shyfemcd="cd $FEMDIR"

# end of routine ----------------------------------------

