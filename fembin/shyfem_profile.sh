#!/bin/sh
#
# shyfem_profile.sh
#
# sets path and aliases
#
# is called from within .bashrc, .bash_profile, .profile
#
#------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
fembin=$FEMDIR/fembin

FEMDIR_INSTALL=${SHYFEM_INSTALL:=$HOME/shyfem}
fembin_install=$FEMDIR_INSTALL/fembin

# set PATH ----------------------------------------

path=`$fembin_install/shyfem_path.pl $PATH`
export PATH=$path:$fembin:$HOME/shyfem/fembin

# set aliases ----------------------------------------

alias shyfemdir=". $fembin_install/shyfem_dir.sh"
alias shyfeminstall=". $fembin_install/shyfem_dir.sh $FEMDIR_INSTALL"
alias shyfemcd="cd $FEMDIR"

# end of routine ----------------------------------------

