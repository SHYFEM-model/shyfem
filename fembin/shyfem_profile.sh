#!/bin/sh
#
# shyfem_profile.sh
#
# sets path and aliases
#
# is called from within .bashrc, .bash_profile, .profile
#
#------------------------------------------------------

femdir=${SHYFEMDIR:=$HOME/shyfem}
fembin=$femdir/fembin

femdir_install=${SHYFEM_INSTALL:=$HOME/shyfem}
fembin_install=$femdir_install/fembin

# set PATH ----------------------------------------

path=`$fembin_install/shyfem_path.pl $PATH`
export PATH=$path:$fembin:$HOME/shyfem/fembin

# set aliases ----------------------------------------

alias shyfemdir=". $fembin_install/shyfem_dir.sh"
alias femdir=". $fembin_install/shyfem_dir.sh"

alias shyfeminstall=". $fembin_install/shyfem_dir.sh $femdir_install"

alias cdshyfem="cd $femdir"
alias cdfem="cd $femdir"

# end of routine ----------------------------------------

