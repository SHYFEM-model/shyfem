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
# shyfem_install_soft.sh
#
# changes symbolic link and dot files (.bashrc, .bash_profile, .profile)
#
# normally called from Makefile
# must be called from within root of shyfem directory
#
# Usage: shyfem_install_soft.sh [-reset]
#
#--------------------------------------------------------

#--------------------------------------------------------
# functions ---------------------------------------------
#--------------------------------------------------------

ChangeDotFiles()
{
  if [ -f $HOME/.bashrc ]; then
    ChangeDot .bashrc
  else
    echo ".bashrc file not exisiting ... created"
    touch $HOME/.bashrc
    ChangeDot .bashrc
  fi

  if [ -f $HOME/.bash_profile ]; then
    ChangeDot .bash_profile
  elif [ -f $HOME/.profile ]; then
    ChangeDot .profile
  else
    echo ".profile file not exisiting ... created"
    touch $HOME/.profile
    ChangeDot .profile
  fi
}

ChangeDot()
{
  afile=$HOME/$1
  if [ -L $afile ]; then	# change file, not link
    afile=`readlink -f $afile`
  fi
  save=$afile.$date.$$

  #echo "$afile -> $save"; return

  mv -f $afile $save
  $femdir/fembin/shyfem_install.pl $option $femdir $save > tmp.tmp
  [ $? != 0 ] && exit 1
  mv -f tmp.tmp $afile

  hfile=`basename $afile`
  hsfile=`basename $save`
  changed="$changed $hfile"

  echo "$hfile has been changed... original saved to $hsfile"
}

CreateSymlink()
{
  link=$HOME/$1
  date=`date "+%Y%m%d"`

  cd $link
  femlink=$( pwd -P )
  cd $femdir

  if [ $femlink = $femdir ]; then
    echo "already in $1... no symbolic link needed"
  elif [ -d $link ]; then                 #directory exists
    if [ -L $link ]; then               #symbolic link
      rm -f $link
    else                                #real directory
      save=$link.$date.$$
      mv -f $link $save
      echo "renaming existing directory shyfem to $save"
    fi
  fi

  if [ $femlink = $femdir ]; then
    :
  elif [ "$reset" = "YES" ]; then
    echo "deleted symbolic link $link"
  else
    echo "creating symbolic link from $femdir to $link"
    ln -s $femdir $link
  fi
}

CheckVersion()
{
  if [ -x ./fembin/shyfem_version.sh ]; then
    version=`./fembin/shyfem_version.sh $dir`
  fi
  if [ -z "$version" -o "$version" = "unknown" ]; then
    echo "cannot get version for $dir ... aborting" 1>&2
    exit 1
  fi
}

#--------------------------------------------------------
# initialize routine ------------------------------------
#--------------------------------------------------------

femdir=`pwd -P`
changed=""
option=$1
reset="NO"
[ "$option" = "-reset" ] && reset="YES"

export SHYFEMDIR=.

CheckVersion

echo "========================================================="
if [ "$reset" = "YES" ]; then
  echo "Uninstalling the SHYFEM model"
else
  echo "Installing the SHYFEM model"
fi
echo "      running shyfem_install_soft.sh $option"
echo "      using directory: $femdir"
echo "========================================================="

#--------------------------------------------------------
# run installation --------------------------------------
#--------------------------------------------------------

CreateSymlink shyfem
ChangeDotFiles

#--------------------------------------------------------
# final message ----------------------------------------
#--------------------------------------------------------

echo ""
if [ "$reset" = "YES" ]; then
  echo "The SHYFEM model has been uninstalled."
else
  echo "The SHYFEM model has been installed."
fi
echo ""
echo "If necessary, please cleanup the following files"
echo "in your home directory: $changed"
echo ""
echo "In order that the new settings are effective,"
echo "please log out and in again, or go to your"
echo "home directory and run one of the following commands:"

for file in $changed
do
  echo "  . $file"
done

[ "$reset" = "YES" ] && exit 0

echo ""
echo "If this is the first time you install SHYFEM,"
echo "please also run 'make check_software' to see if"
echo "all software needed for the SHYFEM model is installed."
echo ""
echo "Please note that you still have to compile the model."
echo "You should first customize the file Rules.make."
echo "After this run"
echo "  make clean"
echo "  make fem"
echo "in order to compile all programs. You have to redo"
echo "this step every time you change Rules.make."
echo ""

#--------------------------------------------------------
# end of routine ----------------------------------------
#--------------------------------------------------------

