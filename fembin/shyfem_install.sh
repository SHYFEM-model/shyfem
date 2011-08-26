#!/bin/sh
#
# shyfem_install.sh
#
# changes symbolic link and dor files (.bashrc, .bash_profile, .profile)
#
# normally called from Makefile
# must be called from within root of shyfem directory
#
# Usage: shyfem_install.sh
#
#------------------------------------------------------

ChangeDot()
{
  file=$1

  afile=$HOME/$1
  if [ -L $afile ]; then	# change file, not link
    afile=`readlink -f $afile`
  fi
  save=$afile.$date.$$

  #echo "$afile -> $save"; return

  mv $afile $save
  $femdir/fembin/shyfem_install.pl $femdir $save > $afile
  [ $? != 0 ] && exit 1

  echo "$afile has been changed... original saved to $save"
}

CreateSymlink()
{
  link=$HOME/$1
  date=`date "+%Y%m%d"`

  if [ -d $link ]; then                 #directory exists
    if [ -L $link ]; then               #symbolic link
      rm -f $link
    else                                #real directory
      save=$link.$date.$$
      mv -f $link $save
      echo "renaming existing directory shyfem to $save"
    fi
  fi

  echo "creating symbolic link from $femdir to $link"
  ln -s $femdir $link
}

# check shyfem directory -----------------------------------

dir=`pwd -P`

#echo "debug message: using dir as $dir"

export SHYFEMDIR=.

if [ -x ./fembin/shyfem_version.sh ]; then
  version=`./fembin/shyfem_version.sh $dir`
fi
if [ -z "$version" -o "$version" = "unknown" ]; then
  echo "cannot get version for $dir ... aborting" 1>&2
  exit 1
fi

femdir=$dir

# make symbolic link -----------------------------------

CreateSymlink shyfem
CreateSymlink fem               #this might be deleted somewhen

# change dot files -----------------------------------

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

# final message ----------------------------------------

echo ""
echo "The SHYFEM model has been installed."
echo ""
echo "If necessary, please cleanup the following files"
echo "in your home directory: .bashrc .bash_profile .profile"
echo ""
echo "In order that the new settings are effective,"
echo "please log out and in again, or go to your"
echo "home directory and run one of the following commands:"
echo "source .bashrc"
echo "source .bash_profile"
echo "source .profile"
echo ""
echo "Please also run 'make check' to see if"
echo "all programs needed for the FEM model are installed."
echo ""

# end of routine ----------------------------------------

