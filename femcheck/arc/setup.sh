#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

SetDir()
{
  pushd $1	> /dev/null
  local tmp=`pwd`
  popd		> /dev/null

  echo $tmp
}

ContinueNo()
{
        echo -n "Do you wish to continue? (y/n) [n] "
        read response
        case "$response" in
        [yY]*)
                echo ""
                ;;
        *)
                echo "Aborting the installation."
                exit 2
                ;;
        esac
}

AddToBashrc()
{
	tmp=bashrc.$$.tmp
	origfile=~/.bashrc
	conffile=$act/bashrc

	cat $origfile				 > $tmp
	echo ''					>> $tmp
	echo '# Definition for FEM model'	>> $tmp
	echo "[ -f $conffile ] && . $conffile"	>> $tmp
	echo ''					>> $tmp

	mv $tmp $origfile

	tmp=profile.$$.tmp
	origfile=~/.bash_profile
	conffile=$act/profile

	cat $origfile				 > $tmp
	echo ''					>> $tmp
	echo '# Definition for FEM model'	>> $tmp
	echo "[ -f $conffile ] && . $conffile"	>> $tmp
	echo ''					>> $tmp

	mv $tmp $origfile
}

#-----------------------------------------------------------------

act=`pwd`
fem=`SetDir $act/..`
home=$HOME
date=`date`

echo "Actual directory: $act"
echo "fem    directory: $fem"
echo "Home   directory: $home"

#-----------------------------------------------------------------

ForceInstall="NO"

cd $fem
if [ -f INSTALLED ]; then
  echo "The set-up routine has already been run..."
  echo "You only have to run this routine once at set-up."
  ContinueNo
  ForceInstall="YES"
  echo "Saving file INSTALLED for later"
  mv --backup=numbered INSTALLED INSTALLED.tmp
fi
echo "latest install of the FEM model: $date" > INSTALLED
cd $act

#-----------------------------------------------------------------

mkdir -p $home/tmp
mkdir -p $home/bin
mkdir -p $home/sbin
mkdir -p $home/lib
mkdir -p $fem/tmp

pushd $fem		> /dev/null
rm -f bin 
rm -f lib 
ln -s fembin bin
ln -s femlib lib
popd			> /dev/null

#-----------------------------------------------------------------

AddToBashrc

cd $home
if [ -f .profile ]; then
  if [ ! -L .profile ]; then
    new=.profile.$$
    echo ".profile existing ... renaming to $new"
    mv .profile $new
  else
    rm -f .profile
  fi
fi
ln -s .bash_profile .profile
cd $act

#-----------------------------------------------------------------

echo ""
echo "The FEM model has been installed."
echo ""
echo "If necessary, please cleanup the files"
echo ".bashrc and .bash_profile in you home directory."
echo ""
echo "In order that the new settings are effective,"
echo "please log out and in again, or go to your"
echo "home directory and run the command"
echo ". .profile"
echo ""
echo "Please also run ./check_setup.sh to see if"
echo "all programs have been installed for the FEM model."
echo ""

#-----------------------------------------------------------------

