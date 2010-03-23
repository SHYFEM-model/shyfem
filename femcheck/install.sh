#!/bin/sh
#
#-----------------------------------------------------------------

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

	cat $origfile | grep -v "FEM model"			 > $tmp
	echo ''							>> $tmp
	echo "[ -f $conffile ] && . $conffile #FEM model"	>> $tmp
	echo ''							>> $tmp

	mv $tmp $origfile

	tmp=profile.$$.tmp
	origfile=~/.bash_profile
	conffile=$act/profile

	cat $origfile | grep -v "FEM model"			 > $tmp
	echo ''							>> $tmp
	echo "[ -f $conffile ] && . $conffile #FEM model"	>> $tmp
	echo ''							>> $tmp

	mv $tmp $origfile
}

IsFemDir()
{
  local vers=$1/fem3d/VERSION
  if [ -f $vers ]; then
    head -1 $vers | sed -e 's/  */ /g' | cut -f4 -d" "
  else
    echo ""
  fi
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

if [ $fem = $home/fem ]; then
  echo ""
  echo "The install routine cannot be run in the directory $home/fem."
  echo "Please go to the original directory where you have unpacked the model."
  echo "If you have unpacked the model in directory fem, please rename the"
  echo "directory and make sure that no directory $home/fem exists."
  echo ""
  exit 1
fi

#-----------------------------------------------------------------

ForceInstall="NO"
cd $fem
if [ -f INSTALLED ]; then
  echo ""
  echo "The install routine has already been run..."
  echo "You only have to run this routine once at model installation."
  ContinueNo
  ForceInstall="YES"
fi
cd $act

#-----------------------------------------------------------------

echo ""
echo "Creating a link $home/fem to the new model directory..."

if [ -a $home/fem ]; then
  version=`IsFemDir $home/fem`
  if [ -n "$version" ]; then
    echo "Old FEM version: $version"
    if [ -L $home/fem ]; then
      rm -f $home/fem
      echo "A link to an old FEM installation has been found."
      echo "It has been adjourned to the new model directory."
    else
      olddir=$home/fem.old.$$
      mv -f $home/fem $olddir
      echo "A directory fem has been found in your home directory."
      echo "This was an old FEM installation."
      echo "It has been renamed to $olddir."
      echo "A link fem has been created to the new model directory."
    fi
  else
    echo "A file or directory fem has been found in your home directory."
    echo "This files does not seem to be an old FEM installation."
    echo "Please rename this file and rerun the installation routine."
    exit 1
  fi
else
  echo "A link fem has been created to the new model directory."
fi

ln -fs $fem $home/fem

#-----------------------------------------------------------------

fem=$home/fem

cd $fem
if [ $ForceInstall = "YES" ]; then
  echo "Saving file INSTALLED for later"
  mv --backup=numbered INSTALLED INSTALLED.tmp
fi
echo "latest install of the FEM model: $date" > INSTALLED
cd $act

#-----------------------------------------------------------------

#mkdir -p $home/tmp
#mkdir -p $home/bin
#mkdir -p $home/sbin
#mkdir -p $home/lib

pushd $fem		> /dev/null
rm -f bin 
rm -f lib 
ln -s fembin bin
ln -s femlib lib
mkdir -p tmp
popd			> /dev/null

#-----------------------------------------------------------------

AddToBashrc

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
echo "Please also run 'make check' to see if"
echo "all programs needed for the FEM model are installed."
echo ""

#-----------------------------------------------------------------

