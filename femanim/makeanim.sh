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
# makes animation from file plot.ps
# uses gifsicle to create animated gif file anim.gif
# uses mencoder to create avi file anim.avi
# view animation with "xanim anim.gif" or "mplayer anim.avi"
#
#---------------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
femanim=$FEMDIR/femanim

#---------------------------------------------------------------------

fps=25			# frames per second for mencoder
delay=10		# delay for gifsicle

dogif="NO"		# produce gifs directly

#---------------------------------------------------------------------

rename_prog=$femanim/rename.pl
rename_prog=$femanim/rename-petras.pl

#---------------------------------------------------------------------

FullUsage()
{
  echo ""
  echo "Usage: makeanim.sh {-gif|-avi} [-h|-help] [-options] [ps-file]"
  echo ""
  echo "  ps-file          ps-file with still images to be animated"
  echo "                   if not given uses plot.ps"
  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -test            tests software availability"
  echo "  -avi             creates avi file (needs mencoder)"
  echo "  -gif             creates animated gif file (needs gifsicle)"
  echo "  -only_anim       recreates only animiation (no image files)"
  echo "  -delay #         sets delay time for gif (default 10)"
  echo "  -fps #           sets frames per second for avi (default 25)"
  echo ""
}

ErrorOption()
{
  echo "No such option : $1"
}

TestSoftware()
{
  if [ $output = "avi" ]; then
    mencoder -list-options > /dev/null 2>&1
    if [ $? -ne 0 ]; then
      echo "*** mencoder is not installed... aborting"
      exit 1
    fi
  elif [ $output = "gif" ]; then
    gifsicle --version > /dev/null 2>&1
    if [ $? -ne 0 ]; then
      echo "*** gifsicle is not installed... aborting"
      exit 1
    fi
  elif [ $output = "test" ]; then
    mencoder -list-options > /dev/null 2>&1
    if [ $? -eq 0 ]; then
      echo "   mencoder is installed..."
    else
      echo "   mencoder is not installed..."
    fi
    gifsicle --version > /dev/null 2>&1
    if [ $? -eq 0 ]; then
      echo "   gifsicle is installed..."
    else
      echo "   gifsicle is not installed..."
    fi
  else
    echo "*** Unknown output format: $output"
    exit 1
  fi

  if [ $rename = "YES" ]; then
    $rename_prog
    if [ $? -ne 0 ]; then
      echo "*** cannot find renaming program (internal error)... aborting"
      exit 1
    elif [ $output = "test" ]; then
      echo "   renaming program is installed..."
    fi
  fi
}

#---------------------------------------------------------------------

delete="YES"
make_eps="YES"
rename="YES"
make_gifs="YES"
make_jpgs="NO"

only_anim="NO"
output="none"

#---------------------------------------------------------------------

while [ -n "$1" ]
do
   case $1 in
        -only_anim)     only_anim="YES";;
        -avi)           output="avi";;
        -gif)           output="gif";;
        -test)          output="test"; shift;;
        -delay)         delay=$2; shift;;
        -fps)           fps=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   [ -n "$1" ] && shift
done

if [ $# -gt 0 ]; then
  plotfile=$1
  echo "filename is $plotfile"
  echo "not yet ready to give file name"
  echo "please use always plot.ps"
  echo "you can do this with the command: mv $plotfile plot.ps"
  exit 1
else
  plotfile=plot.ps
fi

if [ "$output" = "avi" ]; then		#create avi file
  echo "creating anim.avi with fps $fps"
elif [ "$output" = "gif" ]; then	#create gif file
  echo "creating anim.gif with delay $delay"
elif [ "$output" = "test" ]; then	#tests software
  TestSoftware
  exit 0
else
  echo "Unknown output format: $output"
  echo "You must at least specify -avi or -gif"
  FullUsage
  exit 1
fi

TestSoftware

if [ ! -f $plotfile ]; then
  echo "No such file: $plotfile ... aborting"
  exit 1
fi

#---------------------------------------------------------------------

if [ "$output" = "avi" ]; then	#create avi file
  make_gifs="NO"
  make_jpgs="YES"
fi
if [ "$only_anim" = "YES" ]; then	#recreate only animation
  echo "only recreating animation"
  delete="NO"
  make_eps="NO"
  rename="NO"
  make_gifs="NO"
  make_jpgs="NO"
fi

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

if [ $delete = "YES" ]; then
  rm -f plot.*.bak
  rm -f plot.*.eps
  rm -f plot.*.gif
  rm -f plot.*.jpg
  rm -f plot.*.ps
fi

#--------------------------- split and make eps

if [ $make_eps = "YES" ]; then
  gps -split -eps plot.ps
  rm -f plot.*.ps
fi

#--------------------------- rename in numerical order

if [ $rename = "YES" ]; then
  #$femanim/rename.pl plot.*.eps
  $femanim/rename-petras.pl plot.*.eps
fi

#--------------------------- create gifs

if [ $make_gifs = "YES" ]; then
  if [ $dogif = "YES" ]; then
    gps -gif plot.*.eps
  else				#go via jpg
    gps -jpg plot.*.eps
    for file in plot.*.jpg
    do
      name=`basename $file .jpg`
      new="$name.gif"
      echo "$file -> $new"
      convert $file $new
    done
  fi
fi

if [ $make_jpgs = "YES" ]; then
  gps -jpg plot.*.eps
fi

######################################################
#
# delete second frame in gifs (bug in gps -gif)
#
#for file in plot.*.gif
#do
#  echo "deleting frame 1 in $file"
#  gifsicle --batch $file --delete '#1'
#done
#
######################################################

#--------------------------- create animation

if [ $output = "avi" ]; then
  mencoder mf://*.jpg -mf w=800:h=600:fps=$fps:type=jpg -ovc lavc \
       -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o anim.avi
  echo "created anim.avi with fps $fps"
elif [ $output = "gif" ]; then
  gifsicle --delay $delay --optimize --colors 256 plot.*.gif > anim.gif
  echo "created anim.gif with delay $delay"
else
  echo "Unknown output format: $output"
fi

#--------------------------- finish

