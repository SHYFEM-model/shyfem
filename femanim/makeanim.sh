#!/bin/sh
#
# makes gif animation from file plot.ps
# creates file anim.gif
# view animation with "xanim anim.gif"
#
#--------------------------- clean files from last call

delay=10		# delay for gifsicle
dogif="NO"		# produce gifs directly

#---------------------------------------------------------------------

FullUsage()
{
  echo ""
  echo "Usage: makeanim.sh [-h|-help] [-options]"
  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -only_anim       recreates only animiation"
  echo "  -delay #         sets delay time (default 10)"
  echo ""
}

ErrorOption()
{
  echo "No such option : $1"
}

delete="YES"
make_eps="YES"
rename="YES"
make_gifs="YES"

only_anim="NO"

#---------------------------------------------------------------------

while [ -n "$1" ]
do
   case $1 in
        -only_anim)     only_anim="YES";;
        -delay)         delay=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   [ -n "$1" ] && shift
done

if [ -z "$delay" ]; then
  echo "No delay given... aborting"
  exit 1
fi
echo "using delay $delay"

#---------------------------------------------------------------------

if [ "$only_anim" = "YES" ]; then	#recreate only animation
  echo "only recreating animation"
  delete="NO"
  make_eps="NO"
  rename="NO"
  make_gifs="NO"
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
  #./rename.pl plot.*.eps
  ./rename-petras.pl plot.*.eps
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

gifsicle --delay $delay --optimize --colors 256 plot.*.gif > anim.gif

#--------------------------- finish

