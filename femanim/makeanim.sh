#!/bin/sh
#
# makes animation from file plot.ps
# uses gifsicle to create animated gif file anim.gif
# uses mencoder to create avi file anim.avi
# view animation with "xanim anim.gif" or "mplayer anim.avi"
#
#--------------------------- clean files from last call

fps=25			# frames per second for mencoder
delay=10		# delay for gifsicle

dogif="NO"		# produce gifs directly

#---------------------------------------------------------------------

FullUsage()
{
  echo ""
  echo "Usage: makeanim.sh {-gif|-avi} [-h|-help] [-options]"
  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
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
    mencoder --version > /dev/null 2>&1
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
  else
    echo "*** Unknown output format: $output"
    exit 1
  fi
  exit 0
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
        -delay)         delay=$2; shift;;
        -fps)           fps=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   [ -n "$1" ] && shift
done

if [ "$output" = "avi" ]; then		#create avi file
  echo "creating anim.avi with fps $fps"
elif [ "$output" = "gif" ]; then	#create gif file
  echo "creating anim.gif with delay $delay"
else
  echo "Unknown output format: $output"
  echo "You must at least specify -avi or -gif"
  FullUsage
  exit 1
fi

TestSoftware

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

