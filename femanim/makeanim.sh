#!/bin/sh
#
# makes gif animation from file plot.ps
# creates file anim.gif
# view animation with "xanim anim.gif"
#
#--------------------------- clean files from last call

dogif="NO"

delete="YES"
make_eps="YES"
rename="YES"

#delete="NO"
#make_eps="NO"
#rename="NO"

rm -f plot.*.bak
rm -f plot.*.eps
rm -f plot.*.gif
rm -f plot.*.jpg
rm -f plot.*.ps

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

gifsicle --delay 10 --optimize --colors 256 plot.*.gif > anim.gif

#--------------------------- finish

