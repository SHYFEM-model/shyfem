#!/bin/sh
#
# makes gif animation from file plot.ps
# creates file anim.gif
# view animation with "xanim anim.gif"
#
#--------------------------- clean files from last call

rm -f plot.*.bak
rm -f plot.*.eps
rm -f plot.*.gif
rm -f plot.*.ps

#--------------------------- split and make eps

gps -split -eps plot.ps
rm -f plot.*.ps

#--------------------------- rename in numerical order

./rename.pl plot.*.eps

#--------------------------- create gifs

gps -gif plot.*.eps

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

