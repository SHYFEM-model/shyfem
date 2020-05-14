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
#  echo "Usage: ext.sh [simulation [basin]]"

################################# get names of simulation and basin

memory.sh $1 $2

####################################################### run readext

readext << EOI


EOI

##################################### create input file for gnuplot

cat > ext.tmp << EOI

#set terminal postscript portrait
set terminal postscript 
set output "out.ps"
# Total number of points to plot : 5
# Points to plot : 1 2 3 4 5
# Points given   : 1 2 3 4 5
set data style lines
set title "extra points"
se xlabel "time [sec]"
set autoscale

set ylabel "water level [m]"

plot \
        "fort.76" using 1:2 title "5", \
        "fort.76" using 1:3 title "59", \
        "fort.76" using 1:4 title "113", \
        "fort.76" using 1:5 title "167", \
        "fort.76" using 1:6 title "221" 

set ylabel "modulus of water velocity [m/sec]"

plot \
        "fort.77" using 1:2 title "5", \
        "fort.77" using 1:3 title "59", \
        "fort.77" using 1:4 title "113", \
        "fort.77" using 1:5 title "167", \
        "fort.77" using 1:6 title "221" 

EOI

########################################### run gnuplot and ghostview

gnuplot ext.tmp
ghostview -a4 -magstep -1 -landscape out.ps

############################################################ clean up

rm -f ext.tmp

################################################################ end
