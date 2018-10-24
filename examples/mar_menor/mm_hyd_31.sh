#!/bin/sh
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_hyd_31

CleanFiles $sim.ts.shy 

Run $sim

CheckFiles $sim.ts.shy
PlotMapSalt apnbath

#--------------------------------------------------

