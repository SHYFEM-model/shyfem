#!/bin/sh
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_hyd_41

CleanFiles $sim.ts.shy $sim.hydro.shy

Run $sim

CheckFiles $sim.ts.shy $sim.hydro.shy

PlotMapSalt apnbath 1
PlotMapSalt apnbath 5
PlotMapVel apnbath 1
PlotMapVel apnbath 5

#--------------------------------------------------

