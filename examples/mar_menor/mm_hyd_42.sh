#!/bin/sh
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_hyd_42

CleanFiles $sim.ts.shy $sim.hydro.shy

Run $sim

CheckFiles $sim.ts.shy $sim.hydro.shy

PlotMapVel apnbath 1
PlotMapVel apnbath 5

#--------------------------------------------------

