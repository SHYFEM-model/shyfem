#!/bin/sh
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_hyd_13

CleanFiles $sim.hydro.shy

Run $sim

CheckFiles $sim.hydro.shy
PlotMapVel apn_velcolriver

#-------------------

sim=mm_hyd_11

CleanFiles $sim.hydro.shy

Run $sim

CheckFiles $sim.hydro.shy
PlotMapVel apn_velcolriver

#--------------------------------------------------

