#!/bin/sh
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_hyd_21

CleanFiles $sim.hydro.shy 

Run $sim

CheckFiles $sim.hydro.shy
PlotMapVel apn_vel
PlotMapVel apnbath

#---------------

sim=mm_hyd_13

CleanFiles $sim.hydro.shy 

Run $sim

CheckFiles $sim.hydro.shy
PlotMapVel apn_vel

#--------------------------------------------------

