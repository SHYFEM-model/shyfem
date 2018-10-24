#!/bin/sh
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_hyd_22a

CleanFiles $sim.hydro.shy 

Run $sim

CheckFiles $sim.hydro.shy
PlotMapVel apn_vel
PlotMapVel apnbath

#--------------------------------------------------

