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

