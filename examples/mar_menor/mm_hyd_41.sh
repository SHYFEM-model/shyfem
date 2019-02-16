#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
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

