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

sim=mm_hyd_43

CleanFiles $sim.ts.shy $sim.hydro.shy

Run $sim

CheckFiles $sim.ts.shy $sim.hydro.shy

PlotMapVel apnbath 1
PlotMapVel apnbath 5

#--------------------------------------------------

