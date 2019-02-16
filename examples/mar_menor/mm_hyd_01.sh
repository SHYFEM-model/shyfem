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

sim=mm_hyd_01

CleanFiles $sim.ext zeta.2d.* speed.2d.*

Run $sim

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?

CheckFiles zeta.2d.1 speed.2d.1
PlotTsZeta zeta.2d.*
PlotTsVel speed.2d.*

#--------------------------------------------------


