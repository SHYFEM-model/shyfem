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

sim=mm_hyd_32a

CleanFiles $sim.ext salt.2d.*

Run $sim

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?

PlotTsSalt salt.2d.*
cp salt.2d.3 scomp.2d.20

#------------------

sim=mm_hyd_32b

CleanFiles $sim.ext salt.2d.*

Run $sim

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?

PlotTsSalt salt.2d.*
cp salt.2d.3 scomp.2d.50

#------------------

sim=mm_hyd_31

CleanFiles salt.2d.*

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?

PlotTsSalt salt.2d.*
cp salt.2d.3 scomp.2d.35

#------------------

sim=mm_hyd_32
PlotTsSalt scomp.2d.*

#--------------------------------------------------

