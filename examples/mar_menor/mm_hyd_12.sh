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

sim=mm_hyd_12a

CleanFiles $sim.ext zeta.2d.* speed.2d.*

Run $sim

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?
cp zeta.2d.3 zcomp.2d.p30

#----------

sim=mm_hyd_12b

CleanFiles $sim.ext zeta.2d.* speed.2d.*

Run $sim

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?
cp zeta.2d.3 zcomp.2d.m30

#----------

sim=mm_hyd_11

CleanFiles $sim.ext zeta.2d.* speed.2d.*

Run $sim

CheckFiles $sim.ext
shyelab -split $sim.ext
CheckStatus shyelab $?
cp zeta.2d.3 zcomp.2d.p00

#----------

sim=mm_hyd_12ab
PlotTsZeta zcomp.2d.*

#--------------------------------------------------

