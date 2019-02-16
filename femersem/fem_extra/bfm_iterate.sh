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
# changes for compiler errors/warnings throughout directories
# call with -change to actually change
#
# run from this directory
#
#	./bfm_iterate.sh		only for checking
#	./bfm_iterate.sh -change	to change files
#
#--------------------------------------------------

femdir=$HOME/shyfem
export FEMDIR=$femdir

cd ..
iteratedir "$femdir/femersem/fem_extra/bfm_change.sh $1"

