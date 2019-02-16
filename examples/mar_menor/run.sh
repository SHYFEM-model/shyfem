#!/usr/bin/env bash

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# cd to shyfem/examples
pushd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

make clean
make init
./make_sims.sh

#cp -prv . ../Output

