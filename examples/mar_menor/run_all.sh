#!/bin/bash

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# cd to shyfem/examples
pushd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

make clean
make init
./make_all_sims.sh

#cp -prv . ../Output

