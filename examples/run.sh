#!/usr/bin/env bash

# cd to shyfem/examples
pushd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

make clean
make init
./make_sims.sh

cp -prv . ../Output
