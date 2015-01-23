#!/bin/sh

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

fembin=$FEMDIR/fembin
fem3d=$FEMDIR/fem3d

$fem3d/filetype $*


