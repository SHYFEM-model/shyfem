#!/bin/sh

gfortran nudge.f
[ $? -ne 0 ] && exit 1

./a.out

