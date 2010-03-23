#!/bin/sh

make testdash
testdash
ghostview -dsc -scale -2 plot.ps

