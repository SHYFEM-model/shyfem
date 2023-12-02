#!/bin/sh

./regress.pl

gp -t "regression linear" lin_orig.txt lin_reg.txt \
	-style points lin_pert.txt
mv out.ps lin.ps

gp -t "detrended linear" lin_detrended.txt \
	-style points lin_pert.txt
mv out.ps lin_det.ps

gp -t "scatter plot linear" lin_line.txt -style points lin_scatt.txt
mv out.ps lin_scatt.ps

gp -t "regression exponential" exp_orig.txt exp_reg.txt \
	-style points exp_pert.txt
mv out.ps exp.ps

gp -t "detrended exponential" exp_detrended.txt \
	-style points exp_pert.txt
mv out.ps exp_det.ps

gp -t "scatter plot exponential" exp_line.txt -style points exp_scatt.txt
mv out.ps exp_scatt.ps

gpsmerge lin*.ps exp*.ps > all.ps
echo "plots are in all.ps"

gv all.ps

