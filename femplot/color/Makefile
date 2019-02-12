
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


default:

prep: colormap
colormap:
	cd matplotlib; cp matplotlib.dat ..
	cd extra; make prep; cp extra.dat ..
	cd oceanlib; make prep; cp oceanlib.dat ..
	cd oceanlib_txt; make prep; cp oceanlib_txt.dat ..
	cat matplotlib.dat extra.dat \
		oceanlib.dat oceanlib_txt.dat > colormap.dat

clean:
	-rm -f matplotlib.dat extra.dat oceanlib.dat oceanlib_txt.dat
	cd matplotlib; make clean
	cd extra; make clean
	cd oceanlib; make clean
	cd oceanlib_txt; make clean

