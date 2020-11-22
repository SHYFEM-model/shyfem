
#--------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#--------------------------------------------------------------------------

EXES = ncgrid

default: ncgrid

ncgrid: compile
	./ncgrid.sh

compile:
	gfortran -o ncgrid ncgrid.f

clean:
	-rm -f $(EXES)
	-rm -f [xy]data.txt
	-rm -f [xy]header.txt
	-rm -f header.txt

cleanall: clean
	-rm -f ncgrid.grd

