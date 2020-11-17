
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

EXES = mkgeom

default:

mkgeom:
	gfortran -o $@ $@.f

clean:
	-rm -f ggg hhh
	-rm -f *.o
	-rm -f $(EXES)

cleanall: clean
	-rm -f *.grd

