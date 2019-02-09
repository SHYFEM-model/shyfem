
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#------------------------------------------------------

DIR =	spline

FEMDIR = ..
FEMBIN = $(FEMDIR)/fembin

include ../Rules.make

#------------------------------------------------------

LIBFEM  = fem

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

#------------------------------------------------------

OBJS =	hash.o
	
SOBJS =  sgauss.o

EXES =	gridr tsmooth circle

SPECIAL = Makefile README
ZIPFILES   = *.[fFh] $(SPECIAL)

#------------------------------------------------------

default: fem

all: $(EXES) nothing

fem: gridr nothing

nothing:
	@true

gridr:	gridr.o $(OBJS) $(LIBIFS)
	$(LINKER) $(LFLAGS) -o $@ $@.o $(OBJS) $(LIBFS)

tsmooth: tsmooth.o $(SOBJS)
	$(LINKER) $(LFLAGS) -o $@ $@.o $(SOBJS)

circle:	circle.o
	$(LINKER) $(LFLAGS) -o $@ $@.o

#------------------------------------------------------

list:
	pwd; ls

depend: nothing

#------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f $(EXES) a.out
	-rm -f ggg hhh
	-rm -f fort.*
	-rm -f out.ps plot.ps
	-rm -f *.tmp *.bak *.gnu
	-rm -f *.dat
	-rm -f *~
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanall: clean cleangrd
	-rm -f $(DIR).zip

cleangrd:
	-rm -f *.grd

zip: save
save: cleanall
	zip $(DIR) $(ZIPFILES)

#--------------------------------------------- rules

.f.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------

