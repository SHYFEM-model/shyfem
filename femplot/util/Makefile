
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#------------------------------------------------------------

DIR	= femplot

#------------------------------------------------------------ just for default
FEMDIR = ../..
FEMBIN = $(FEMDIR)/fembin
DIRLIB = $(FEMDIR)/femlib
F77    = g77
FFLAGS =  -O
LINKER = $(F77)
#------------------------------------------------------------

include $(FEMDIR)/Rules.make

#------------------------------------------------------------

LIBFEM  = fem
LIBPOST = post

COLOR  = supcol.o col.o suputi.o

RGFOBJS = subrgf.o

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#------------------------------------------------------------

SPECIAL = Makefile VERSION COMMIT

EXES    = 
EXESAUX = plotcolb windrose testdash plotbug lagxitest
EXESALL = $(EXES) $(EXESAUX)
EXESLINKS = $(EXES)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------------------------------------------

default: all

fem: all

all: $(EXESALL) links

help:
	@echo "help         this screen"
	@echo "all          compile utility programs"

plotcolb: plotcolb.o $(COLOR) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(COLOR) $(LIBFS) $(LIBPS) -o $@

windrose: windrose.o $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBPS) -o $@

testdash: testdash.o $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBFS) $(LIBPS) -o $@

plotbug: plotbug.o $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBFS) $(LIBPS) -o $@

lagxitest: lagxitest.o $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBFS) $(LIBPS) -o $@

links:

#------------------------------------------------------------

cleanall: clean cleanzip cleandiff

clean: cleanobj cleanexe cleanvar

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanexe:
	-rm -f $(EXESALL)
	-rm -f *.exp
	-rm -f *.exe

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out gmon.out
	-rm -f fort.*
	-rm -f *.0
	-rm -f *.ps
	-rm -f .memory
	-rm -f *~
	-rm -rf cvsdiff_tmp

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXESALL)

#------------------------------------------------------------

param:
	
diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f
	
list:
	pwd; ls

rcs:
	ci -l *.f *.h *.c $(SPECIAL)

save:   cleanall
	@echo $(DIR)
	-zip $(DIR) *.[fhc]
	-zip $(DIR) $(SPECIAL)

zip:	save

strip:	$(EXES)
	-strip $(EXES)

#------------------------------------------------------------
# Rules -----------------------------------------------------
#------------------------------------------------------------

.c.o:
	$(CC)  -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../femlib/mod/mod_testdash.mod: testdash.o

