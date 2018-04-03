#
# $Id: Makefile,v 1.9 2009-04-21 10:20:50 georg Exp $
#
#--------------------------------------------- 

DIR	= femadj

#--------------------------------------------- just for default
FEMDIR = ..
FEMBIN = $(FEMDIR)/fembin
DIRLIB = $(FEMDIR)/femlib
F77    = g77
FFLAGS =  -O -Wall
LINKER = $(F77)
#---------------------------------------------

include ../Rules.make

#--------------------------------------------- libraries

LIBFEM  = fem
LIBPOST = post

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#--------------------------------------------- macros

OBJS	=	subsrt.o subssv.o subsss.o subfil.o subgeo.o \
		subscn.o

ADJOBJS	=	adj4el.o adj5el.o adj7el.o adj575.o \
		adjplo.o adjgrd.o adjneu.o adjvar.o \
		mod_adj_grade.o

EXES	=	shyadj
EXELINKS = 	$(FEMBIN) $(EXES)

SPECIAL =       Makefile TODO COMMIT VERSION
SAVEFILES =     *.[fh] $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	shyadj links
all:	$(EXES)

links:
	@$(FEMBIN)/symlink $(EXELINKS)

shyadj:  shyadj.o $(ADJOBJS) $(OBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(ADJOBJS) $(OBJS) $(LIBFS) $(LIBPS) -o $@

#--------------------------------------------- clean

cleanall: clean cleanzip cleangrd cleandiff

clean: cleanobj cleanexe cleanvar cleantex cleanlinks

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(EXELINKS)

cleanexe:
	-rm -f $(EXES)
	-rm -f *.exe

cleantex:
	-rm -f *.dvi *.aux

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f errout.dat
	-rm -f gmon.out
	-rm -f out.ps plot.ps
	-rm -f .memory
	-rm -fr cvsdiff_tmp

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleangrd:
	-rm -f *.grd

#--------------------------------------------- various

param:

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f

list:
	pwd; ls

save:   zip
zip:	cleanall
	-zip $(DIR) $(SAVEFILES)

strip:
	-strip $(EXES)

#--------------------------------------------- rules

.SUFFIXES: .F

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.f.o:
	$(F77) -c $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

#--------------------------------------------- makedepend

# DO NOT DELETE THIS LINE -- make depend depends on it.

../femlib/mod/mod_adj_grade.mod: mod_adj_grade.o
../femlib/mod/scan_string.mod: subscn.o
adj4el.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod
adj575.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod
adj5el.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod \
		param.h 
adj7el.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod
adjgrd.o: nbstatic.h
adjneu.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod \
		../femlib/mod/mod_depth.mod param.h 
adjplo.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod \
		param.h 
adjvar.o: ../femlib/mod/basin.mod ../femlib/mod/mod_adj_grade.mod \
		nbstatic.h param.h 
shyadj.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/mod_adj_grade.mod \
		../femlib/mod/mod_depth.mod param.h 

