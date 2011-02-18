#
# $Id: Makefile,v 1.20 2010-03-11 15:32:16 georg Exp $
#
#------------------------------------------------------------

DIR	= femplot

#------------------------------------------------------------ just for default
FEMDIR = ..
DIRLIB = $(FEMDIR)/femlib
F77    = g77
FFLAGS =  -O
LINKER = $(F77)
#------------------------------------------------------------

include ../Rules.make

#------------------------------------------------------------

FEMBIN    = $(FEMDIR)/fembin

LIBFEM  = fem
LIBPOST = post

NPOBJS  = supbas.o supiso.o supdep.o supano.o supcol.o \
		suptim.o supout.o suputi.o supsim.o \
		supint.o supver.o col.o suplag.o \
		suplin.o suplin1.o

LPOBJS = subfil.o subbas.o links.o sublnku.o subnsu.o	# not used anymore
LPOBJS = 

RGFOBJS = subrgf.o

PSOBJS	= supsim.o $(NPOBJS)
PSOBJS	= $(NPOBJS)

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#------------------------------------------------------------

SPECIAL = Makefile VERSION COMMIT

EXES    = plotsim plotcolb ousinfo windrose turb testdash scalar \
		reg2fem line_nodes

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------------------------------------------

default: fem

fem: param plotsim line_nodes

all: param $(EXES)

plotsim:  plotsim.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(PSOBJS) $(LIBFS) $(LIBPS) -o $@

plotcolb:  plotcolb.o $(NPOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(NPOBJS) $(LIBFS) $(LIBPS) -o $@

ouvinfo:  ouvinfo.o $(NPOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(NPOBJS) $(LIBIFS) $(LIBIPS) -o $@

ousinfo:  ousinfo.o $(NPOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(NPOBJS) $(LIBIFS) $(LIBIPS) -o $@

readout:  readout.o $(NPOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(NPOBJS) $(LIBIFS) $(LIBIPS) -o $@

windrose: windrose.o $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBPS) -o $@

testdash: testdash.o $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBFS) $(LIBPS) -o $@

turb:  turb.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(PSOBJS) $(LIBFS) $(LIBPS) -o $@

scalar: scalar.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(PSOBJS) $(LIBFS) $(LIBPS) -o $@

reg2fem: reg2fem.o $(RGFOBJS) $(LIBIFS)
	$(LINKER) $(LFLAGS) $@.o $(RGFOBJS) $(LIBFS) -o $@

line_nodes: line_nodes.o $(LPOBJS) $(LIBIFS)
	$(LINKER) $(LFLAGS) $@.o $(LPOBJS) $(LIBFS) -o $@

#------------------------------------------------------------

cleanall: clean cleanzip

clean: cleanobj cleanexe cleanvar

cleanobj:
	-rm -f *.o
	-rm -f *.stb

cleanexe:
	-rm -f $(EXES)
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

#------------------------------------------------------------

param:
	-@$(FEMBIN)/make_param.pl $(FEMDIR) param.h
	-@$(FEMBIN)/param.sh $(FEMDIR) param.h
	
diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	$(MAKEDEPEND) *.f
	
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

col.o: color.h
line_nodes.o: param.h basin.h
ousinfo.o: param.h basin.h
plotsim.o: param.h evmain.h basin.h
supano.o: color.h legend.h
supcol.o: color.h
supiso.o: color.h
suplin.o: color.h param.h
suplin1.o: param.h ev.h
supout.o: param.h
supsim.o: param.h
supver.o: param.h ev.h

