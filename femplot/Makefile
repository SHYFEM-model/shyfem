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
		suplin.o suplin1.o \
		mod_plot2d.o mod_plot3d.o mod_hydro_plot.o \
		supdum.o 

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

EXES    = plotsim plotcolb windrose testdash \
		reg2fem line_nodes

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------------------------------------------

default: fem

fem: param plotsim line_nodes

all: param $(EXES)

plotsim:  plotsim.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(PSOBJS) $(LIBFS) $(LIBPS) -o $@

plotbug:  plotbug.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
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

lagxitest:  lagxitest.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(PSOBJS) $(LIBFS) $(LIBPS) -o $@

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
	-rm -f *.mod
	-rm -f *_genmod.f90

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
	$(MAKEDEPEND) -moddir=$(MODDIR) *.f
	
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
lagxitest.o: param.h evmain.h basin.h
line_nodes.o: param.h basin.h geom.h
../femlib/mod/mod_hydro_plot.mod: mod_hydro_plot.o
../femlib/mod/mod_plot2d.mod: mod_plot2d.o
../femlib/mod/mod_plot3d.mod: mod_plot3d.o
plotsim.o: depth.h plot_aux_3d.h nlevel.h mkonst.h evmain.h \
		nbasin.h plot_aux.h geom.h hydro_vel.h param.h \
		hydro.h pkonst.h levels.h basin.h hydro_print.h \
		hydro_plot.h 
supano.o: color.h param.h basin.h legend.h simul.h
supbas.o: nbasin.h param.h supout.h basin.h geom.h
supcol.o: color.h
supdep.o: depth.h hydro.h param.h ev.h nlevel.h nbasin.h levels.h \
		basin.h plot_aux.h 
supdum.o: param.h
supint.o: param.h levels.h
supiso.o: color.h param.h basin.h
suplag.o: color.h param.h basin.h plot_aux.h
suplin.o: depth.h color.h param.h levels.h nlevel.h basin.h \
		hydro_print.h 
suplin1.o: hydro_vel.h hydro.h param.h plot_aux_3d.h ev.h nlevel.h \
		nbasin.h levels.h basin.h hydro_print.h plot_aux.h 
supout.o: depth.h hydro.h param.h plot_aux_3d.h nlevel.h simul.h \
		nbasin.h levels.h supout.h basin.h hydro_print.h \
		plot_aux.h hydro_plot.h 
supsim.o: hydro_vel.h depth.h hydro.h param.h plot_aux_3d.h \
		nlevel.h nbasin.h levels.h basin.h geom.h \
		plot_aux.h hydro_plot.h 
suptim.o: timlim.h
supver.o: hydro_vel.h hydro.h param.h plot_aux_3d.h levels.h \
		basin.h nlevel.h ev.h plot_aux.h 

