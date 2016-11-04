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

NPOBJS  = supbas.o supiso.o supdep.o supano.o \
		supcol.o col.o \
		suptim.o supout.o suputi.o supsim.o \
		supint.o supver.o suplag.o \
		suplin.o suplin1.o \
		mod_plot2d.o mod_plot3d.o mod_hydro_plot.o \
		plotutil.o

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

EXE_AUX = plotcolb windrose testdash reg2fem
EXES    = plotsim line_nodes shyplot

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------------------------------------------

default: fem

fem: $(EXES) links

all: $(EXES)

shyplot:  shyplot.o $(PSOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(PSOBJS) $(LIBFS) $(LIBPS) -o $@

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

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

#------------------------------------------------------------

cleanall: clean cleanzip cleandiff

clean: cleanobj cleanexe cleanvar

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanexe:
	-rm -f $(EXES)
	-rm -f $(EXE_AUX)
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
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

#------------------------------------------------------------

param:
	
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

../femlib/mod/color.mod: col.o
../femlib/mod/mod_hydro_plot.mod: mod_hydro_plot.o
../femlib/mod/mod_plot.mod: supint.o
../femlib/mod/mod_plot2d.mod: mod_plot2d.o
../femlib/mod/mod_plot3d.mod: mod_plot3d.o
../femlib/mod/plotutil.mod: plotutil.o
lagxitest.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod
line_nodes.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod
plotsim.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod mkonst.h pkonst.h 
plotutil.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/nls.mod ../femlib/mod/shyfile.mod 
shyplot.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/plotutil.mod \
		../femlib/mod/shyfile.mod \
		../femlib/mod/shyutil.mod 
supano.o: ../femlib/mod/basin.mod ../femlib/mod/color.mod legend.h \
		simul.h 
supbas.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod \
		supout.h 
supcol.o: ../femlib/mod/color.mod
supdep.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod 
supint.o: ../femlib/mod/levels.mod ../femlib/mod/para.mod
supiso.o: ../femlib/mod/basin.mod ../femlib/mod/color.mod
suplag.o: ../femlib/mod/basin.mod ../femlib/mod/color.mod \
		../femlib/mod/mod_hydro_plot.mod 
suplin.o: ../femlib/mod/basin.mod ../femlib/mod/color.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro_print.mod 
suplin1.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod 
supout.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod simul.h supout.h 
supsim.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod \
		../femlib/mod/mod_hydro_vel.mod 
suptim.o: timlim.h
supver.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_plot.mod \
		../femlib/mod/mod_hydro_vel.mod 

