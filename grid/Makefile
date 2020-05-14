
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#---------------------------------------------

DIR    = grid

#---------------------------------------------	just for default
FEMDIR  = ..
FEMBIN = $(FEMDIR)/fembin
LIBX = -L/usr/X11R6/lib -L/usr/X11/lib -L/usr/lib/X11  -lX11
CC     = gcc
CFLAGS =  -O -Wall -pedantic
#---------------------------------------------

include ../Rules.make

#---------------------------------------------

OBJECTS =	gridma1.o gridma2.o \
		gridfi.o gridge.o \
		gridhs.o gridlo.o gridop.o gridpl.o  \
		gridut.o gridwi.o gridps.o \
		gridky.o gridco.o gridvs.o \
		gridnl.o \
		menu.o stack.o screen.o

OBUTIL =	general.o gustd.o debug.o color.o options.o \
		 hash.o list.o args.o queue.o psgraph.o 

OBJS   =	$(OBJECTS) $(OBUTIL)

GOBDEVICE = 	gmouse.o ggraph.o gevents.o gkeybd.o

XOBDEVICE = 	xmouse.o xgraph.o xevents.o xkeybd.o

TOBDEVICE = 	mouse.o graph.o 

TOBJECTS =	gridma1.obj gridma2.obj \
		gridbl.obj gridbu.obj griddi.obj gridfi.obj gridge.obj \
		gridhs.obj gridlo.obj gridop.obj gridpl.obj gridss.obj \
		gridut.obj gridwi.obj gridff.obj gridps.obj gridmu.obj \
		gridky.obj gridco.obj

TOBUTIL =	general.obj gustd.obj keybd.obj \
		hash.obj list.obj args.obj queue.obj psgraph.obj \
		mouse.obj graph.obj tevents.obj

EXES	= 	grid tstgr tstcol tstev tstgrx short
TEXES	=	grid.exe

MAKEDEPEND = $(FEMBIN)/mkdp.pl -nomake
MAKEDEPEND = $(FEMBIN)/mkdp.pl

OBJECTS_SOURCE = ${OBJECTS:.o=.c}
OBUTIL_SOURCE = ${OBUTIL:.o=.c}
GOBDEVICE_SOURCE = ${GOBDEVICE:.o=.c}
XOBDEVICE_SOURCE = ${XOBDEVICE:.o=.c}
TOBDEVICE_SOURCE = ${TOBDEVICE:.o=.c}

#---------------------------------------------

LIBSM = $(LIBX) -lXt -lm

LIBGRX = -lgrx -lm
LIBGRX = -lgrx20 -lm
LIBGR = d:\tc\lib\graphics.lib
LIBG = grid.lib
TCC    = tcc
TCCOPT = -O -G -Z -w -d -mh

#----------------------

default: fem
fem: grid links

links:
	@$(FEMBIN)/symlink $(FEMBIN) grid

#----------------------

grid: grid.o $(OBJS) $(XOBDEVICE)
	$(CC) $(LCFLAGS) grid.o $(OBJS) $(XOBDEVICE) $(LIBSM) -o grid
	chmod a+rx grid

gcc: grid.o $(OBJS) $(GOBDEVICE)
	$(CC) $(LCFLAGS) grid.o $(OBJS) $(GOBDEVICE) $(LIBGRX) -o grid

tcc: grid.obj $(OBJECTS) $(OBUTIL)
	$(TCC) $(TCCOPT) -egrid grid.obj $(LIBG) $(LIBGR)

itstgr: tstgr.o 
	$(CC) $(LCFLAGS) tstgr.o  $(LIBGRX) -o tstgr

tstcol: tstcol.o color.o $(OBUTIL)
	$(CC) $(LCFLAGS) tstcol.o color.o $(OBUTIL) $(LIBSM) -o tstcol

tstgr: tstgr.o ggraph.o gevents.o event32.o nextevnt.o debug.o general.o
	$(CC) $(LCFLAGS) tstgr.o ggraph.o gevents.o event32.o nextevnt.o \
		debug.o general.o  $(LIBGRX) -o tstgr

tstev: tstev.o gevents.o gmouse.o event32.o nextevnt.o
	$(CC) $(LCFLAGS) tstev.o gevents.o gmouse.o \
		event32.o nextevnt.o -o tstev

tstgrx: tstgrx.o
	$(CC) $(LCFLAGS) tstgrx.o $(LIBGRX) -o tstgrx

readgrd: readgrd.o $(OBJECTS) $(OBUTIL)
	$(CC) $(LCFLAGS) readgrd.o $(OBJECTS) $(OBUTIL) $(LIBSM) -o readgrd
	chmod a+rx readgrd

short: short.o general.o $(XOBDEVICE)
	$(CC) $(LCFLAGS) short.o general.o $(XOBDEVICE) $(LIBSM) -o short
	chmod a+rx short

#----------------------

clean: cleanobj cleanexe cleanvar cleanlinks

cleanall: clean cleanzip cleandiff

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) grid

cleanobj:
	-rm -f *.o
	-rm -f *.obj

cleanexe:
	-rm -f $(EXES)
	-rm -f $(TEXES)

cleanvar:
	-rm -f a.out core
	-rm -f *.bak
	-rm -f ggg hhh
	-rm -f status.me
	-rm -fr cvsdiff_tmp

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f *.zip

cleangrd:
	-rm -f *.grd

#----------------------

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend: dependgeneral

dependgeneral:
	@$(MAKEDEPEND) $(OBJECTS_SOURCE) $(OBUTIL_SOURCE) $(XOBDEVICE_SOURCE)

dependx11:
	$(MAKEDEPEND) $(OBJECTS_SOURCE) $(OBUTIL_SOURCE) $(XOBDEVICE_SOURCE)

dependgcc:
	$(MAKEDEPEND) $(OBJECTS_SOURCE) $(OBUTIL_SOURCE) $(GOBDEVICE_SOURCE)

dependtcc:
	$(MAKEDEPEND) $(OBJECTS_SOURCE) $(OBUTIL_SOURCE) $(TOBDEVICE_SOURCE)

#----------------------

zip: save
save:   cleanall cleangrd cleanzip
	@echo $(DIR)
	-zip $(DIR) *

#----------------------

install: $(FEMBIN)/grid

$(FEMBIN)/grid: grid
	cp -f grid $(FEMBIN)

strip:
	-strip grid

stub:
	strip grid
	coff2exe grid

list:
	@ls

.PHONY: links

# Rules ------------------

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.c.obj:
	$(TCC) -c $(TCCOPT) -o$*.obj $*.c
	tlib /c grid.lib -+$*.obj

#----------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

args.o: general.h
color.o: general.h graph.h
debug.o: debug.h general.h
general.o: general.h
gridco.o: color.h fund.h general.h graph.h grid.h grid_df.h \
		grid_ex.h grid_fp.h grid_ty.h gridhs.h gustd.h \
		hash.h list.h menu.h queue.h screen.h 
gridfi.o: args.h fund.h general.h grid.h grid_df.h grid_ex.h \
		grid_fp.h grid_ty.h gridhs.h gustd.h hash.h \
		keybd.h list.h menu.h psgraph.h queue.h screen.h 
gridge.o: fund.h general.h grid.h grid_df.h grid_ex.h grid_fp.h \
		grid_ty.h gridhs.h gustd.h hash.h list.h menu.h \
		queue.h screen.h 
gridhs.o: fund.h general.h grid_df.h grid_ex.h grid_ty.h gridhs.h \
		gridnl.h hash.h keybd.h list.h queue.h 
gridky.o: fund.h general.h graph.h grid.h grid_df.h grid_ex.h \
		grid_fp.h grid_ty.h gridhs.h gustd.h hash.h list.h \
		menu.h queue.h screen.h 
gridlo.o: events.h fund.h general.h graph.h grid.h grid_df.h \
		grid_ex.h grid_fp.h grid_ty.h gridhs.h gustd.h \
		hash.h keybd.h list.h menu.h mouse.h queue.h \
		screen.h 
gridma1.o: events.h fund.h general.h graph.h grid.h grid_df.h \
		grid_ex.h grid_fp.h grid_ty.h gridhs.h gridky.h \
		gustd.h hash.h keybd.h list.h menu.h mouse.h \
		queue.h screen.h 
gridma2.o: fund.h general.h graph.h grid.h grid_df.h grid_ex.h \
		grid_fp.h grid_ty.h gridhs.h gridky.h gridnl.h \
		gustd.h hash.h keybd.h list.h menu.h mouse.h \
		queue.h screen.h 
gridnl.o: fund.h general.h gridnl.h list.h
gridop.o: fund.h general.h generalx.h graph.h grid.h grid_df.h \
		grid_ex.h grid_fp.h grid_ty.h gridhs.h gustd.h \
		hash.h list.h menu.h options.h queue.h screen.h \
		xgraph.h 
gridpl.o: fund.h general.h graph.h grid.h grid_df.h grid_ex.h \
		grid_fp.h grid_ty.h gridhs.h hash.h list.h menu.h \
		mouse.h queue.h screen.h 
gridps.o: args.h fund.h general.h grid.h grid_df.h grid_ex.h \
		grid_fp.h grid_ty.h gridhs.h gustd.h hash.h \
		keybd.h list.h menu.h psgraph.h queue.h screen.h 
gridut.o: fund.h general.h grid.h grid_df.h grid_ex.h grid_fp.h \
		grid_ty.h gridhs.h hash.h list.h menu.h queue.h \
		screen.h 
gridwi.o: fund.h general.h graph.h grid.h grid_df.h grid_ex.h \
		grid_fp.h grid_ty.h gridhs.h hash.h list.h menu.h \
		mouse.h queue.h screen.h 
gustd.o: general.h gustd.h
hash.o: general.h hash.h
list.o: general.h list.h
menu.o: assert.h events.h general.h graph.h menu.h mouse.h \
		screen.h stack.h 
psgraph.o: general.h gustd.h psgraph.h
queue.o: general.h queue.h
screen.o: general.h graph.h
stack.o: general.h stack.h
xevents.o: debug.h events.h general.h generalx.h xgraph.h
xgraph.o: general.h generalx.h graph.h xgraph.h
xkeybd.o: general.h keybd.h
xmouse.o: mouse.h

