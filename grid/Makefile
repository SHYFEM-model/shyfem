#
# $Id: Makefile,v 1.7 2009-04-21 10:23:05 georg Exp $
#
#---------------------------------------------

DIR    = grid

#---------------------------------------------	just for default
FEMDIR  = ..
LIBX = -L/usr/X11R6/lib -L/usr/X11/lib -L/usr/lib/X11  -lX11
CC     = gcc
CFLAGS =  -O -Wall -pedantic
#---------------------------------------------

include ../Rules.make

#---------------------------------------------

FEMBIN = $(FEMDIR)/bin

OBJECTS =	gridma1.o gridma2.o \
		gridfi.o gridge.o \
		gridhs.o gridlo.o gridop.o gridpl.o  \
		gridut.o gridwi.o gridps.o \
		gridky.o gridco.o gridvs.o \
		gridnl.o \
		menu.o stack.o screen.o

OBUTIL =	general.o gustd.o debug.o color.o \
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
	$(MAKEDEPEND) $(OBJECTS_SOURCE) $(OBUTIL_SOURCE) $(XOBDEVICE_SOURCE)

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

# Rules ------------------

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.c.obj:
	$(TCC) -c $(TCCOPT) -o$*.obj $*.c
	tlib /c grid.lib -+$*.obj

#----------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

gridma1.o: gridky.h gridhs.h grid_fp.h events.h hash.h grid.h \
		mouse.h screen.h general.h gustd.h grid_df.h \
		grid_ty.h queue.h list.h grid_ex.h fund.h keybd.h \
		menu.h graph.h 
gridma2.o: general.h screen.h hash.h mouse.h grid.h grid_fp.h \
		gridhs.h gridnl.h gridky.h graph.h keybd.h menu.h \
		queue.h grid_ty.h grid_ex.h list.h fund.h gustd.h \
		grid_df.h 
gridfi.o: general.h screen.h grid.h hash.h grid_fp.h gridhs.h \
		psgraph.h args.h menu.h keybd.h grid_ex.h fund.h \
		list.h grid_ty.h queue.h grid_df.h gustd.h 
gridge.o: screen.h hash.h grid.h general.h grid_fp.h gustd.h \
		grid_df.h gridhs.h menu.h grid_ty.h queue.h fund.h \
		grid_ex.h list.h 
gridhs.o: general.h hash.h keybd.h grid_ex.h list.h fund.h \
		grid_ty.h queue.h grid_df.h gridhs.h gridnl.h 
gridlo.o: graph.h menu.h keybd.h fund.h grid_ex.h list.h queue.h \
		grid_ty.h grid_df.h gustd.h general.h screen.h \
		grid.h mouse.h hash.h events.h grid_fp.h gridhs.h 
gridop.o: gridhs.h grid_fp.h general.h hash.h grid.h screen.h \
		grid_ty.h queue.h grid_ex.h list.h fund.h \
		generalx.h menu.h gustd.h grid_df.h graph.h \
		xgraph.h 
gridpl.o: graph.h general.h screen.h hash.h grid.h mouse.h menu.h \
		queue.h grid_ty.h grid_ex.h fund.h list.h \
		grid_fp.h gridhs.h grid_df.h 
gridut.o: general.h grid.h hash.h screen.h grid_ex.h fund.h list.h \
		queue.h grid_ty.h menu.h grid_df.h gridhs.h \
		grid_fp.h 
gridwi.o: grid_fp.h grid_df.h gridhs.h menu.h grid_ty.h queue.h \
		grid_ex.h list.h fund.h screen.h hash.h grid.h \
		mouse.h graph.h general.h 
gridps.o: grid_fp.h gridhs.h general.h screen.h grid.h hash.h \
		menu.h keybd.h grid_ex.h fund.h list.h grid_ty.h \
		queue.h grid_df.h gustd.h psgraph.h args.h 
gridky.o: screen.h hash.h grid.h graph.h general.h grid_fp.h \
		gustd.h grid_df.h gridhs.h menu.h grid_ty.h \
		queue.h grid_ex.h list.h fund.h 
gridco.o: queue.h grid_ty.h grid_ex.h list.h fund.h menu.h gustd.h \
		color.h grid_df.h graph.h gridhs.h grid_fp.h \
		general.h hash.h grid.h screen.h 
gridnl.o: general.h fund.h list.h gridnl.h
menu.o: screen.h mouse.h general.h graph.h stack.h menu.h assert.h \
		events.h 
stack.o: general.h stack.h
screen.o: graph.h general.h
general.o: general.h
gustd.o: gustd.h general.h
debug.o: general.h debug.h
color.o: general.h graph.h
hash.o: general.h hash.h
list.o: list.h general.h
args.o: general.h
queue.o: general.h queue.h
psgraph.o: general.h psgraph.h gustd.h
xmouse.o: mouse.h
xgraph.o: xgraph.h graph.h generalx.h general.h
xevents.o: events.h generalx.h general.h xgraph.h debug.h
xkeybd.o: keybd.h general.h

