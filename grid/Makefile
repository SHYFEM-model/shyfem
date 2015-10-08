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

cleanall: clean cleanzip

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

gridma1.o: graph.h general.h menu.h grid.h list.h grid_fp.h \
		gridky.h queue.h screen.h gustd.h fund.h grid_ex.h \
		grid_df.h mouse.h hash.h gridhs.h events.h keybd.h \
		grid_ty.h 
gridma2.o: screen.h gustd.h fund.h gridky.h queue.h list.h \
		grid_fp.h grid.h menu.h general.h graph.h keybd.h \
		grid_ty.h gridhs.h hash.h mouse.h gridnl.h \
		grid_df.h grid_ex.h 
gridfi.o: list.h grid_fp.h psgraph.h grid.h gustd.h fund.h \
		screen.h queue.h args.h menu.h general.h hash.h \
		grid_ty.h keybd.h gridhs.h grid_ex.h grid_df.h 
gridge.o: grid_ex.h grid_df.h general.h menu.h grid.h hash.h \
		grid_fp.h list.h gridhs.h queue.h screen.h gustd.h \
		grid_ty.h fund.h 
gridhs.o: list.h hash.h grid_ty.h fund.h keybd.h queue.h gridhs.h \
		grid_ex.h grid_df.h general.h gridnl.h 
gridlo.o: menu.h general.h graph.h screen.h gustd.h fund.h queue.h \
		grid_fp.h list.h grid.h mouse.h grid_df.h \
		grid_ex.h keybd.h grid_ty.h gridhs.h events.h \
		hash.h 
gridop.o: generalx.h xgraph.h hash.h gridhs.h grid_ty.h grid_ex.h \
		grid_df.h grid.h grid_fp.h list.h queue.h screen.h \
		gustd.h fund.h graph.h general.h menu.h 
gridpl.o: grid.h hash.h grid_fp.h list.h gridhs.h queue.h screen.h \
		fund.h grid_ty.h graph.h grid_ex.h general.h \
		grid_df.h menu.h mouse.h 
gridut.o: fund.h grid_ty.h screen.h queue.h gridhs.h list.h \
		grid_fp.h hash.h grid.h menu.h general.h grid_df.h \
		grid_ex.h 
gridwi.o: list.h grid_fp.h hash.h grid.h grid_ty.h fund.h screen.h \
		queue.h gridhs.h graph.h grid_ex.h menu.h mouse.h \
		grid_df.h general.h 
gridps.o: args.h queue.h screen.h gustd.h fund.h grid.h psgraph.h \
		list.h grid_fp.h general.h menu.h gridhs.h keybd.h \
		grid_ty.h hash.h grid_df.h grid_ex.h 
gridky.o: grid_ex.h graph.h menu.h general.h grid_df.h hash.h \
		list.h grid_fp.h grid.h screen.h grid_ty.h gustd.h \
		fund.h gridhs.h queue.h 
gridco.o: color.h grid_ex.h grid_df.h hash.h gridhs.h grid_ty.h \
		graph.h general.h menu.h grid.h list.h grid_fp.h \
		queue.h fund.h gustd.h screen.h 
gridnl.o: general.h gridnl.h fund.h list.h
menu.o: graph.h assert.h menu.h stack.h mouse.h general.h screen.h \
		events.h 
stack.o: stack.h general.h
screen.o: graph.h general.h
general.o: general.h
gustd.o: general.h gustd.h
debug.o: general.h debug.h
color.o: general.h graph.h
hash.o: hash.h general.h
list.o: general.h list.h
args.o: general.h
queue.o: queue.h general.h
psgraph.o: psgraph.h general.h gustd.h
xmouse.o: mouse.h
xgraph.o: general.h generalx.h graph.h xgraph.h
xevents.o: general.h events.h xgraph.h debug.h generalx.h
xkeybd.o: keybd.h general.h

