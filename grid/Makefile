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

gridma1.o: grid_ex.h gridky.h grid_df.h screen.h grid_fp.h keybd.h \
		list.h mouse.h graph.h queue.h gustd.h hash.h \
		gridhs.h grid_ty.h fund.h events.h grid.h menu.h \
		general.h 
gridma2.o: grid_ex.h gridky.h grid_df.h screen.h grid_fp.h keybd.h \
		list.h mouse.h graph.h queue.h gustd.h hash.h \
		gridhs.h grid_ty.h fund.h gridnl.h grid.h menu.h \
		general.h 
gridfi.o: grid_ex.h grid_df.h screen.h grid_fp.h keybd.h list.h \
		psgraph.h queue.h gustd.h hash.h gridhs.h \
		grid_ty.h fund.h grid.h general.h menu.h args.h 
gridge.o: grid_ex.h gustd.h queue.h gridhs.h hash.h grid_df.h \
		screen.h grid_ty.h grid_fp.h fund.h list.h grid.h \
		general.h menu.h 
gridhs.o: grid_ex.h queue.h grid_df.h hash.h gridhs.h grid_ty.h \
		fund.h gridnl.h list.h keybd.h general.h 
gridlo.o: grid_ex.h grid_df.h screen.h grid_fp.h keybd.h list.h \
		mouse.h graph.h queue.h gustd.h hash.h gridhs.h \
		grid_ty.h fund.h events.h grid.h menu.h general.h 
gridop.o: grid_ex.h xgraph.h grid_df.h screen.h grid_fp.h list.h \
		generalx.h graph.h queue.h gustd.h hash.h gridhs.h \
		grid_ty.h fund.h grid.h general.h menu.h 
gridpl.o: graph.h grid_ex.h queue.h gridhs.h hash.h grid_df.h \
		screen.h grid_ty.h grid_fp.h fund.h list.h mouse.h \
		grid.h general.h menu.h 
gridut.o: grid_ex.h queue.h gridhs.h hash.h grid_df.h screen.h \
		grid_ty.h grid_fp.h fund.h list.h grid.h general.h \
		menu.h 
gridwi.o: graph.h grid_ex.h queue.h gridhs.h hash.h grid_df.h \
		screen.h grid_ty.h grid_fp.h fund.h list.h mouse.h \
		grid.h general.h menu.h 
gridps.o: grid_ex.h grid_df.h screen.h grid_fp.h keybd.h list.h \
		psgraph.h queue.h gustd.h hash.h gridhs.h \
		grid_ty.h fund.h grid.h general.h menu.h args.h 
gridky.o: graph.h grid_ex.h gustd.h queue.h gridhs.h hash.h \
		grid_df.h screen.h grid_ty.h grid_fp.h fund.h \
		list.h grid.h general.h menu.h 
gridco.o: color.h grid_ex.h grid_df.h screen.h grid_fp.h list.h \
		graph.h queue.h gustd.h hash.h gridhs.h grid_ty.h \
		fund.h grid.h general.h menu.h 
gridnl.o: fund.h gridnl.h list.h general.h
menu.o: stack.h graph.h assert.h screen.h events.h mouse.h \
		general.h menu.h 
stack.o: stack.h general.h
screen.o: graph.h general.h
general.o: general.h
gustd.o: gustd.h general.h
debug.o: debug.h general.h
color.o: graph.h general.h
hash.o: hash.h general.h
list.o: list.h general.h
args.o: general.h
queue.o: queue.h general.h
psgraph.o: psgraph.h gustd.h general.h
xmouse.o: mouse.h
xgraph.o: graph.h xgraph.h generalx.h general.h
xevents.o: debug.h events.h xgraph.h general.h generalx.h
xkeybd.o: keybd.h general.h

