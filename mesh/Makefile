#
# $Id: Makefile,v 1.8 2009-04-21 10:25:35 georg Exp $
#
#---------------------------------------------

DIR    = mesh

#---------------------------------------------  just for default
FEMDIR  = ..
CC     = gcc
CFLAGS =  -O -Wall -pedantic
#---------------------------------------------

include ../Rules.make

#---------------------------------------------

FEMBIN    = $(FEMDIR)/fembin

OBJECTS =        meshut.o meshcv.o meshfi.o meshck.o meshin.o \
		 meshbd.o meshge.o meshhs.o meshvs.o meshop.o \
		 meshgd.o meshty.o

OBUTIL =         general.o gustd.o debug.o options.o \
		 hash.o list.o args.o queue.o \
		 heap.o stack.o nlist.o assert.o

EXOBJECTS = 	 exgrdvs.o exgrdop.o

SPOBJECTS = 	 spgrdvs.o spgrdop.o spgrdut.o spgrdnl.o

BCKOBJS	=	 backg.o

GRDOBJECTS = 	 grdio.o grdhs.o grdut.o \
		 gustd.o args.o hash.o queue.o stack.o general.o options.o

GRDFILES = 	 final.grd M_*.grd

EXES =		 mesh elips exgrd maskgrd intgrd spgrd

MAKEDEPEND = $(FEMBIN)/mkdp.pl

OBJECTS_SOURCE = ${OBJECTS:.o=.c}
OBUTIL_SOURCE = ${OBUTIL:.o=.c}
EXOBJECTS_SOURCE = ${EXOBJECTS:.o=.c}
GRDOBJECTS_SOURCE = ${GRDOBJECTS:.o=.c}

LIBM = -lm 
LIBMALLOC = -L/usr/local/lib -ldmalloc 
LIBEXT = $(LIBMALLOC)
LIBEXT =

TCC    = tcc
TCCOPT = -O -G -Z -w -d -mh

#----------------------

default: fem
all: depend $(EXES)
fem: mesh exgrd links

links:
	@$(FEMBIN)/symlink $(FEMBIN) mesh exgrd

#----------------------

mesh: mesh.o $(OBJECTS) $(OBUTIL)
	$(CC) $(LCFLAGS) mesh.o $(OBJECTS) $(OBUTIL) $(LIBEXT) $(LIBM) -o mesh

elips: elips.o
	$(CC) $(LCFLAGS) elips.o $(LIBM) -o elips

exgrd: exgrd.o $(EXOBJECTS) $(GRDOBJECTS)
	$(CC) $(LCFLAGS) exgrd.o $(EXOBJECTS) $(GRDOBJECTS) -o exgrd

spgrd: spgrd.o $(SPOBJECTS) $(GRDOBJECTS)
	$(CC) $(LCFLAGS) spgrd.o $(SPOBJECTS) $(GRDOBJECTS) -lm -o spgrd

maskgrd: maskgrd.o $(GRDOBJECTS)
	$(CC) $(LCFLAGS) maskgrd.o $(GRDOBJECTS) $(LIBM) -o maskgrd

intgrd: intgrd.o $(GRDOBJECTS) $(BCKOBJS)
	$(CC) $(LCFLAGS) intgrd.o $(GRDOBJECTS) $(BCKOBJS) $(LIBM) -o intgrd

#----------------------

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	$(MAKEDEPEND) mesh.c exgrd.c 		\
			$(OBJECTS_SOURCE)	\
			$(OBUTIL_SOURCE)	\
			$(EXOBJECTS_SOURCE)	\
			$(GRDOBJECTS_SOURCE)

#----------------------

clean: cleanobj cleanexe cleanvar cleanlinks

cleanall: clean cleangrd cleanzip cleandiff

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) mesh exgrd

cleanobj:
	-rm -f *.o
	-rm -f *.obj

cleanexe:
	-rm -f $(EXES)
	-rm -f *.exe

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -fr cvsdiff_tmp

cleanarc:
	-rm -f *.~*~

cleangrd:
	-rm -f $(GRDFILES)

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f mesh.zip
	-rm -f save.zip

#----------------------

zip: save
save:   cleanall
	-zip $(DIR) *

#----------------------

install: $(FEMBIN)/mesh $(FEMBIN)/exgrd

$(FEMBIN)/mesh: mesh
	cp -f mesh $(FEMBIN)

$(FEMBIN)/exgrd: exgrd
	cp -f exgrd $(FEMBIN)

strip:
	-strip mesh exgrd

stub:
	strip mesh exgrd
	coff2exe mesh
	coff2exe exgrd

list:
	@ls

# Rules ------------------

.SUFFIXES:	.o

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.c.obj:
	$(TCC) -c $(TCCOPT) -o$*.obj $*.c
	tlib /c grid.lib -+$*.obj

#----------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

mesh.o: meshty.h general.h meshut.h meshhs.h meshin.h stack.h \
		list.h hash.h meshcv.h meshop.h meshck.h queue.h \
		fund.h meshbd.h meshge.h nlist.h gustd.h meshfi.h \
		mesh.h heap.h 
exgrd.o: fund.h queue.h grdio.h general.h grd.h hash.h grdut.h \
		gustd.h stack.h grdhs.h exgrdop.h 
meshut.o: fund.h meshty.h queue.h meshhs.h meshge.h nlist.h \
		meshut.h general.h assert.h hash.h mesh.h list.h \
		heap.h 
meshcv.o: queue.h meshck.h meshty.h fund.h meshut.h general.h \
		nlist.h meshhs.h stack.h list.h hash.h mesh.h \
		heap.h meshcv.h 
meshfi.o: meshfi.h hash.h mesh.h list.h gustd.h fund.h queue.h \
		meshty.h meshop.h nlist.h meshhs.h meshgd.h \
		general.h meshut.h 
meshck.o: hash.h list.h stack.h meshty.h meshhs.h meshin.h \
		meshut.h general.h mesh.h gustd.h heap.h fund.h \
		meshop.h queue.h meshck.h meshge.h nlist.h 
meshin.o: queue.h meshck.h meshop.h fund.h nlist.h meshge.h mesh.h \
		heap.h meshty.h meshut.h general.h meshin.h \
		meshhs.h stack.h list.h hash.h 
meshbd.o: mesh.h heap.h fund.h queue.h meshop.h meshge.h nlist.h \
		meshbd.h hash.h stack.h list.h meshty.h meshin.h \
		meshhs.h meshut.h general.h 
meshge.o: fund.h queue.h meshhs.h meshge.h nlist.h general.h \
		meshfi.h hash.h mesh.h gustd.h list.h heap.h 
meshhs.o: queue.h list.h fund.h hash.h mesh.h
meshop.o: options.h hash.h mesh.h gustd.h list.h general.h fund.h \
		meshop.h queue.h 
meshgd.o: list.h args.h gustd.h hash.h mesh.h meshgd.h meshut.h \
		general.h nlist.h meshhs.h queue.h fund.h 
meshty.o: fund.h meshty.h queue.h general.h assert.h mesh.h hash.h \
		list.h 
general.o: general.h
gustd.o: gustd.h general.h
debug.o: debug.h general.h
hash.o: general.h hash.h
list.o: general.h list.h
args.o: general.h
queue.o: queue.h general.h
heap.o: heap.h general.h
stack.o: general.h stack.h
nlist.o: general.h nlist.h
assert.o: assert.h
exgrdop.o: grdut.h gustd.h hash.h options.h exgrdop.h queue.h \
		fund.h grd.h general.h 
grdio.o: hash.h grdut.h gustd.h args.h grdhs.h fund.h queue.h \
		general.h grd.h 
grdhs.o: grd.h queue.h hash.h fund.h
grdut.o: grdhs.h grd.h general.h hash.h fund.h queue.h grdut.h
gustd.o: gustd.h general.h
args.o: general.h
hash.o: hash.h general.h
queue.o: general.h queue.h
stack.o: stack.h general.h
general.o: general.h

