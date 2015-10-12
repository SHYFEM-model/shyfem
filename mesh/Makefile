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

cleanall: clean cleangrd cleanzip

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

mesh.o: meshfi.h stack.h general.h meshop.h mesh.h queue.h \
		meshhs.h hash.h meshge.h meshin.h meshcv.h \
		meshut.h fund.h nlist.h list.h meshck.h heap.h \
		meshty.h meshbd.h gustd.h 
exgrd.o: grdhs.h fund.h grdio.h grdut.h grd.h hash.h gustd.h \
		queue.h exgrdop.h stack.h general.h 
meshut.o: meshhs.h queue.h meshty.h mesh.h heap.h general.h list.h \
		nlist.h assert.h fund.h meshut.h meshge.h hash.h 
meshcv.o: list.h meshck.h general.h stack.h meshhs.h queue.h \
		heap.h meshty.h mesh.h meshcv.h meshut.h hash.h \
		nlist.h fund.h 
meshfi.o: list.h general.h meshop.h queue.h meshty.h mesh.h \
		meshgd.h meshhs.h gustd.h hash.h meshut.h meshfi.h \
		nlist.h fund.h 
meshck.o: gustd.h heap.h meshty.h list.h meshck.h nlist.h fund.h \
		meshut.h meshin.h meshge.h hash.h meshhs.h queue.h \
		mesh.h meshop.h general.h stack.h 
meshin.o: meshty.h heap.h list.h meshck.h fund.h nlist.h meshge.h \
		meshin.h meshut.h hash.h meshhs.h mesh.h queue.h \
		meshop.h general.h stack.h 
meshbd.o: meshop.h stack.h general.h meshhs.h mesh.h queue.h \
		meshin.h meshge.h meshut.h hash.h fund.h nlist.h \
		list.h meshbd.h meshty.h heap.h 
meshge.o: meshhs.h gustd.h queue.h mesh.h heap.h general.h list.h \
		meshfi.h nlist.h fund.h meshge.h hash.h 
meshhs.o: queue.h mesh.h fund.h hash.h list.h
meshop.o: fund.h hash.h options.h mesh.h queue.h gustd.h list.h \
		general.h meshop.h 
meshgd.o: gustd.h meshhs.h meshgd.h mesh.h queue.h general.h \
		args.h list.h fund.h nlist.h meshut.h hash.h 
meshty.o: general.h list.h queue.h mesh.h meshty.h hash.h assert.h \
		fund.h 
general.o: general.h
gustd.o: gustd.h general.h
debug.o: general.h debug.h
hash.o: general.h hash.h
list.o: general.h list.h
args.o: general.h
queue.o: queue.h general.h
heap.o: heap.h general.h
stack.o: stack.h general.h
nlist.o: general.h nlist.h
assert.o: assert.h
exgrdop.o: fund.h grd.h hash.h options.h grdut.h queue.h gustd.h \
		general.h exgrdop.h 
grdio.o: grdut.h grd.h hash.h grdhs.h fund.h args.h general.h \
		gustd.h queue.h 
grdhs.o: hash.h grd.h queue.h fund.h
grdut.o: grdut.h grd.h general.h hash.h grdhs.h queue.h fund.h
gustd.o: general.h gustd.h
args.o: general.h
hash.o: hash.h general.h
queue.o: general.h queue.h
stack.o: general.h stack.h
general.o: general.h

