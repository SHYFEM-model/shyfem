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

mesh.o: meshut.h meshop.h meshcv.h gustd.h hash.h meshge.h \
		meshck.h heap.h meshty.h meshbd.h mesh.h meshhs.h \
		stack.h meshfi.h nlist.h fund.h queue.h general.h \
		list.h meshin.h 
exgrd.o: gustd.h hash.h grdio.h grdhs.h fund.h queue.h general.h \
		exgrdop.h grd.h grdut.h stack.h 
meshut.o: fund.h meshge.h hash.h nlist.h meshut.h meshhs.h \
		meshty.h mesh.h list.h assert.h heap.h queue.h \
		general.h 
meshcv.o: meshut.h meshcv.h fund.h hash.h nlist.h heap.h meshck.h \
		queue.h general.h meshhs.h stack.h meshty.h mesh.h \
		list.h 
meshfi.o: meshgd.h list.h mesh.h meshty.h meshhs.h general.h \
		queue.h meshfi.h nlist.h hash.h gustd.h fund.h \
		meshop.h meshut.h 
meshck.o: meshut.h meshop.h gustd.h meshge.h hash.h meshck.h \
		heap.h mesh.h meshty.h stack.h meshhs.h nlist.h \
		fund.h general.h queue.h list.h meshin.h 
meshin.o: fund.h nlist.h general.h queue.h meshin.h list.h \
		meshop.h meshut.h hash.h meshge.h meshck.h heap.h \
		stack.h meshhs.h meshty.h mesh.h 
meshbd.o: queue.h general.h list.h meshin.h nlist.h fund.h heap.h \
		mesh.h meshty.h meshbd.h meshhs.h stack.h meshut.h \
		meshop.h meshge.h hash.h 
meshge.o: meshge.h gustd.h hash.h nlist.h meshfi.h fund.h queue.h \
		general.h heap.h mesh.h list.h meshhs.h 
meshhs.o: queue.h hash.h mesh.h list.h fund.h
meshop.o: queue.h general.h mesh.h list.h meshop.h hash.h gustd.h \
		options.h fund.h 
meshgd.o: meshhs.h mesh.h list.h meshgd.h queue.h general.h fund.h \
		hash.h gustd.h nlist.h meshut.h args.h 
meshty.o: queue.h general.h assert.h meshty.h mesh.h list.h hash.h \
		fund.h 
general.o: general.h
gustd.o: general.h gustd.h
debug.o: debug.h general.h
hash.o: hash.h general.h
list.o: general.h list.h
args.o: general.h
queue.o: queue.h general.h
heap.o: general.h heap.h
stack.o: general.h stack.h
nlist.o: nlist.h general.h
assert.o: assert.h
exgrdop.o: fund.h options.h gustd.h hash.h grdut.h grd.h exgrdop.h \
		general.h queue.h 
grdio.o: queue.h general.h grd.h grdut.h args.h grdhs.h fund.h \
		hash.h gustd.h 
grdhs.o: queue.h hash.h grd.h fund.h
grdut.o: queue.h general.h grdhs.h fund.h grd.h hash.h grdut.h
gustd.o: general.h gustd.h
args.o: general.h
hash.o: hash.h general.h
queue.o: queue.h general.h
stack.o: general.h stack.h
general.o: general.h

