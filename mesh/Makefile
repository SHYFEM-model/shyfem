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

mesh.o: meshcv.h general.h meshop.h hash.h meshge.h meshbd.h \
		gustd.h heap.h meshfi.h list.h stack.h meshut.h \
		fund.h meshhs.h nlist.h meshin.h meshck.h mesh.h \
		meshty.h queue.h 
exgrd.o: exgrdop.h grdio.h gustd.h stack.h queue.h grdhs.h grd.h \
		grdut.h general.h fund.h hash.h 
meshut.o: fund.h general.h assert.h meshge.h meshhs.h hash.h \
		nlist.h heap.h list.h mesh.h meshty.h meshut.h \
		queue.h 
meshcv.o: stack.h meshty.h queue.h meshut.h mesh.h list.h meshck.h \
		heap.h meshhs.h hash.h nlist.h meshcv.h general.h \
		fund.h 
meshfi.o: general.h meshop.h fund.h meshhs.h hash.h nlist.h \
		gustd.h meshgd.h queue.h meshty.h meshut.h mesh.h \
		meshfi.h list.h 
meshck.o: queue.h meshty.h mesh.h meshck.h meshin.h nlist.h \
		meshhs.h fund.h meshut.h stack.h list.h gustd.h \
		heap.h meshge.h hash.h general.h meshop.h 
meshin.o: meshck.h queue.h meshty.h mesh.h fund.h meshin.h \
		meshhs.h nlist.h heap.h stack.h meshut.h list.h \
		meshop.h general.h hash.h meshge.h 
meshbd.o: queue.h meshty.h mesh.h fund.h meshin.h meshhs.h nlist.h \
		heap.h stack.h meshut.h list.h meshop.h general.h \
		meshbd.h hash.h meshge.h 
meshge.o: queue.h meshfi.h list.h mesh.h gustd.h heap.h meshge.h \
		nlist.h meshhs.h hash.h general.h fund.h 
meshhs.o: fund.h queue.h list.h mesh.h hash.h
meshop.o: gustd.h queue.h options.h list.h mesh.h general.h \
		meshop.h fund.h hash.h 
meshgd.o: fund.h general.h args.h nlist.h meshhs.h hash.h gustd.h \
		list.h mesh.h meshut.h queue.h meshgd.h 
meshty.o: general.h assert.h fund.h hash.h queue.h meshty.h list.h \
		mesh.h 
general.o: general.h
gustd.o: gustd.h general.h
debug.o: debug.h general.h
hash.o: hash.h general.h
list.o: general.h list.h
args.o: general.h
queue.o: queue.h general.h
heap.o: general.h heap.h
stack.o: stack.h general.h
nlist.o: nlist.h general.h
assert.o: assert.h
exgrdop.o: hash.h fund.h general.h grdut.h grd.h queue.h options.h \
		gustd.h exgrdop.h 
grdio.o: hash.h args.h general.h grdut.h fund.h queue.h grd.h \
		grdhs.h gustd.h 
grdhs.o: fund.h hash.h grd.h queue.h
grdut.o: fund.h grdut.h general.h hash.h grdhs.h grd.h queue.h
gustd.o: general.h gustd.h
args.o: general.h
hash.o: hash.h general.h
queue.o: general.h queue.h
stack.o: stack.h general.h
general.o: general.h

