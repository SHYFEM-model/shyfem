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

mesh.o: meshut.h meshty.h queue.h fund.h meshck.h gustd.h meshop.h \
		meshin.h meshcv.h hash.h stack.h list.h nlist.h \
		meshbd.h mesh.h general.h meshge.h meshfi.h \
		meshhs.h heap.h 
exgrd.o: hash.h stack.h grdio.h grd.h queue.h fund.h grdut.h \
		general.h exgrdop.h gustd.h grdhs.h 
meshut.o: mesh.h meshty.h nlist.h list.h meshut.h hash.h heap.h \
		meshhs.h meshge.h assert.h general.h fund.h \
		queue.h 
meshcv.o: meshhs.h meshcv.h heap.h queue.h fund.h general.h \
		meshck.h list.h nlist.h mesh.h meshty.h hash.h \
		stack.h meshut.h 
meshfi.o: meshty.h mesh.h nlist.h list.h meshut.h hash.h meshgd.h \
		meshop.h meshhs.h meshfi.h gustd.h general.h \
		fund.h queue.h 
meshck.o: meshut.h meshty.h gustd.h meshck.h fund.h queue.h \
		meshin.h meshop.h stack.h hash.h mesh.h nlist.h \
		list.h meshge.h general.h heap.h meshhs.h 
meshin.o: meshhs.h heap.h general.h meshge.h list.h nlist.h mesh.h \
		hash.h stack.h meshop.h meshin.h queue.h fund.h \
		meshck.h meshty.h meshut.h 
meshbd.o: hash.h stack.h nlist.h meshbd.h mesh.h list.h general.h \
		meshge.h meshhs.h heap.h meshut.h meshty.h queue.h \
		fund.h meshop.h meshin.h 
meshge.o: list.h nlist.h mesh.h hash.h meshfi.h meshhs.h heap.h \
		queue.h fund.h general.h gustd.h meshge.h 
meshhs.o: fund.h queue.h hash.h mesh.h list.h
meshop.o: queue.h fund.h general.h gustd.h options.h meshop.h \
		hash.h list.h mesh.h 
meshgd.o: meshut.h meshgd.h hash.h nlist.h mesh.h list.h general.h \
		gustd.h queue.h fund.h meshhs.h args.h 
meshty.o: hash.h list.h mesh.h meshty.h queue.h fund.h general.h \
		assert.h 
general.o: general.h
gustd.o: general.h gustd.h
debug.o: general.h debug.h
hash.o: hash.h general.h
list.o: list.h general.h
args.o: general.h
queue.o: queue.h general.h
heap.o: heap.h general.h
stack.o: general.h stack.h
nlist.o: general.h nlist.h
assert.o: assert.h
exgrdop.o: options.h queue.h fund.h grdut.h exgrdop.h general.h \
		gustd.h grd.h hash.h 
grdio.o: gustd.h general.h grdut.h fund.h queue.h args.h grdhs.h \
		hash.h grd.h 
grdhs.o: grd.h queue.h hash.h fund.h
grdut.o: grdhs.h grd.h general.h queue.h hash.h grdut.h fund.h
gustd.o: gustd.h general.h
args.o: general.h
hash.o: hash.h general.h
queue.o: queue.h general.h
stack.o: general.h stack.h
general.o: general.h

