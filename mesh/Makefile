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

mesh.o: list.h heap.h meshin.h nlist.h meshop.h meshhs.h meshge.h \
		general.h meshbd.h queue.h meshty.h mesh.h \
		meshut.h stack.h fund.h meshcv.h meshck.h meshfi.h \
		hash.h gustd.h 
exgrd.o: general.h fund.h exgrdop.h grdut.h hash.h queue.h gustd.h \
		grdio.h grd.h stack.h grdhs.h 
meshut.o: heap.h list.h assert.h meshhs.h nlist.h meshut.h \
		meshge.h general.h fund.h mesh.h queue.h meshty.h \
		hash.h 
meshcv.o: list.h heap.h stack.h meshut.h nlist.h meshhs.h fund.h \
		meshck.h meshcv.h general.h queue.h meshty.h \
		hash.h mesh.h 
meshfi.o: meshgd.h list.h meshut.h meshop.h meshhs.h nlist.h \
		fund.h general.h mesh.h gustd.h queue.h meshty.h \
		hash.h meshfi.h 
meshck.o: general.h mesh.h queue.h meshty.h heap.h meshin.h list.h \
		meshge.h meshop.h meshhs.h nlist.h meshck.h fund.h \
		gustd.h hash.h meshut.h stack.h 
meshin.o: heap.h meshin.h list.h meshge.h meshhs.h meshop.h \
		nlist.h general.h mesh.h queue.h meshty.h meshut.h \
		stack.h meshck.h fund.h hash.h 
meshbd.o: meshbd.h general.h queue.h meshty.h mesh.h list.h \
		meshin.h heap.h meshge.h nlist.h meshhs.h meshop.h \
		fund.h hash.h stack.h meshut.h 
meshge.o: hash.h queue.h mesh.h gustd.h meshfi.h fund.h general.h \
		meshge.h nlist.h meshhs.h list.h heap.h 
meshhs.o: list.h fund.h queue.h hash.h mesh.h
meshop.o: meshop.h options.h list.h mesh.h gustd.h queue.h hash.h \
		general.h fund.h 
meshgd.o: fund.h general.h hash.h queue.h mesh.h gustd.h meshgd.h \
		list.h meshut.h nlist.h meshhs.h args.h 
meshty.o: fund.h general.h mesh.h hash.h queue.h meshty.h assert.h \
		list.h 
general.o: general.h
gustd.o: gustd.h general.h
debug.o: general.h debug.h
hash.o: hash.h general.h
list.o: general.h list.h
args.o: general.h
queue.o: general.h queue.h
heap.o: heap.h general.h
stack.o: general.h stack.h
nlist.o: general.h nlist.h
assert.o: assert.h
exgrdop.o: grdut.h gustd.h queue.h hash.h general.h fund.h \
		exgrdop.h grd.h options.h 
grdio.o: grd.h grdhs.h args.h queue.h hash.h gustd.h grdut.h \
		fund.h general.h 
grdhs.o: hash.h queue.h grd.h fund.h
grdut.o: hash.h queue.h grd.h grdhs.h grdut.h fund.h general.h
gustd.o: general.h gustd.h
args.o: general.h
hash.o: hash.h general.h
queue.o: queue.h general.h
stack.o: stack.h general.h
general.o: general.h

