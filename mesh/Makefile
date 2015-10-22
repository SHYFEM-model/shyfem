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

mesh.o: meshge.h meshhs.h meshbd.h queue.h heap.h hash.h meshut.h \
		meshop.h meshfi.h meshty.h meshin.h list.h \
		general.h meshcv.h nlist.h meshck.h gustd.h \
		stack.h mesh.h fund.h 
exgrd.o: queue.h exgrdop.h grd.h hash.h gustd.h stack.h grdut.h \
		grdio.h grdhs.h general.h fund.h 
meshut.o: queue.h heap.h nlist.h meshhs.h meshge.h meshty.h mesh.h \
		list.h general.h fund.h meshut.h hash.h assert.h 
meshcv.o: nlist.h heap.h meshcv.h queue.h meshhs.h fund.h list.h \
		general.h meshty.h stack.h mesh.h meshck.h \
		meshut.h hash.h 
meshfi.o: meshhs.h nlist.h queue.h meshfi.h gustd.h meshop.h \
		meshut.h hash.h general.h list.h fund.h meshgd.h \
		meshty.h mesh.h 
meshck.o: list.h general.h meshin.h meshty.h meshut.h meshop.h \
		hash.h heap.h queue.h meshhs.h meshge.h fund.h \
		stack.h mesh.h gustd.h meshck.h nlist.h 
meshin.o: meshge.h meshhs.h queue.h heap.h hash.h meshut.h \
		meshop.h meshty.h list.h meshin.h general.h \
		nlist.h meshck.h stack.h mesh.h fund.h 
meshbd.o: fund.h stack.h mesh.h nlist.h hash.h meshop.h meshut.h \
		meshin.h list.h general.h meshty.h meshbd.h \
		meshge.h meshhs.h heap.h queue.h 
meshge.o: meshhs.h meshge.h nlist.h heap.h queue.h meshfi.h \
		gustd.h hash.h list.h fund.h general.h mesh.h 
meshhs.o: hash.h mesh.h queue.h fund.h list.h
meshop.o: hash.h meshop.h gustd.h mesh.h list.h fund.h general.h \
		options.h queue.h 
meshgd.o: general.h list.h fund.h meshgd.h mesh.h gustd.h args.h \
		meshut.h hash.h nlist.h queue.h meshhs.h 
meshty.o: assert.h hash.h fund.h list.h general.h mesh.h meshty.h \
		queue.h 
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
exgrdop.o: options.h queue.h exgrdop.h grd.h hash.h gustd.h \
		grdut.h fund.h general.h 
grdio.o: queue.h fund.h general.h grdhs.h grdut.h gustd.h args.h \
		hash.h grd.h 
grdhs.o: queue.h fund.h grd.h hash.h
grdut.o: hash.h grd.h general.h fund.h grdhs.h queue.h grdut.h
gustd.o: gustd.h general.h
args.o: general.h
hash.o: general.h hash.h
queue.o: queue.h general.h
stack.o: stack.h general.h
general.o: general.h

