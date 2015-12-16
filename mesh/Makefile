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

mesh.o: queue.h gustd.h meshty.h meshfi.h meshhs.h heap.h meshbd.h \
		meshop.h meshge.h list.h meshin.h hash.h meshck.h \
		mesh.h meshcv.h nlist.h meshut.h general.h fund.h \
		stack.h 
exgrd.o: grdut.h queue.h hash.h gustd.h grdhs.h grd.h fund.h \
		stack.h exgrdop.h grdio.h general.h 
meshut.o: meshut.h general.h meshge.h nlist.h fund.h hash.h \
		assert.h queue.h mesh.h list.h meshhs.h heap.h \
		meshty.h 
meshcv.o: nlist.h meshut.h general.h fund.h stack.h list.h hash.h \
		meshck.h queue.h meshcv.h mesh.h meshty.h heap.h \
		meshhs.h 
meshfi.o: general.h meshop.h meshut.h nlist.h fund.h meshgd.h \
		mesh.h gustd.h queue.h hash.h list.h meshhs.h \
		meshfi.h meshty.h 
meshck.o: list.h meshin.h mesh.h hash.h meshck.h nlist.h general.h \
		meshut.h fund.h stack.h gustd.h queue.h meshty.h \
		meshhs.h heap.h meshge.h meshop.h 
meshin.o: meshin.h list.h mesh.h meshck.h hash.h nlist.h general.h \
		meshut.h fund.h stack.h queue.h meshty.h heap.h \
		meshhs.h meshop.h meshge.h 
meshbd.o: meshbd.h meshop.h meshge.h meshty.h heap.h meshhs.h \
		queue.h fund.h stack.h nlist.h meshut.h general.h \
		list.h meshin.h hash.h mesh.h 
meshge.o: meshfi.h meshhs.h heap.h list.h hash.h queue.h mesh.h \
		gustd.h fund.h nlist.h meshge.h general.h 
meshhs.o: fund.h list.h queue.h hash.h mesh.h
meshop.o: queue.h hash.h gustd.h mesh.h list.h options.h general.h \
		meshop.h fund.h 
meshgd.o: mesh.h meshgd.h gustd.h queue.h hash.h list.h meshhs.h \
		args.h general.h meshut.h nlist.h fund.h 
meshty.o: fund.h general.h meshty.h mesh.h queue.h hash.h assert.h \
		list.h 
general.o: general.h
gustd.o: general.h gustd.h
debug.o: debug.h general.h
hash.o: general.h hash.h
list.o: general.h list.h
args.o: general.h
queue.o: queue.h general.h
heap.o: heap.h general.h
stack.o: general.h stack.h
nlist.o: nlist.h general.h
assert.o: assert.h
exgrdop.o: general.h exgrdop.h fund.h grd.h gustd.h queue.h hash.h \
		options.h grdut.h 
grdio.o: gustd.h grd.h grdhs.h queue.h hash.h grdut.h args.h \
		general.h fund.h 
grdhs.o: fund.h queue.h hash.h grd.h
grdut.o: fund.h grdut.h hash.h queue.h grd.h grdhs.h general.h
gustd.o: gustd.h general.h
args.o: general.h
hash.o: general.h hash.h
queue.o: queue.h general.h
stack.o: stack.h general.h
general.o: general.h

