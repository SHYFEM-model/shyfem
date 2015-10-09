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

mesh.o: meshge.h meshfi.h gustd.h meshbd.h heap.h meshhs.h mesh.h \
		meshty.h nlist.h meshcv.h meshck.h meshop.h \
		meshin.h hash.h list.h meshut.h stack.h fund.h \
		general.h queue.h 
exgrd.o: stack.h fund.h grdhs.h gustd.h hash.h exgrdop.h grd.h \
		queue.h grdut.h grdio.h general.h 
meshut.o: fund.h meshhs.h meshut.h assert.h heap.h list.h meshge.h \
		hash.h queue.h general.h nlist.h meshty.h mesh.h 
meshcv.o: meshut.h stack.h meshhs.h fund.h hash.h list.h heap.h \
		queue.h mesh.h meshty.h nlist.h meshck.h meshcv.h \
		general.h 
meshfi.o: meshop.h queue.h general.h mesh.h meshty.h nlist.h \
		fund.h meshut.h meshhs.h list.h hash.h meshfi.h \
		meshgd.h gustd.h 
meshck.o: general.h queue.h list.h hash.h fund.h meshut.h stack.h \
		meshck.h mesh.h meshty.h nlist.h meshin.h meshop.h \
		heap.h meshge.h gustd.h meshhs.h 
meshin.o: meshop.h meshin.h nlist.h meshty.h mesh.h meshck.h \
		meshhs.h meshge.h heap.h queue.h general.h stack.h \
		meshut.h fund.h hash.h list.h 
meshbd.o: nlist.h meshty.h mesh.h meshop.h meshin.h meshge.h \
		heap.h meshbd.h meshhs.h general.h queue.h hash.h \
		list.h stack.h meshut.h fund.h 
meshge.o: general.h mesh.h nlist.h queue.h list.h heap.h meshge.h \
		hash.h meshfi.h gustd.h fund.h meshhs.h 
meshhs.o: list.h hash.h mesh.h fund.h queue.h
meshop.o: fund.h list.h hash.h gustd.h options.h meshop.h queue.h \
		general.h mesh.h 
meshgd.o: queue.h general.h mesh.h nlist.h fund.h meshut.h \
		meshhs.h args.h list.h hash.h meshgd.h gustd.h 
meshty.o: list.h assert.h hash.h fund.h general.h mesh.h meshty.h \
		queue.h 
general.o: general.h
gustd.o: gustd.h general.h
debug.o: debug.h general.h
hash.o: general.h hash.h
list.o: list.h general.h
args.o: general.h
queue.o: general.h queue.h
heap.o: heap.h general.h
stack.o: stack.h general.h
nlist.o: general.h nlist.h
assert.o: assert.h
exgrdop.o: options.h queue.h grdut.h general.h fund.h grd.h \
		gustd.h exgrdop.h hash.h 
grdio.o: general.h grdut.h queue.h hash.h gustd.h args.h grd.h \
		grdhs.h fund.h 
grdhs.o: fund.h queue.h grd.h hash.h
grdut.o: hash.h grd.h general.h queue.h grdut.h grdhs.h fund.h
gustd.o: general.h gustd.h
args.o: general.h
hash.o: general.h hash.h
queue.o: general.h queue.h
stack.o: stack.h general.h
general.o: general.h

