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

args.o: general.h
args.o: general.h
assert.o: assert.h
debug.o: debug.h general.h
exgrd.o: exgrdop.h fund.h general.h grd.h grdhs.h grdio.h grdut.h \
		gustd.h hash.h queue.h stack.h 
exgrdop.o: exgrdop.h fund.h general.h grd.h grdut.h gustd.h hash.h \
		options.h queue.h 
general.o: general.h
general.o: general.h
grdhs.o: fund.h grd.h hash.h queue.h
grdio.o: args.h fund.h general.h grd.h grdhs.h grdut.h gustd.h \
		hash.h queue.h 
grdut.o: fund.h general.h grd.h grdhs.h grdut.h hash.h queue.h
gustd.o: general.h gustd.h
gustd.o: general.h gustd.h
hash.o: general.h hash.h
hash.o: general.h hash.h
heap.o: general.h heap.h
list.o: general.h list.h
mesh.o: fund.h general.h gustd.h hash.h heap.h list.h mesh.h \
		meshbd.h meshck.h meshcv.h meshfi.h meshge.h \
		meshhs.h meshin.h meshop.h meshty.h meshut.h \
		nlist.h queue.h stack.h 
meshbd.o: fund.h general.h hash.h heap.h list.h mesh.h meshbd.h \
		meshge.h meshhs.h meshin.h meshop.h meshty.h \
		meshut.h nlist.h queue.h stack.h 
meshck.o: fund.h general.h gustd.h hash.h heap.h list.h mesh.h \
		meshck.h meshge.h meshhs.h meshin.h meshop.h \
		meshty.h meshut.h nlist.h queue.h stack.h 
meshcv.o: fund.h general.h hash.h heap.h list.h mesh.h meshck.h \
		meshcv.h meshhs.h meshty.h meshut.h nlist.h \
		queue.h stack.h 
meshfi.o: fund.h general.h gustd.h hash.h list.h mesh.h meshfi.h \
		meshgd.h meshhs.h meshop.h meshty.h meshut.h \
		nlist.h queue.h 
meshgd.o: args.h fund.h general.h gustd.h hash.h list.h mesh.h \
		meshgd.h meshhs.h meshut.h nlist.h queue.h 
meshge.o: fund.h general.h gustd.h hash.h heap.h list.h mesh.h \
		meshfi.h meshge.h meshhs.h nlist.h queue.h 
meshhs.o: fund.h hash.h list.h mesh.h queue.h
meshin.o: fund.h general.h hash.h heap.h list.h mesh.h meshck.h \
		meshge.h meshhs.h meshin.h meshop.h meshty.h \
		meshut.h nlist.h queue.h stack.h 
meshop.o: fund.h general.h gustd.h hash.h list.h mesh.h meshop.h \
		options.h queue.h 
meshty.o: assert.h fund.h general.h hash.h list.h mesh.h meshty.h \
		queue.h 
meshut.o: assert.h fund.h general.h hash.h heap.h list.h mesh.h \
		meshge.h meshhs.h meshty.h meshut.h nlist.h \
		queue.h 
nlist.o: general.h nlist.h
queue.o: general.h queue.h
queue.o: general.h queue.h
stack.o: general.h stack.h
stack.o: general.h stack.h

