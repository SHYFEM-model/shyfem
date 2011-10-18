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

mesh.o: mesh.h meshin.h meshfi.h meshge.h list.h meshbd.h meshut.h \
		meshty.h stack.h meshcv.h gustd.h queue.h meshop.h \
		meshck.h hash.h fund.h heap.h nlist.h meshhs.h \
		general.h 
exgrd.o: stack.h grdio.h queue.h gustd.h grdut.h hash.h exgrdop.h \
		fund.h grdhs.h grd.h general.h 
meshut.o: meshty.h queue.h hash.h mesh.h assert.h heap.h fund.h \
		nlist.h meshge.h meshhs.h list.h meshut.h \
		general.h 
meshcv.o: meshty.h stack.h meshcv.h queue.h meshck.h hash.h mesh.h \
		heap.h fund.h nlist.h meshhs.h list.h meshut.h \
		general.h 
meshfi.o: meshty.h meshgd.h queue.h gustd.h meshop.h hash.h mesh.h \
		fund.h meshfi.h nlist.h meshhs.h list.h general.h \
		meshut.h 
meshck.o: mesh.h meshin.h meshge.h list.h meshut.h stack.h \
		meshty.h gustd.h queue.h meshck.h meshop.h hash.h \
		fund.h heap.h nlist.h meshhs.h general.h 
meshin.o: mesh.h meshin.h meshge.h list.h meshut.h stack.h \
		meshty.h queue.h meshck.h meshop.h hash.h heap.h \
		fund.h nlist.h meshhs.h general.h 
meshbd.o: mesh.h meshin.h meshge.h list.h meshbd.h meshut.h \
		meshty.h stack.h queue.h meshop.h hash.h fund.h \
		heap.h nlist.h meshhs.h general.h 
meshge.o: queue.h gustd.h hash.h mesh.h heap.h fund.h meshfi.h \
		nlist.h meshhs.h meshge.h list.h general.h 
meshhs.o: queue.h fund.h list.h hash.h mesh.h
meshop.o: queue.h gustd.h options.h meshop.h hash.h mesh.h fund.h \
		list.h general.h 
meshgd.o: meshgd.h queue.h gustd.h hash.h mesh.h fund.h nlist.h \
		meshhs.h list.h meshut.h general.h args.h 
meshty.o: meshty.h assert.h queue.h fund.h hash.h list.h mesh.h \
		general.h 
general.o: general.h
gustd.o: gustd.h general.h
debug.o: debug.h general.h
hash.o: hash.h general.h
list.o: list.h general.h
args.o: general.h
queue.o: queue.h general.h
heap.o: heap.h general.h
stack.o: stack.h general.h
nlist.o: nlist.h general.h
assert.o: assert.h
exgrdop.o: queue.h gustd.h grdut.h options.h hash.h exgrdop.h \
		fund.h grd.h general.h 
grdio.o: queue.h gustd.h grdut.h hash.h fund.h grd.h grdhs.h \
		general.h args.h 
grdhs.o: queue.h fund.h grd.h hash.h
grdut.o: queue.h fund.h grdut.h grdhs.h grd.h hash.h general.h
gustd.o: gustd.h general.h
args.o: general.h
hash.o: hash.h general.h
queue.o: queue.h general.h
stack.o: stack.h general.h
general.o: general.h

