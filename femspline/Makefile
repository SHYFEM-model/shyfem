#
# $Id:$
#
#------------------------------------------------------

DIR =	spline

FEMDIR = ..
include ../Rules.make

#------------------------------------------------------

FEMBIN    = $(FEMDIR)/fembin

LIBFEM  = fem

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

#------------------------------------------------------

OBJS =	subscn.o subfil.o subsss.o hash.o spline.o sort.o subver.o 
OBJS =	subscn.o subfil.o subsss.o hash.o spline.o sort.o
	
SOBJS =  sgauss.o
VOBJS =  subrnd.o

EXES =	gridr sss ven circle

SPECIAL = Makefile README
ZIPFILES   = *.[fFh] $(SPECIAL)

#------------------------------------------------------

default: fem

all: gridr gridf grids

fem: gridr

gridr:	gridr.o $(OBJS) $(LIBIFS)
	$(LINKER) $(LFLAGS) -o $@ $@.o $(OBJS) $(LIBFS)

sss:	sss.o $(SOBJS)
	$(LINKER) $(LFLAGS) -o $@ $@.o $(SOBJS)

ven:	ven.o $(VOBJS)
	$(LINKER) $(LFLAGS) -o $@ $@.o $(VOBJS)

circle:	circle.o $(VOBJS)
	$(LINKER) $(LFLAGS) -o $@ $@.o $(VOBJS)

#------------------------------------------------------

list:
	pwd; ls

depend:

#------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f $(EXES) a.out
	-rm -f ggg hhh
	-rm -f fort.*
	-rm -f out.ps plot.ps
	-rm -f *.tmp *.bak *.gnu
	-rm -f *.dat
	-rm -f *~
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanall: clean cleangrd
	-rm -f $(DIR).zip

cleangrd:
	-rm -f *.grd

zip: save
save: cleanall
	zip $(DIR) $(ZIPFILES)

#--------------------------------------------- rules

.f.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------

