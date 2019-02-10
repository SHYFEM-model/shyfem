
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#-----------------------------------------------------------

DIR    = hcbs

#-----------------------------------------------------------

FEMDIR = ..
DIRLIB = $(FEMDIR)/femlib
FEMBIN = $(FEMDIR)/fembin
F77    = g77
FFLAGS =  -O
LINKER = $(F77)
LFLAGS = -O
CC     = gcc
CFLAGS = -O -Wall -pedantic

#------------------------------------------------------------

include ../Rules.make

#------------------------------------------------------------

LIBX = -L/usr/X11R6/lib -L/usr/X11/lib -L/usr/lib/X11  -lX11

LIBGRAPPA = grappa

LIBGS   = -L$(DIRLIB) -l$(LIBGRAPPA)
LIBIGS  = $(DIRLIB)/lib$(LIBGRAPPA).a

OBJS = xgraph.o xgraphf.o gustd.o general.o

EXES = xhcbs tvis

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------

default: fem

fem: grappa nothing

nothing:
	@true

all: fem $(EXES)

libs: grappa

xhcbs: xhcbs.o $(OBJS)
	$(F77) -o $@ $@.o $(OBJS) $(LIBX)

tvis:   tvis.o $(LIBIGS)
	$(F77) $(LFLAGS) -o $@ $@.o $(LIBGS) $(LIBX)

list:
	@ls

diff:

param:

depend:
	@$(MAKEDEPEND) *.f *.c

#------------------------------------------------------------

grappa: $(LIBIGS)

$(LIBIGS): general.o gustd.o xgraph.o xgraphf.o
	ar rv $(LIBIGS) $?
	ranlib $(LIBIGS)

#------------------------------------------------------------

cleanall: clean cleanlib cleanzip

clean: cleanobj cleanexe cleanvar

cleanlib:
	-rm -f $(LIBIGS)

cleanexe:
	-rm -f $(EXES)

cleanobj:
	-rm -f *.o

cleanvar:
	-rm -f a.out core
	-rm -f *.bak
	-rm -f ggg hhh
	-rm -f fort.*
	-rm -f gmon.out

cleanzip:
	-rm -f *.zip

#----------------------------

save:   zip
zip: cleanall
	zip $(DIR) *

#----------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

#----------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

general.o: general.h
gustd.o: general.h gustd.h
xgraph.o: general.h generalx.h graph.h xgraph.h
xgraphf.o: general.h generalx.h graph.h gustd.h xgraph.h

