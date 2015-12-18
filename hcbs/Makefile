#
# $Id:$
#
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

EXES = xhcbs ttt

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------

default: fem

fem: grappa

all: grappa $(EXES)

libs: grappa

xhcbs: xhcbs.o $(OBJS)
	$(F77) -o $@ $@.o $(OBJS) $(LIBX)

ttt:   ttt.o $(LIBIGS)
	$(F77) $(LFLAGS) -o $@ $@.o $(LIBGS) $(LIBX)

list:
	@ls

diff:

param:

depend:
	$(MAKEDEPEND) *.f *.c

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
gustd.o: gustd.h general.h
xgraph.o: general.h graph.h xgraph.h generalx.h
xgraphf.o: generalx.h xgraph.h gustd.h graph.h general.h

