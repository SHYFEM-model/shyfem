#
# $Id: Makefile,v 1.11 2009-04-21 10:21:45 georg Exp $
#
#-----------------------------------------------------------

DIR    = post

#----------------------------------------------------------- just for default
FEMDIR = ..
DIRLIB = $(FEMDIR)/femlib
CC     = gcc
CFLAGS =  -O -Wall -pedantic
F77    = g77
FFLAGS =  -O -Wall
#-----------------------------------------------------------

include ../Rules.make

#-----------------------------------------------------------

FEMBIN    = $(FEMDIR)/fembin

LIBPOST = $(DIRLIB)/libpost
LIBCALP = $(DIRLIB)/libcalp

LIBPS =	-L$(DIRLIB) -lpost
LIBCS =	-L$(DIRLIB) -lcalp

#--------------------------------

EXES = democalp demopost psdemo colorbar pf pc layout

SPECIAL = Makefile CALCOMP COMMIT EPS LOG VERSION pspost.tex ps/cnr.ps

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#---------------------------

default: fem

fem: libs demopost

all: $(EXES) libs

libs: post calp

democalp: democalp.o pgraph.o general.o pcalp.o
	$(F77) $@.o pgraph.o general.o pcalp.o -o $@

demopost: demopost.o pgraph.o general.o
	$(CC) $@.o pgraph.o general.o -o $@

psdemo: psdemo.o post
	$(F77) $@.o $(LIBPS) -o $@

colorbar: colorbar.o post
	$(F77) $@.o $(LIBPS) -o $@

pf: pf.o post
	$(F77) $@.o $(LIBPS) -o $@

pc: pc.o post
	$(CC) $@.o $(LIBPS) -o $@

layout: layout.o
	$(CC) $@.o -o $@

#---------------------------

calp: $(LIBCALP).a

$(LIBCALP).a: general.o pgraph.o pcalp.o pcalpf.o
	echo $(HOME) $(DIRLIB)
	ar rvs $(LIBCALP).a $?

post: $(LIBPOST).a

$(LIBPOST).a: general.o gustd.o psgraph.o psgraphf.o colorutil.o
	ar rvs $(LIBPOST).a $?

#----------------------------- documentation

docs: tex
latex: tex
tex:
	latex pspost
	latex pspost
	latex pspost
	pdflatex pspost

ps:
	dvips pspost

#----------------------------- dependencies

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff -quiet * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	$(MAKEDEPEND) *.[cf]

#---------------------------

cleanall : clean cleanlib cleanzip cleandiff

clean: cleanobj cleanexe cleanvar cleantex

cleanlib:
	-rm -f $(LIBCALP).a
	-rm -f $(LIBPOST).a

cleanexe:
	-rm -f $(EXES)
	-rm -f *.exe

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanvar:
	-rm -f a.out core
	-rm -f *.bak
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.* *.tmp
	-rm -f plot.ps plot.*.ps demo.ps
	-rm -f plot.pdf plot.*.pdf
	-rm -f out.ps out.*.ps
	-rm -f *.stackdump
	-rm -fr cvsdiff_tmp

cleantex:
	-rm -f pspost.aux pspost.log pspost.dvi pspost.ps pspost.toc pspost.pdf

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f *.zip

#----------------------------

cvschange:
	cvschange -d

list:
	pwd; ls

zip: save
save:   cleanall
	@echo $(DIR)
	-zip $(DIR) *.[fch] $(SPECIAL)

#----------------------------

.c.o:
	$(CC)  -c -g $(CFLAGS) $<
		 
.f.o:
	$(F77) -c -g $(FFLAGS) $<

#----------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

colorutil.o: general.h
demopost.o: general.h pgraph.h
general.o: general.h
gustd.o: general.h gustd.h
pc.o: pspost.h
pcalp.o: pgraph.h general.h
pgraph.o: general.h pgraph.h
psgraph.o: general.h colorutil.h psgraph.h gustd.h
psgraphf.o: general.h psgraph.h

