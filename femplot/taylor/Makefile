#
# $Id: Makefile,v 1.1 2008-04-16 10:19:24 georg Exp $
#
#------------------------------------------------------------

DIR	= taylor

include ../../Rules.make

FEMDIR  = ../..

#------------------------------------------------------------

LIBFEM  = fem
LIBPOST = post

PSOBJS	= supsim.o $(NPOBJS)
PSOBJS	= $(NPOBJS)

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#----------------

SPECIAL = Makefile DATA taylor.txt

EXES    = taylor_test
ALLEXES    = taylor_test taylor

#----------------

fem: all

all: $(EXES) data

taylor: taylor.o $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o $(LIBFS) $(LIBPS) -o $@

taylor_test: taylor_test.o taylor.o $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $@.o taylor.o $(LIBFS) $(LIBPS) -o $@

#----------------------------

data:
	grep ' C1 ' DATA > DATA_C1.dat
	grep ' A '  DATA > DATA_A.dat
	grep ' B '  DATA > DATA_B.dat

#----------------------------

cleanall: clean cleanzip

clean: cleanobj cleanexe cleanvar

cleanobj:
	-rm -f *.o

cleanexe:
	-rm -f $(ALLEXES)
	-rm -f *.exp
	-rm -f *.exe

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f *.0
	-rm -f *.ps
	-rm -f .memory _memory
	-rm -f apnpar.tmp
	-rm -f gmon.out
	-rm -f DATA_*.dat

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

#----------------------

list:
	pwd; ls

rcs:
	ci -l *.f *.h *.c $(SPECIAL)

save:   cleanall
	@echo $(DIR)
	-zip $(DIR) *.[fhc] $(SPECIAL)

zip:	save

strip:	$(EXES)
	-strip $(EXES)

# Rules ------------------

.c.o:
	$(CC)  -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

#----------------------

