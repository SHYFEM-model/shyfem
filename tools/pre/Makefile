#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

include ../../Rules.make

#----------------------	INCLUDE	----------------------

INCLUDE = -I$(FEMDIR)/../include

#--------------------   PATH   ---------------------------

DIRLIB=$(FEMDIR)/../libs

#---------------------  LIBRARY -----------------------

LIBS = -L$(DIRLIB) -lshy

#------------------------------------------------------------------------------------------------
EXE            = shypre
#------------------------------------------------------------

#OOBJ = ../2nc/strings.o ../2nc/apn_util.o ../2nc/apnfile.o 
#OBJ = subclo.o grd.o grdutil.o cuthill_mckee.o rosen.o shypre.o
OBJ = subclo.o strings.o apn_util.o apnfile.o grd.o grdutil.o cuthill_mckee.o rosen.o shypre.o

#-------------------------------------------------------------------------------------------------

all:	pre

pre:	$(EXE)

#------------------------------------------------

%.o: %.f 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

%.o: %.f90 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

#------------------------------------------------

shypre:  $(MODULES) $(OOBJ) $(OBJ) 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBS)
	cp $@ $(FEMDIR)/../bin/
clean:
	rm *.o
	rm $(EXE)

