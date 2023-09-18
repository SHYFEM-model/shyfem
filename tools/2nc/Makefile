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

#---------------------	NETCDF	-----------------------------

LIBNETCDF  = netcdff
DIRLIB_NETCDF = $(NETCDFDIR)/lib/
DIRINC_NETCDF = $(NETCDFDIR)/include
LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)

#----------------------	LIBRARY	-----------------------------

LIBGS      = $(LIBG_NETCDF) $(LIBS)

#------------------------------------------------------------
EXENC            = nos2nc ous2nc ts2nc
#----------------------------------------------------------------------

OBJNOS = evgeom_2nc.o regular_2nc.o strings.o apn_util.o apnfile.o depth_util_2nc.o netcdf_util.o \
	 output_util.o netcdf.o ionos_2nc.o nosutil.o nos2nc.o

OBJOUS = evgeom_2nc.o regular_2nc.o strings.o apn_util.o apnfile.o depth_util_2nc.o netcdf_util.o \
	 output_util.o netcdf.o ousutil.o ous2nc.o

OBJTS = evgeom_2nc.o netcdf.o ts2nc.o

#----------------------------------------------------------------------

all:	nc

nc:	$(EXENC)

ous: ous2nc

nos: nos2nc

ts: ts2nc

#------------------------------------------------

%.o: %.f 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

%.o: %.f90 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

#------------------------------------------------

nos2nc: $(MODULES) $(OBJNOS) 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBGS)
	cp $@ $(FEMDIR)/../bin/

ous2nc: $(MODULES) $(OBJOUS) 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBGS)
	cp $@ $(FEMDIR)/../bin/

ts2nc:  $(MODULES) $(OBJTS)
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBGS)
	cp $@ $(FEMDIR)/../bin/
clean:
	rm *.o
	rm $(EXENC)

