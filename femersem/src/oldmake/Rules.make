## Rules.make for BFM
## Modified after the GOTM Makefile (Karsten Bolding & Hans Burchard)

#SHELL   = /bin/sh



#F77=ifort
#DEFINES += -DFORTRAN95
#can_do_F90=true
#MODULES=-module $(MODDIR)
#EXTRAS  = -w95
#PROD_FLAGS  = -O3
#REAL_4B = real\(4\)


## Local settings
## Required environment variables (see bfm_env.sh)
#BFMDIR=../
#FORTRAN_COMPILER=IFORT
#NETCDFINC=/usr/include
#NETCDFLIBDIR=/usr/lib

## The compilation mode is obtained from $COMPILATION_MODE -
#### default production - else debug or profiling
#	compilation=production

## Netcdf specifications
#INCDIRS		+= -I$(NETCDFINC)
## NETCDFLIB	= $(NETCDFLIBNAME)
#NETCDFLIB	= -lnetcdf
#LDFLAGS		+= -L$(NETCDFLIBDIR)

##
## phony targets
##
#.PHONY: clean realclean distclean dummy


#CPP = /lib/cpp

#EXTRA_LIBS += $(NETCDFLIB)

## Directory related settings.

#ifndef BINDIR
#BINDIR	= $(BFMDIR)/bin
#endif

#ifndef LIBDIR
#LIBDIR	= $(BFMDIR)/lib/$(FORTRAN_COMPILER)
#endif

#ifndef MODDIR
#MODDIR	= $(BFMDIR)/modules/$(FORTRAN_COMPILER)
#endif

#INCDIRS	+= -I/usr/local/include -I$(BFMDIR)/include -I$(MODDIR)

## BFM include files and the library
#BFMINCDIR = $(BFMDIR)/src/BFM/include
#INCDIRS		+= -I$(BFMINCDIR)
##DEFINES += -DBFM_NOPOINTERS -DNOT_STANDALONE
## DEFINES += -DNOT_STANDALONE
## FEATURE_LIBS += -lbio$(buildtype)

## Normally this should not be changed - unless you want something very specific.

## The Fortran compiler is determined from the EV FORTRAN_COMPILER - options 
## so far NAG(linux), FUJITSU(Linux), DECF90 (OSF1 and likely Linux on alpha),
## SunOS, PGF90 - Portland Group Fortran Compiler (on Intel Linux), Intel (ifc, ifort)

#include $(BFMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

#DEFINES += -DREAL_4B=$(REAL_4B)

## Sets options for production compilation
#ifeq ($(compilation),production)
#buildtype = _prod
#DEFINES += -DPRODUCTION $(STATIC)
#FLAGS   = $(PROD_FLAGS) -assume 2underscores 
#endif

## For making the source code documentation.
#PROTEX	= protex -b -n -s

#.SUFFIXES:
#.SUFFIXES: .F90 .f90

#LINKDIR	= -L$(LIBDIR)

#CPPFLAGS	= $(DEFINES) $(INCDIRS)
#FFLAGS          = $(DEFINES) $(FLAGS) 
#MDFLAGS         = $(MODULES) $(INCDIRS) 
#F90FLAGS  	= $(FFLAGS)
#LDFLAGS		+= $(FFLAGS) $(LINKDIR)

##
## Common rules
##
ifeq  ($(can_do_F90),true)
%.o: %.F90
	$(F77) $(F90FLAGS) $(EXTRAS) $(MDFLAGS) $(EXTRA_FFLAGS) -c $< -o $@
%.o: %.f90 
	$(F77) $(F90FLAGS) $(EXTRAS) $(MDFLAGS) $(EXTRA_FFLAGS) -c $< -o $@
else
%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< -o $@
	$(F90_to_f90)
%.o: %.f90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
endif
