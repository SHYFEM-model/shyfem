## Rules.make for BFM
## Modified after the GOTM Makefile (Karsten Bolding & Hans Burchard)

SHELL   = /bin/sh

## Local settings
## Required environment variables (see bfm_env.sh)
BFMDIR=../
#FORTRAN_COMPILER=IFORT
NETCDFINC=/usr/local/include
NETCDFLIBDIR=/usr/local/lib

## The compilation mode is obtained from $COMPILATION_MODE -
#### default production - else debug or profiling
ifndef COMPILATION_MODE
	compilation=production
else
	compilation=$(COMPILATION_MODE)
endif

FEATURES	=
FEATURE_LIBS	=
EXTRA_LIBS	=
INCDIRS		=
LDFLAGS		=

## Netcdf specifications
INCDIRS		+= -I$(NETCDFINC)
## NETCDFLIB	= $(NETCDFLIBNAME)
NETCDFLIB	= -lnetcdf
LDFLAGS		+= -L$(NETCDFLIBDIR)

##
## phony targets
##
.PHONY: clean realclean distclean dummy

## Default root directory BFM
ifndef BFMDIR
BFMDIR  := $(HOME)
endif

CPP = /lib/cpp

EXTRA_LIBS += $(NETCDFLIB)

## Directory related settings.

ifndef BINDIR
BINDIR	= $(BFMDIR)/bin
endif

ifndef LIBDIR
LIBDIR	= $(BFMDIR)/lib/$(FORTRAN_COMPILER)
endif

ifndef MODDIR
MODDIR	= $(BFMDIR)/modules/$(FORTRAN_COMPILER)
endif

INCDIRS	+= -I/usr/local/include -I$(BFMDIR)/include -I$(MODDIR)

## BFM include files and the library
BFMINCDIR = $(BFMDIR)/src/BFM/include
INCDIRS		+= -I$(BFMINCDIR)
##DEFINES += -DBFM_NOPOINTERS -DNOT_STANDALONE
## DEFINES += -DNOT_STANDALONE
## FEATURE_LIBS += -lbio$(buildtype)

## Normally this should not be changed - unless you want something very specific.

## The Fortran compiler is determined from the EV FORTRAN_COMPILER - options 
## so far NAG(linux), FUJITSU(Linux), DECF90 (OSF1 and likely Linux on alpha),
## SunOS, PGF90 - Portland Group Fortran Compiler (on Intel Linux), Intel (ifc, ifort)

include $(BFMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

DEFINES += -DREAL_4B=$(REAL_4B)

## Sets options for debug compilation
ifeq ($(compilation),debug)
buildtype = _debug
DEFINES += -DDEBUG $(STATIC)
FLAGS   = $(DEBUG_FLAGS) 
endif

## Sets options for profiling compilation
ifeq ($(compilation),profiling)
buildtype = _prof
DEFINES += -DPROFILING $(STATIC)
FLAGS   = $(PROF_FLAGS) 
endif

## Sets options for production compilation
ifeq ($(compilation),production)
buildtype = _prod
DEFINES += -DPRODUCTION $(STATIC)
FLAGS   = $(PROD_FLAGS) 
endif

## For making the source code documentation.
PROTEX	= protex -b -n -s

.SUFFIXES:
.SUFFIXES: .F90 .f90

LINKDIR	= -L$(LIBDIR)

CPPFLAGS	= $(DEFINES) $(INCDIRS)
FFLAGS  	= $(DEFINES) $(FLAGS) $(MODULES) $(INCDIRS) $(EXTRAS)
F90FLAGS  	= $(FFLAGS)
LDFLAGS		+= $(FFLAGS) $(LINKDIR)

##
## Common rules
##
ifeq  ($(can_do_F90),true)
%.o: %.F90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
%.o: %.f90 
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
else
%.f90: %.F90
##	$(CPP) $(CPPFLAGS) $< -o $@
	$(F90_to_f90)
%.o: %.f90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
endif
