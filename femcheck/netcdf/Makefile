
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

SHYFEMDIR=../..

FEMBIN = $(SHYFEMDIR)/fembin

include $(SHYFEMDIR)/Rules.make

NETCDF = true

ifeq ($(NETCDF),true)
  NCD := $(strip $(NETCDFDIR))
  FC := $(strip $(FORTRAN_COMPILER))
  NCSETPROG := $(SHYFEMDIR)/femcheck/netcdf/set_netcdf.sh
  DIRLIB_NETCDF := $(shell $(NCSETPROG) lib $(FC) $(NCD))
  DIRINC_NETCDF := $(shell $(NCSETPROG) include $(FC) $(NCD))
endif
$(info DIRLIB_NETCDF = $(DIRLIB_NETCDF) )
$(info DIRINC_NETCDF = $(DIRINC_NETCDF) )
ifeq ($(DIRLIB_NETCDF),unknown)
  $(info DIRLIB_NETCDF = $(DIRLIB_NETCDF) )
  $(info DIRINC_NETCDF = $(DIRINC_NETCDF) )
  $(error cannot find DIRLIB_NETCDF )
endif

EXENC =
NETCDF_OBJ =
NETCDF_ELAB = elab_nc_dummy.o
NETCDF_FEM =
LIBF_NETCDF =
NETCDF_2FEM = netcdf.o nc_util.o nc_tutil.o nc_domain.o nc_dim_coords.o \
              subbas.o subdts.o subscn.o subfil.o subclo.o \
              substrings.o subsss.o subround.o \
              subfemfile.o subiso8601.o mod_depth.o subnev.o newlevels.o \
              subreg.o elab_expand.o $(MPI_OBJ)
NETCDF_EXE = ts2nc nc2fem

ifeq ($(NETCDF),true)
  NETCDF_OBJ       = netcdf.o netcdf_util.o
  NETCDF_ELAB      = netcdf.o netcdf_util.o elab_nc.o
  NETCDF_FEM       = netcdf.o nc_util.o nc_tutil.o
  LIBG_NETCDF      = $(DIRLIB_NETCDF)
  LIBF_NETCDF      =
  #LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDFF) -l$(LIBNETCDFC)
  #LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDFF).a \
#                       $(DIRLIB_NETCDF)/lib$(LIBNETCDFC).a
  EXENC            = $(NETCDF_EXE)
endif

FFLAGS = $(DIRINC_NETCDF)

EXES = nc_include nc_mod

#--------------------------------------------------

default: $(EXES)

test: test_netcdf
test_netcdf:
	./test_netcdf.sh

nc_include: nc_include.o
	$(F77)  $@.o $(LIBG_NETCDF) -o $@

nc_mod: nc_mod.o
	$(F77)  $@.o $(LIBG_NETCDF) -o $@

clean:
	-rm -f $(EXES)
	-rm -f a.out
	-rm -f *.o

cleanall: clean

# Rules ----------------------------------------------

.SUFFIXES: .F90 .f90

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.f.o:
	$(F77) -c $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

.f90.o:
	$(F77) -c $(FFLAGS) $<

#--------------------------------------------------

