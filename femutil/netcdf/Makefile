
#--------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#--------------------------------------------------------------------------

#--------------------------------------------------------------

DIR = ncf

#--------------------------------------------------------------

F77 = gfortran
FFLAGS = -g

NETCDFDDD := $(strip $(NETCDFDIR))
NETCDFDDD := /usr
LIBNETCDF  = netcdff

  ifneq ("$(wildcard $(NETCDFDDD)/lib/libnetcdff.a)","")
    DIRLIB_NETCDF = $(NETCDFDDD)/lib/
  else ifneq ("$(wildcard $(NETCDFDDD)/lib/x86_64-linux-gnu/libnetcdff.a)","")
    DIRLIB_NETCDF = $(NETCDFDDD)/lib/x86_64-linux-gnu
  else ifneq ("$(wildcard /opt/sw/netcdf/lib/libnetcdff.a)","")
    NETCDFDDD = /opt/sw/netcdf
    DIRLIB_NETCDF = $(NETCDFDDD)/lib/
  else ifneq ("$(wildcard /usr/local/netcdf/lib/libnetcdff.a)","")
    NETCDFDDD = /usr/local/netcdf
    DIRLIB_NETCDF = $(NETCDFDDD)/lib/
  else
    $(info "cannot find netcdf library libnetcdff.a ... aborting")
  endif
  #$(info DIRLIB_NETCDF = $(DIRLIB_NETCDF))

LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)
LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a

LIBGS      = $(LIBG_NETCDF)
LIBFS      = $(LIBF_NETCDF)

OBJS = ncf_util.o subclo.o

EXES = gr gw ncinfo
EXES = gr ncinfo

DATA = $(HOME)/work/nc2fem/data

#--------------------------------------------------------------

default: all

all: link $(EXES)

ncinfo: $(OBJS) $(LIBFS) ncinfo.o
	gfortran -o $@ $@.o $(OBJS) $(LIBGS)

gr: $(OBJS) $(LIBFS) gr.o
	gfortran -o $@ $@.o $(OBJS) $(LIBGS)

link: data
data:
	@if [ -e $(DATA) ]; then ln -s $(DATA) data; else true; fi

#--------------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f $(EXES)
	-rm -f out.nc
	-rm -f ncf.zip
	-rm -f data

zip: save
save: clean
	zip $DIR).zip Makefile *.f *.inc

#--------------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

#--------------------------------------------------------------

