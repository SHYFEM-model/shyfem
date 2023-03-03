
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
FFLAGS = -fcheck=all -g
FFLAGS = -g -O

MAKEDEPEND = ../../bin/mkdp.pl

#--------------------------------------------------------------
# please adjust directory where to look for libnetcdff.a
#--------------------------------------------------------------

NETCDFDIR := /usr

  ifneq ("$(wildcard $(NETCDFDIR)/lib/libnetcdff.a)","")
    DIRLIB_NETCDF = $(NETCDFDIR)/lib/
  else ifneq ("$(wildcard $(NETCDFDIR)/lib/x86_64-linux-gnu/libnetcdff.a)","")
    DIRLIB_NETCDF = $(NETCDFDIR)/lib/x86_64-linux-gnu
  else ifneq ("$(wildcard /opt/sw/netcdf/lib/libnetcdff.a)","")
    NETCDFDIR = /opt/sw/netcdf
    DIRLIB_NETCDF = $(NETCDFDIR)/lib/
  else ifneq ("$(wildcard /usr/local/netcdf/lib/libnetcdff.a)","")
    NETCDFDIR = /usr/local/netcdf
    DIRLIB_NETCDF = $(NETCDFDIR)/lib/
  else
    $(info "cannot find netcdf library libnetcdff.a ... aborting")
  endif
  #$(info DIRLIB_NETCDF = $(DIRLIB_NETCDF))

#--------------------------------------------------------------
# please do not change anything beyond here
#--------------------------------------------------------------

LIBNETCDF := netcdff

LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)
LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a

LIBGS      = $(LIBG_NETCDF)
LIBFS      = $(LIBF_NETCDF)

#--------------------------------------------------------------

OBJS = ncf_util.o subclo.o ncf_dim_coords.o subscn.o
OBJS = ncf_util.o subclo.o
OBJSBGC = ncf_util.o subclo.o ncf_dim_coords.o subscn.o bgc_util.o subgeo.o
OBJSTIME = ncf_tutil.o subiso8601.o subdts.o
OBJSBGC += $(OBJSTIME)

EXES = gr gw ncinfo
EXES = ncinfo ncbgc

DATA = $(HOME)/work/nc2fem/data

#--------------------------------------------------------------

default: $(EXES)

all: link $(EXES)

ncinfo: $(OBJS) $(LIBFS) ncinfo.o
	gfortran -o $@ $@.o $(OBJS) $(LIBGS)

ncbgc: $(OBJSBGC) $(LIBFS) ncbgc.o
	gfortran -o $@ $@.o $(OBJSBGC) $(LIBGS)

gr: $(OBJS) $(LIBFS) gr.o
	gfortran -o $@ $@.o $(OBJS) $(LIBGS)

link: data
data:
	@if [ -e $(DATA) ]; then ln -s $(DATA) data; else true; fi

run: ncbgc
	./run.sh

compile_test:
	./compile_test.sh

depend:
	@$(MAKEDEPEND) *.f

#--------------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f $(EXES)
	-rm -f out.nc
	-rm -f fort.*
	-rm -f *.bak *.tmp
	-rm -f fake.f
	-rm -f data

cleanzip:
	-rm -f $(DIR).zip

cleanall: clean cleanzip

zip: save
save: clean
	zip $(DIR).zip Makefile *.f *.inc

#--------------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

#--------------------------------------------------------------


# DO NOT DELETE THIS LINE -- make depend depends on it.

bgc_util.o: ncf.mod
clo.mod: subclo.o
dts.mod: subdts.o
gr.o: ncf.mod
iso8601.mod: subiso8601.o
mod_ncinfo.mod: ncbgc.o
mod_ncinfo.mod: ncinfo.o
nc_time.mod: ncf_tutil.o
ncbgc.o: clo.mod ncf.mod
ncf.mod: ncf_util.o
ncf_dim_coords.o: ncf.mod
ncf_tutil.o: iso8601.mod ncf.mod
ncf_util.o: netcdf.inc
ncinfo.o: clo.mod ncf.mod
ncnames.mod: ncf_dim_coords.o
scan_string.mod: subscn.o

