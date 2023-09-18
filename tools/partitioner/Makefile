#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

include ../../Rules.make

#----------------------	INCLUDE	----------------------

INCLUDE = -I$(FEMDIR)/../include

#--------------------   PATH   ---------------------------

DIRLIB=$(FEMDIR)/../libs

#---------------------	LIBRARY	-----------------------

LIBS = -L$(DIRLIB) -lshy -L/zeus/opt/intel19.5/metis/5.1.0/lib -lmetis

#----------------------------------------------------------------------------------------------------------------

OBJ = initialize_par.o util_part.o partitioning.o bnd_admin_par.o iostr_par.o constants_par.o shypart.o

#OBJ = levels.o fil.o basin.o geom.o lnku.o mpi_communication_struct.o mpi_common_struct.o \
      convert.o utility.o para.o def_para.o defnames.o sigma.o sigma_admin.o links.o bnd.o shypart_elem.o \
      bnd_geom.o fem_util.o nls.o bnd_admin_par.o initialize_par.o mpi_global_graph_elem.o shypart_util.o \
      mpi_partitioning_elem.o line_admin.o flux_par.o geometrical.o volume_par.o iostr_par.o constants_par.o shypart.o

#----------------------------------------------------------------------------------------------------------------

%.o: %.f 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

%.o: %.f90 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

all: shypart

shypart: $(MODULES) $(OBJ)
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBS)
	cp $@ ../../bin/
clean:
	rm *.o
	rm shypart
#----------------------------------------------------------------------------------------------------------------
