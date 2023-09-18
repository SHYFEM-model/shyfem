#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

include ../Rules.make

#--------------------	VARS	---------------------

EXEC = shympi

#----------------------	INCLUDE	----------------------

INCLUDE = -I$(FEMDIR)/include -I$(PETSCDIR)/include

#--------------------	GOTM	-----------------------
ifeq ($(GOTM),true)
  LIBF_GOTM      = $(DIRLIB)/libgotm.a
endif

#---------------------	PETSC	-----------------------
ifeq ($(GOTM),true)
  LIBG_PETSC = -L$(PETSCDIR)/lib -lpetsc
endif

#----------------------	LIBRARY	------------------------
LIBGS= $(LIBG_GOTM) $(LIBG_PETSC)
LIBFS= $(LIBF_GOTM)

#----------------------------------------------------------------------------------------------------------------

OBJ = convert.o utility.o para.o time_util.o levels.o fil.o basin.o system_matrix.o mpi_communication_struct.o \
      mpi_common_struct.o timing.o geom.o lnku.o tvd.o mpi_communication.o shympi_elem_internal.o shympi_elem.o \
      links.o hydro_admin.o evgeom.o regular.o openmp_admin.o mpi_graph_elem.o bnd_geom.o \
      meteo.o waves.o turbulence_util.o sinking.o nudging.o internal.o geom_dynamic.o depth.o \
      bnd_aux.o gotm_aux.o diff_aux.o bnd_dynamic.o area.o ts.o roughness.o diffusion.o hydro_baro.o \
      hydro_print.o hydro_vel.o dts.o femfile.o tsfile.o sigma.o intp.o defnames.o mpi_io_elem.o output.o \
      version.o ionos.o conz_util.o vert_intp.o nls.o bnd.o fem_util.o bnd_admin.o intp_fem_file.o \
      bnd_scalar.o projection.o coordinates.o subset.o petsc_admin.o tidef.o layer_thickness.o bndo.o \
      conz_common.o bclfix_util.o nohyd.o outputd.o elems_dealing.o bndo_admin.o heat_util.o meteo_forcing.o \
      heat_util2.o heat_admin.o topological.o heat_temp.o heat_mfsbulk.o heat_coare.o heat_gotm.o \
      heat_areg.o heat_gill.o heat_default.o chk_NaN.o meteo_admin.o heat_admin2.o sigma_admin.o \
      bit10.o ecological.o restart.o transforms.o tsfile_admin.o fluidmud.o mud_admin.o chezy.o line_admin.o \
      discharge_flux.o flux_util.o ioflux.o flux.o geometrical.o volume.o ioext.o random.o depth_util.o \
      vec_util.o extra.o ets.o ioets.o waves_admin.o ets_admin.o modls.o ioous.o ioshy.o shy_util.o \
      ous_admin.o flux_box.o residual.o tvd_admin.o check.o wetdry.o stability.o time_admin.o sedim_util.o \
      model3d_util.o trace.o keps_util.o kepsilon.o eq_state.o offline_data.o shy_turbulence.o intp_tst.o \
      closing.o nohydro.o nos_util.o diffusion_admin.o bcvel.o lagrange_data.o lagrange_sedim.o \
      lagrange_flux.o lagrange_util.o random_gen.o int_to_cart.o lagrange_vertical.o lagrange_connect.o \
      lagrange_inout.o lagrange_cont.o lagrange_larve.o lagrange_decay.o lagrange_init.o lagrange_util_tr.o \
      lagrange_track.o lagrange_dif.o ioeos.o lagrange_main.o lagrange_back.o nudge.o topological_admin.o \
      explicit.o sort_util.o coo_matrix.o spk_blas1_mod.o matrix_inv.o spk_matvec.o \
      spk_ilut.o spk_itaux.o spk_iters.o sparskit_admin.o system_spk.o hydrodynamic.o submir.o initialize.o \
      bnd_routines.o concentration_omp.o histogram.o concentration.o baroclinic.o def_para.o sedim_admin.o \
      iostr.o constants.o f_vol.o mkdistance.o conz_admin.o water_ren_time.o toxical.o custom_admin.o \
      befor_after.o #shympi.o

#----------------------------------------------------------------------------------------------------------------

%.o: %.f 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

%.o: %.f90 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

$(EXEC):$(MODULES) $(OBJ) $(LIBFS) shympi.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBGS) 
	cp $(EXEC) ../bin/	
	ar r $(DIRLIB)/libshy.a $(OBJ)
clean:
	rm $(MODDIR)/*.mod *.o
	rm $(EXEC)
	rm $(DIRLIB)/libshy.a 
#----------------------------------------------------------------------------------------------------------------
