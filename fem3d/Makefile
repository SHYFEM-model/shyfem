#
# $Id: Makefile,v 1.81 2010-03-11 15:36:38 georg Exp $
#
#------------------------------------------------------------

DIR	= fem3d

#------------------------------------------------------------ just for default
FEMDIR = ..
FEMBIN = $(FEMDIR)/fembin
DIRLIB = $(FEMDIR)/femlib
F77    = gfortran
FFLAGS =  -O
LINKER = $(F77)
#------------------------------------------------------------

include ../Rules.make

#------------------------------------------------------------

LIBFEM     = fem
LIBGOTM    = gotm
LIBNETCDF  = netcdff
LIBERSEM   = bfm

LIBFFEM    = $(DIRLIB)/lib$(LIBFEM).a

LIBLS      = -L$(DIRLIB)

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#------------------------------------------------------------
# fluid mud files
#------------------------------------------------------------

ifeq ($(FLUID_MUD),true)
  MUD_FILES = submud.o subgotm_mud.o
else
  MUD_FILES = submud_dummy.o subgotm.o
endif

#------------------------------------------------------------
# solver files
#------------------------------------------------------------

LIBG_MKL =
SOLVER_FILES = simsys_lp.o
ifeq ($(SOLVER),PARDISO)
  LIBG_MKL = `./bin/find_pardlib.sh`
  SOLVER_FILES = simsys_pard.o subpard.o subcoo.o
endif
ifeq ($(SOLVER),SPARSKIT)
  SOLVER_FILES = simsys_spk.o subspk.o subcoo.o subcooutil.o \
		spk_ilut.o spk_iters.o spk_itaux.o \
		spk_itpack.o spk_matvec.o spk_blas1_mod.o
endif

#------------------------------------------------------------
# netcdf files
#------------------------------------------------------------

NETCDFDDD := $(strip $(NETCDFDIR))
ifeq ($(NETCDF),true)
  ifneq ("$(wildcard $(NETCDFDDD)/lib/libnetcdff.a)","")
    DIRLIB_NETCDF = $(NETCDFDDD)/lib/
  else ifneq ("$(wildcard $(NETCDFDDD)/lib/x86_64-linux-gnu/libnetcdff.a)","")
    DIRLIB_NETCDF = $(NETCDFDDD)/lib/x86_64-linux-gnu
  else
    $(info "cannot find netcdf library libnetcdff.a ... aborting")
  endif
  #$(info DIRLIB_NETCDF = $(DIRLIB_NETCDF))
endif

EXENC = 
NETCDF_OBJ = 
NETCDF_ELAB = elab_nc_dummy.o
NETCDF_FEM = 
LIBF_NETCDF =
NETCDF_2FEM = netcdf.o nc_util.o nc_tutil.o nc_domain.o nc_dim_coords.o \
	      subbas.o subdts.o subscn.o subfil.o subclo.o \
	      substrings.o subsss.o \
              subfemfile.o subiso8601.o mod_depth.o subnev.o newlevels.o \
	      subreg.o elab_expand.o $(MPI_OBJ)
NETCDF_EXE = ts2nc nc2fem
DIRINC_NETCDF = $(NETCDFDIR)/include
ifeq ($(NETCDF),true)
  NETCDF_OBJ       = netcdf.o netcdf_util.o
  NETCDF_ELAB      = netcdf.o netcdf_util.o elab_nc.o
  NETCDF_FEM       = netcdf.o nc_util.o nc_tutil.o
  LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)
  LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a
  EXENC            = $(NETCDF_EXE)
endif

#------------------------------------------------------------
# mpi files
#------------------------------------------------------------

#PHGOBJ = mod_shypart_dummy.o
#MPI_OBJ = mod_shympi_node.o mod_shympi_dummy.o basin_mpi.o
MPI_DUMMY_OBJ = mod_shympi_dummy.o
MPI_NODE_OBJ = mod_shympi_node.o mod_shympi_node_internal.o
MPI_GENERAL_OBJ = submpi_ghost.o submpi_basin.o submpi_custom.o \
		submpi_buffer.o submpi_time.o
MPI_SHYFEM_OBJ = submpi_assert.o

MPI_OBJ = $(MPI_DUMMY_OBJ) $(MPI_GENERAL_OBJ)
ifeq ($(PARALLEL_MPI),NODE)
  MPI_OBJ = $(MPI_NODE_OBJ) $(MPI_GENERAL_OBJ)
  FFLAGS += $(MPI_COMPILE_FLAGS)
  LFLAGS += $(MPI_LINK_FLAGS)
endif


#------------------------------------------------------------
# omp files
#------------------------------------------------------------

OMP_OBJ = subomp_dummy.o
ifeq ($(PARALLEL_OMP),true)
  OMP_OBJ = subomp.o
endif

#------------------------------------------------------------
# gotm files
#------------------------------------------------------------

GOTM_OBJ = gotm3d.o
ifeq ($(GOTM),true)
  GOTM_OBJ =
  LIBG_GOTM      = -l$(LIBGOTM)
  LIBF_GOTM      = $(DIRLIB)/lib$(LIBGOTM).a
endif

#------------------------------------------------------------
# ecological module files
#------------------------------------------------------------

LIBF_ERSEM =
ECOLOGICAL_OBJ = ecological_dummy.o
ifeq ($(ECOLOGICAL),EUTRO)
  ECOLOGICAL_OBJ = bio3d.o bio3d_util.o \
		weutro.o weutro_sedim.o weutro_light.o \
		weutro_seed.o weutro_shell.o weutro_ulva.o
endif
ifeq ($(ECOLOGICAL),ERSEM)
  ECOLOGICAL_OBJ = subbfm.o
  LIBG_ERSEM     = -l$(LIBERSEM)
  LIBF_ERSEM     = $(DIRLIB)/lib$(LIBERSEM).a
endif
ifeq ($(ECOLOGICAL),AQUABC)
  ECOLOGICAL_OBJ = aquabc_II_pelagic_lib.o aquabc_II_fem_interface.o \
		aquabc_II_co2sys.o aquabc_II_apara.o aquabc_II_util.o \
		aquabc_II_sediment_model_1.o aquabc_II_sediment_lib.o \
		aquabc_II_pelagic_model.o aquabc_II.o \
		aquabc_II_ased_ini.o
endif

#------------------------------------------------------------
# all libraries
#------------------------------------------------------------

LIBGS      = $(LIBLS) $(LIBG_NETCDF) $(LIBG_GOTM) $(LIBG_ERSEM) $(LIBG_MKL)
LIBFS      = $(LIBF_NETCDF) $(LIBF_GOTM) $(LIBF_ERSEM)
INCLUDE    = $(INCL_NETCDF)

#------------------------------------------------------------

SPECIAL =	Makefile \
		COMMIT README TESTLOG TODO VERSION WARNINGS \
		ACTUAL LOG

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------------------------------------------

# 3D objects for 3D model

GENERAL_OBJS = subver.o subclo.o

MODULES = newlevels.o \
		mod_hydro.o mod_hydro_vel.o \
		mod_hydro_print.o mod_hydro_baro.o \
		mod_diff_visc_fric.o mod_roughness.o \
		mod_ts.o mod_area.o mod_bound_dynamic.o \
		mod_diff_aux.o mod_gotm_aux.o \
		mod_bnd_aux.o mod_depth.o mod_geom_dynamic.o \
		mod_internal.o mod_nohyd.o mod_nudging.o \
		mod_bclfix.o mod_turbulence.o \
		mod_fluidmud.o mod_sinking.o mod_waves.o \
		mod_meteo.o mod_geom.o mod_conz.o \
		mod_bound_geom.o mod_bnd.o \
		mod_tvd.o mod_bndo.o mod_system.o \
		mod_layer_thickness.o \
		mod_subset.o mod_lagrange.o

NEWOBJS = \
	  new3di.o \
	  newrog.o newcra.o \
	  newbcl.o newbcl0.o \
	  newini.o newtra.o \
	  newchao.o \
	  newcon.o newconz.o newcon_omp.o \
	  newchk.o newfix.o \
	  newtvd.o newstab.o newsig.o \
	  newnudge.o newpoi.o newbfm.o \
	  $(MODULES) \
	  $(MPI_OBJ) \
	  $(MPI_SHYFEM_OBJ)

# objects for library

SOBJS   = \
	  subver.o \
	  subnsv.o subsst.o \
          subcmv.o subdry.o subssv.o \
          $(PARMOD) \
          sublst.o subsys.o \
          subcst.o submir.o subros.o \
          subn11.o subsrt.o \
          subdef.o subext.o subssc.o \
          subbas.o subnev.o \
          subnls.o subwat.o \
          subfil.o subn35.o subout.o subous.o subssg.o \
          subssm.o \
          subnsu.o subsss.o subscn.o subres.o \
          subapn.o subnsa.o \
          subreg.o subfind.o submask.o \
          subbnd.o \
          subexta.o subdep.o subflxa.o subflx3d.o subflx.o \
	  subdum.o subcus.o \
	  $(BCTSOBJS) \
	  $(QFLUX) \
	  $(GRDOBJS) \
	  subnos.o subcon1.o subrnd.o sublin.o \
	  submod.o subuti.o subele.o subgeo.o subdts.o \
	  links.o sublnka.o sublnks.o sublnku.o \
	  debug.o subhisto.o subnosa.o subfvl.o \
	  subeos.o sigmautil.o \
	  nosutil.o ousutil.o subfemfile.o subvintp.o \
	  genutil.o subfemintp.o $(TSOBJS) subxi.o \
	  $(MODULES) \
	  subshy.o shyutil.o subshyutil.o subshyutil2.o \
	  shyelab_util.o \
	  subbit10.o shyutil.o subclo.o substrings.o subpar3.o \
	  elabutil.o elabtime.o elab_gis.o elab_check.o subiso8601.o \
	  unitutils.o \
	  $(MPI_OBJ)

# general objects for 3D model shyfem

HTOBJS = \
	subver.o \
	subcst.o subsys.o subnsh.o subtime.o suboutput.o suboutputd.o \
	subbas.o subnev.o \
	subnsu.o subn35.o subdep.o \
	subdef.o subfil.o subnls.o \
	$(PARMOD) \
	$(BCTSOBJS) \
	subsss.o subscn.o subssv.o subssm.o \
	submir.o \
	subext.o subexta.o \
	subous.o subousa.o \
	subbnd.o subnos.o \
	subwat.o subn11.o subdry.o \
	subcus.o subcon1.o \
	subtidef.o subres.o subuti.o \
	sublin.o subvola.o subgeo.o debug.o \
	submod.o subele.o \
	$(SOLVER_FILES) \
	$(QFLUX) \
	$(MUD_FILES) $(GOTM_OBJ) \
	subbndo.o subbnds.o \
	new36.o intp.o subssf.o \
	bdist.o links.o sublnka.o sublnkd.o sublnks.o sublnku.o \
	sedi3d.o sedtrans05.o subwaves.o \
	subhisto.o \
	$(LAGRANGE) \
	subdts.o loading.o \
	subrst.o subnosa.o subdif.o subreg.o subfind.o \
	atoxi3d.o atoxi.o kepsd.o keps_util.o \
	subinv.o newexpl.o subrnd.o \
	submet.o subrgf.o subcoord.o subsrt.o subfvl.o \
	subtrace.o subproj.o subeos.o sigmautil.o \
	subwrt.o subvintp.o \
	subfemfile.o subtsuvfile.o \
	subflxa.o subflx3d.o subflx.o subflxu.o \
	subboxa.o \
	newnohydro.o suboff.o \
	subets.o subetsa.o subetsm.o \
	$(IFF_FILES) \
	$(ECOLOGICAL_OBJ) \
	$(OMP_OBJ) \
	$(GRDOBJS) \
	subxi.o subshy.o subshyutil.o substrings.o \
	subbit10.o shyutil.o \
	subfloodfill.o mercury.o merc_water.o subssed.o \
	subclo.o

	#$(NETCDF_OBJ) \

# objects for pre-processor shypre

VPOBJS = \
	subver.o \
        subros.o subcmv.o \
        subdef.o subbas.o \
        subfil.o \
        $(PARMOD) \
        subsrt.o \
        subsss.o subscn.o subnsu.o \
	subapn.o subnsa.o subnls.o subsys.o subdum.o \
	$(GRDOBJS) subclo.o substrings.o

# objects for output routines

OUTOBJS = $(PARMOD) subnls.o subsss.o subscn.o \
	  subfil.o subdef.o subnsu.o subssv.o \
	  subbas.o subsys.o \
	  subapn.o subnsa.o subdum.o \
	  subdts.o subdep.o \
	  debug.o subreg.o subnev.o \
	  subfemfile.o subfemutil.o \
	  $(TSOBJS) \
	  subvintp.o subver.o sigmautil.o suboutput.o suboutputd.o \
	  subclo.o \
	  mod_depth.o mod_hydro.o mod_geom.o mod_hydro_print.o \
	  mod_bound_geom.o mod_bnd.o mod_hydro_baro.o newlevels.o \
	  $(GRDOBJS) subsrt.o subgeo.o substrings.o subiso8601.o \
	  $(MPI_OBJ)

PARTOBJS = subfloodfill.o
TSOBJS = subtsfile1.o
TSOBJS = subtsfile.o

RSTOBJS = subrst.o ecological_dummy.o \
	  mod_ts.o mod_conz.o mod_hydro_vel.o mod_geom_dynamic.o \
	  subbit10.o

BASOBJS = basutil.o basutil_resol.o basutil_hsigma.o basutil_bathy.o \
	  basutil_interp.o basutil_smooth.o basutil_box.o basutil_part.o \
	  subnos.o \
	  sublnka.o sublnku.o links.o
BASEXTRA = sublnka.o links.o sublnku.o \
		subinterpol.o subinterpol1.o subnos.o

GENOBJS = newlevels.o
NOSOBJS = subnos.o nosutil.o genutil.o subeos.o elabutil.o \
		elabtime.o elab_gis.o subdates.o
SHYOBJS = genutil.o \
		subous.o ousutil.o \
		subnos.o nosutil.o subeos.o \
		subflx.o flxelab1.o \
		subext.o extelab1.o \
		tselab1.o femelab1.o \
		elabutil.o elabtime.o elab_gis.o elab_off.o \
		subshy.o subshyutil.o subshyutil2.o \
		shyutil.o subdates.o
SHYCOMPAT = shyelab_nodes.o shyelab_average.o elab_regular.o elab_expand.o \
		shyelab_util.o unitutils.o elab_check.o
SHYELAB = shyelab1.o shyelab_output.o shyelab_nodes.o \
		shyelab_area.o shyelab_average.o \
		shyelab_util.o unitutils.o \
		line_util.o \
		elab_regular.o elab_expand.o $(NETCDF_ELAB)
SHYNOBJS = genutil.o elabutil.o subshy.o subshyutil.o subshyutil2.o shyutil.o
ETSOBJS = subets.o etsutil.o genutil.o subetsm.o
NOSOUSOBJS = subnos.o nosutil.o subous.o ousutil.o subets.o etsutil.o \
		genutil.o
GRDOBJS = subgrd.o grdutil.o
PROJOBJS = subproj.o
UTILOBJS = subsrt.o
UTILOBJS = 

BATHYOBJS = $(BASEXTRA)
#BATHYOBJS = subgeo.o subgrd.o subsrt.o subinterpol.o subinterpol1.o \
#		subnos.o

LAPLACE = suboi.o submir.o subnos.o links.o sublnku.o \
	  subnosa.o subssm.o subintp.o

LAGRANGE =	lagrange_main.o \
		lagrange_back.o lagrange_dif.o \
		lagrange_cont.o lagrange_init.o \
		lagrange_inout.o lagrange_util.o \
		lagrange_larve.o lagrange_flux.o \
		lagrange_track.o lagrange_util_tr.o \
		lagrange_decay.o lagrange_sedim.o \
		lagrange_connect.o lagrange_vertical.o \
		line_util.o

# heat flux module

QFLUX = subqfxt.o subqfxf.o \
	subqfxm1.o subqfxm2.o subqfxm3.o subqfxm4.o subqfxm5.o \
        subqfxm6.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o

QFLUX_TEST = subqfxf.o \
	subqfxm1.o subqfxm2.o subqfxm3.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o \
	subfil.o

# experimental stuff

#PARMOD = subpar.o subfnm.o		# old parameter modules
#PARMOD = subpar1.o			# newer parameter modules
PARMOD = subpar3.o
IFF_FILES = 
IFF_FILES = subfemintp.o $(TSOBJS) submeteo2.o subiso8601.o

BCTSOBJS = subsse.o
BCTSOBJS = subintp.o subbc.o
BCTSOBJS = subintp.o subbc99.o
BCTSOBJS = subintp.o

# objects for interpolation

INTP    = subreg.o

# objects for lagrangian interpolation

LAGINTP = submir.o

#------------------------------------------------

EXEMAIN  = shyfem shypre shybas

EXEINFO  = rstinf tselab \
		offinf lgrinf \
		filetype \
		extelab flxelab \
		shyelab \
		femelab femintp femcombine
EXESPLIT = 
EXELAP   =  
EXEBAS   = 
EXEVAR   = scal_intp
EXEUTIL  = shyproj
EXECOMPAT = bc2fem nos2shy ous2shy

NOSVAR   = 
NOSAVER  = 
NOSEXTR  = 
NOSALL   = $(NOSEXTR) $(NOSAVER) $(NOSVAR)

OUSEXTR  = 

EXEEXTRA = \
		basproj baswork \
		ts2nc \
		ouswork ouswork1 ouswrite prep_input \
		basbox
EXEAUX   = baschk bascomp check_debug laplap mkgeom \
		optintp ousaver ousvdiff

EXECLEAN = $(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXELAP) $(EXEBAS) \
		$(EXEVAR) $(EXEUTIL) $(NOSEXTR) $(OUSEXTR) $(NOSAVER) \
		$(EXEEXTRA) $(EXENC) $(NOSVAR) $(EXEAUX) $(NETCDF_EXE) \
		$(EXECOMPAT)

EXELINKS = $(FEMBIN) $(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXEUTIL) $(EXENC) \
		$(EXECOMPAT)

#------------------------------------------------

default:	fem

all:	fem
fem:	libfem main info split bas var util nc links
main:	$(EXEMAIN)
info:	$(EXEINFO)
split:	$(EXESPLIT)
bas:	$(EXEBAS)
var:	$(EXEVAR)
nc:	$(EXENC)
util:	$(NOSEXTR) $(OUSEXTR) $(NOSAVER) $(EXEUTIL)
compat: $(EXECOMPAT)
compatibility: compat

links:
	@$(FEMBIN)/symlink $(EXELINKS)

#--------------------------------------------- executables

bedload:  bedload.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS) -o $@

shyfem:  shyfem.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS) -o $@

shypre:  shypre.o $(VPOBJS)
	$(LINKER) $(LFLAGS) $@.o $(VPOBJS) -o $@

#------------------------------------------------

zlevel:  zlevel.o 
	$(LINKER) $(LFLAGS) $@.o -o $@

lagint:  lagint.o sublag.o $(LAGINTP) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o sublag.o $(LAGINTP) $(OUTOBJS) -o $@

lgrinf:  lgrinf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

zinit:  zinit.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

mkgeom:   mkgeom.o
	$(LINKER) $(LFLAGS) $@.o -o $@

cubic:   cubic.o $(OUTOBJS) $(BCTSOBJS) subssm.o
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BCTSOBJS) subssm.o -o $@

#------------------------------------------------

shybas:  shybas.o $(BASOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(BASOBJS) $(OUTOBJS) -o $@

shypart:  shypart.o $(BASOBJS) $(OUTOBJS) $(PARTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(BASOBJS) $(OUTOBJS) $(PARTOBJS) -o $@

basinf:  basinf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

baswork:  baswork.o $(OUTOBJS) $(BASEXTRA)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BASEXTRA) -o $@

basbox:  basbox.o $(OUTOBJS) $(BASEXTRA)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BASEXTRA) -o $@

bashsigma:  bashsigma.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

baschk:  baschk.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

basbathy:  basbathy.o $(BATHYOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(BATHYOBJS) $(OUTOBJS) -o $@

basres:  basres.o $(OUTOBJS) $(NOSOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(NOSOBJS) -o $@

basproj:  basproj.o $(OUTOBJS) $(PROJOBJS) $(UTILOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(PROJOBJS) $(UTILOBJS) -o $@

shyproj:  shyproj.o $(OUTOBJS) $(PROJOBJS) $(UTILOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(PROJOBJS) $(UTILOBJS) -o $@

bastreat:  bastreat.o $(BASEXTRA) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(BASEXTRA) $(OUTOBJS) -o $@

scalintp:  scalintp.o $(OUTOBJS) $(BATHYOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BATHYOBJS) -o $@

scal_intp:  scal_intp.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@

bascomp:  bascomp.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

filetype:  filetype.o $(NOSOUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOUSOBJS) $(OUTOBJS) -o $@

femintp:  femintp.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT)
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT) -o $@

femcombine:  femcombine.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT)
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT) -o $@

bc2fem:  bc2fem.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT) -o $@

check_debug: check_debug.o
	$(LINKER) $(LFLAGS) $@.o -o $@

2nc: ts2nc nc2fem

OBJS_NC = $(SHYOBJS) $(OUTOBJS) $(SHYCOMPAT) $(NETCDF_OBJ) $(LIBGS)

ts2nc:  ts2nc.o $(SHYOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) $(LFLAGS) $@.o $(OBJS_NC) -o $@

nc2fem:  nc2fem.o $(NETCDF_2FEM)
	$(LINKER) $(LFLAGS) $@.o $(NETCDF_2FEM) $(LIBGS) -o $@

prep_input:  prep_input.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

offinf: offinf.o $(GENERAL_OBJS)
	$(LINKER) $(LFLAGS) $@.o $(GENERAL_OBJS) -o $@

#TS_OBJS = subscn.o subdts.o
#ts2ts: ts2ts.o $(TS_OBJS)
#	$(LINKER) $(LFLAGS) $@.o $(TS_OBJS) -o $@

#-------------------------------------------------------

.PHONY: ousextr nosaver nos

nos: $(NOSALL)

NOSOUS_ALL = $(SHYOBJS) $(OUTOBJS) $(SHYCOMPAT)

nos2shy: nos2shy.o $(NOSOUS_ALL)
	$(LINKER) $(LFLAGS) $@.o $(NOSOUS_ALL) -o $@

ous2shy: ous2shy.o $(NOSOUS_ALL)
	$(LINKER) $(LFLAGS) $@.o $(NOSOUS_ALL) -o $@

#femelab:  femelab.o femelab1.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT)
#	$(LINKER) $(LFLAGS) $@.o femelab1.o $(NOSOBJS) $(OUTOBJS) $(SHYCOMPAT) -o $@

femelab:  femelab.o $(NOSOUS_ALL)
	$(LINKER) $(LFLAGS) $@.o $(NOSOUS_ALL) -o $@

extelab: extelab.o $(NOSOUS_ALL)
	$(LINKER) $(LFLAGS) $@.o $(NOSOUS_ALL) -o $@

flxelab: flxelab.o $(NOSOUS_ALL)
	$(LINKER) $(LFLAGS) $@.o $(NOSOUS_ALL) -o $@

tselab:  tselab.o $(NOSOUS_ALL)
	$(LINKER) $(LFLAGS) $@.o $(NOSOUS_ALL) -o $@

SHYELAB_EXTRA = shyelab1.o shyelab_output.o shyelab_area.o \
		line_util.o $(NETCDF_ELAB)
SHYELAB_ALL = $(SHYOBJS) $(OUTOBJS) $(SHYCOMPAT) $(SHYELAB_EXTRA)

shyelab: shyelab.o $(SHYELAB_ALL) $(LIBFS)
	$(LINKER) $(LFLAGS) $@.o $(SHYELAB_ALL) $(LIBGS) -o $@

#shydiff: shydiff.o $(SHYNOBJS) $(OUTOBJS)
#	$(LINKER) $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) -o $@

#-------------------------------------------------------

lgrelab:  lgrelab.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

optintp:   optintp.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@ 

rstinf:   rstinf.o $(OUTOBJS) $(RSTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(RSTOBJS) -o $@ 

qflux_test:   qflux_test.o $(QFLUX_TEST)
	$(LINKER) $(LFLAGS) $@.o $(QFLUX_TEST) -o $@ 

#------------------------------------------------

# to test interpolation routines in intp.f

INTPOBJS = $(BCTSOBJS) subfil.o subsss.o subscn.o subssf.o subssm.o

intp:	intp.o $(INTPOBJS)
	$(LINKER) $(LFLAGS) $@.o $(INTPOBJS) -o $@ 

#------------------------------------------------

# extra treatment of GOTM routines because of compiler bug

gotm3d.o: gotm3d.f gotmturb.i
	$(F77) -c $(FFLAGS_NOOPT) gotm3d.f

#------------------------------------------------

libs: libfem
femlib: libfem
libfem: $(LIBFFEM)

$(LIBFFEM): $(SOBJS)
	ar rvs $(LIBFFEM) $?

#------------------------------------------------

headers:
	header *.[fh]

#------------------------------------------------

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f *.f90

#---------------------------------------------------

options:
	@echo "compiler  $(COMPILER) $(F77)"
	@echo "solver    $(SOLVER)"

#---------------------------------------------------

clean: cleanobj cleanexe cleanvar cleantex cleanlahey cleanexport cleandoc \
	cleanifc cleanlinks

cleanall: clean cleanzip cleanlib cleanextra cleandiff cleanmod

cleanlinks:
	@$(FEMBIN)/symlink -remove $(EXELINKS)

cleanifc:
	-rm -f ifc*

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanexe:
	-rm -f $(EXECLEAN)
	-rm -f *.exp
	-rm -f *.exe

cleandoc:
	-rm -f *.html
	-rm -f head.txt revlog.txt

cleantex:
	-rm -f *.dvi *.aux

cleanlahey:
	-rm -fr lahey

cleanexport:
	-rm -fr export

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp *.gnu
	-rm -f *.lst *.map
	-rm -f *~
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f errout.dat
	-rm -f out.ps out.*.ps
	-rm -f .memory
	-rm -f html/cvs*
	-rm -fr cvsdiff_tmp
	-rm -f *.rev *.dcl
	-rm -f CVSSTATUS CVSRAWSTATUS
	-rm -f COMMONF COMMONL
	-rm -f gmon.out tags
	-rm -f main_dummy.f

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleanlib:
	-rm -f $(LIBFFEM)

cleanmod:
	-rm -f $(MODDIR)/*.mod

cleanextra:
	-rm -f LASTZIP

#--------------------------------------------------

lahey:	cleanall
	bin/mklahey.sh
	cp $(SPECIAL) lahey
	bin/mklmake.pl Makefile > lahey/Makefile
	cp Makefile lahey/Makefile.org
	bin/unix2other.pl -lahey lahey/*.for
	cd lahey; rm -f *.for.bak
	cd lahey; zip -l lahey *

export:	cleanall
	bin/mkexport.sh
	cp $(SPECIAL) export
	cd export; zip export *

zip: save

save:   cleanall
	date > LASTZIP
	zip $(DIR) *.[fhi]
	zip $(DIR) $(SPECIAL)
	zip $(DIR) bin/*

param:

list:
	pwd; ls

last:
	ls -latrd [A-Z]* Makefile *.[fh]

#
# Rules ----------------------------------------------
#

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

# DO NOT DELETE THIS LINE -- make depend depends on it.

../femlib/mod/AQUABC_II_GLOBAL.mod: aquabc_II_co2sys.o
../femlib/mod/CO2SYS_CDIAC.mod: aquabc_II_co2sys.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_dummy.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/PROCEDURE.mod: mod_shympi_node.o
../femlib/mod/VECTOR_MATRIX_UTILS.mod: aquabc_II_co2sys.o
../femlib/mod/aquabc_II_sed_ini.mod: aquabc_II_ased_ini.o
../femlib/mod/basin.mod: subbas.o
../femlib/mod/basutil.mod: basutil.o
../femlib/mod/box.mod: subboxa.o
../femlib/mod/box_arrays.mod: subboxa.o
../femlib/mod/check_unit.mod: newchk.o
../femlib/mod/chezy.mod: subn35.o
../femlib/mod/clo.mod: subclo.o
../femlib/mod/coordinates.mod: subcoord.o
../femlib/mod/custom_dates.mod: subdates.o
../femlib/mod/default_names.mod: subdef.o
../femlib/mod/dts.mod: subdts.o
../femlib/mod/elabtime.mod: elabtime.o
../femlib/mod/elabutil.mod: elabutil.o
../femlib/mod/ets.mod: subetsm.o
../femlib/mod/eutro.mod: bio3d.o
../femlib/mod/evgeom.mod: subnev.o
../femlib/mod/extfile.mod: subext.o
../femlib/mod/extra.mod: subexta.o
../femlib/mod/fem_file.mod: subfemfile.o
../femlib/mod/fem_util.mod: subfemutil.o
../femlib/mod/flux.mod: subflxa.o
../femlib/mod/flxfile.mod: subflx.o
../femlib/mod/grd.mod: subgrd.o
../femlib/mod/intp_fem_file.mod: subfemintp.o
../femlib/mod/iso8601.mod: subiso8601.o
../femlib/mod/levels.mod: newlevels.o
../femlib/mod/lgr_sedim_module.mod: lagrange_sedim.o
../femlib/mod/lines.mod: sublin.o
../femlib/mod/mercury.mod: mercury.o
../femlib/mod/meteo_forcing_module.mod: submeteo2.o
../femlib/mod/mod_area.mod: mod_area.o
../femlib/mod/mod_bclfix.mod: mod_bclfix.o
../femlib/mod/mod_bfm.mod: newbfm.o
../femlib/mod/mod_bnd.mod: mod_bnd.o
../femlib/mod/mod_bnd_aux.mod: mod_bnd_aux.o
../femlib/mod/mod_bndo.mod: mod_bndo.o
../femlib/mod/mod_bound_dynamic.mod: mod_bound_dynamic.o
../femlib/mod/mod_bound_geom.mod: mod_bound_geom.o
../femlib/mod/mod_conz.mod: mod_conz.o
../femlib/mod/mod_debug.mod: debug.o
../femlib/mod/mod_depth.mod: mod_depth.o
../femlib/mod/mod_diff_aux.mod: mod_diff_aux.o
../femlib/mod/mod_diff_visc_fric.mod: mod_diff_visc_fric.o
../femlib/mod/mod_fluidmud.mod: mod_fluidmud.o
../femlib/mod/mod_geom.mod: mod_geom.o
../femlib/mod/mod_geom_dynamic.mod: mod_geom_dynamic.o
../femlib/mod/mod_gotm_aux.mod: mod_gotm_aux.o
../femlib/mod/mod_hydro.mod: mod_hydro.o
../femlib/mod/mod_hydro_baro.mod: mod_hydro_baro.o
../femlib/mod/mod_hydro_print.mod: mod_hydro_print.o
../femlib/mod/mod_hydro_vel.mod: mod_hydro_vel.o
../femlib/mod/mod_internal.mod: mod_internal.o
../femlib/mod/mod_lagrange.mod: mod_lagrange.o
../femlib/mod/mod_layer_thickness.mod: mod_layer_thickness.o
../femlib/mod/mod_meteo.mod: mod_meteo.o
../femlib/mod/mod_nohyd.mod: mod_nohyd.o
../femlib/mod/mod_nudge.mod: newnudge.o
../femlib/mod/mod_nudging.mod: mod_nudging.o
../femlib/mod/mod_offline.mod: suboff.o
../femlib/mod/mod_parwaves.mod: subwaves.o
../femlib/mod/mod_random.mod: mod_random.o
../femlib/mod/mod_renewal_time.mod: subwrt.o
../femlib/mod/mod_restart.mod: subrst.o
../femlib/mod/mod_roughness.mod: mod_roughness.o
../femlib/mod/mod_sediment.mod: sedi3d.o
../femlib/mod/mod_sediment_para.mod: sedi3d.o
../femlib/mod/mod_sedtrans05.mod: sedtrans05.o
../femlib/mod/mod_shypre.mod: shypre.o
../femlib/mod/mod_sinking.mod: mod_sinking.o
../femlib/mod/mod_subset.mod: mod_subset.o
../femlib/mod/mod_system.mod: mod_system.o
../femlib/mod/mod_system_global.mod: simsys_spk.o
../femlib/mod/mod_system_interface.mod: mod_system.o
../femlib/mod/mod_tides.mod: mod_tides.o
../femlib/mod/mod_ts.mod: mod_ts.o
../femlib/mod/mod_turbulence.mod: mod_turbulence.o
../femlib/mod/mod_tvd.mod: mod_tvd.o
../femlib/mod/mod_waves.mod: mod_waves.o
../femlib/mod/nc_domain.mod: nc_domain.o
../femlib/mod/nc_interpol.mod: nc_domain.o
../femlib/mod/nc_time.mod: nc_tutil.o
../femlib/mod/ncnames.mod: nc_dim_coords.o
../femlib/mod/netcdf_params.mod: netcdf.o
../femlib/mod/nls.mod: subnls.o
../femlib/mod/para.mod: subpar3.o
../femlib/mod/para_aqua.mod: aquabc_II_apara.o
../femlib/mod/poisson.mod: newpoi.o
../femlib/mod/projection.mod: subproj.o
../femlib/mod/regular.mod: subreg.o
../femlib/mod/scan_string.mod: subscn.o
../femlib/mod/semi_implicit.mod: subnsh.o
../femlib/mod/shyelab_out.mod: shyelab_output.o
../femlib/mod/shyelab_unit.mod: unitutils.o
../femlib/mod/shyfem_strings.mod: substrings.o
../femlib/mod/shyfem_version.mod: subver.o
../femlib/mod/shyfile.mod: subshy.o
../femlib/mod/shympi.mod: mod_shympi_dummy.o
../femlib/mod/shympi.mod: mod_shympi_node.o
../femlib/mod/shympi_aux.mod: mod_shympi_node_internal.o
../femlib/mod/shympi_time.mod: submpi_time.o
../femlib/mod/shyutil.mod: shyutil.o
../femlib/mod/sigma.mod: sigmautil.o
../femlib/mod/simple_sediments.mod: subssed.o
../femlib/mod/stab.mod: newstab.o
../femlib/mod/tidef.mod: subtidef.o
aquabc_II.o: ../femlib/mod/aquabc_II_sed_ini.mod \
		../femlib/mod/para_aqua.mod param.h 
aquabc_II_ased_ini.o: ../femlib/mod/para_aqua.mod aquabc_II.h \
		nbasin.h 
aquabc_II_fem_interface.o: ../femlib/mod/aquabc_II_sed_ini.mod \
		../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/para_aqua.mod aquabc_II.h \
		aquabc_II_aout.h femtime.h mkonst.h param.h 
aquabc_II_pelagic_lib.o: ../femlib/mod/AQUABC_II_GLOBAL.mod
aquabc_II_pelagic_model.o: ../femlib/mod/AQUABC_II_GLOBAL.mod \
		../femlib/mod/CO2SYS_CDIAC.mod \
		../femlib/mod/para_aqua.mod param.h 
aquabc_II_sediment_lib.o: ../femlib/mod/AQUABC_II_GLOBAL.mod
aquabc_II_sediment_model_1.o: ../femlib/mod/AQUABC_II_GLOBAL.mod \
		../femlib/mod/CO2SYS_CDIAC.mod \
		../femlib/mod/para_aqua.mod param.h 
aquabc_II_util.o: param.h
atoxi3d.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod mkonst.h \
		param.h 
basutil.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/mod_depth.mod 
basutil_bathy.o: ../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/grd.mod \
		../femlib/mod/mod_depth.mod 
basutil_box.o: ../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod 
basutil_hsigma.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod 
basutil_interp.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod param.h 
basutil_part.o: ../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/grd.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod 
basutil_resol.o: ../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/mod_depth.mod 
basutil_smooth.o: ../femlib/mod/basin.mod \
		../femlib/mod/basutil.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/grd.mod ../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod 
bdist.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod \
		../femlib/mod/shympi.mod 
bio3d.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod mkonst.h 
bio3d_util.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		donata.h param.h 
elab_check.o: ../femlib/mod/iso8601.mod
elab_gis.o: ../femlib/mod/basin.mod \
		../femlib/mod/shyfem_strings.mod 
elab_nc.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyelab_out.mod 
elab_nc_dummy.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/shyelab_out.mod 
elab_off.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod 
elab_regular.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod
elabutil.o: ../femlib/mod/clo.mod ../femlib/mod/shyfem_strings.mod
extelab1.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyfem_strings.mod 
femcombine.o: ../femlib/mod/clo.mod ../femlib/mod/fem_util.mod
femelab1.o: ../femlib/mod/clo.mod ../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/shyfem_strings.mod 
femintp.o: ../femlib/mod/clo.mod ../femlib/mod/fem_util.mod
filetype.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod
flxelab1.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyfem_strings.mod 
gotm3d.o: gotmturb.i
grdutil.o: ../femlib/mod/basin.mod ../femlib/mod/grd.mod
intp.o: ../femlib/mod/basin.mod ../femlib/mod/mod_hydro.mod \
		param.h 
lagrange_back.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_lagrange.mod param.h 
lagrange_cont.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_lagrange.mod param.h 
lagrange_decay.o: ../femlib/mod/mod_lagrange.mod femtime.h param.h
lagrange_dif.o: ../femlib/mod/mod_lagrange.mod femtime.h param.h
lagrange_flux.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_lagrange.mod \
		../femlib/mod/mod_layer_thickness.mod 
lagrange_init.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_lagrange.mod 
lagrange_inout.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/lgr_sedim_module.mod \
		../femlib/mod/mod_lagrange.mod \
		../femlib/mod/mod_layer_thickness.mod femtime.h 
lagrange_larve.o: ../femlib/mod/mod_lagrange.mod param.h
lagrange_main.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/lgr_sedim_module.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_lagrange.mod femtime.h 
lagrange_oil.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_lagrange.mod femtime.h param.h 
lagrange_sedim.o: ../femlib/mod/levels.mod \
		../femlib/mod/mod_lagrange.mod param.h 
lagrange_track.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_lagrange.mod param.h 
lagrange_util.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		param.h 
lagrange_util_tr.o: ../femlib/mod/basin.mod param.h
lagrange_vertical.o: ../femlib/mod/basin.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_lagrange.mod 
lgrelab.o: ../femlib/mod/basin.mod ../femlib/mod/mod_depth.mod \
		param.h 
line_util.o: ../femlib/mod/basin.mod param.h
loading.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		param.h 
mercury.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod femtime.h 
mod_shympi_dummy.o: ../femlib/mod/basin.mod \
		../femlib/mod/levels.mod 
mod_shympi_node.o: ../femlib/mod/basin.mod \
		../femlib/mod/levels.mod 
mod_shympi_node_internal.o: ../femlib/mod/shympi.mod
mod_subset.o: ../femlib/mod/basin.mod
nc2fem.o: ../femlib/mod/clo.mod ../femlib/mod/iso8601.mod \
		../femlib/mod/ncnames.mod \
		../femlib/mod/shyfem_strings.mod 
nc_domain.o: ../femlib/mod/basin.mod
nc_tutil.o: ../femlib/mod/iso8601.mod
netcdf.o: netcdf.inc
netcdf_util.o: ../femlib/mod/basin.mod
new36.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/shympi.mod close.h femtime.h \
		mkonst.h 
new3di.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_nudging.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_waves.mod \
		../femlib/mod/shympi.mod ../femlib/mod/tidef.mod \
		mkonst.h pkonst.h 
newbcl.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_ts.mod mkonst.h pkonst.h 
newbfm.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod mkonst.h 
newchao.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod femtime.h param.h 
newchk.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_nohyd.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod \
		femtime.h mkonst.h param.h 
newcon.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_aux.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod \
		mkonst.h 
newcon_omp.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_aux.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_subset.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod 
newconz.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/para.mod \
		../femlib/mod/shympi.mod mkonst.h 
newcra.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_bnd_aux.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod femtime.h \
		mkonst.h param.h pkonst.h 
newexpl.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_ts.mod pkonst.h 
newfix.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_bclfix.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod param.h 
newini.o: ../femlib/mod/basin.mod ../femlib/mod/coordinates.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/shympi.mod mkonst.h pkonst.h 
newnudge.o: ../femlib/mod/basin.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_nudging.mod \
		../femlib/mod/shympi.mod 
newpoi.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_system.mod femtime.h mkonst.h \
		pkonst.h 
newrog.o: ../femlib/mod/basin.mod param.h
newstab.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro.mod femtime.h pkonst.h 
newtra.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/shympi.mod 
newtvd.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_tvd.mod 
nos2shy.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/iso8601.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyfile.mod 
offinf.o: ../femlib/mod/clo.mod
optintp.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod 
ous2shy.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/iso8601.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/shyfile.mod \
		../femlib/mod/shyutil.mod 
ousutil.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod
rstinf.o: ../femlib/mod/clo.mod
scal_intp.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		param.h 
scalintp.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod param.h pkonst.h \
		simul.h 
sedi3d.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_sedtrans05.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_waves.mod ../femlib/mod/para.mod \
		femtime.h mkonst.h 
shybas.o: ../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/clo.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/shympi.mod 
shydiff.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyfile.mod \
		../femlib/mod/shyutil.mod param.h 
shyelab.o: ../femlib/mod/clo.mod ../femlib/mod/elabutil.mod \
		../femlib/mod/shyfile.mod 
shyelab1.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/custom_dates.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyfile.mod ../femlib/mod/shympi.mod \
		../femlib/mod/shyutil.mod 
shyelab_area.o: ../femlib/mod/elabutil.mod
shyelab_nodes.o: ../femlib/mod/basin.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shyfem_strings.mod 
shyelab_output.o: ../femlib/mod/basin.mod \
		../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_offline.mod \
		../femlib/mod/shyfem_strings.mod \
		../femlib/mod/shyfile.mod 
shyelab_util.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod 
shyfem.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/coordinates.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_bclfix.mod \
		../femlib/mod/mod_bfm.mod \
		../femlib/mod/mod_bnd.mod \
		../femlib/mod/mod_bnd_aux.mod \
		../femlib/mod/mod_bndo.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_aux.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_gotm_aux.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_nohyd.mod \
		../femlib/mod/mod_nudging.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_sediment.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_subset.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_turbulence.mod \
		../femlib/mod/mod_tvd.mod \
		../femlib/mod/mod_waves.mod \
		../femlib/mod/projection.mod \
		../femlib/mod/shympi.mod ../femlib/mod/tidef.mod \
		femtime.h mkonst.h pkonst.h 
shypart.o: ../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/clo.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/shympi.mod 
shypre.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/grd.mod 
shyproj.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/grd.mod ../femlib/mod/projection.mod 
shyutil.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod 
simsys_lp.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_system.mod 
simsys_pard.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_system.mod 
simsys_spk.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_system.mod \
		../femlib/mod/mod_system_interface.mod \
		../femlib/mod/shympi.mod 
simsys_spkmpi.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_system.mod \
		../femlib/mod/shympi.mod 
subapn.o: ../femlib/mod/basin.mod
subbfm.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_ts.mod bfm_common.h param.h \
		pkonst.h 
subbnd.o: ../femlib/mod/mod_bnd.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/nls.mod ../femlib/mod/para.mod \
		../femlib/mod/shympi.mod 
subbndo.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_bndo.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/shympi.mod 
subbnds.o: ../femlib/mod/intp_fem_file.mod
subboxa.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_ts.mod simul.h 
subcon1.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod femtime.h simul.h 
subcoo.o: ../femlib/mod/mod_system.mod
subcoord.o: ../femlib/mod/basin.mod ../femlib/mod/shympi.mod
subcst.o: ../femlib/mod/basin.mod ../femlib/mod/mod_bnd.mod \
		../femlib/mod/mod_bound_geom.mod mkonst.h \
		modules.h param.h pkonst.h 
subcus.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod \
		femtime.h param.h 
subdep.o: ../femlib/mod/basin.mod ../femlib/mod/mod_depth.mod \
		../femlib/mod/shympi.mod 
subdif.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_aux.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/shympi.mod 
subdry.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/shympi.mod femtime.h mkonst.h 
subele.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/shympi.mod param.h 
subets.o: etsinf.h
subetsa.o: ../femlib/mod/basin.mod ../femlib/mod/coordinates.mod \
		../femlib/mod/ets.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_waves.mod ../femlib/mod/nls.mod \
		femtime.h modules.h simul.h 
subexta.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_sediment.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_waves.mod ../femlib/mod/nls.mod \
		../femlib/mod/shympi.mod modules.h simul.h 
subfind.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom.mod 
subfloodfill.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod
subflx3d.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_vel.mod femtime.h 
subflxa.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_sediment.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/nls.mod \
		modules.h simul.h 
subflxu.o: ../femlib/mod/levels.mod ../femlib/mod/mod_geom.mod
subfvl.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_layer_thickness.mod 
subgotm.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_gotm_aux.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_turbulence.mod \
		../femlib/mod/shympi.mod femtime.h pkonst.h 
subgotm_mud.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_gotm_aux.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_turbulence.mod femtime.h param.h \
		pkonst.h 
subhisto.o: histo.h
sublin.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod
sublnka.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod
sublnkd.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/shympi.mod 
sublnks.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/shympi.mod femtime.h 
sublnku.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod
sublpl.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod param.h
submask.o: ../femlib/mod/basin.mod
submet.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/meteo_forcing_module.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_ts.mod 
submeteo2.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/shyfem_strings.mod 
submpi_assert.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod 
submpi_basin.o: ../femlib/mod/basin.mod ../femlib/mod/shympi.mod
submpi_custom.o: ../femlib/mod/basin.mod ../femlib/mod/shympi.mod
submpi_ghost.o: ../femlib/mod/basin.mod ../femlib/mod/shympi.mod
submud.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_gotm_aux.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_ts.mod param.h pkonst.h 
submud_dummy.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_fluidmud.mod 
subn11.o: ../femlib/mod/basin.mod ../femlib/mod/intp_fem_file.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_bnd_aux.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/shympi.mod mkonst.h pkonst.h 
subn35.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/nls.mod \
		../femlib/mod/shympi.mod pkonst.h 
subnev.o: ../femlib/mod/basin.mod ../femlib/mod/shympi.mod
subnls.o: ../femlib/mod/para.mod
subnos.o: nosinf.h
subnosa.o: ../femlib/mod/basin.mod ../femlib/mod/mod_depth.mod \
		param.h 
subnsa.o: param.h simul.h
subnsh.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/nls.mod ../femlib/mod/shympi.mod \
		femtime.h modules.h simul.h 
subnsu.o: ../femlib/mod/basin.mod
suboff.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod \
		femtime.h param.h 
subous.o: ousinf.h
subousa.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/shyfile.mod ../femlib/mod/shympi.mod \
		simul.h 
suboutput.o: femtime.h
suboutputd.o: femtime.h
subpard.o: ../femlib/mod/basin.mod ../femlib/mod/mod_system.mod
subqfxf.o: subqfx.h
subqfxm1.o: subqfxm.h
subqfxm2.o: subqfxm.h
subqfxm3.o: subqfxm.h
subqfxm4.o: subqfxm.h
subqfxm5.o: subqfxm.h
subqfxt.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_ts.mod subqfxm.h 
subqfxu4.o: subqfxm.h
subreg.o: ../femlib/mod/basin.mod
subres.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_ts.mod femtime.h simul.h 
subrst.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_ts.mod femtime.h 
subshyutil.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/shyfile.mod ../femlib/mod/shympi.mod \
		simul.h 
subshyutil2.o: ../femlib/mod/clo.mod ../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/shyfile.mod \
		../femlib/mod/shyutil.mod 
subspk.o: ../femlib/mod/mod_system.mod \
		../femlib/mod/mod_system_interface.mod 
subssed.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod 
subsys.o: ../femlib/mod/para.mod
subtidef.o: ../femlib/mod/basin.mod ../femlib/mod/coordinates.mod \
		../femlib/mod/iso8601.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod 
subtime.o: ../femlib/mod/shympi.mod femtime.h
subtrace.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_vel.mod femtime.h param.h 
subtsfile.o: ../femlib/mod/iso8601.mod
subtsuvfile.o: ../femlib/mod/intp_fem_file.mod
subuti.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shympi.mod \
		pkonst.h 
subvola.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod \
		femtime.h modules.h param.h volcomp.h 
subwat.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod 
subwaves.o: ../femlib/mod/basin.mod ../femlib/mod/coordinates.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_waves.mod param.h pkonst.h 
subwrt.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/shympi.mod simul.h 
ts2nc.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro.mod param.h simul.h 
tselab1.o: ../femlib/mod/clo.mod ../femlib/mod/elabtime.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/shyfem_strings.mod 
weutro.o: donata.h weutro.h
weutro_sedim.o: donata.h weutro.h
weutro_seed.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod
weutro_shell.o: donata.h
weutro_ulva.o: donata.h
weutro_wdia.o: ../femlib/mod/eutro.mod weutro.h

