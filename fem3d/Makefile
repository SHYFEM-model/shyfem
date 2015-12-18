#
# $Id: Makefile,v 1.81 2010-03-11 15:36:38 georg Exp $
#
#------------------------------------------------------------

DIR	= fem3d

#------------------------------------------------------------ just for default
FEMDIR = ..
DIRLIB = $(FEMDIR)/femlib
F77    = gfortran
FFLAGS =  -O
LINKER = $(F77)
#------------------------------------------------------------

include ../Rules.make

#------------------------------------------------------------

FEMBIN    = $(FEMDIR)/fembin

LIBFEM     = fem
LIBGOTM    = gotm
LIBNETCDF  = netcdff
LIBERSEM   = bfm

LIBFFEM    = $(DIRLIB)/lib$(LIBFEM).a

LIBLS      = -L$(DIRLIB)

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#------------------------------------------------------------

ifeq ($(FLUID_MUD),true)
  MUD_FILES = submud.o subgotm_mud.o
else
  MUD_FILES = submud_dummy.o subgotm.o
endif

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

#ifneq ("$(wildcard /usr/lib/libnetcdff.a","")
#FILE_EXISTS = 1
#else
#FILE_EXISTS = 0
#endif

EXENC = 
NETCDF_OBJ = 
DIRLIB_NETCDF = $(NETCDFDIR)/lib/
DIRLIB_NETCDF = $(NETCDFDIR)/lib/x86_64-linux-gnu
DIRINC_NETCDF = $(NETCDFDIR)/include
ifeq ($(NETCDF),true)
  NETCDF_OBJ       = netcdf.o netcdf_util.o 
  LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)
  LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a
  #INCL_NETCDF      = -I$(DIRINC_NETCDF)
  EXENC            = nos2nc ous2nc ts2nc
endif

#  NETCDF_OBJ       = netcdf.o netcdf_util.o \
#			/usr/lib/x86_64-linux-gnu/libnetcdff.so

PARALLEL_OBJ = subomp_dummy.o
ifeq ($(PARALLEL),true)
  PARALLEL_OBJ = subomp.o
endif

GOTM_OBJ = gotm3d.o
ifeq ($(GOTM),true)
  GOTM_OBJ =
  LIBG_GOTM      = -l$(LIBGOTM)
  LIBF_GOTM      = $(DIRLIB)/lib$(LIBGOTM).a
endif

ECOLOGICAL_OBJ = ecological_dummy.o
ifeq ($(ECOLOGICAL),EUTRO)
  ECOLOGICAL_OBJ = bio3d.o bio3d_util.o \
		weutro.o weutro_sedim.o weutro_light.o
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
#  ECOLOGICAL_OBJ = aquabc_fem_interface.o alukas.o aquabc.o \
#		aquabc_sediment_model_1.o aquabc_util.o

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
	  newnudge.o newpoi.o \
	  $(MODULES) \

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
          subreg.o \
          subbnd.o \
          subexta.o subouta.o subdep.o subflxa.o subflx3d.o subflx.o \
	  subdum.o subcus.o \
	  $(BCTSOBJS) \
	  $(QFLUX) \
	  $(GRDOBJS) \
	  subnos.o subcon1.o subrnd.o sublin.o \
	  submod.o subuti.o subele.o subgeo.o subdts.o \
	  links.o sublnka.o sublnks.o sublnku.o subisd.o \
	  debug.o subhisto.o subnosa.o subfvl.o \
	  subeos.o sigmautil.o \
	  nosutil.o ousutil.o subfemfile.o subvintp.o \
	  genutil.o subfemintp.o subtsfile.o subxi.o \
	  $(MODULES) \
	  subshy.o subshyutil.o subbit10.o

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
	subrst.o subnosa.o subdif.o subreg.o \
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
	$(NETCDF_OBJ) $(PARALLEL_OBJ) \
	$(GRDOBJS) \
	subxi.o subshy.o subshyutil.o subbit10.o

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
	$(GRDOBJS) subclo.o

# objects for output routines

OUTOBJS = $(PARMOD) subnls.o subsss.o subscn.o \
	  subfil.o subdef.o subnsu.o subssv.o \
	  subbas.o subsys.o \
	  subapn.o subnsa.o subdum.o \
	  subdts.o subdep.o \
	  debug.o subreg.o subnev.o \
	  subfemfile.o subtsfile.o \
	  subvintp.o subver.o sigmautil.o suboutput.o suboutputd.o \
	  subclo.o \
	  mod_depth.o mod_hydro.o mod_geom.o mod_hydro_print.o \
	  mod_bound_geom.o mod_bnd.o mod_hydro_baro.o newlevels.o \
	  $(GRDOBJS) subsrt.o subgeo.o

RSTOBJS = subrst.o ecological_dummy.o \
	  mod_ts.o mod_conz.o mod_hydro_vel.o mod_geom_dynamic.o \
	  subbit10.o

BASOBJS = basutil.o basutil_resol.o basutil_hsigma.o basutil_bathy.o \
	  basutil_interp.o basutil_smooth.o basutil_box.o \
	  subnos.o \
	  sublnka.o sublnku.o links.o
BASEXTRA = sublnka.o links.o sublnku.o \
		subinterpol.o subinterpol1.o subnos.o

GENOBJS = newlevels.o
NOSOBJS = subnos.o nosutil.o genutil.o subeos.o elabutil.o noselab1.o
OUSOBJS = subous.o ousutil.o genutil.o elabutil.o ouselab1.o
SHYOBJS = subous.o ousutil.o subnos.o nosutil.o subeos.o genutil.o elabutil.o \
		noselab1.o ouselab1.o subext.o extelab1.o \
		shyelab1.o subshy.o subshyutil.o \
		subflx.o flxelab1.o
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
	  subnosa.o subssm.o

LAGRANGE =	lagrange_main.o \
		lagrange_back.o lagrange_dif.o \
		lagrange_cont.o lagrange_init.o \
		lagrange_inout.o lagrange_util.o \
		lagrange_larve.o lagrange_flux.o \
		lagrange_track.o lagrange_util_tr.o \
		lagrange_decay.o lagrange_sedim.o \
		lagrange_connect.o lagrange_vertical.o

# heat flux module

QFLUX = subqfxt.o subqfxf.o \
	subqfxm1.o subqfxm2.o subqfxm3.o subqfxm4.o subqfxm5.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o

QFLUX_TEST = subqfxf.o \
	subqfxm1.o subqfxm2.o subqfxm3.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o \
	subfil.o

# experimental stuff

#PARMOD = subpar.o subfnm.o		# old parameter modules
#PARMOD = subpar1.o			# new parameter modules
PARMOD = subpar3.o			# newest parameter modules
IFF_FILES = 
IFF_FILES = subfemintp.o subtsfile.o submeteo2.o

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

EXEINFO  = flxinf extinf rstinf \
		feminf tsinf filetype \
		offinf lgrinf etsinf \
		noselab ouselab shyelab nosextr_nodes extelab flxelab \
		femelab
EXESPLIT = splitext splitflx splitets
EXELAP   =  
EXEBAS   = 
EXEVAR   = nosresidence scal_intp
EXEUTIL  = bc2fem

NOSVAR   = 
NOSAVER  = 
NOSEXTR  = 
NOSALL   = $(NOSEXTR) $(NOSAVER) $(NOSVAR)

OUSEXTR  = 

EXEEXTRA = \
		basproj baswork eosinf \
		nos2nc ous2nc ts2nc \
		ouswork ouswork1 ouswrite prep_input \
		basbox
EXEAUX   = baschk bascomp check_debug femelab laplap mkgeom \
		optintp ousaver ousvdiff

EXES     = $(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXELAP) $(EXEBAS) \
		$(EXEVAR) $(EXEUTIL) $(NOSEXTR) $(OUSEXTR) $(NOSAVER) \
		$(EXEEXTRA) $(EXENC) $(NOSVAR) $(EXEAUX)

#------------------------------------------------

default:	fem

all:	fem
fem:	libfem main info split bas var util links nc
main:	$(EXEMAIN)
info:	$(EXEINFO)
split:	$(EXESPLIT)
bas:	$(EXEBAS)
var:	$(EXEVAR)
nc:	$(EXENC)
util:	$(NOSEXTR) $(OUSEXTR) $(NOSAVER) $(EXEUTIL)

links:
	@$(FEMBIN)/symlink $(FEMBIN) \
		$(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXEUTIL)

#--------------------------------------------- EXES

bedload:  bedload.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS)

shyfem:  shyfem.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS)

shypre:  shypre.o $(VPOBJS)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(VPOBJS)

#------------------------------------------------

zlevel:  zlevel.o 
	$(LINKER) $(LFLAGS) $@.o -o $@

lagint:  lagint.o sublag.o $(LAGINTP) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o sublag.o $(LAGINTP) $(OUTOBJS) -o $@

readext:  readext.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subext.o $(OUTOBJS) -o $@

splitext:  splitext.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subext.o $(OUTOBJS) -o $@

splitext_big:  splitext_big.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subext.o $(OUTOBJS) -o $@

splitflx:  splitflx.o subflx.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subflx.o $(OUTOBJS) -o $@

splitets:  splitets.o $(ETSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(ETSOBJS) $(OUTOBJS) -o $@

extinf:  extinf.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subext.o $(OUTOBJS) -o $@

flxinf:  flxinf.o subflx.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subflx.o $(OUTOBJS) -o $@

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

bastreat:  bastreat.o $(BASEXTRA) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(BASEXTRA) $(OUTOBJS) -o $@

scalintp:  scalintp.o $(OUTOBJS) $(BATHYOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BATHYOBJS) -o $@

scal_intp:  scal_intp.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@

bascomp:  bascomp.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

etsinf:  etsinf.o $(ETSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(ETSOBJS) $(OUTOBJS) -o $@

filetype:  filetype.o $(NOSOUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOUSOBJS) $(OUTOBJS) -o $@

feminf:  feminf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

tsinf:  tsinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

femelab:  femelab.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

bc2fem:  bc2fem.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

check_debug: check_debug.o
	$(LINKER) $(LFLAGS) $@.o -o $@

2nc: ous2nc nos2nc ts2nc

ous2nc:  ous2nc.o $(SHYOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) \
			$(NETCDF_OBJ) $(LIBGS)

nos2nc:  nos2nc.o $(SHYOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) \
			$(NETCDF_OBJ) $(LIBGS)

ts2nc:  ts2nc.o $(SHYOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) \
			$(NETCDF_OBJ) $(LIBGS)

eosinf:  eosinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

prep_input:  prep_input.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

offinf: offinf.o $(GENERAL_OBJS)
	$(LINKER) $(LFLAGS) $@.o $(GENERAL_OBJS) -o $@

#-------------------------------------------------------

.PHONY: nosextr ousextr nosaver nos

nos: $(NOSALL)

ouselab: ouselab.o $(SHYOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) -o $@

noselab: noselab.o $(SHYOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) -o $@

extelab: extelab.o $(SHYOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) -o $@

flxelab: flxelab.o $(SHYOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) -o $@

shyelab: shyelab.o $(SHYOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(SHYOBJS) $(OUTOBJS) -o $@

nosextr_nodes:  nosextr_nodes.o $(NOSOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

#-------------------------------------------------------

map_influence:  map_influence.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosresidence:  nosresidence.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

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
	$(MAKEDEPEND) -moddir=$(MODDIR) *.f *.f90

#---------------------------------------------------

options:
	@echo "compiler  $(COMPILER) $(F77)"
	@echo "solver    $(SOLVER)"

#---------------------------------------------------

clean: cleanobj cleanexe cleanvar cleantex cleanlahey cleanexport cleandoc \
	cleanifc cleanlinks

cleanall: clean cleanzip cleanlib cleanextra cleandiff cleanmod

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) \
		$(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXEUTIL)

cleanifc:
	-rm -f ifc*

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanexe:
	-rm -f $(EXES)
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

strip:
	-strip $(EXES)

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

../femlib/mod/para_aqua.mod: aquabc_II_apara.o
aquabc_II_fem_interface.o: ../femlib/mod/levels.mod \
		../femlib/mod/evgeom.mod param.h femtime.h \
		../femlib/mod/aquabc_II_sed_ini.mod \
		../femlib/mod/basin.mod aquabc_II_aout.h \
		../femlib/mod/para_aqua.mod \
		../femlib/mod/mod_diff_visc_fric.mod mkonst.h \
		aquabc_II.h 
aquabc_II_util.o: param.h
atoxi3d.o: param.h ../femlib/mod/levels.mod mkonst.h \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_diff_visc_fric.mod 
basproj.o: ../femlib/mod/projection.mod ../femlib/mod/basin.mod
../femlib/mod/basutil.mod: basutil.o
basutil.o: ../femlib/mod/mod_depth.mod ../femlib/mod/clo.mod \
		../femlib/mod/basin.mod 
basutil_bathy.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod ../femlib/mod/basutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/grd.mod 
basutil_box.o: ../femlib/mod/evgeom.mod ../femlib/mod/mod_geom.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod 
basutil_hsigma.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
basutil_interp.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod param.h 
basutil_resol.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod ../femlib/mod/basutil.mod 
basutil_smooth.o: ../femlib/mod/mod_geom.mod \
		../femlib/mod/basutil.mod ../femlib/mod/grd.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod 
bdist.o: ../femlib/mod/mod_geom.mod ../femlib/mod/basin.mod
bio3d.o: param.h ../femlib/mod/levels.mod ../femlib/mod/basin.mod \
		mkonst.h ../femlib/mod/mod_diff_visc_fric.mod 
bio3d_util.o: donata.h ../femlib/mod/basin.mod param.h \
		../femlib/mod/levels.mod 
debug.o: debug_aux1.h
../femlib/mod/elabutil.mod: elabutil.o
elabutil.o: ../femlib/mod/levels.mod ../femlib/mod/clo.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
eosinf.o: param.h
etsinf.o: ../femlib/mod/ets.mod param.h
extelab1.o: ../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/elabutil.mod ../femlib/mod/basin.mod \
		../femlib/mod/clo.mod ../femlib/mod/mod_depth.mod 
femelab.o: ../femlib/mod/clo.mod
feminf.o: ../femlib/mod/clo.mod
filetype.o: ../femlib/mod/clo.mod
flxelab1.o: ../femlib/mod/mod_depth.mod ../femlib/mod/clo.mod \
		../femlib/mod/basin.mod ../femlib/mod/elabutil.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod 
flxinf.o: param.h
gotm3d.o: gotmturb.i debug_aux1.h
grdutil.o: ../femlib/mod/grd.mod ../femlib/mod/basin.mod
intp.o: ../femlib/mod/basin.mod param.h \
		../femlib/mod/mod_hydro.mod 
lagrange_back.o: param.h ../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_lagrange.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod 
lagrange_cont.o: ../femlib/mod/mod_lagrange.mod param.h \
		../femlib/mod/mod_geom.mod ../femlib/mod/basin.mod 
lagrange_decay.o: param.h ../femlib/mod/mod_lagrange.mod femtime.h
lagrange_dif.o: ../femlib/mod/mod_lagrange.mod femtime.h param.h
lagrange_flux.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_lagrange.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod 
lagrange_init.o: ../femlib/mod/basin.mod param.h \
		../femlib/mod/mod_lagrange.mod 
lagrange_inout.o: ../femlib/mod/levels.mod \
		../femlib/mod/evgeom.mod param.h femtime.h \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/lgr_sedim_module.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_lagrange.mod 
lagrange_larve.o: ../femlib/mod/mod_lagrange.mod param.h
lagrange_main.o: param.h ../femlib/mod/levels.mod femtime.h \
		../femlib/mod/mod_lagrange.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/lgr_sedim_module.mod 
lagrange_oil.o: param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_lagrange.mod femtime.h \
		../femlib/mod/basin.mod 
../femlib/mod/lgr_sedim_module.mod: lagrange_sedim.o
lagrange_sedim.o: ../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_lagrange.mod 
lagrange_track.o: param.h ../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_lagrange.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro_vel.mod 
lagrange_util.o: ../femlib/mod/evgeom.mod param.h \
		../femlib/mod/basin.mod 
lagrange_util_tr.o: param.h ../femlib/mod/basin.mod
lagrange_vertical.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_lagrange.mod 
lgrelab.o: param.h ../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod 
lgrinf.o: param.h
loading.o: param.h ../femlib/mod/levels.mod \
		../femlib/mod/basin.mod 
map_influence.o: param.h ../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
../femlib/mod/mod_area.mod: mod_area.o
../femlib/mod/mod_bclfix.mod: mod_bclfix.o
../femlib/mod/mod_bnd.mod: mod_bnd.o
../femlib/mod/mod_bnd_aux.mod: mod_bnd_aux.o
../femlib/mod/mod_bndo.mod: mod_bndo.o
../femlib/mod/mod_bound_dynamic.mod: mod_bound_dynamic.o
../femlib/mod/mod_bound_geom.mod: mod_bound_geom.o
../femlib/mod/mod_conz.mod: mod_conz.o
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
../femlib/mod/mod_nudging.mod: mod_nudging.o
../femlib/mod/mod_roughness.mod: mod_roughness.o
../femlib/mod/mod_sinking.mod: mod_sinking.o
../femlib/mod/mod_subset.mod: mod_subset.o
mod_subset.o: ../femlib/mod/basin.mod
../femlib/mod/mod_system.mod: mod_system.o
../femlib/mod/mod_tides.mod: mod_tides.o
../femlib/mod/mod_ts.mod: mod_ts.o
../femlib/mod/mod_turbulence.mod: mod_turbulence.o
../femlib/mod/mod_tvd.mod: mod_tvd.o
../femlib/mod/mod_waves.mod: mod_waves.o
netcdf.o: ../femlib/mod/levels.mod param.h ../femlib/mod/basin.mod \
		netcdf.inc netcdf.h ../femlib/mod/mod_depth.mod 
netcdf_util.o: ../femlib/mod/basin.mod param.h
new36.o: femtime.h param.h ../femlib/mod/mod_bound_dynamic.mod \
		close.h ../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_diff_visc_fric.mod mkonst.h 
new3di.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_bound_geom.mod mkonst.h \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_nudging.mod \
		../femlib/mod/mod_waves.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/tidef.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod pkonst.h \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod femtime.h param.h 
newbcl.o: pkonst.h ../femlib/mod/mod_diff_visc_fric.mod mkonst.h \
		../femlib/mod/mod_layer_thickness.mod femtime.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod 
newchao.o: ../femlib/mod/basin.mod param.h \
		../femlib/mod/levels.mod femtime.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod 
newchk.o: ../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro_vel.mod mkonst.h \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		femtime.h ../femlib/mod/mod_hydro.mod param.h 
newcon.o: ../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_hydro_vel.mod mkonst.h \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_diff_aux.mod 
newcon_omp.o: ../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_diff_aux.mod \
		../femlib/mod/mod_subset.mod \
		../femlib/mod/mod_hydro_vel.mod 
newconz.o: ../femlib/mod/mod_diff_visc_fric.mod mkonst.h \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_conz.mod femtime.h \
		../femlib/mod/levels.mod param.h 
newcra.o: pkonst.h ../femlib/mod/mod_bnd_aux.mod \
		../femlib/mod/mod_depth.mod mkonst.h \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/mod_hydro_baro.mod param.h \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/basin.mod 
newexpl.o: ../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_diff_visc_fric.mod pkonst.h \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro.mod 
newfix.o: ../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_bclfix.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/basin.mod 
newini.o: ../femlib/mod/basin.mod ../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_hydro_print.mod femtime.h \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_hydro.mod param.h \
		../femlib/mod/levels.mod \
		../femlib/mod/intp_fem_file.mod mkonst.h \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro_vel.mod pkonst.h \
		../femlib/mod/coordinates.mod 
../femlib/mod/levels.mod: newlevels.o
../femlib/mod/mod_nudge.mod: newnudge.o
newnudge.o: ../femlib/mod/mod_internal.mod ../femlib/mod/basin.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/mod_nudging.mod 
../femlib/mod/poisson.mod: newpoi.o
newpoi.o: ../femlib/mod/mod_internal.mod ../femlib/mod/basin.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_system.mod femtime.h \
		../femlib/mod/mod_layer_thickness.mod mkonst.h \
		../femlib/mod/mod_depth.mod pkonst.h 
newrog.o: ../femlib/mod/basin.mod param.h
newstab.o: param.h ../femlib/mod/levels.mod femtime.h \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_diff_visc_fric.mod stab.h 
newtra.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod 
newtvd.o: ../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/evgeom.mod param.h \
		../femlib/mod/levels.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_tvd.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod 
nos2nc.o: ../femlib/mod/mod_hydro.mod ../femlib/mod/levels.mod \
		param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
noselab1.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		param.h ../femlib/mod/elabutil.mod 
nosextr_nodes.o: param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
nosresidence.o: ../femlib/mod/levels.mod param.h \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
offinf.o: ../femlib/mod/clo.mod
optintp.o: ../femlib/mod/mod_depth.mod ../femlib/mod/clo.mod \
		../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		param.h 
ous2nc.o: ../femlib/mod/basin.mod ../femlib/mod/mod_depth.mod \
		../femlib/mod/evgeom.mod param.h \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod 
ouselab1.o: ../femlib/mod/clo.mod ../femlib/mod/mod_depth.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/elabutil.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/basin.mod 
ousutil.o: ../femlib/mod/evgeom.mod
readext.o: param.h ../femlib/mod/basin.mod
rstinf.o: ../femlib/mod/clo.mod
scal_intp.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod \
		param.h 
scalintp.o: pkonst.h param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_depth.mod simul.h \
		../femlib/mod/basin.mod 
sedi3d.o: ../femlib/mod/mod_waves.mod sed_param.h \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_geom.mod param.h \
		../femlib/mod/basin.mod ../femlib/mod/mod_ts.mod 
shybas.o: ../femlib/mod/basutil.mod ../femlib/mod/mod_geom.mod \
		../femlib/mod/evgeom.mod param.h \
		../femlib/mod/basin.mod ../femlib/mod/clo.mod \
		../femlib/mod/mod_depth.mod 
shyelab.o: ../femlib/mod/elabutil.mod ../femlib/mod/shyfile.mod \
		../femlib/mod/clo.mod 
shyelab1.o: ../femlib/mod/mod_depth.mod ../femlib/mod/clo.mod \
		../femlib/mod/basin.mod param.h \
		../femlib/mod/evgeom.mod ../femlib/mod/shyfile.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/elabutil.mod 
shyfem.o: ../femlib/mod/mod_diff_aux.mod pkonst.h \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_turbulence.mod \
		../femlib/mod/mod_nohyd.mod \
		../femlib/mod/mod_gotm_aux.mod \
		../femlib/mod/projection.mod \
		../femlib/mod/mod_bclfix.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/mod_bndo.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_bnd.mod femtime.h param.h \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/intp_fem_file.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_subset.mod mkonst.h \
		../femlib/mod/coordinates.mod \
		../femlib/mod/mod_bnd_aux.mod \
		../femlib/mod/mod_waves.mod \
		../femlib/mod/mod_nudging.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_tvd.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/tidef.mod 
shypre.o: ../femlib/mod/basin.mod ../femlib/mod/clo.mod param.h
sigmautil.o: sigma.h
simsys_lp.o: ../femlib/mod/basin.mod ../femlib/mod/mod_system.mod \
		../femlib/mod/levels.mod param.h 
simsys_pard.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_system.mod 
simsys_spk.o: ../femlib/mod/basin.mod ../femlib/mod/mod_system.mod \
		../femlib/mod/levels.mod 
splitets.o: param.h
splitflx.o: param.h
subapn.o: ../femlib/mod/basin.mod
../femlib/mod/basin.mod: subbas.o
subbas.o: param.h
subbfm.o: ../femlib/mod/basin.mod ../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_hydro_print.mod bfm_common.h \
		../femlib/mod/mod_hydro.mod param.h \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_diff_visc_fric.mod pkonst.h 
subbnd.o: param.h bound_names.h ../femlib/mod/mod_bnd.mod \
		../femlib/mod/mod_bound_geom.mod 
subbndo.o: ../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/basin.mod ../femlib/mod/mod_bndo.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/levels.mod 
subbnds.o: ../femlib/mod/intp_fem_file.mod
subboxa.o: ../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom_dynamic.mod subboxa.h \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_hydro.mod femtime.h param.h \
		modules.h ../femlib/mod/mod_ts.mod \
		../femlib/mod/basin.mod 
../femlib/mod/clo.mod: subclo.o
subcon1.o: ../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_depth.mod simul.h mkonst.h \
		../femlib/mod/basin.mod femtime.h \
		../femlib/mod/levels.mod 
subcoo.o: ../femlib/mod/basin.mod ../femlib/mod/mod_system.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod param.h 
../femlib/mod/coordinates.mod: subcoord.o
subcoord.o: ../femlib/mod/basin.mod
subcst.o: pkonst.h ../femlib/mod/mod_bnd.mod param.h \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/basin.mod mkonst.h modules.h 
subcus.o: femtime.h ../femlib/mod/mod_hydro.mod param.h \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_area.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_internal.mod 
subdep.o: ../femlib/mod/mod_depth.mod ../femlib/mod/basin.mod
subdif.o: ../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_diff_aux.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/evgeom.mod param.h \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro_print.mod 
subdry.o: ../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/basin.mod femtime.h \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro.mod param.h \
		../femlib/mod/evgeom.mod mkonst.h 
subdts.o: subdts.h
subele.o: ../femlib/mod/mod_depth.mod param.h \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/basin.mod ../femlib/mod/mod_area.mod 
subets.o: etsinf.h
subetsa.o: ../femlib/mod/mod_waves.mod \
		../femlib/mod/coordinates.mod \
		../femlib/mod/nls.mod ../femlib/mod/ets.mod \
		simul.h ../femlib/mod/mod_depth.mod \
		../femlib/mod/levels.mod param.h femtime.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_print.mod modules.h \
		../femlib/mod/basin.mod ../femlib/mod/mod_ts.mod 
../femlib/mod/ets.mod: subetsm.o
../femlib/mod/extra.mod: subexta.o
subexta.o: simul.h ../femlib/mod/mod_hydro_print.mod modules.h \
		../femlib/mod/nls.mod femtime.h 
../femlib/mod/intp_fem_file.mod: subfemintp.o
subflx3d.o: femtime.h ../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_hydro_vel.mod 
../femlib/mod/flux.mod: subflxa.o
subflxa.o: ../femlib/mod/levels.mod param.h ../femlib/mod/nls.mod \
		../femlib/mod/mod_conz.mod femtime.h modules.h \
		../femlib/mod/mod_ts.mod 
subflxu.o: ../femlib/mod/levels.mod ../femlib/mod/mod_geom.mod
subfvl.o: ../femlib/mod/levels.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/basin.mod ../femlib/mod/mod_area.mod 
subgotm.o: pkonst.h ../femlib/mod/mod_turbulence.mod femtime.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_gotm_aux.mod debug_aux2.h \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_hydro_print.mod 
subgotm_mud.o: pkonst.h param.h ../femlib/mod/mod_turbulence.mod \
		femtime.h ../femlib/mod/mod_hydro.mod \
		../femlib/mod/basin.mod ../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_gotm_aux.mod \
		../femlib/mod/mod_diff_visc_fric.mod debug_aux2.h \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_hydro_print.mod 
../femlib/mod/grd.mod: subgrd.o
subhisto.o: histo.h
sublin.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod \
		param.h 
sublnka.o: ../femlib/mod/mod_geom.mod ../femlib/mod/basin.mod
sublnkd.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/evgeom.mod param.h \
		../femlib/mod/mod_geom.mod 
sublnks.o: femtime.h ../femlib/mod/mod_geom.mod \
		../femlib/mod/evgeom.mod param.h \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/basin.mod 
sublnku.o: param.h ../femlib/mod/mod_geom.mod \
		../femlib/mod/basin.mod 
sublpl.o: ../femlib/mod/basin.mod ../femlib/mod/evgeom.mod param.h
submet.o: param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod femtime.h \
		../femlib/mod/basin.mod ../femlib/mod/mod_ts.mod \
		../femlib/mod/meteo_forcing_module.mod \
		../femlib/mod/mod_bound_dynamic.mod 
../femlib/mod/meteo_forcing_module.mod: submeteo2.o
submeteo2.o: femtime.h param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/intp_fem_file.mod 
submud.o: ../femlib/mod/mod_gotm_aux.mod pkonst.h \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_sinking.mod \
		../femlib/mod/mod_bound_geom.mod \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_bound_dynamic.mod 
submud_dummy.o: ../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/basin.mod ../femlib/mod/levels.mod 
subn11.o: mkonst.h ../femlib/mod/intp_fem_file.mod \
		../femlib/mod/mod_bound_geom.mod pkonst.h \
		../femlib/mod/mod_bnd_aux.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_bound_dynamic.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod femtime.h 
../femlib/mod/chezy.mod: subn35.o
subn35.o: ../femlib/mod/nls.mod pkonst.h \
		../femlib/mod/mod_fluidmud.mod \
		../femlib/mod/mod_diff_visc_fric.mod \
		../femlib/mod/mod_roughness.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/basin.mod 
../femlib/mod/evgeom.mod: subnev.o
subnev.o: ../femlib/mod/basin.mod param.h
../femlib/mod/nls.mod: subnls.o
subnos.o: nosinf.h
subnosa.o: param.h ../femlib/mod/mod_depth.mod \
		../femlib/mod/basin.mod 
subnsa.o: simul.h param.h
subnsh.o: modules.h ../femlib/mod/basin.mod femtime.h \
		../femlib/mod/levels.mod semi.h param.h simul.h \
		../femlib/mod/nls.mod 
subnsu.o: param.h ../femlib/mod/basin.mod
../femlib/mod/mod_offline.mod: suboff.o
suboff.o: param.h ../femlib/mod/levels.mod femtime.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro_print.mod \
		../femlib/mod/mod_hydro_vel.mod 
subous.o: ousinf.h
subousa.o: ../femlib/mod/mod_depth.mod simul.h \
		../femlib/mod/basin.mod femtime.h \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/levels.mod 
subouta.o: ../femlib/mod/mod_hydro_baro.mod femtime.h \
		../femlib/mod/mod_hydro.mod param.h \
		../femlib/mod/basin.mod simul.h \
		../femlib/mod/mod_hydro_print.mod 
suboutput.o: femtime.h
suboutputd.o: femtime.h
../femlib/mod/para.mod: subpar3.o
subpard.o: ../femlib/mod/basin.mod ../femlib/mod/mod_system.mod
../femlib/mod/projection.mod: subproj.o
subqfxf.o: subqfx.h
subqfxm1.o: subqfxm.h
subqfxm2.o: subqfxm.h
subqfxm3.o: subqfxm.h
subqfxm4.o: subqfxm.h
subqfxm5.o: subqfxm.h
subqfxt.o: ../femlib/mod/mod_meteo.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_ts.mod subqfxm.h 
subqfxu4.o: subqfxm.h
subreg.o: reg.h ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_geom.mod param.h \
		../femlib/mod/evgeom.mod 
subres.o: simul.h ../femlib/mod/mod_depth.mod param.h \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod 
../femlib/mod/mod_restart.mod: subrst.o
subrst.o: ../femlib/mod/basin.mod ../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_conz.mod femtime.h \
		../femlib/mod/mod_hydro_vel.mod 
../femlib/mod/shyfile.mod: subshy.o
subshyutil.o: simul.h ../femlib/mod/basin.mod \
		../femlib/mod/levels.mod ../femlib/mod/shyfile.mod 
subspk.o: ../femlib/mod/mod_system.mod ../femlib/mod/basin.mod
../femlib/mod/tidef.mod: subtidef.o
subtidef.o: ../femlib/mod/basin.mod ../femlib/mod/mod_depth.mod \
		../femlib/mod/coordinates.mod \
		../femlib/mod/mod_hydro.mod 
subtime.o: femtime.h
subtrace.o: ../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_hydro.mod femtime.h \
		../femlib/mod/levels.mod param.h \
		../femlib/mod/mod_hydro_vel.mod 
subtsuvfile.o: ../femlib/mod/intp_fem_file.mod param.h
subuti.o: ../femlib/mod/mod_depth.mod pkonst.h \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/basin.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/evgeom.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_geom.mod 
subvola.o: ../femlib/mod/basin.mod volcomp.h modules.h femtime.h \
		param.h ../femlib/mod/mod_geom.mod 
subwat.o: ../femlib/mod/basin.mod ../femlib/mod/mod_geom.mod \
		../femlib/mod/levels.mod ../femlib/mod/evgeom.mod \
		param.h ../femlib/mod/mod_hydro.mod 
subwaves.o: ../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_depth.mod \
		../femlib/mod/coordinates.mod \
		../femlib/mod/mod_waves.mod \
		../femlib/mod/mod_internal.mod \
		../femlib/mod/evgeom.mod \
		../femlib/mod/mod_meteo.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_layer_thickness.mod \
		../femlib/mod/mod_hydro_baro.mod \
		../femlib/mod/mod_roughness.mod pkonst.h \
		../femlib/mod/basin.mod param.h \
		../femlib/mod/mod_geom.mod \
		../femlib/mod/mod_hydro.mod femtime.h 
../femlib/mod/mod_renewal_time.mod: subwrt.o
subwrt.o: ../femlib/mod/evgeom.mod param.h \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod femtime.h \
		../femlib/mod/basin.mod simul.h \
		../femlib/mod/mod_depth.mod 
ts2nc.o: ../femlib/mod/levels.mod param.h ../femlib/mod/evgeom.mod \
		../femlib/mod/mod_hydro.mod simul.h \
		../femlib/mod/basin.mod \
		../femlib/mod/mod_depth.mod 
tsinf.o: ../femlib/mod/clo.mod
weutro.o: weutro.h donata.h
weutro_sedim.o: weutro.h donata.h
zinit.o: param.h ../femlib/mod/basin.mod simul.h
aquabc_II.o: param.h ../femlib/mod/para_aqua.mod \
		../femlib/mod/aquabc_II_sed_ini.mod 
../femlib/mod/aquabc_II_sed_ini.mod: aquabc_II_ased_ini.o
aquabc_II_ased_ini.o: aquabc_II.h nbasin.h \
		../femlib/mod/para_aqua.mod 
../femlib/mod/AQUABC_II_GLOBAL.mod: aquabc_II_co2sys.o
../femlib/mod/VECTOR_MATRIX_UTILS.mod: aquabc_II_co2sys.o
../femlib/mod/CO2SYS_CDIAC.mod: aquabc_II_co2sys.o
aquabc_II_pelagic_lib.o: ../femlib/mod/AQUABC_II_GLOBAL.mod
aquabc_II_pelagic_model.o: ../femlib/mod/AQUABC_II_GLOBAL.mod \
		../femlib/mod/para_aqua.mod \
		../femlib/mod/CO2SYS_CDIAC.mod param.h 
aquabc_II_sediment_lib.o: ../femlib/mod/AQUABC_II_GLOBAL.mod
aquabc_II_sediment_model_1.o: ../femlib/mod/AQUABC_II_GLOBAL.mod \
		param.h ../femlib/mod/CO2SYS_CDIAC.mod \
		../femlib/mod/para_aqua.mod 

