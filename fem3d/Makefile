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
SOLVER_FILES = simsys_lp.o nonhydro_lp.o
ifeq ($(SOLVER),PARDISO)
  LIBG_MKL = `./bin/find_pardlib.sh`
  SOLVER_FILES = simsys_pard.o subpard.o subcoo.o nonhydro_pard.o
endif
ifeq ($(SOLVER),SPARSKIT)
  SOLVER_FILES = simsys_spk.o subspk.o subcoo.o nonhydro_spk.o \
		spk_ilut.o spk_iters.o spk_itaux.o \
		spk_itpack.o spk_matvec.o spk_blas1_mod.o
endif

EXENC = 
NETCDF_OBJ = 
DIRLIB_NETCDF = $(NETCDFDIR)/lib/
DIRINC_NETCDF = $(NETCDFDIR)/include
ifeq ($(NETCDF),true)
  NETCDF_OBJ       = netcdf.o netcdf_util.o
  LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)
  LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a
  INCL_NETCDF      = -I$(DIRINC_NETCDF)
  EXENC            = nos2nc ous2nc ts2nc
endif

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

MODULES = mod_tides.o newlevels.o \
		mod_hydro.o mod_hydro_vel.o \
		mod_hydro_print.o mod_hydro_baro.o \
		mod_diff_visc_fric.o mod_roughness.o \
		mod_ts.o mod_area.o mod_bound_dynamic.o \
		mod_aux_array.o mod_diff_aux.o mod_gotm_aux.o \
		mod_bnd_aux.o mod_depth.o mod_geom_dynamic.o \
		mod_internal.o mod_nohyd.o mod_nudging.o \
		mod_bclfix.o mod_turbulence.o \
		mod_fluidmud.o mod_sinking.o mod_waves.o \
		mod_meteo.o mod_geom.o mod_conz.o \
		mod_bound_geom.o mod_bnd.o \
		mod_tvd.o mod_bndo.o

NEWOBJS = \
	  new3di.o \
	  newrog.o newcra.o \
	  newbcl.o newbcl0.o \
	  newini.o newtra.o \
	  newchao.o \
	  newcon.o newconz.o \
	  newchk.o newdum.o newdum1.o newfix.o \
	  newtvd.o newstab.o newsig.o \
	  newnudge.o \
	  $(MODULES)

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
          subbas.o subnev.o subdim.o \
          subnls.o subwat.o \
          subfil.o subn35.o subout.o subous.o subssg.o \
          subssm.o \
          subnsu.o subsss.o subscn.o subres.o \
          subapn.o subnsa.o \
          subgrd.o subreg.o \
          subbnd.o \
          subexta.o subouta.o subdep.o subflxa.o subflx3d.o subflx.o \
	  subdum.o subcus.o subdwq.o \
	  $(BCTSOBJS) \
	  $(QFLUX) \
	  subnos.o subcon1.o subrnd.o sublin.o \
	  submod.o subuti.o subele.o subgeo.o subdts.o \
	  links.o sublnka.o sublnks.o sublnku.o subisd.o \
	  debug.o subhisto.o subnosa.o subfvl.o \
	  subeos.o sigmautil.o \
	  nosutil.o ousutil.o subfemfile.o subvintp.o \
	  genutil.o subfemintp.o subtsfile.o subxi.o \
	  $(MODULES)

	  #mod_bnd.o

# general objects for 3D model ht

HTOBJS = \
	subver.o \
	subcst.o subsys.o subnsh.o subtime.o suboutput.o \
	subbas.o subdim.o subnev.o \
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
	tideg.o tide1.o subres.o subuti.o \
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
	submet.o subrgf.o subcoord.o subgrd.o subsrt.o subfvl.o \
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
	subxi.o

# objects for pre-processor vp

VPOBJS = \
	subver.o \
        subros.o subcmv.o \
        subdef.o subbas.o subdim.o \
        subfil.o \
        $(PARMOD) \
        subsrt.o subgrd.o \
        subsss.o subscn.o subnsu.o \
	subapn.o subnsa.o subnls.o subsys.o subdum.o \
	grdutil.o \
	newdum.o

# objects for output routines

OUTOBJS = $(PARMOD) subnls.o subsss.o subscn.o \
	  subfil.o subdef.o subnsu.o subssv.o \
	  subbas.o subsys.o subdim.o \
	  subapn.o subnsa.o subdum.o \
	  subdts.o subdep.o \
	  debug.o subreg.o subnev.o \
	  subfemfile.o subtsfile.o \
	  subvintp.o subver.o sigmautil.o suboutput.o \
	  subclo.o newdum.o \
	  mod_depth.o mod_hydro.o mod_geom.o mod_hydro_print.o

RSTOBJS = subrst.o ecological_dummy.o newlevels.o \
	  mod_ts.o mod_conz.o mod_hydro_vel.o mod_geom_dynamic.o

BASEXTRA = subgeo.o subgrd.o subsrt.o sublnka.o links.o sublnku.o \
		subinterpol.o subinterpol1.o subnos.o

GENOBJS = newlevels.o
NOSOBJS = subnos.o nosutil.o genutil.o subeos.o newlevels.o
OUSOBJS = subous.o ousutil.o genutil.o newlevels.o
ETSOBJS = subets.o etsutil.o genutil.o subetsm.o
NOSOUSOBJS = subnos.o nosutil.o subous.o ousutil.o subets.o etsutil.o \
		genutil.o
GRDOBJS = subgrd.o grdutil.o subcoord.o
UTILOBJS = subsrt.o

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

EXEMAIN  = ht vp

EXEINFO  = ousinf flxinf volinf basinf nosinf wininf extinf rstinf \
		feminf tsinf filetype \
		offinf lgrinf etsinf
EXESPLIT = nossplit splitext splitflx splitets
EXELAP   =  
EXEBAS   = bastreat basbathy scalintp
EXEVAR   = nosinit nosresidence scal_intp bashsigma
EXEUTIL  = bc2fem

NOSVAR   = nos_elab nossum nosdiff noswrite nosadmin noselab
NOSAVER  = nosaver_basin nosaver_vert nosaver_time
NOSEXTR  = nosextr_nodes nosextr_gis nosextr_records nosextr_custom \
		nosextr_cat
NOSALL   = $(NOSEXTR) $(NOSAVER) $(NOSVAR)

OUSEXTR  = ousextr_nodes ousextr_records ousextr_gis \
		ousextr_cat

EXEEXTRA = wincrea windcrea \
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

all:	param fem
fem:	param libfem main info split bas var util links nc
main:	param $(EXEMAIN)
info:	param $(EXEINFO)
split:	param $(EXESPLIT)
bas:	param $(EXEBAS)
var:	param $(EXEVAR)
nc:	param $(EXENC)
util:	param $(NOSEXTR) $(OUSEXTR) $(NOSAVER) $(EXEUTIL)

links:
	@$(FEMBIN)/symlink $(FEMBIN) \
		$(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXEUTIL)

#--------------------------------------------- EXES

bedload:  bedload.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS)

ht:  ht.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS)

vp:  vp.o $(VPOBJS)
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

volinf:  volinf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

lgrinf:  lgrinf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

zinit:  zinit.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

wincrea:   wincrea.o
	$(LINKER) $(LFLAGS) $@.o -o $@

windcrea:  windcrea.o $(OUTOBJS) subgeo.o
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) subgeo.o -o $@

wininf:   wininf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

winshift:   winshift.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

winmerge:   winmerge.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

mkgeom:   mkgeom.o
	$(LINKER) $(LFLAGS) $@.o -o $@

cubic:   cubic.o $(OUTOBJS) $(BCTSOBJS) subssm.o
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BCTSOBJS) subssm.o -o $@

#------------------------------------------------

basinf:  basinf.o $(OUTOBJS) $(UTILOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(UTILOBJS) -o $@

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

basproj:  basproj.o $(OUTOBJS) $(GRDOBJS) $(UTILOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(GRDOBJS) $(UTILOBJS) -o $@

bastreat:  bastreat.o $(BASEXTRA) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(BASEXTRA) $(OUTOBJS) -o $@

scalintp:  scalintp.o $(OUTOBJS) $(BATHYOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BATHYOBJS) -o $@

scal_intp:  scal_intp.o subgeo.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) $(LAPLACE) -o $@

bascomp:  bascomp.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

ousinf:  ousinf.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ouswork:  ouswork.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ouswork1:  ouswork1.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousintp:  ousintp.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousaver:  ousaver.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ouswrite:  ouswrite.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousvdiff:  ousvdiff.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

etsinf:  etsinf.o $(ETSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(ETSOBJS) $(OUTOBJS) -o $@

nosinf:  nosinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosinf2:  nosinf2.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

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

nossum:  nossum.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

noselab_time:  noselab_time.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

check_debug: check_debug.o
	$(LINKER) $(LFLAGS) $@.o -o $@

2nc: ous2nc nos2nc ts2nc

ous2nc:  ous2nc.o $(OUSOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) \
			$(NETCDF_OBJ) $(LIBGS)

nos2nc:  nos2nc.o $(NOSOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) \
			$(NETCDF_OBJ) $(LIBGS)

ts2nc:  ts2nc.o $(NOSOBJS) $(OUTOBJS) $(NETCDF_OBJ)
	$(LINKER) -o $@ $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) \
			$(NETCDF_OBJ) $(LIBGS)

eosinf:  eosinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

prep_input:  prep_input.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

offinf: offinf.o
	$(LINKER) $(LFLAGS) $@.o -o $@

#-------------------------------------------------------

.PHONY: nosextr ousextr nosaver nos

ousextr: $(OUSEXTR)

ousextr_nodes:  ousextr_nodes.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousextr_records:  ousextr_records.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousextr_gis:  ousextr_gis.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousextr_cat:  ousextr_cat.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

nos: $(NOSALL)

nosextr: $(NOSEXTR)

nosextr_nodes:  nosextr_nodes.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_gis:  nosextr_gis.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_records:  nosextr_records.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_custom:  nosextr_custom.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_cat:  nosextr_cat.o $(NOSOUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOUSOBJS) $(OUTOBJS) -o $@

nosextr_fluxes:  nosextr_fluxes.o $(NOSOUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOUSOBJS) $(OUTOBJS) -o $@

noselab: noselab.o $(NOSOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosaver: $(NOSAVER)

nosaver_basin:  nosaver_basin.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosaver_vert:  nosaver_vert.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosaver_time:  nosaver_time.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

noswrite:  noswrite.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosadmin:  nosadmin.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

#-------------------------------------------------------

map_influence:  map_influence.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nossplit:  nossplit.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosresidence:  nosresidence.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosmaniq:  nosmaniq.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nostransf:  nostransf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosdiff:  nosdiff.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosinit:  nosinit.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nos_elab:  nos_elab.o $(NOSOBJS) $(OUTOBJS) 
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
	-@$(FEMBIN)/make_param.pl $(FEMDIR) param.h
	-@$(FEMBIN)/param.sh $(FEMDIR) param.h

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

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
aquabc_II_fem_interface.o: aquabc_II.h param.h \
		../femlib/mod/aquabc_II_sed_ini.mod ev.h femtime.h \
		nlevel.h mkonst.h nbasin.h \
		../femlib/mod/para_aqua.mod aquabc_II_aout.h \
		levels.h basin.h diff_visc_fric.h 
aquabc_II_util.o: param.h
atoxi3d.o: nbasin.h param.h levels.h nlevel.h mkonst.h \
		diff_visc_fric.h 
basbathy.o: depth.h pkonst.h param.h evmain.h basin.h simul.h
basbox.o: geom_aux.h depth.h param.h links.h evmain.h basin.h \
		geom.h 
baschk.o: pkonst.h param.h evmain.h basin.h simul.h
bascomp.o: param.h evmain.h basin.h
bashsigma.o: depth.h param.h evmain.h basin.h
basinf.o: depth.h param.h evmain.h basin.h
basproj.o: param.h basin.h
basres.o: param.h basin.h
bastreat.o: geom_aux.h depth.h param.h evmain.h nbasin.h basin.h \
		links.h geom.h 
baswork.o: depth.h param.h evmain.h basin.h geom.h
bdist.o: param.h links.h basin.h
bio3d.o: param.h nlevel.h mkonst.h nbasin.h levels.h basin.h \
		diff_visc_fric.h aux_array.h 
bio3d_util.o: nbasin.h param.h levels.h donata.h
debug.o: debug_aux1.h
eosinf.o: param.h
etsinf.o: param.h ../femlib/mod/ets.mod
femelab.o: ../femlib/mod/clo.mod
feminf.o: ../femlib/mod/clo.mod
filetype.o: ../femlib/mod/clo.mod
flux.o: regflux.h hydro.h param.h ev.h femtime.h nbasin.h \
		hydro_baro.h basin.h links.h 
flxinf.o: param.h
gotm3d.o: gotmturb.i debug_aux1.h
grdutil.o: grd_extra.h param.h basin.h
ht.o: internal.h ts.h evmain.h nohyd.h nudging.h tides.h nbasin.h \
		volcomp.h turbulence.h diff_visc_fric.h gotm_aux.h \
		bnd_aux.h param.h bclfix.h diff_aux.h pkonst.h \
		fluidmud.h geom_dynamic.h levels.h area.h extra.h \
		depth.h bound_dynamic.h mkonst.h nlevel.h \
		bound_geom.h tvd.h simul.h \
		../femlib/mod/mod_bndo.mod meteo.h roughness.h \
		aux_array.h waves.h geom.h hydro_vel.h conz.h \
		hydro.h bound_names.h subbndo.h femtime.h \
		sinking.h nbound.h hydro_baro.h meteo_aux.h \
		basin.h ../femlib/mod/intp_fem_file.mod \
		hydro_print.h 
intp.o: nbasin.h hydro.h param.h
lagrange_back.o: depth.h nbasin.h hydro.h param.h basin.h \
		hydro_print.h lagrange.h 
lagrange_cont.o: param.h links.h basin.h lagrange.h
lagrange_decay.o: param.h femtime.h lagrange.h
lagrange_dif.o: param.h femtime.h lagrange.h
lagrange_flux.o: depth.h hydro.h param.h nlevel.h ev.h femtime.h \
		lagrange.h nbasin.h levels.h basin.h links.h 
lagrange_init.o: nbasin.h param.h basin.h lagrange.h aux_array.h
lagrange_inout.o: depth.h param.h ev.h \
		../femlib/mod/lgr_sedim_module.mod nlevel.h \
		femtime.h lagrange.h levels.h basin.h 
lagrange_larve.o: param.h lagrange.h
lagrange_main.o: nbasin.h param.h basin.h femtime.h \
		../femlib/mod/lgr_sedim_module.mod lagrange.h 
lagrange_oil.o: nbasin.h param.h femtime.h ev.h lagrange.h \
		aux_array.h 
../femlib/mod/lgr_sedim_module.mod: lagrange_sedim.o
lagrange_sedim.o: param.h levels.h lagrange.h
lagrange_track.o: hydro_vel.h param.h levels.h basin.h nlevel.h \
		lagrange.h geom.h 
lagrange_util.o: nbasin.h param.h ev.h basin.h
lagrange_util_tr.o: param.h basin.h
lagrange_vertical.o: depth.h hydro_vel.h hydro.h param.h nlevel.h \
		lagrange.h levels.h basin.h geom.h 
lgrelab.o: nbasin.h depth.h param.h basin.h
lgrinf.o: param.h
loading.o: nbasin.h param.h levels.h
map_influence.o: depth.h param.h basin.h
matrix49.o: regflux.h
../femlib/mod/mod_area.mod: mod_area.o
../femlib/mod/mod_aux_array.mod: mod_aux_array.o
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
../femlib/mod/mod_meteo.mod: mod_meteo.o
../femlib/mod/mod_nohyd.mod: mod_nohyd.o
../femlib/mod/mod_nudging.mod: mod_nudging.o
../femlib/mod/mod_roughness.mod: mod_roughness.o
../femlib/mod/mod_sinking.mod: mod_sinking.o
../femlib/mod/mod_system.mod: mod_system.o
../femlib/mod/mod_tides.mod: mod_tides.o
../femlib/mod/mod_ts.mod: mod_ts.o
../femlib/mod/mod_turbulence.mod: mod_turbulence.o
../femlib/mod/mod_tvd.mod: mod_tvd.o
../femlib/mod/mod_waves.mod: mod_waves.o
netcdf.o: depth.h param.h netcdf.h levels.h netcdf.inc basin.h
netcdf_util.o: param.h basin.h
new36.o: param.h bound_dynamic.h mkonst.h femtime.h close.h \
		nbasin.h basin.h hydro_print.h diff_visc_fric.h 
new3di.o: testbndo.h depth.h internal.h bound_dynamic.h mkonst.h \
		ev.h nlevel.h tides.h nudging.h nbasin.h \
		roughness.h diff_visc_fric.h waves.h aux_array.h \
		hydro_vel.h param.h hydro.h femtime.h pkonst.h \
		fluidmud.h meteo_aux.h geom_dynamic.h levels.h \
		hydro_baro.h area.h basin.h hydro_print.h 
newbcl.o: depth.h hydro.h param.h bound_names.h ts.h nlevel.h \
		mkonst.h femtime.h nbasin.h pkonst.h levels.h \
		basin.h diff_visc_fric.h hydro_print.h aux_array.h 
newchao.o: hydro.h param.h hydro_baro.h basin.h femtime.h nlevel.h
newchk.o: testbndo.h depth.h internal.h bound_dynamic.h ev.h \
		mkonst.h ts.h nlevel.h nbasin.h diff_visc_fric.h \
		hydro_vel.h param.h hydro.h femtime.h levels.h \
		geom_dynamic.h meteo_aux.h hydro_baro.h area.h \
		basin.h hydro_print.h 
newcon.o: testbndo.h depth.h bound_dynamic.h ts.h ev.h mkonst.h \
		nlevel.h nbasin.h aux_array.h geom.h hydro_vel.h \
		param.h hydro.h femtime.h diff_aux.h levels.h \
		area.h basin.h hydro_print.h 
newconz.o: depth.h conz.h param.h bound_names.h ts.h nlevel.h \
		mkonst.h femtime.h nbasin.h levels.h \
		diff_visc_fric.h aux_array.h 
newcra.o: depth.h hydro.h param.h bound_dynamic.h mkonst.h \
		femtime.h ev.h pkonst.h hydro_baro.h basin.h \
		bnd_aux.h 
newdum.o: nbasin.h param.h nlevel.h bound_geom.h nbound.h
newdum1.o: param.h ../femlib/mod/levels.mod levels.h nlevel.h
newexpl.o: depth.h internal.h ts.h ev.h nlevel.h bound_geom.h \
		nbasin.h diff_visc_fric.h aux_array.h geom.h \
		hydro_vel.h param.h hydro.h pkonst.h levels.h \
		geom_dynamic.h basin.h hydro_print.h 
newfix.o: depth.h hydro_vel.h hydro.h internal.h param.h \
		bound_names.h nlevel.h femtime.h bclfix.h nbasin.h \
		levels.h geom_dynamic.h basin.h aux_array.h 
newini.o: hydro.h internal.h param.h mkonst.h nlevel.h nbasin.h \
		pkonst.h levels.h roughness.h basin.h \
		diff_visc_fric.h 
../femlib/mod/levels.mod: newlevels.o
newlevels.o: ../femlib/mod/nls.mod
../femlib/mod/mod_nudge.mod: newnudge.o
newnudge.o: hydro.h param.h nudge.h femtime.h basin.h nudging.h
newrog.o: param.h basin.h
newstab.o: conz.h param.h nlevel.h nbasin.h levels.h basin.h \
		diff_visc_fric.h stab.h aux_array.h 
newtra.o: depth.h hydro_vel.h hydro.h param.h nlevel.h ev.h \
		nbasin.h geom_dynamic.h levels.h hydro_baro.h \
		basin.h hydro_print.h aux_array.h 
newtvd.o: depth.h hydro_vel.h param.h ev.h tvd.h levels.h basin.h \
		hydro_print.h 
nos2nc.o: depth.h hydro.h param.h levels.h evmain.h basin.h
nos_elab.o: param.h basin.h
nosadmin.o: param.h
nosaver_basin.o: depth.h param.h basin.h
nosaver_time.o: param.h
nosaver_vert.o: depth.h param.h basin.h
nosdiff.o: param.h
noselab.o: param.h ../femlib/mod/clo.mod basin.h
nosextr_cat.o: param.h
nosextr_custom.o: param.h basin.h
nosextr_fluxes.o: param.h basin.h simul.h
nosextr_gis.o: param.h basin.h
nosextr_nodes.o: depth.h param.h levels.h nlevel.h evmain.h \
		basin.h 
nosextr_records.o: param.h
nosinf.o: param.h ../femlib/mod/clo.mod
nosinit.o: param.h basin.h
nosintp.o: param.h basin.h
nosmaniq.o: param.h basin.h
nosresidence.o: nbasin.h depth.h param.h levels.h nlevel.h basin.h \
		evmain.h 
nossplit.o: param.h
nossum.o: param.h basin.h
noswrite.o: hydro.h param.h levels.h basin.h hydro_print.h simul.h
optintp.o: depth.h param.h ../femlib/mod/clo.mod evmain.h basin.h
ous2nc.o: depth.h hydro.h param.h levels.h evmain.h basin.h
ousaver.o: hydro.h param.h levels.h basin.h hydro_print.h simul.h
ousextr_cat.o: hydro_vel.h depth.h hydro.h param.h levels.h \
		nlevel.h basin.h 
ousextr_gis.o: depth.h hydro.h param.h levels.h nlevel.h basin.h
ousextr_nodes.o: depth.h hydro.h param.h levels.h nlevel.h \
		evmain.h basin.h hydro_print.h 
ousextr_records.o: hydro_vel.h depth.h hydro.h param.h levels.h \
		nlevel.h basin.h 
ousinf.o: depth.h hydro.h param.h nlevel.h evmain.h simul.h \
		levels.h basin.h 
ousintp.o: hydro.h param.h levels.h basin.h hydro_print.h simul.h
ousutil.o: ev.h
ousvdiff.o: depth.h param.h levels.h basin.h simul.h
ouswork.o: depth.h hydro.h param.h levels.h basin.h evmain.h \
		simul.h 
ouswork1.o: depth.h hydro.h param.h levels.h basin.h evmain.h \
		simul.h 
ouswrite.o: hydro.h param.h levels.h basin.h hydro_print.h simul.h
readext.o: param.h basin.h
scal_intp.o: param.h evmain.h basin.h
scalintp.o: depth.h pkonst.h param.h evmain.h basin.h simul.h
sedi3d.o: testbndo.h depth.h ts.h evmain.h nlevel.h nbasin.h \
		roughness.h diff_visc_fric.h geom.h waves.h \
		aux_array.h hydro_vel.h param.h hydro.h \
		bound_names.h sed_param.h geom_dynamic.h levels.h \
		area.h basin.h 
sigmautil.o: sigma.h
simsys_lp.o: nbasin.h param.h common.h
simsys_pard.o: param.h common.h
simsys_spk.o: param.h common.h
splitets.o: param.h
splitflx.o: param.h
subapn.o: ../femlib/mod/basin.mod
../femlib/mod/basin.mod: subbas.o
subbas.o: param.h basin.h
subbfm.o: depth.h hydro.h param.h bound_names.h sinking.h ts.h \
		nlevel.h nbasin.h pkonst.h levels.h basin.h \
		diff_visc_fric.h hydro_print.h bfm_common.h 
subbnd.o: param.h bound_names.h bnd.h nbvdim.h bound_geom.h \
		nbound.h 
subbndo.o: hydro.h param.h subbndo.h bound_geom.h nbasin.h \
		levels.h basin.h geom.h 
subbnds.o: ../femlib/mod/intp_fem_file.mod
subboxa.o: depth.h ev.h ts.h nlevel.h subboxa.h meteo.h \
		diff_visc_fric.h hydro_vel.h param.h conz.h \
		hydro.h femtime.h modules.h levels.h \
		geom_dynamic.h basin.h 
../femlib/mod/clo.mod: subclo.o
subcon1.o: depth.h param.h bound_dynamic.h nlevel.h femtime.h \
		mkonst.h simul.h nbasin.h levels.h 
subcoo.o: param.h common.h
subcoord.o: coords_cpp.h coords.h coords_utm.h coords_gb.h
subcst.o: param.h mkonst.h bound_geom.h nbound.h pkonst.h nbasin.h \
		modules.h basin.h 
subcus.o: depth.h internal.h ev.h nlevel.h ts.h nbasin.h \
		diff_visc_fric.h aux_array.h hydro_vel.h param.h \
		hydro.h conz.h femtime.h fluidmud.h meteo_aux.h \
		geom_dynamic.h levels.h hydro_baro.h basin.h \
		area.h hydro_print.h 
subdep.o: depth.h nbasin.h param.h basin.h
subdif.o: depth.h hydro.h param.h nlevel.h femtime.h ev.h \
		diff_aux.h nbasin.h geom_dynamic.h levels.h \
		basin.h hydro_print.h diff_visc_fric.h geom.h \
		aux_array.h 
subdim.o: nbasin.h
subdry.o: hydro.h param.h mkonst.h femtime.h ev.h geom_dynamic.h \
		hydro_baro.h basin.h aux_array.h 
subdts.o: subdts.h
subdwq.o: nbasin.h param.h geom_dynamic.h hydro_baro.h basin.h \
		femtime.h 
subele.o: depth.h hydro.h param.h nlevel.h ev.h nbasin.h levels.h \
		basin.h area.h aux_array.h 
subets.o: etsinf.h
subetsa.o: depth.h hydro.h param.h ../femlib/mod/ets.mod ts.h \
		nlevel.h femtime.h tides.h simul.h nbasin.h \
		levels.h modules.h hydro_print.h \
		../femlib/mod/nls.mod waves.h 
../femlib/mod/ets.mod: subetsm.o
../femlib/mod/extra.mod: subexta.o
subexta.o: param.h modules.h femtime.h hydro_print.h aux_array.h \
		simul.h ../femlib/mod/nls.mod 
../femlib/mod/intp_fem_file.mod: subfemintp.o
subflx3d.o: hydro_vel.h testbndo.h hydro.h param.h bound_dynamic.h \
		femtime.h nlevel.h ev.h geom_dynamic.h levels.h \
		hydro_baro.h basin.h links.h geom.h 
../femlib/mod/flux.mod: subflxa.o
subflxa.o: conz.h param.h modules.h nlevel.h ts.h femtime.h \
		../femlib/mod/nls.mod 
subflxu.o: param.h levels.h links.h nlevel.h
subfvl.o: depth.h nbasin.h param.h levels.h area.h nlevel.h \
		aux_array.h 
subgotm.o: depth.h ev.h ts.h nlevel.h nbasin.h turbulence.h \
		roughness.h diff_visc_fric.h gotm_aux.h \
		hydro_vel.h param.h hydro.h debug_aux2.h femtime.h \
		pkonst.h meteo_aux.h levels.h basin.h \
		hydro_print.h 
subgotm_mud.o: depth.h ev.h ts.h nlevel.h nbasin.h turbulence.h \
		roughness.h diff_visc_fric.h gotm_aux.h \
		hydro_vel.h hydro.h param.h debug_aux2.h femtime.h \
		sinking.h pkonst.h fluidmud.h levels.h meteo_aux.h \
		basin.h hydro_print.h 
../femlib/mod/grd.mod: subgrd.o
subgrd.o: pkonst.h param.h subgrd.h basin.h
subhisto.o: histo.h
subinterpol.o: depth.h param.h basin.h
subinterpol1.o: nbasin.h depth.h param.h basin.h
sublin.o: param.h links.h basin.h
sublnka.o: nbasin.h geom_aux.h param.h links.h basin.h geom.h
sublnkd.o: nbasin.h param.h geom_dynamic.h ev.h basin.h geom.h
sublnks.o: param.h ev.h femtime.h nbasin.h geom_dynamic.h links.h \
		basin.h geom.h 
sublnku.o: nbasin.h param.h links.h basin.h
sublpl.o: nbasin.h param.h ev.h basin.h
submet.o: param.h bound_dynamic.h ts.h nlevel.h femtime.h ev.h \
		nbasin.h ../femlib/mod/meteo_forcing_module.mod \
		meteo.h basin.h aux_array.h 
../femlib/mod/meteo_forcing_module.mod: submeteo2.o
submeteo2.o: depth.h param.h ev.h nlevel.h femtime.h nbasin.h \
		meteo_aux.h levels.h meteo.h basin.h \
		../femlib/mod/intp_fem_file.mod 
submud.o: depth.h param.h bound_names.h bound_dynamic.h sinking.h \
		ts.h nlevel.h bound_geom.h pkonst.h nbasin.h \
		fluidmud.h levels.h roughness.h basin.h gotm_aux.h \
		diff_visc_fric.h hydro_print.h 
submud_dummy.o: param.h fluidmud.h
subn11.o: hydro.h param.h bound_dynamic.h nlevel.h mkonst.h \
		femtime.h bound_geom.h pkonst.h nbasin.h \
		geom_dynamic.h levels.h basin.h hydro_print.h \
		../femlib/mod/intp_fem_file.mod bnd_aux.h 
../femlib/mod/chezy.mod: subn35.o
subn35.o: depth.h hydro.h param.h pkonst.h fluidmud.h levels.h \
		basin.h roughness.h hydro_print.h diff_visc_fric.h \
		../femlib/mod/nls.mod 
subnev.o: nbasin.h param.h basin.h ev.h
../femlib/mod/evgeom.mod: subnev_new.o
subnev_new.o: nbasin.h param.h basin.h
../femlib/mod/nls.mod: subnls.o
subnos.o: nosinf.h
subnosa.o: depth.h nbasin.h param.h
subnsa.o: param.h simul.h
subnsh.o: param.h semi.h nlevel.h femtime.h simul.h modules.h \
		levels.h basin.h 
subnsu.o: param.h basin.h
suboff.o: hydro_vel.h nbasin.h hydro.h param.h levels.h ts.h \
		femtime.h hydro_print.h 
subous.o: ousinf.h
subousa.o: depth.h hydro.h param.h nlevel.h femtime.h simul.h \
		nbasin.h levels.h 
subouta.o: nbasin.h hydro.h param.h hydro_baro.h femtime.h \
		hydro_print.h simul.h 
suboutput.o: femtime.h
../femlib/mod/para.mod: subpar3.o
subpard.o: ../femlib/mod/omp_lib.mod param.h common.h
subproj.o: param.h basin.h tides.h
subqfxf.o: subqfx.h
subqfxm1.o: subqfxm.h
subqfxm2.o: subqfxm.h
subqfxm3.o: subqfxm.h
subqfxm4.o: subqfxm.h
subqfxm5.o: subqfxm.h
subqfxt.o: param.h levels.h meteo.h ts.h subqfxm.h
subqfxu4.o: subqfxm.h
subreg.o: nbasin.h reg.h hydro.h param.h ev.h basin.h geom.h
subres.o: depth.h hydro.h param.h ts.h nlevel.h femtime.h simul.h \
		nbasin.h hydro_baro.h levels.h basin.h aux_array.h 
subrst.o: hydro_vel.h conz.h hydro.h param.h ts.h nlevel.h \
		femtime.h geom_dynamic.h basin.h 
subspk.o: param.h common.h basin.h
subtime.o: femtime.h
subtrace.o: hydro_vel.h hydro.h param.h femtime.h geom_dynamic.h \
		hydro_baro.h levels.h basin.h 
subtsuvfile.o: param.h ../femlib/mod/intp_fem_file.mod
subuti.o: depth.h hydro.h param.h ts.h ev.h pkonst.h \
		geom_dynamic.h levels.h hydro_baro.h basin.h \
		links.h 
subvola.o: param.h femtime.h nbasin.h volcomp.h modules.h basin.h \
		geom.h aux_array.h 
subwat.o: nbasin.h hydro.h param.h levels.h links.h ev.h basin.h
subwaves.o: depth.h internal.h evmain.h nlevel.h nbasin.h meteo.h \
		roughness.h geom.h aux_array.h waves.h hydro_vel.h \
		param.h hydro.h femtime.h pkonst.h hydro_baro.h \
		levels.h meteo_aux.h basin.h 
subwrt.o: depth.h conz.h param.h ev.h nlevel.h femtime.h simul.h \
		nbasin.h levels.h basin.h 
tideg.o: depth.h nbasin.h pkonst.h hydro.h param.h basin.h tides.h
ts2nc.o: depth.h hydro.h param.h levels.h basin.h evmain.h simul.h
tsinf.o: ../femlib/mod/clo.mod
vp.o: param.h basin.h
weutro.o: donata.h weutro.h
weutro_sedim.o: donata.h weutro.h
windcrea.o: depth.h param.h evmain.h basin.h
winmerge.o: param.h
winshift.o: param.h
zinit.o: param.h basin.h simul.h
aquabc_II.o: param.h ../femlib/mod/aquabc_II_sed_ini.mod \
		../femlib/mod/para_aqua.mod 
../femlib/mod/aquabc_II_sed_ini.mod: aquabc_II_ased_ini.o
aquabc_II_ased_ini.o: nbasin.h aquabc_II.h \
		../femlib/mod/para_aqua.mod 
../femlib/mod/AQUABC_II_GLOBAL.mod: aquabc_II_co2sys.o
../femlib/mod/VECTOR_MATRIX_UTILS.mod: aquabc_II_co2sys.o
../femlib/mod/CO2SYS_CDIAC.mod: aquabc_II_co2sys.o
aquabc_II_pelagic_lib.o: ../femlib/mod/AQUABC_II_GLOBAL.mod
aquabc_II_pelagic_model.o: ../femlib/mod/CO2SYS_CDIAC.mod param.h \
		../femlib/mod/para_aqua.mod \
		../femlib/mod/AQUABC_II_GLOBAL.mod 
aquabc_II_sediment_lib.o: ../femlib/mod/AQUABC_II_GLOBAL.mod
aquabc_II_sediment_model_1.o: ../femlib/mod/CO2SYS_CDIAC.mod \
		param.h ../femlib/mod/para_aqua.mod \
		../femlib/mod/AQUABC_II_GLOBAL.mod 

