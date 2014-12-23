#
# $Id: Makefile,v 1.81 2010-03-11 15:36:38 georg Exp $
#
#------------------------------------------------------------

DIR	= fem3d

#------------------------------------------------------------ just for default
FEMDIR = ..
DIRLIB = $(FEMDIR)/femlib
F77    = g77
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
  ECOLOGICAL_OBJ = aquabc_fem_interface.o alukas.o aquabc.o \
		aquabc_sediment_model_1.o aquabc_util.o
endif

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

NEWOBJS = \
	  new3di.o \
	  newrog.o newcra.o \
	  newbcl.o newbcl0.o \
	  newini.o newtra.o \
	  newchao.o \
	  newcon.o newconz.o \
	  newchk.o newdum.o newfix.o \
	  newtvd.o newstab.o newsig.o \
	  newnudge.o

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
          subssm.o subzzz.o \
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
	  genutil.o subfemintp.o subtsfile.o

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
	tideg.o tide1.o subzzz.o subres.o subuti.o \
	sublin.o subvola.o subgeo.o debug.o \
	submod.o subele.o \
	$(SOLVER_FILES) \
	$(QFLUX) \
	$(MUD_FILES) $(GOTM_OBJ) \
	subbndo.o subbnds.o \
	new36.o intp.o subssf.o \
	bdist.o links.o sublnka.o sublnks.o sublnku.o \
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
	subets.o subetsa.o \
	$(IFF_FILES) \
	$(ECOLOGICAL_OBJ) \
	$(NETCDF_OBJ) $(PARALLEL_OBJ)

# objects for pre-processor vp

VPOBJS = \
	subver.o \
        subros.o subcmv.o \
        subdef.o subbas.o subdim.o \
        subfil.o \
        $(PARMOD) \
        subsrt.o subgrd.o \
        subsss.o subscn.o subnsu.o \
	subapn.o subnsa.o subnls.o subsys.o subdum.o subbnd.o

# objects for output routines

OUTOBJS = $(PARMOD) subnls.o subsss.o subscn.o \
	  subfil.o subdef.o subnsu.o subssv.o \
	  subbas.o subsys.o subdim.o \
	  subapn.o subnsa.o subdum.o \
	  subbnd.o subdts.o subdep.o \
	  debug.o subreg.o subnev.o \
	  subfemfile.o subtsfile.o \
	  subvintp.o subver.o sigmautil.o suboutput.o \
	  subclo.o

RSTOBJS = subrst.o ecological_dummy.o

BASEXTRA = subgeo.o subgrd.o subsrt.o sublnka.o links.o sublnku.o \
		subinterpol.o subinterpol1.o subnos.o

NOSOBJS = subnos.o nosutil.o genutil.o subeos.o
OUSOBJS = subous.o ousutil.o genutil.o 
ETSOBJS = subets.o etsutil.o genutil.o 
NOSOUSOBJS = subnos.o nosutil.o subous.o ousutil.o subets.o etsutil.o \
		genutil.o
GRDOBJS = subgrd.o grdutil.o subsrt.o subcoord.o

BATHYOBJS = $(BASEXTRA)
#BATHYOBJS = subgeo.o subgrd.o subsrt.o subinterpol.o subinterpol1.o \
#		subnos.o

LAPLACE = sublpl.o submir.o subnos.o links.o sublnku.o \
	  subnosa.o subssm.o

LAGRANGE =	lagrange_main.o \
		lagrange_back.o lagrange_dif.o \
		lagrange_cont.o lagrange_init.o \
		lagrange_inout.o lagrange_util.o \
		lagrange_larve.o lagrange_flux.o \
		lagrange_track.o lagrange_util_tr.o \
		lagrange_decay.o \
		lagrange_connect.o

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

INTP    = subzzz.o subreg.o

# objects for lagrangian interpolation

LAGINTP = submir.o

#------------------------------------------------

EXEMAIN  = ht vp

EXEINFO  = ousinf flxinf volinf basinf nosinf wininf extinf rstinf feminf \
		offinf lgrinf etsinf
EXESPLIT = nossplit splitext splitflx splitets
EXELAP   =  
EXEBAS   = bastreat basbathy scalintp
EXEVAR   = nosinit nosresidence scal_intp bashsigma
EXEUTIL  = bc2fem

NOSVAR   = nos_elab nossum nosdiff noswrite nosadmin
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

EXES     = $(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXELAP) $(EXEBAS) \
		$(EXEVAR) $(EXEUTIL) $(NOSEXTR) $(OUSEXTR) $(NOSAVER) \
		$(EXEEXTRA) $(EXENC) $(NOSVAR)

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

basproj:  basproj.o $(OUTOBJS) $(GRDOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(GRDOBJS) -o $@

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

fileinf:  fileinf.o $(NOSOUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOUSOBJS) $(OUTOBJS) -o $@

feminf:  feminf.o $(NOSOBJS) $(OUTOBJS) 
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
	$(MAKEDEPEND) *.f

#---------------------------------------------------

options:
	@echo "compiler  $(COMPILER) $(F77)"
	@echo "solver    $(SOLVER)"

#---------------------------------------------------

clean: cleanobj cleanexe cleanvar cleantex cleanlahey cleanexport cleandoc \
	cleanifc cleanlinks

cleanall: clean cleanzip cleanlib cleanextra cleandiff

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

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.f.o:
	$(F77) -c $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

#--------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

aquabc.o: param.h
aquabc_fem_interface.o: param.h aquabc_aout.h aquabc.h femtime.h \
		ev.h 
atoxi3d.o: param.h
basbathy.o: param.h evmain.h
basbox.o: param.h links.h evmain.h basin.h
baschk.o: param.h evmain.h basin.h
bascomp.o: param.h evmain.h basin.h
bashsigma.o: param.h evmain.h basin.h
basinf.o: param.h evmain.h basin.h
basproj.o: param.h basin.h
basres.o: param.h basin.h
bastreat.o: param.h links.h evmain.h basin.h
baswork.o: param.h evmain.h basin.h
bdist.o: param.h nkonst.h links.h
bio3d.o: param.h
bio3d_util.o: param.h donata.h
eosinf.o: param.h
etsinf.o: param.h
feminf.o: clo.mod
flux.o: regflux.h param.h links.h ev.h femtime.h
flxinf.o: param.h
gotm3d.o: gotmturb.i
grdutil.o: grd_extra.h param.h basin.h
ht.o: intp_fem_file.mod param.h meteo.h evmain.h femtime.h
lagrange_back.o: param.h lagrange.h
lagrange_connect.o: param.h femtime.h lagrange.h \
		lagrange_connect.h 
lagrange_cont.o: param.h links.h lagrange.h
lagrange_decay.o: param.h femtime.h lagrange.h
lagrange_dif.o: param.h femtime.h lagrange.h
lagrange_flux.o: param.h ev.h femtime.h links.h lagrange.h
lagrange_init.o: param.h
lagrange_inout.o: param.h ev.h femtime.h lagrange.h
lagrange_larve.o: param.h lagrange.h
lagrange_main.o: param.h femtime.h lagrange.h
lagrange_oil.o: param.h femtime.h ev.h lagrange.h
lagrange_track.o: param.h lagrange.h
lagrange_util.o: ev.h
lgrelab.o: param.h
lgrinf.o: param.h
loading.o: param.h
map_influence.o: param.h
matrix49.o: regflux.h
netcdf.o: param.h netcdf.h netcdf.inc
netcdf_util.o: param.h
new36.o: nkonst.h femtime.h close.h
new3di.o: testbndo.h param.h ev.h femtime.h
newbcl.o: param.h femtime.h
newchao.o: param.h femtime.h
newchk.o: testbndo.h param.h ev.h femtime.h
newcon.o: testbndo.h param.h ev.h femtime.h
newconz.o: param.h femtime.h
newcra.o: femtime.h ev.h
newexpl.o: param.h ev.h
newfix.o: param.h nkonst.h femtime.h bclfix.h
newini.o: param.h
newnohydro.o: param.h nohydlinks.h links.h femtime.h ev.h
newnudge.o: param.h nudge.h femtime.h basin.h
newstab.o: param.h
newtra.o: param.h ev.h
newtvd.o: param.h ev.h tvd.h
noh_ggu.o: param.h femtime.h ev.h
nonhydro_pard.o: param.h nohydlinks.h common.h nkonst.h
nonhydro_spk.o: param.h nohydlinks.h
nos2nc.o: param.h evmain.h basin.h
nos_elab.o: param.h
nosadmin.o: param.h
nosaver_basin.o: param.h basin.h
nosaver_time.o: param.h
nosaver_vert.o: param.h basin.h
nosdiff.o: param.h
nosextr_cat.o: param.h
nosextr_custom.o: param.h
nosextr_fluxes.o: param.h
nosextr_gis.o: param.h
nosextr_nodes.o: param.h evmain.h basin.h
nosextr_records.o: param.h
nosinf.o: param.h
nosinit.o: param.h
nosintp.o: param.h
nosmaniq.o: param.h
nosresidence.o: param.h basin.h evmain.h
nossplit.o: param.h
nossum.o: param.h
nosutil.o: ev.h
noswrite.o: param.h
optintp.o: param.h evmain.h basin.h
ous2nc.o: param.h evmain.h basin.h
ousaver.o: param.h
ousextr_cat.o: param.h
ousextr_gis.o: param.h
ousextr_nodes.o: param.h evmain.h basin.h
ousextr_records.o: param.h
ousinf.o: param.h evmain.h
ousintp.o: param.h
ousutil.o: param.h ev.h
ousvdiff.o: param.h
ouswork.o: param.h evmain.h
ouswork1.o: param.h evmain.h
ouswrite.o: param.h
readext.o: param.h
scal_intp.o: param.h evmain.h basin.h
scalintp.o: param.h evmain.h
sedi3d.o: testbndo.h param.h evmain.h sed_param.h waves.h
simsys_lp.o: param.h common.h nkonst.h
simsys_pard.o: param.h common.h nkonst.h
simsys_spk.o: param.h common.h nkonst.h
splitets.o: param.h
splitflx.o: param.h
subbas.o: param.h basin.h
subbfm.o: param.h bfm_common.h
subbnd.o: param.h nkonst.h bnd.h nbound.h
subbndo.o: param.h subbndo.h nkonst.h
subbnds.o: intp_fem_file.mod
subboxa.o: param.h meteo.h modules.h ev.h basin.h subboxa.h \
		femtime.h 
clo.mod: subclo.o
subcon1.o: femtime.h
subcoo.o: param.h common.h nkonst.h
subcst.o: modules.h
subcus.o: param.h nkonst.h ev.h femtime.h
subdep.o: ev.h
subdif.o: param.h femtime.h ev.h
subdim.o: nbound.h
subdry.o: femtime.h ev.h
subdts.o: subdts.h
subdwq.o: param.h femtime.h
subele.o: param.h ev.h
subets.o: etsinf.h
subetsa.o: param.h modules.h femtime.h waves.h subetsa.h
subexta.o: modules.h femtime.h
intp_fem_file.mod: subfemintp.o
subflx3d.o: testbndo.h param.h femtime.h links.h ev.h
subflxa.o: param.h subflxa.h modules.h femtime.h
subflxu.o: param.h links.h
subfvl.o: param.h
subgotm.o: param.h ev.h femtime.h
subgotm_mud.o: param.h ev.h femtime.h
subinterpol.o: param.h
subinterpol1.o: param.h
sublin.o: links.h
sublnka.o: nkonst.h ev.h links.h
sublnks.o: ev.h femtime.h links.h
sublnku.o: links.h
sublpl.o: ev.h
submet.o: param.h meteo_forcing_module.mod meteo.h femtime.h ev.h
meteo_forcing_module.mod: submeteo2.o
submeteo2.o: intp_fem_file.mod meteo.h femtime.h
submud.o: param.h
submud_dummy.o: param.h
subn11.o: intp_fem_file.mod param.h nkonst.h femtime.h
subn35.o: param.h chezy.h
subnev.o: ev.h
subnos.o: nosinf.h
subnsh.o: param_dummy.h modules.h basin.h femtime.h nbound.h
suboff.o: param.h femtime.h
subous.o: ousinf.h
subousa.o: param.h femtime.h
subouta.o: femtime.h
suboutput.o: femtime.h
para.mod: subpar3.o
subpard.o: param.h omp_lib.mod common.h nkonst.h
subproj.o: param.h
subqfxf.o: subqfx.h
subqfxm1.o: subqfxm.h
subqfxm2.o: subqfxm.h
subqfxm3.o: subqfxm.h
subqfxm4.o: subqfxm.h
subqfxm5.o: subqfxm.h
subqfxt.o: param.h meteo.h subqfxm.h
subqfxu4.o: subqfxm.h
subreg.o: ev.h
subres.o: param.h femtime.h
subrst.o: param.h femtime.h
subspk.o: param.h common.h nkonst.h
subtime.o: femtime.h
subtrace.o: param.h femtime.h
subtsuvfile.o: intp_fem_file.mod param.h
subuti.o: param.h ev.h links.h
subvola.o: modules.h femtime.h
subwat.o: links.h ev.h
subwaves.o: param.h meteo.h evmain.h femtime.h waves.h
subwrt.o: param.h femtime.h
tideg.o: param.h
ts2nc.o: param.h evmain.h
vp.o: param.h basin.h
weutro.o: donata.h weutro.h
weutro_sedim.o: donata.h weutro.h
windcrea.o: param.h evmain.h basin.h
winmerge.o: param.h
winshift.o: param.h
zinit.o: param.h

