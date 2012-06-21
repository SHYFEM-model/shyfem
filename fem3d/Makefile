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

LIBG_MKL =
SOLVER_FILES = simsys_lp.o
ifeq ($(SOLVER),PARDISO)
  LIBG_MKL = `./bin/find_pardlib.sh`
  SOLVER_FILES = simsys_pard.o subpard.o subcoo.o
endif
ifeq ($(SOLVER),SPARSKIT)
  SOLVER_FILES = simsys_spk.o subspk.o subcoo.o \
		spk_ilut.o spk_iters.o spk_itaux.o \
		spk_itpack.o spk_matvec.o spk_blas1_mod.o
endif

NETCDF_OBJ = subnetcdf_dummy.o
DIRLIB_NETCDF = $(NETCDFDIR)/lib/
DIRINC_NETCDF = $(NETCDFDIR)/include
ifeq ($(NETCDF),true)
  NETCDF_OBJ       = subnetcdf_admin.o subnetcdf_hydro.o
  NETCDF_OBJ       = netcdf.o
  LIBG_NETCDF      = -L$(DIRLIB_NETCDF) -l$(LIBNETCDF)
  LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a
  INCL_NETCDF      = -I$(DIRINC_NETCDF)
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
	  new3di.o newpr.o \
	  newrog.o newcra.o \
	  newbcl.o newbcl0.o \
	  newini.o newtra.o \
	  newchao.o newparam.o \
	  newcon.o newcon_1.o newconz.o \
	  newchk.o newdum.o newfix.o \
	  newbsig.o newtvd.o newstab.o newsig.o

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
          subnls.o subsse.o subwat.o \
          subfil.o subn35.o subout.o subous.o subssg.o \
          subnsh.o subssm.o subzzz.o \
          subnsu.o subsss.o subscn.o subres.o \
          subapn.o subnsa.o \
          subgrd.o subreg.o \
          subwin.o subbnd.o \
          subexta.o subouta.o subdep.o subflxa.o subflx3d.o subflx.o \
	  subdum.o subcus.o subdwq.o \
	  $(QFLUX) \
	  subnos.o subcon1.o subrnd.o sublin.o \
	  submod.o subuti.o subele.o subgeo.o subdts.o \
	  links.o sublnka.o sublnks.o sublnku.o subisd.o \
	  debug.o subhisto.o subnosa.o subfvl.o \
	  subeos.o sigmautil.o \
	  nosutil.o ousutil.o

# general objects for 3D model ht

HTOBJS = \
	subver.o \
	subcst.o subsys.o subnsh.o \
	subbas.o subdim.o subnev.o \
	subnsu.o subn35.o subdep.o \
	subdef.o subfil.o subnls.o \
	$(PARMOD) \
	subsss.o subscn.o subssv.o subsse.o subssm.o \
	submir.o \
	subext.o subexta.o \
	subous.o subousa.o \
	subbnd.o subnos.o \
	subwin.o \
	subflxa.o subflx3d.o subflx.o subwat.o subn11.o subdry.o \
	subcus.o subcon1.o \
	tideg.o tide1.o subzzz.o subres.o subuti.o \
	sublin.o subvola.o subgeo.o debug.o \
	submod.o subele.o \
	$(SOLVER_FILES) \
	$(QFLUX) \
	subgotm.o $(GOTM_OBJ) \
	subbndo.o subbnds.o \
	new36.o intp.o subssf.o \
	bdist.o links.o sublnka.o sublnks.o sublnku.o \
	sedi3d.o sedtrans05.o subwaves.o \
	subhisto.o \
	$(LAGRANGE) \
	subdts.o loading.o \
	subrst.o subnosa.o subdif.o subreg.o \
	atoxi3d.o atoxi.o subwwm.o kepsd.o keps_util.o \
	subinv.o newexpl.o subrnd.o \
	submet.o subrgf.o subcoord.o subgrd.o sort.o subfvl.o \
	subtrace.o subproj.o subeos.o sigmautil.o \
	submeteo.o submeteo1.o submeteo_inout.o \
	subwrt.o \
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
	  debug.o subreg.o subnev.o

RSTOBJS = subrst.o ecological_dummy.o

BASEXTRA = subgeo.o subgrd.o subsrt.o sublnka.o links.o sublnku.o

NOSOBJS = subnos.o nosutil.o genutil.o sigmautil.o subeos.o
OUSOBJS = subous.o ousutil.o genutil.o sigmautil.o
NOSOUSOBJS = subnos.o nosutil.o subous.o ousutil.o genutil.o sigmautil.o
GRDOBJS = subgrd.o grdutil.o subsrt.o subcoord.o

LAPLACE = sublpl.o submir.o subnos.o links.o sublnku.o \
	  subwin.o subnosa.o

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
	subqfxm1.o subqfxm2.o subqfxm3.o subqfxm4.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o

QFLUX_TEST = subqfxf.o \
	subqfxm1.o subqfxm2.o subqfxm3.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o \
	subfil.o

PARMOD = subpar1.o
PARMOD = subpar.o subfnm.o

# objects for interpolation

INTP    = subzzz.o subreg.o

# objects for lagrangian interpolation

LAGINTP = submir.o

#------------------------------------------------

EXEMAIN  = ht vp

EXEINFO  = ousinf flxinf volinf basinf nosinf wininf extinf rstinf
EXESPLIT = nossplit splitext splitflx
EXELAP   = laplap lapwin
EXEBAS   = bastreat basbathy
EXEVAR   = nosinit nosresidence scal_intp bashsigma

NOSAVER  = nosaver_basin nosaver_vert nosaver_time
NOSEXTR  = nosextr_nodes nosextr_gis nosextr_records nosextr_custom \
		nosextr_cat
OUSEXTR  = ousextr_nodes ousextr_records ousextr_gis \
		ousextr_cat

EXEEXTRA = bedload wincrea windcrea

EXES     = $(EXEMAIN) $(EXEINFO) $(EXESPLIT) $(EXELAP) $(EXEBAS) \
		$(EXEVAR) $(NOSEXTR) $(OUSEXTR) $(NOSAVER) \
		$(EXEEXTRA)

#------------------------------------------------

default:	fem

all:	param fem
fem:	param libfem main info split lap bas var util links
main:	param $(EXEMAIN)
info:	param $(EXEINFO)
split:	param $(EXESPLIT)
lap:	param $(EXELAP)
bas:	param $(EXEBAS)
var:	param $(EXEVAR)
util:	param $(NOSEXTR) $(OUSEXTR) $(NOSAVER)

links:
	@$(FEMBIN)/symlink $(FEMBIN) \
		$(EXEMAIN) $(EXEINFO) $(EXESPLIT)

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

cubic:   cubic.o $(OUTOBJS) subsse.o subssm.o
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) subsse.o subssm.o -o $@

#------------------------------------------------

basinf:  basinf.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

baswork:  baswork.o $(OUTOBJS) $(BASEXTRA)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(BASEXTRA) -o $@

bashsigma:  bashsigma.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

baschk:  baschk.o subnev.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnev.o $(OUTOBJS) -o $@

basbathy:  basbathy.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) -o $@

basres:  basres.o $(OUTOBJS) $(NOSOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(NOSOBJS) -o $@

basproj:  basproj.o $(OUTOBJS) $(GRDOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(GRDOBJS) -o $@

bastreat:  bastreat.o $(BASEXTRA) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(BASEXTRA) $(OUTOBJS) -o $@

scal_intp:  scal_intp.o subgeo.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) $(LAPLACE) -o $@

bascomp:  bascomp.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

ousinf:  ousinf.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ouswork:  ouswork.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousintp:  ousintp.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousaver:  ousaver.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousvdiff:  ousvdiff.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

nosinf:  nosinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosinf2:  nosinf2.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nossum:  nossum.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

ous2nc:  ous2nc.o $(OUSOBJS) $(OUTOBJS) netcdf.o
	$(LINKER) -o $@ $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) netcdf.o $(LIBGS)

eosinf:  eosinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

#-------------------------------------------------------

.PHONY: nosextr ousextr nosaver

ousextr: $(OUSEXTR)

ousextr_nodes:  ousextr_nodes.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousextr_records:  ousextr_records.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousextr_gis:  ousextr_gis.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

ousextr_cat:  ousextr_cat.o $(OUSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUSOBJS) $(OUTOBJS) -o $@

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

laplap:   laplap.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@ 

lapwin:   lapwin.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@ 

rstinf:   rstinf.o $(OUTOBJS) $(RSTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(RSTOBJS) -o $@ 

qflux_test:   qflux_test.o $(QFLUX_TEST)
	$(LINKER) $(LFLAGS) $@.o $(QFLUX_TEST) -o $@ 

#------------------------------------------------

# to test interpolation routines in intp.f

INTPOBJS = subsse.o subfil.o subsss.o subscn.o subssf.o subssm.o

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
		$(EXEMAIN) $(EXEINFO) $(EXESPLIT)

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
aquabc_fem_interface.o: param.h aquabc_aout.h aquabc.h ev.h
atoxi3d.o: param.h
basbathy.o: param.h evmain.h
baschk.o: param.h evmain.h
bascomp.o: param.h evmain.h
bashsigma.o: param.h evmain.h basin.h
basinf.o: param.h evmain.h basin.h
basproj.o: param.h basin.h
basres.o: param.h basin.h
bastreat.o: param.h links.h evmain.h basin.h
baswork.o: param.h evmain.h basin.h
bdist.o: param.h links.h
bedload.o: param.h evmain.h
bio3d.o: param.h
bio3d_util.o: param.h donata.h
eosinf.o: param.h
flux.o: regflux.h param.h links.h ev.h
flxinf.o: param.h
gotm3d.o: gotmturb.i
grdutil.o: grd_extra.h param.h basin.h
ht.o: param.h evmain.h
lagrange_back.o: param.h lagrange.h
lagrange_connect.o: param.h lagrange.h lagrange_connect.h
lagrange_cont.o: param.h links.h lagrange.h
lagrange_decay.o: param.h lagrange.h
lagrange_dif.o: param.h lagrange.h
lagrange_flux.o: param.h ev.h links.h lagrange.h
lagrange_init.o: param.h
lagrange_inout.o: param.h ev.h lagrange.h
lagrange_larve.o: param.h lagrange.h
lagrange_main.o: param.h lagrange.h
lagrange_oil.o: param.h ev.h lagrange.h
lagrange_track.o: param.h lagrange.h
lagrange_util.o: ev.h
laplap.o: param.h evmain.h
lapwin.o: param.h evmain.h
lgrelab.o: param.h
lgrinf.o: param.h
loading.o: param.h
map_influence.o: param.h
matrix49.o: regflux.h
netcdf.o: param.h netcdf.h netcdf.inc
new36.o: close.h
new3di.o: testbndo.h param.h ev.h
newbcl.o: param.h
newbsig.o: bsig.h param.h
newchao.o: param.h
newchk.o: testbndo.h param.h ev.h
newcon.o: testbndo.h param.h ev.h
newcon_1.o: testbndo.h param.h ev.h
newconz.o: param.h
newcra.o: ev.h
newexpl.o: param.h ev.h
newfix.o: param.h bclfix.h
newini.o: param.h
newstab.o: param.h
newtra.o: param.h ev.h
newtvd.o: param.h ev.h tvd.h
nos_elab.o: param.h
nosaver_basin.o: param.h basin.h
nosaver_time.o: param.h
nosaver_vert.o: param.h basin.h
nosdiff.o: param.h
nosextr_cat.o: param.h
nosextr_custom.o: param.h
nosextr_fluxes.o: param.h
nosextr_gis.o: param.h
nosextr_nodes.o: param.h
nosextr_records.o: param.h
nosinf.o: param.h
nosinit.o: param.h
nosmaniq.o: param.h
nosresidence.o: param.h basin.h evmain.h
nossplit.o: param.h
nossum.o: param.h
nosutil.o: ev.h
ous2nc.o: param.h evmain.h
ousaver.o: param.h
ousextr_cat.o: param.h
ousextr_gis.o: param.h
ousextr_nodes.o: param.h
ousextr_records.o: param.h
ousinf.o: param.h evmain.h
ousintp.o: param.h
ousvdiff.o: param.h
readext.o: param.h
scal_intp.o: param.h evmain.h
sedi3d.o: param.h links.h ev.h sed_param.h
simsys_lp.o: param.h common.h
simsys_pard.o: param.h common.h
simsys_spk.o: param.h common.h
subbfm.o: param.h bfm_common.h
subbndo.o: param.h subbndo.h
subcoo.o: param.h common.h
subcst.o: modules.h
subcus.o: param.h ev.h
subdep.o: ev.h
subdif.o: param.h ev.h
subdry.o: ev.h
subdwq.o: param.h
subele.o: param.h ev.h
subexta.o: modules.h
subflx3d.o: testbndo.h param.h links.h ev.h
subflxa.o: param.h modules.h links.h
subfvl.o: param.h
subgotm.o: param.h
sublin.o: links.h
sublnka.o: ev.h links.h
sublnks.o: ev.h links.h
sublnku.o: links.h
sublpl.o: ev.h
submet.o: ev.h
subn11.o: param.h
subn35.o: param.h chezy.h
subnetcdf_admin.o: param.h
subnetcdf_hydro.o: netcdf.inc
subnev.o: ev.h
subnsh.o: modules.h
subousa.o: param.h
subpar1.o: subpar.h
subpard.o: param.h common.h
subproj.o: param.h
subqfxf.o: subqfx.h
subqfxm4.o: subqfxm4.h
subreg.o: ev.h
subres.o: param.h
subrst.o: param.h
subspk.o: param.h common.h
subtrace.o: param.h
subuti.o: param.h ev.h links.h
subvola.o: modules.h
subwat.o: links.h ev.h
subwaves.o: param.h
subwrt.o: param.h
subwwm.o: param.h evmain.h
tideg.o: param.h
vp.o: param.h
weutro.o: donata.h weutro.h
weutro_sedim.o: donata.h weutro.h
windcrea.o: param.h evmain.h
wininf.o: param.h
winmerge.o: param.h
winshift.o: param.h
zinit.o: param.h

