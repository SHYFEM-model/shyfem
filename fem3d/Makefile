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
LIBNETCDF  = netcdf
LIBERSEM   = bfm

LIBFFEM    = $(DIRLIB)/lib$(LIBFEM).a

LIBLS      = -L$(DIRLIB) -L$(DIRLIB_NETCDF)

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a

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
DIRLIB_NETCDF = $(NETCDFDIR)/lib
DIRINC_NETCDF = $(NETCDFDIR)/include
ifeq ($(NETCDF),true)
  NETCDF_OBJ       = subnetcdf_admin.o subnetcdf_hydro.o
  LIBG_NETCDF      = -l$(LIBNETCDF)
  LIBF_NETCDF      = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a
  INCL_NETCDF      = -I$(DIRINC_NETCDF)
endif

GOTM_OBJ = gotm3d.o
ifeq ($(GOTM),true)
  GOTM_OBJ =
  LIBG_GOTM      = -l$(LIBGOTM)
  LIBF_GOTM      = $(DIRLIB)/lib$(LIBGOTM).a
endif

ERSEM_OBJ = subbfm_dummy.o
ifeq ($(ERSEM),true)
  ERSEM_OBJ	 = subbfm.o
  LIBG_ERSEM     = -l$(LIBERSEM)
  LIBF_ERSEM     = $(DIRLIB)/lib$(LIBERSEM).a
endif

#------------------------------------------------------------

LIBGS      = $(LIBLS) $(LIBG_NETCDF) $(LIBG_GOTM) $(LIBG_ERSEM) $(LIBG_MKL)
LIBFS      = $(LIBF_NETCDF) $(LIBF_GOTM) $(LIBF_ERSEM)
INCLUDE    = $(INCL_NETCDF)

#------------------------------------------------------------

SPECIAL =	Makefile \
		COMMIT README TESTLOG TODO VERSION WARNINGS \
		ACTUAL LOG

NOSEXTR = nosextr_nodes nosextr_gis nosextr_records nosextr_custom
OUSEXTR = ousextr_records ousextr_gis

EXES	= ht hp vp \
	  ousinf whdf whdfs \
	  zlevel ttt readext zinit wincrea outinf \
	  basinf flxinf wininf nosinf nossplit nosintp \
	  nosmaniq $(NOSAVER) \
	  splitext outintp volinf extinf bastreat basbathy baschk \
	  nosdiff nostransf mkgeom lagint splitflx \
	  laplap lapwin nosinit ousintp nosexport \
	  bascomp rstinf nosresidence ousaver ousvdiff \
	  lgrelab lgrinf qflux_test bedload scal_intp \
	  $(NOSEXTR) $(OUSEXTR) map_influence \
	  splitext_big windcrea winshift winmerge

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
	  newsigma.o newtvd.o newstab.o

# objects for library

SOBJS   = \
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
          subexta.o subouta.o subdep.o subflxa.o \
	  subdum.o subcus.o subdwq.o \
	  $(QFLUX) \
	  subnos.o subcon1.o subrnd.o sublin.o \
	  submod.o subuti.o subele.o subgeo.o subdts.o \
	  links.o sublnka.o sublnks.o sublnku.o subisd.o \
	  debug.o subhisto.o subnosa.o subfvl.o

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
	subflxa.o subwat.o subn11.o subdry.o \
	subcus.o subcon1.o \
	tideg.o tide1.o subzzz.o subres.o subuti.o \
	sublin.o subvola.o subgeo.o debug.o \
	submod.o subele.o \
	$(SOLVER_FILES) \
	$(QFLUX) \
	subgotm.o $(GOTM_OBJ) \
	$(NETCDF_OBJ) \
	$(ERSEM_OBJ) \
	subbndo.o subbnds.o \
	weutro.o bio3d.o bio3d_util.o new36.o intp.o subssf.o \
	bdist.o links.o sublnka.o sublnks.o sublnku.o \
	sedi3d.o sedtrans05.o subwaves.o \
	subhisto.o \
	$(LAGRANGE) \
	weutro_sedim.o weutro_light.o subdts.o loading.o \
	subrst.o subnosa.o subdif.o subreg.o \
	atoxi3d.o atoxi.o subwwm.o kepsd.o keps_util.o \
	subinv.o newexpl.o subrnd.o \
	submeteo.o subrgf.o subcoord.o subgrd.o sort.o subfvl.o \
	subtrace.o

# objects for pre-processor vp

VPOBJS = \
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
	  debug.o subreg.o subnev.o subrst.o

NOSOBJS = subnos.o nosutil.o

LAPLACE = sublpl.o submir.o subnos.o links.o sublnku.o \
	  subwin.o subnosa.o

LAGRANGE =	lagrange_main.o \
		lagrange_back.o lagrange_dif.o \
		lagrange_cont.o lagrange_init.o \
		lagrange_inout.o lagrange_util.o \
		lagrange_larve.o lagrange_flux.o \
		lagrange_track.o lagrange_util_tr.o \
		lagrange_decay.o

# heat flux module

QFLUX = subqfxt.o subqfxf.o \
	subqfxm1.o subqfxm2.o subqfxm3.o \
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

default:	fem

all:	param fem
fem:	param libfem main info util 
main:	param ht vp
info:	param ousinf flxinf volinf basinf nosinf wininf extinf
util:	param nossplit splitext ousinf splitflx laplap lapwin \
	nosinit bastreat basbathy scal_intp nosresidence rstinf \
	$(NOSEXTR) $(OUSEXTR) $(NOSAVER)

#--------------------------------------------- EXES

bedload:  bedload.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS) -o $@

ht:  ht.o $(NEWOBJS) $(HTOBJS) $(LIBFS)
	$(LINKER) $(LFLAGS) $@.o $(NEWOBJS) $(HTOBJS) $(LIBGS) -o $@

vp:  vp.o $(VPOBJS)
	$(LINKER) $(LFLAGS) $@.o $(VPOBJS) -o $@

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

splitflx:  splitflx.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

extinf:  extinf.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subext.o $(OUTOBJS) -o $@

flxinf:  flxinf.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subext.o $(OUTOBJS) -o $@

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

baschk:  baschk.o subnev.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnev.o $(OUTOBJS) -o $@

basbathy:  basbathy.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) -o $@

bastreat:  bastreat.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) -o $@

scal_intp:  scal_intp.o subgeo.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) $(LAPLACE) -o $@

bascomp:  bascomp.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

ousinf:  ousinf.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousintp:  ousintp.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousaver:  ousaver.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousvdiff:  ousvdiff.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

nosinf:  nosinf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

concinf:  concinf.o $(NOSOBJS) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosinf2:  nosinf2.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

.PHONY: nosextr ousextr 

ousextr: $(OUSEXTR)

ousextr_records:  ousextr_records.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousextr_gis:  ousextr_gis.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

nosextr: $(NOSEXTR)

nosextr_nodes:  nosextr_nodes.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_records:  nosextr_records.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_gis:  nosextr_gis.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosextr_custom:  nosextr_custom.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

map_influence:  map_influence.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nossplit:  nossplit.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

#-------------------------------------------------------

.PHONY: nosaver
NOSAVER = nosaver_basin nosaver_vert nosaver_time
NOSOBJS = subnos.o nosutil.o

nosaver: $(NOSAVER)

nosaver_basin:  nosaver_basin.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosaver_vert:  nosaver_vert.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosaver_time:  nosaver_time.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosresidence:  nosresidence.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

#-------------------------------------------------------

nosmaniq:  nosmaniq.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nostransf:  nostransf.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosdiff:  nosdiff.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosinit:  nosinit.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosexport:  nosexport.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nos_elab:  nos_elab.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

lgrelab:  lgrelab.o $(NOSOBJS) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(OUTOBJS) -o $@

nosintp:  nosintp.o $(NOSOBJS) $(INTP) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(NOSOBJS) $(INTP) $(OUTOBJS) -o $@

outintp:  outintp.o subout.o $(INTP) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subout.o $(INTP) $(OUTOBJS) -o $@

outinf:   outinf.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

outinfs:   outinfs.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

laplap:   laplap.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@ 

lapwin:   lapwin.o $(OUTOBJS) $(LAPLACE)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) $(LAPLACE) -o $@ 

rstinf:   rstinf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@ 

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

cleanall: clean cleanzip cleanlib cleanextra

clean: cleanobj cleanexe cleanvar cleantex cleanlahey cleanexport cleandoc \
	cleanifc

cleanifc:
	-rm -f ifc*

cleanobj:
	-rm -f *.o
	-rm -f *.stb

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
	-rm -f *.rev
	-rm -f CVSSTATUS CVSRAWSTATUS
	-rm -f COMMONF COMMONL
	-rm -f gmon.out

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

atoxi3d.o: param.h
basbathy.o: param.h evmain.h
baschk.o: param.h evmain.h
bascomp.o: param.h evmain.h
basinf.o: param.h evmain.h
bastreat.o: param.h evmain.h
bdist.o: param.h links.h
bedload.o: param.h evmain.h
bio3d.o: param.h
bio3d_util.o: param.h donata.h
concinf.o: param.h
flux.o: regflux.h param.h links.h ev.h
gotm3d.o: gotmturb.i
ht.o: param.h evmain.h
lagrange_back.o: param.h lagrange.h
lagrange_decay.o: param.h lagrange.h
lagrange_dif.o: param.h lagrange.h
lagrange_flux.o: param.h ev.h links.h lagrange.h
lagrange_in.o: param.h lagrange.h
lagrange_init.o: param.h
lagrange_inout.o: param.h ev.h lagrange.h
lagrange_larve.o: param.h lagrange.h
lagrange_main.o: param.h lagrange.h
lagrange_oil.o: param.h lagrange.h
lagrange_track.o: param.h lagrange.h
lagrange_util.o: ev.h
laplap.o: param.h evmain.h
lapwin.o: param.h evmain.h
lgrelab.o: param.h
lgrinf.o: param.h
loading.o: param.h
map_influence.o: param.h
matrix49.o: regflux.h
new36.o: close.h
new3di.o: testbndo.h param.h ev.h
newbcl.o: param.h
newchao.o: param.h
newchk.o: testbndo.h param.h ev.h
newcon.o: testbndo.h param.h ev.h
newcon_1.o: testbndo.h param.h ev.h
newconz.o: param.h
newcra.o: ev.h
newexpl.o: param.h ev.h
newfix.o: param.h bclfix.h
newini.o: param.h
newsigma.o: param.h sigma.h
newstab.o: param.h
newtra.o: param.h ev.h
newtvd.o: param.h ev.h tvd.h
nos_elab.o: param.h
nosaver.o: param.h
nosaver_basin.o: param.h basin.h
nosaver_time.o: param.h
nosaver_vert.o: param.h basin.h
nosaverf.o: param.h
nosdiff.o: param.h
nosexport.o: param.h
nosextr_custom.o: param.h
nosextr_gis.o: param.h
nosextr_nodes.o: param.h
nosextr_records.o: param.h
nosinf.o: param.h
nosinit.o: param.h
nosintp.o: param.h
nosmaniq.o: param.h
nosresidence.o: param.h ev.h basin.h evmain.h
nossplit.o: param.h
ousaver.o: param.h
ousextr_nodes.o: param.h
ousextr_records.o: param.h
ousinf.o: param.h evmain.h
ousintp.o: param.h
ousvdiff.o: param.h
outinf.o: param.h
outintp.o: param.h
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
subcus.o: param.h
subdep.o: ev.h
subdif.o: param.h ev.h
subdry.o: ev.h
subdwq.o: param.h
subele.o: param.h ev.h
subexta.o: modules.h
subflxa.o: testbndo.h modules.h links.h ev.h
subfvl.o: param.h
subgotm.o: param.h
sublin.o: links.h
sublnka.o: ev.h links.h
sublnks.o: ev.h links.h
sublnku.o: links.h
sublpl.o: ev.h
subn11.o: param.h ev.h
subn35.o: param.h chezy.h
subnetcdf_admin.o: param.h
subnetcdf_hydro.o: netcdf.inc
subnev.o: param.h ev.h
subnsh.o: modules.h
subousa.o: param.h
subpar1.o: subpar.h
subpard.o: param.h common.h
subqfxf.o: subqfx.h
subreg.o: ev.h
subres.o: param.h
subrst.o: param.h
subspk.o: param.h common.h
subuti.o: ev.h links.h
subvola.o: modules.h
subwat.o: links.h ev.h
subwaves.o: param.h
subwwm.o: param.h
tideg.o: param.h
vp.o: param.h
weutro.o: donata.h weutro.h
weutro_sedim.o: donata.h weutro.h
wininf.o: param.h
zinit.o: param.h

