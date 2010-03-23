#
# $Id: Makefile,v 1.49 2007-06-08 11:22:00 georg Exp $
#

#------------------------------------------------------------

DIR	= fem3d

include ../Rules.make

#------------------------------------------------------------

LIBFEM     = fem
LIBGOTM    = gotm
LIBNETCDF  = netcdf

DIRLIB_NETCDF = $(NETCDFDIR)/lib

LIBF_FEM    = $(DIRLIB)/lib$(LIBFEM).a
LIBF_GOTM   = $(DIRLIB)/lib$(LIBGOTM).a
LIBF_NETCDF = $(DIRLIB_NETCDF)/lib$(LIBNETCDF).a

LIBLS      = -L$(DIRLIB)
LIBFS      = 
LIBGS      =

INCLUDE    =

#------------------------------------------------------------

NetCDF_OBJ = subnetcdf_dummy.o
ifeq ($(NetCDF),true)
  NetCDF_OBJ = subnetcdf_admin.o subnetcdf_hydro.o
  LIBLS      = $(LIBLS) -L$(DIRLIB_NETCDF)
  LIBGS      = $(LIBGS) -l$(LIBNETCDF)
  LIBFS      = $(LIBFS) $(LIBF_NETCDF)
  INCLUDE    = $(INCLUDE) -I$(NETCDFDIR)/include
endif

GOTM_OBJ = gotm3d.o
ifeq ($(GOTM),true)
  GOTM_OBJ =
  LIBGS      = $(LIBGS) -l$(LIBGOTM)
  LIBFS      = $(LIBFS) $(LIBF_GOTM)
endif

#------------------------------------------------------------

SPECIAL =	Makefile apnpar.str \
		COMMIT README TESTLOG TODO VERSION WARNINGS \
		ACTUAL LASTZIP LOG

EXES	= ht hp vp \
	  ousinf whdf whdfs \
	  zlevel ttt readext zinit wincrea outinf \
	  basinf flxinf wininf nosinf nossplit nosextr nosintp \
	  nosaver nosaverf nosmaniq nosextract \
	  splitext outintp volinf extinf bastreat basbathy baschk \
	  nosdiff nostransf mkgeom lagint splitflx \
	  laplap lapwin nosinit ousintp nosexport \
	  bascomp ousextr rstinf nosresidence ousaver ousvdiff \
	  nosextr_inter

#------------------------------------------------------------

# 3D objects for 3D model

NEWOBJS = \
	  new3di.o newpr.o \
	  newrog.o newcra.o \
	  newcon.o newbcl.o \
	  newini.o newchk.o newtra.o \
	  newbcl0.o \
	  newchao.o newparam.o \
	  newdum.o

# objects for library

SOBJS   = \
	  subnsv.o subsst.o \
          subcmv.o subdry.o subpar.o subssv.o \
          sublst.o subsys.o \
          subcst.o submir.o subros.o \
          subn11.o subsrt.o \
          subdef.o subext.o subssc.o \
          subbas.o subnev.o subdim.o \
          subnls.o subsse.o subwat.o \
          subfil.o subn35.o subout.o subous.o subssg.o \
          subfnm.o subnsh.o subssm.o subzzz.o \
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
	  debug.o subhisto.o subnosa.o

# general objects for 3D model ht

HTOBJS = \
	subver.o \
	subcst.o subsys.o subnsh.o \
	subbas.o subdim.o subnev.o \
	subnsu.o subn35.o subdep.o \
	subdef.o subfil.o subnls.o \
	subpar.o subfnm.o \
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
	$(QFLUX) \
	subgotm.o $(GOTM_OBJ) \
	$(NETCDF_OBJ) \
	subbndo.o subbnds.o \
	weutro.o bio3d.o new36.o intp.o subssf.o \
	bdist.o links.o sublnka.o sublnks.o sublnku.o \
	sedi3d.o sedtrans05.o subwaves.o \
	subhisto.o \
	lagrange.o lagrange_in.o lagrback.o lagrange_dif.o lagrange_cont.o \
	weutro_sedim.o weutro_light.o subdts.o loading.o \
	subrst.o subnosa.o subdif.o subreg.o \
	atoxi3d.o atoxi.o subwwm.o kepsd.o keps_util.o \
	subinv.o newexpl.o subrnd.o

# objects for pre-processor vp

VPOBJS = \
        subros.o subcmv.o \
        subdef.o subbas.o subdim.o \
        subfil.o \
        subpar.o subfnm.o \
        subsrt.o subgrd.o \
        subsss.o subscn.o subnsu.o \
	subapn.o subnsa.o subnls.o subsys.o subdum.o subbnd.o

# objects for output routines

OUTOBJS = subpar.o subfnm.o subnls.o subsss.o subscn.o \
	  subfil.o subdef.o subnsu.o subssv.o \
	  subbas.o subsys.o subdim.o \
	  subapn.o subnsa.o subdum.o \
	  subbnd.o subdts.o subdep.o \
	  debug.o subreg.o subnev.o subrst.o

LAPLACE = sublpl.o submir.o subnos.o links.o sublnku.o \
	  subwin.o subnosa.o

# heat flux module

QFLUX = subqfxt.o \
	subqfxf.o \
	subqfxm1.o subqfxm2.o \
	subqfxu1.o subqfxu2.o subqfxu3.o subqfxu4.o

# objects for interpolation

INTP    = subzzz.o subreg.o

# objects for lagrangian interpolation

LAGINTP = submir.o

#------------------------------------------------

default:	fem

all:	fem
fem:	libfem main info util 
main:	ht vp
info:	outinf flxinf volinf basinf nosinf wininf extinf
util:	nosextr nossplit splitext ousinf splitflx laplap lapwin \
	nosinit bastreat basbathy nosextract nosextr_inter

#--------------------------------------------- EXES

ht:  ht.o $(NEWOBJS) $(HTOBJS) $(LIBFGOTM)
	$(LINKER) $(LFLAGS) ht.o $(NEWOBJS) $(HTOBJS) $(LIBGS) -o $@

vp:  vp.o $(VPOBJS)
	$(LINKER) $(LFLAGS) vp.o $(VPOBJS) -o $@

#------------------------------------------------

zlevel:  zlevel.o 
	$(LINKER) $(LFLAGS) zlevel.o -o $@

lagint:  lagint.o sublag.o $(LAGINTP) $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o sublag.o $(LAGINTP) $(OUTOBJS) -o $@

readext:  readext.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) readext.o subext.o $(OUTOBJS) -o $@

splitext:  splitext.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) splitext.o subext.o $(OUTOBJS) -o $@

splitflx:  splitflx.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) splitflx.o $(OUTOBJS) -o $@

extinf:  extinf.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) extinf.o subext.o $(OUTOBJS) -o $@

flxinf:  flxinf.o subext.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) flxinf.o subext.o $(OUTOBJS) -o $@

volinf:  volinf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) volinf.o $(OUTOBJS) -o $@

zinit:  zinit.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) zinit.o $(OUTOBJS) -o $@

wincrea:   wincrea.o
	$(LINKER) $(LFLAGS) wincrea.o -o $@

wininf:   wininf.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) wininf.o $(OUTOBJS) -o $@

mkgeom:   mkgeom.o
	$(LINKER) $(LFLAGS) mkgeom.o -o $@

cubic:   cubic.o $(OUTOBJS) subsse.o subssm.o
	$(LINKER) $(LFLAGS) cubic.o $(OUTOBJS) subsse.o subssm.o -o $@

#------------------------------------------------

basinf:  basinf.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o $(OUTOBJS) -o $@

baschk:  baschk.o subnev.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnev.o $(OUTOBJS) -o $@

bastreat:  bastreat.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

bascomp:  bascomp.o subgeo.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o $(OUTOBJS) -o $@

basbathy:  basbathy.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subgeo.o subgrd.o subsrt.o $(OUTOBJS) -o $@

ousinf:  ousinf.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousextr:  ousextr.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousintp:  ousintp.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousaver:  ousaver.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

ousvdiff:  ousvdiff.o subous.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subous.o $(OUTOBJS) -o $@

nosinf:  nosinf.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

concinf:  concinf.o subnos.o $(OUTOBJS)
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosinf2:  nosinf2.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosextr:  nosextr.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosextr_inter:  nosextr_inter.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nossplit:  nossplit.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosaver:  nosaver.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosaverf:  nosaverf.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosresidence:  nosresidence.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosmaniq:  nosmaniq.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosextract:  nosextract.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nostransf:  nostransf.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosdiff:  nosdiff.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosinit:  nosinit.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosexport:  nosexport.o subnos.o $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(OUTOBJS) -o $@

nosintp:  nosintp.o subnos.o $(INTP) $(OUTOBJS) 
	$(LINKER) $(LFLAGS) $@.o subnos.o $(INTP) $(OUTOBJS) -o $@

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
	makedepend $(FSOURCE)		#this only deletes last part of Mfile
	bin/include.pl -make *.f *.h

#---------------------------------------------------

cleanall: clean cleanzip cleanlib cleanextra

clean: cleanobj cleanexe cleanvar cleantex cleanlahey cleanexport cleandoc \
	cleanifc

cleanifc:
	-rm -f ifc*

cleanobj:
	-rm -f *.o

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
	-rm -f .memory _memory
	-rm -f html/cvs*
	-rm -fr cvsdiff_tmp
	-rm -f CVSSTATUS CVSRAWSTATUS
	-rm -f COMMONF COMMONL
	-rm -f gmon.out

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleanlib:
	-rm -f $(LIBFFEM) $(DLIBFFEM)

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
	$(F77) $(FFOLD) -c $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)


# DO NOT DELETE THIS LINE -- make depend depends on it.

newbcl.o:  param.h 
subqfxf.o:  subqfx.h 
subwaves.o:  param.h 
sedtrans05.o:  sed_param.h 
nosinf.o:  param.h 
subres.o:  param.h 
subcus.o:  param.h 
nosmaniq.o:  param.h 
newchao.o:  param.h 
lagrange_in.o:  lagrange.h param.h 
basinf.o:  basinf.h 
atoxi3d.o:  param.h 
nosintp.o:  param.h 
bascomp.o:  basinf.h 
sedi3d.o:  param.h sed_param.h 
loading.o:  param.h 
nosaver.o:  param.h 
laplap.o:  basinf.h 
nosextract.o:  param.h 
subousa.o:  param.h 
bastreat.o:  basinf.h 
nosresidence.o:  param.h 
newcon.o:  param.h 
nosinit.o:  param.h 
gotm3d.o:  gotmturb.i 
concinf.o:  param.h 
subele.o:  param.h 
newchk.o:  param.h 
basbathy.o:  basinf.h 
subbndo.o:  subbndo.h 
newpr.o:  param.h 
subn11.o:  param.h 
subn35.o:  chezy.h param.h 
subvola.o:  modules.h 
baschk.o:  basinf.h 
ht.o:  param.h 
flux.o:  param.h regflux.h 
new3di.o:  param.h 
lagrange.o:  lagrange.h param.h 
lagrback.o:  lagrange.h param.h 
bdist.o:  param.h 
lagrange_dif.o:  lagrange.h param.h 
newtra.o:  param.h 
subexta.o:  modules.h 
weutro.o:  donata.h weutro.h 
subflxa.o:  modules.h 
subrst.o:  param.h 
matrix49.o:  regflux.h 
newini.o:  param.h 
weutro_sedim.o:  donata.h param.h weutro.h 
subgotm.o:  param.h 
nossplit.o:  param.h 
nosexport.o:  param.h 
nosextr.o:  param.h 
lapwin.o:  basinf.h 
subdif.o:  param.h 
subnsh.o:  modules.h 
newexpl.o: param.h
bio3d.o:  param.h 
outintp.o:  param.h 
subcst.o:  modules.h 
new36.o:  close.h 
subwwm.o: param.h
subnetcdf.o: param.h
Hydro_netcdfs_fem.o:
subnetcdf_dummy.o:
