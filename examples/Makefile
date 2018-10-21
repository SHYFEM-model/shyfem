
DIR = shyfem_tutorial_hydro

BASEDIR = ..
GRID  = MM_coarse
GRID1  = MM_coarse_smooth

#-----------------------------------------------------------

default:

init:
	cd GRID; shypre $(GRID).grd
	cd GRID; shypre $(GRID1).grd

#-----------------------------------------------------------

clean:
	-rm -f rkdiff.nos dist.nos bndo_info.dat
	-rm -f *.tmp *.bak
	-rm -f *.[0123]d.[0-9]*
	-rm -f *.[23]d.txt
	-rm -f p.[1-9]*
	-rm -f out*.ps plot*.ps
	-rm -f out.fem out.ext out.flx
	-rm -f fort.*
	-rm -f errout.dat gmon.out
	-rm -f .memory

cleanall: clean cleansims cleanzip cleanfigs cleanbas

cleanbas:
	-rm -f GRID/*.bas
	-rm -f GRID/errout.dat

cleanfigs:
	-rm -f *.ps
	-rm -f *.eps
	-rm -f *.jpg

cleansims:
	-rm -f mm_hyd*.shy
	-rm -f mm_hyd*.inf
	-rm -f mm_hyd*.ext
	-rm -f mm_hyd*.flx
	-rm -f out.shy

cleanzip:
	-rm -f $(DIR).zip

#-----------------------------------------------------------

zip:
	#-zip -r $(DIR).zip *
	-zip  $(DIR).zip Makefile *.sh *.str GRID/* INPUT/*

#-----------------------------------------------------------

