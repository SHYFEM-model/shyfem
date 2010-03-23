
#--------------------------------------------------------------

TEST=rain
GRD=canale_vert

BASEDIR = $(HOME)/fem/femregres
GRDDIR = $(BASEDIR)/grds
DATADIR = $(BASEDIR)/skel/data/venlag

LINKS = ht splitflx splitext basinf flxinf
DATA = 

FEMDIR = $(HOME)/fem
BINDIR = $(FEMDIR)/fem3d

#--------------------------------------------------------------

default:

all: pre run post

pre: links grd bas data memory

#----------------------

links: $(LINKS)

$(LINKS):
	ln -fs $(BINDIR)/$@

grd:
	ln -fs $(GRDDIR)/$(GRD).grd

data: $(DATA)

$(DATA):
	ln -fs $(DATADIR)/$@

#----------------------

run:
	./ht < $(TEST).str

#----------------------

bas: grd
	vpgrd $(GRD)

memory:
	memory -s $(TEST) -b $(GRD)

post: memory
	./regress.sh

#----------------------

clean: cleanvar

cleanvar:
	-rm -f [uvzm].[0-9]*
	-rm -f errout.dat
	-rm -f *.bak *.tmp
	-rm -f *.ps
	-rm -f bndo_info.dat rkdiff.nos
	-rm -f PASSED

cleanbas:
	-rm -f *.bas

cleanlinks:
	-rm -f *.grd
	-rm -f $(LINKS)
	-rm -f $(DATA)

cleanall: clean cleanlinks cleanbas
	-rm -f .memory
	-rm -f fort.*
	-rm -f $(TEST).log $(TEST).inf
	-rm -f $(TEST).ext $(TEST).flx
	-rm -f $(TEST).cc1 $(TEST).ccc

#--------------------------------------------------------------

