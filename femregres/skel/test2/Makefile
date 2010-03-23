
#--------------------------------------------------------------

TEST=test2
GRD=venlag62

FEMDIR = $(HOME)/fem/fem3d
GRDDIR = $(DATADIR)
DATADIR = $(HOME)/fem/femregres/skel/data/venlag

LINKS = ht splitflx splitext basinf flxinf
DATA = wind.dat rain.dat ds2000.dat prev.dat

#--------------------------------------------------------------

default:

all: links grd bas data memory

#----------------------

links: $(LINKS)

$(LINKS):
	ln -fs $(FEMDIR)/$@

grd:
	ln -fs $(GRDDIR)/$(GRD).grd

data: $(DATA)

$(DATA):
	ln -fs $(DATADIR)/$@

#----------------------

bas: grd
	vpgrd $(GRD)

memory:
	memory -s $(TEST) -b $(GRD)

post:
	splitext < ~/bin/CR

#----------------------

clean: cleanlinks cleanbas
	-rm -f [uvzm].[0-9]*
	-rm -f .memory errout.dat

cleanbas:
	-rm -f *.bas

cleanlinks:
	-rm -f *.grd
	-rm -f $(LINKS)
	-rm -f $(DATA)

cleanall: clean
	-rm -f fort.*
	-rm -f $(TEST).log $(TEST).inf
	-rm -f $(TEST).ext $(TEST).flx
	-rm -f $(TEST).cc1 $(TEST).ccc

#--------------------------------------------------------------

