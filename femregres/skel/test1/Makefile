
#--------------------------------------------------------------

TEST=test1
GRD=venlag62

FEMDIR = $(HOME)/fem/fem3d
GRDDIR = $(DATADIR)
DATADIR = $(HOME)/fem/femregres/skel/data/venlag

LINKS = ht splitflx splitext basinf flxinf vp laplap

#--------------------------------------------------------------

default:

all: links grd bas memory

#----------------------

links: $(LINKS)

$(LINKS):
	ln -fs $(FEMDIR)/$@

grd:
	ln -fs $(GRDDIR)/$(GRD).grd

data:
	ln -fs $(DATADIR) data

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
	-rm -f apnstd_.str plot.ps

cleanbas:
	-rm -f *.bas

cleanlinks:
	-rm -f *.grd
	-rm -f $(LINKS)

cleanall: clean
	-rm -f fort.*
	-rm -f $(TEST).log $(TEST).inf
	-rm -f $(TEST).ext $(TEST).flx
	-rm -f $(TEST).nos
	-rm -f laplace.nos

#--------------------------------------------------------------

