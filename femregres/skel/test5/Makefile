
#--------------------------------------------------------------

TEST=test5
GRD=hakata22

FEMDIR = $(HOME)/fem/fem3d
GRDDIR = $(DATADIR)
DATADIR = $(HOME)/fem/femregres/skel/data/haka

LINKS = ht splitflx splitext basinf flxinf vp nosinf ousinf
DATA = 1998.qfx 1998.win arm4ex.mar 9899.salt 9899.temp

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

clean: cleanvar cleanlinks cleanbas

cleanvar:
	-rm -f [uvzm].[0-9]*
	-rm -f .memory errout.dat
	-rm -f *.tmp
	-rm -f ggg hhh

cleanlinks:
	-rm -f *.grd
	-rm -f $(LINKS)
	-rm -f $(DATA)

cleanbas:
	-rm -f *.bas

cleanall: clean cleansims

cleansims:
	-rm -f fort.*
	-rm -f *.nos
	-rm -f $(TEST).log $(TEST).inf
	-rm -f $(TEST).ext $(TEST).flx
	-rm -f $(TEST).out $(TEST).ouv
	-rm -f $(TEST).con $(TEST).hak
	-rm -f $(TEST).93  $(TEST).o3d
	-rm -f $(TEST).sss $(TEST).ttt $(TEST).tsa
	-rm -f $(TEST).ous 

#--------------------------------------------------------------

