
DIR = examples

FEMDIR = ..
FEMBIN = ../fembin

SUBDIRS = mar_menor

#-----------------------------------------------------------

default: nothing

fem: nothing

nothing:
	@true

list:
	@pwd; ls

depend: nothing

clean:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanall:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

#-----------------------------------------------------------

