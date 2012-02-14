# 
# general makefile for dir fem
#
#---------------------------------------------------------------
#
# targets that should be in all Makefiles of SUBDIRS:
#
# fem all 
# clean cleanall 
# list
# save zip
#
# really necessary are only: fem clean cleanall
#
#---------------------------------------------------------------

DIR = fem
RULES_MAKE_VERSION = 0.0	# default if Rules.make cannot be read

#---------------------------------------------------------------

include ./Rules.make

#---------------------------------------------------------------

RULES_MAKE_EXPECTED = 1.0
RULES_MAKE_COMPATIBILITY = RULES_MAKE_OK
ifneq ($(RULES_MAKE_VERSION),"0.0")
  ifneq ($(RULES_MAKE_VERSION),$(RULES_MAKE_EXPECTED))
    RULES_MAKE_COMPATIBILITY = RULES_MAKE_NOT_COMPATIBLE
  endif
endif

#---------------------------------------------------------------

FEMDIR    = .
DIRLIB    = $(FEMDIR)/femlib
FEMSRC    = $(FEMDIR)/fem3d
FEMBIN    = $(FEMDIR)/fembin
TMPDIR    = $(HOME)/fem/tmp

SUBDIRS   = `ls -dF * | grep  '/' | sed -e 's/\///'`
FEMLIBS   = femcheck post hcbs
FEMC      = grid mesh
FEMPROG   = fem3d femplot femadj femspline
FEMUTIL   = femregres femdoc fembin femlib femanim
FEMOPT    = femgotm femersem
FEMEXTRA  = 
PARAMDIRS = fem3d femplot femadj #femspline

SPECIAL   = Makefile Rules.make param.h README LASTTAR CHANGES
SPECIAL   = Makefile Rules.make param.h \
		BUG COMMIT FAQ LASTTAR LOG README VERSION

VERSION = `head -1 $(FEMDIR)/VERSION | sed -e 's/  */ /g' | cut -f4 -d" "`
VERSNAME = `head -1 $(FEMDIR)/VERSION | sed -e 's/  */ /g' | cut -f4 -d" " | sed -e 's/VERS_//'`
DATE = `date "+%F"`
TARNAME  = $(VERSION)
BETANAME = $(VERSNAME)_beta_$(DATE)
TARDIR   = $(TMPDIR)/fem_$(TARNAME)
BETADIR  = $(TMPDIR)/shyfem-$(BETANAME)

#---------------------------------------------------------------

FEMTOTS   = $(FEMLIBS) $(FEMC) $(FEMPROG) $(FEMUTIL) $(FEMOPT)

ifeq ($(GOTM),true)
  SUBDIRS += 
  FEMEXTRA += femgotm
endif

ifeq ($(ECOLOGICAL),ERSEM)
  SUBDIRS += femersem/src
  FEMEXTRA += femersem/src
endif

FEMDIRS   = $(FEMLIBS) $(FEMEXTRA) $(FEMC) $(FEMPROG) $(FEMUTIL)

#---------------------------------------------------------------

default: fem

.IGNORE: clean
.PHONY: publish version

all: fem doc

fem: checkv directories links
	$(FEMBIN)/recursivemake $@ $(FEMDIRS)

doc:
	cd femdoc; make doc

dirs:
	@echo "listing subdirectories (first level)"
	@echo $(SUBDIRS)

dirf:
	@echo "listing femdirectories (first level)"
	@echo $(FEMDIRS)

list:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

diff:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

depend:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

param:
	$(FEMBIN)/recursivemake $@ $(PARAMDIRS)

#---------------------------------------------------------------

clean: cleanlocal
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanall: cleanlocal
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanlocal:
	rm -f fem.tar fem.tar.gz
	rm -f changed_zip.zip
	rm -f *~
	rm -f *.tmp *.bak

cleandiff:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanbck:
	-rm -rf *.bck

#---------------------------------------------------------------

directories:
	-mkdir -p tmp

directories_old:
	-mkdir -p $(HOME)/lib

check:
	$(FEMBIN)/recursivemake $@ femcheck

install: checkv
	$(FEMBIN)/shyfem_install.sh

install_hard: checkv
	$(FEMBIN)/shyfem_install_hard.sh

install_hard_reset: checkv
	$(FEMBIN)/shyfem_install_hard.sh -reset

links:
	-rm -f bin lib
	-ln -sf fembin bin
	-ln -sf femlib lib

tar:
	@echo "Please use targets vers or beta"

beta: cleanall
	date > LASTTAR
	echo "$(BETANAME)"
	rm -rf $(TMPDIR)/fem.tar $(TMPDIR)/fem_VERS_* $(TMPDIR)/shyfem-*
	mkdir -p $(BETADIR)
	cp -al $(FEMTOTS) $(SPECIAL) $(BETADIR)
	cd $(BETADIR); ./fembin/shyfem_beta.sh
	cd $(TMPDIR); tar cvf fem.tar shyfem-*
	mv -f $(TMPDIR)/fem.tar .
	gzip -f fem.tar
	rm -rf $(TMPDIR)/fem.tar $(BETADIR)
	mv fem.tar.gz shyfem-$(BETANAME).tar.gz

vers: 
	@echo "please do not use this command - use gittar instead"

versold: cleanall
	date > LASTTAR
	rm -rf $(TMPDIR)/fem.tar $(TMPDIR)/fem_VERS_*; mkdir -p $(TARDIR)
	cp -al $(FEMTOTS) $(SPECIAL) $(TARDIR)
	cd $(TMPDIR); tar cvf fem.tar fem_VERS_*
	mv -f $(TMPDIR)/fem.tar .
	gzip -f fem.tar
	rm -rf $(TMPDIR)/fem.tar $(TARDIR)
	mv fem.tar.gz shyfem_$(TARNAME).tar.gz

publish:
	@fembin/fempub.sh shyfem_$(TARNAME).tar.gz

cvschange:
	 $(FEMBIN)/iterate1dir "$@ -d -t" $(FEMDIRS)

revision:
	 $(FEMBIN)/revision_last

version:
	@echo $(VERSION)

changed: modified
modified:
	@find . -newer VERSION -type f

changed_zip:
	zip changed_zip.zip `find . -newer VERSION -type f`

#---------------------------------------------------------------

checkv: $(RULES_MAKE_COMPATIBILITY)

RULES_MAKE_OK:

RULES_MAKE_NOT_COMPATIBLE:
	@echo "*** incompatible version of Rules.make file"
	@echo "provided: $(RULES_MAKE_VERSION)"
	@echo "expected: $(RULES_MAKE_EXPECTED)"
	@echo "please adapt your Rules.make file to the latest version"
	@exit 1

#---------------------------------------------------------------

