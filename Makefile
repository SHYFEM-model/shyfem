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

RULES_MAKE_EXPECTED = 1.1
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

.IGNORE: clean
.PHONY: publish version public

#---------------------------------------------------------------
# compiling and recursive targets
#---------------------------------------------------------------

default: fem

all: fem doc

fem: checkv directories links
	$(FEMBIN)/recursivemake $@ $(FEMDIRS)
	@femcheck/check_compilation.sh -quiet

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

directories:
	-mkdir -p tmp

links:
	-rm -f bin lib
	-ln -sf fembin bin
	-ln -sf femlib lib

#---------------------------------------------------------------
# cleaning
#---------------------------------------------------------------

cleanlocal:
	-rm -f fem.tar fem.tar.gz
	-rm -f changed_zip.zip
	-rm -f *~
	-rm -f *.tmp *.bak *.out
	-rm -f ggg hhh
	-rm -f errout.dat a.out plot.ps
	-rm -f .memory
	-rm -f CHECKLOG

clean: cleanlocal
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanall: cleanlocal
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleandist: cleanall

cleandiff:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanbck:
	-rm -rf *.bck

#---------------------------------------------------------------
# public commands
#---------------------------------------------------------------

version:
	@echo $(VERSION)

check: check_software
check_software:
	@cd femcheck; ./check_software.sh

check_compilation:
	@femcheck/check_compilation.sh

modified: changed
changed:
	@find . -newer VERSION -type f

changed_zip:
	zip changed_zip.zip `find . -newer VERSION -type f`

#--------------------------------------------------------
# installing
#--------------------------------------------------------

install: install_hard install_soft

install_soft: checkv
	$(FEMBIN)/shyfem_install_soft.sh

install_hard: checkv
	$(FEMBIN)/shyfem_install_hard.sh

uninstall: install_hard_reset install_soft_reset

install_soft_reset: checkv
	$(FEMBIN)/shyfem_install_soft.sh -reset

install_hard_reset: checkv
	$(FEMBIN)/shyfem_install_hard.sh -reset

#--------------------------------------------------------
# private and admin commands
#--------------------------------------------------------

test_compile:
	@femcheck/test_compile.sh

revision:
	 $(FEMBIN)/revision_last

ggu:
	cp -f rules/Rules.save ./Rules.make

dist: cleandist
	mv --backup=numbered ./Rules.make rules/Rules.save
	cp -f rules/Rules.skel ./Rules.make
	make doc; make clean

public:
	@public/make_public.sh shyfem-$(VERSNAME).tar.gz

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

publish:
	@fembin/fempub.sh shyfem_$(TARNAME).tar.gz

#---------------------------------------------------------------
# compatibility checks
#---------------------------------------------------------------

checkv: $(RULES_MAKE_COMPATIBILITY) $(RULES_MAKE_PARAMETERS)

RULES_MAKE_OK:

RULES_MAKE_NOT_COMPATIBLE:
	@echo "*** incompatible version of Rules.make file"
	@echo "provided: $(RULES_MAKE_VERSION)"
	@echo "expected: $(RULES_MAKE_EXPECTED)"
	@echo "please adapt your Rules.make file to the latest version"
	@exit 1

RULES_MAKE_PARAMETER_ERROR:
	@echo "*** incompatible parameters in Rules.make file"
	@echo "$(RULES_MAKE_MESSAGE)"
	@echo "please adapt your Rules.make file"
	@exit 1

#---------------------------------------------------------------
# end of Makefile
#---------------------------------------------------------------

