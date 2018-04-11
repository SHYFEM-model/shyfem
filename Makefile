# 
# general makefile for dir fem
#
#---------------------------------------------------------------
#
# targets that should be in all Makefiles of SUBDIRS:
#
# clean cleanall 
#
# really necessary are only: fem clean cleanall
#
#---------------------------------------------------------------

DIR = fem
RULES_MAKE_VERSION = 0.0	# default if Rules.make cannot be read
FEMDIR    = .
FEMBIN    = $(FEMDIR)/fembin

#---------------------------------------------------------------

include ./Rules.make

#---------------------------------------------------------------

RULES_MAKE_EXPECTED = 1.6
RULES_MAKE_COMPATIBILITY = RULES_MAKE_OK
ifneq ($(RULES_MAKE_VERSION),"0.0")
  ifneq ($(RULES_MAKE_VERSION),$(RULES_MAKE_EXPECTED))
    RULES_MAKE_COMPATIBILITY = RULES_MAKE_NOT_COMPATIBLE
  endif
endif

#---------------------------------------------------------------

ifndef ($(SHYFEMDIR))
        SHYFEMDIR := $(HOME)/shyfem
endif
export SHYFEMDIR

#---------------------------------------------------------------

FEMDIR    = .
DIRLIB    = $(FEMDIR)/femlib
FEMSRC    = $(FEMDIR)/fem3d
FEMBIN    = $(FEMDIR)/fembin
TMPDIR    = $(HOME)/fem/tmp
ACTFEMDIR = `pwd`

REGRESSDIR = femregress

SUBDIRS   = `ls -dF * | grep  '/' | sed -e 's/\///'`
FEMLIBS   = femcheck post hcbs
FEMC      = grid mesh
FEMPROG   = fem3d femplot femadj femspline
FEMUTIL   = $(REGRESSDIR) femdoc fembin femlib femanim
FEMOPT    = femgotm femersem
FEMEXTRA  = 
PARAMDIRS = fem3d femplot femadj #femspline

SPECIAL   = Makefile Rules.make param.h README LASTTAR CHANGES
SPECIAL   = Makefile Rules.make param.h \
		BUG COMMIT FAQ LASTTAR LOG README VERSION

VERSION = `head -1 $(FEMDIR)/VERSION | sed -e 's/  */ /g' | cut -f4 -d" "`
COMMIT = `head -1 $(FEMDIR)/VERSION | sed -e 's/  */ /g' | cut -f5 -d" "`
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
.PHONY: publish version stable

#---------------------------------------------------------------
# compiling and recursive targets
#---------------------------------------------------------------

default:
	@echo "   shyfem - version $(VERSION) $(COMMIT)"
	@echo '   run "make help" for more information'
	@echo '   if you are new to shyfem run "make first_time"'

all: fem doc
	@cd fem3d; make compatibility

fem: checkv directories links test_executable
	@$(FEMBIN)/recursivemake $@ $(FEMDIRS)
	@femcheck/check_compilation.sh -quiet

docs: doc
doc:
	@cd femdoc; make doc

dirs:
	@echo "listing subdirectories (first level)"
	@echo $(SUBDIRS)

dirf:
	@echo "listing femdirectories (first level)"
	@echo $(FEMDIRS)

list:
	@$(FEMBIN)/recursivemake $@ $(FEMDIRS)

depend:
	@$(FEMBIN)/recursivemake $@ $(FEMDIRS)

param:
	@$(FEMBIN)/recursivemake $@ $(PARAMDIRS)

directories:
	@-mkdir -p tmp
	@-mkdir -p femlib/mod
	@if [ ! -f ./tmp/Makefile ]; then cp ./femdummy/Makefile ./tmp; fi
	@if [ ! -f ./arc/Makefile ]; then cp ./femdummy/Makefile ./arc; fi

links:
	@-rm -f bin lib
	@-ln -sf fembin bin
	@-ln -sf femlib lib
	@if [ ! -d ./femregress ]; then ln -fs femdummy femregress; fi

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

cleanall: cleanlocal cleanregress
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleandist: cleanall

cleanregress:
	if [ -d $(REGRESSDIR) ]; then cd $(REGRESSDIR); make cleanall; fi
	
cleanbck:
	-rm -rf *.bck

#---------------------------------------------------------------
# public commands
#---------------------------------------------------------------

help:
	@echo "help                this screen"
	@echo "help_dev            more help for developers"
	@echo "first_time          what to do for the first time"
	@echo "install             installs the model"
	@echo "configure           configures the model"
	@echo "fem                 compiles everything"
	@echo "doc                 makes documentation"
	@echo "all                 compiles everything and makes documentation"
	@echo "clean, cleanall     cleans installation from tmp and obj files"
	@echo "version             returns version of this distribution"
	@echo "info                general options from Rules.make"
	@echo "check_software      checks needed software"
	@echo "check_compilation   checks if all executables have compiled"
	@echo "changed             finds files changed after installation"
	@echo "changed_zip         zips files changed after installation"

first_time:
	@echo 'Recommended use if you see shyfem for the first time:'
	@echo '   make help            gives overview of possible commands'
	@echo 'Commands to run only for the first time:'
	@echo '   make check_software  checks availability of software'
	@echo '   make configure       configures the use of the model'
	@echo '   make install         installs the model'
	@echo '   make fem             compiles the model'
	@echo 'Commands to run everytime you change something in the model:'
	@echo '   make cleanall        cleans directories'
	@echo '   make fem             compiles the model'

version:
	@echo $(VERSION) $(COMMIT)

info: version
	@echo "general:"
	@echo "  version           = $(VERSNAME)"
	@echo "  DISTRIBUTION_TYPE = $(DISTRIBUTION_TYPE)"
	@echo "  SHYFEM directory  = $(SHYFEMDIR)"
	@echo "macros:"
	@echo "  COMPILER     = $(COMPILER)"
	@echo "  PARALLEL_OMP = $(PARALLEL_OMP)"
	@echo "  PARALLEL_MPI = $(PARALLEL_MPI)"
	@echo "  SOLVER       = $(SOLVER)"
	@echo "  NETCDF       = $(NETCDF)"
	@echo "  GOTM         = $(GOTM)"
	@echo "  ECOLOGICAL   = $(ECOLOGICAL)"
	@echo "parameters:"
	@echo "  NKNDIM = $(NKNDIM)"
	@echo "  NELDIM = $(NELDIM)"
	@echo "  NLVDIM = $(NLVDIM)"
	@echo "compiler:"
	@echo "  FORTRAN  = $(F77)"
	@echo "  CC       = $(CC)"
	@echo "compiler options:"
	@echo "  PROFILE  = $(PROFILE)"
	@echo "  DEBUG    = $(DEBUG)"
	@echo "  OPTIMIZE = $(OPTIMIZE)"
	@echo "  WARNING  = $(WARNING)"
	@echo "compiler version:"
	@make compiler_version

check: check_software
check_software:
	@cd femcheck; ./check_software.sh

configure:
	@cd femcheck; ./configure.sh

check_compilation:
	@femcheck/check_compilation.sh

modified: changed
changed:
	@femcheck/find_changed.sh

changed_zip:
	zip changed_zip.zip `femcheck/find_changed.sh`
	@echo "changed files have bin zipped to changed_zip.zip"

advance_time:
	touch -d "`date -R -r VERSION` + 5 seconds" VERSION

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

ggu_help: help_ggu
help_ggu:
	@echo "rules_ggu_save     saves my Rules.make file"
	@echo "rules_ggu_restore  restores my Rules.make file"
	@echo "nemon              set special treatment nemunas server on"
	@echo "nemoff             set special treatment nemunas server off"
	@echo "rules_nemunas      copy Rules.make for nemunas server"
	@echo "rules_lagoon       copy Rules.make for lagoon"
	@echo "rules_carbonium    copy Rules.make for carbonium"
	@echo "git_nemunas        enable git for push on nemunas server"

help_dev:
	@echo "help_dev           this screen"
	@echo "test_compile       compiles model with different configs"
	@#echo "test_stable        compiles stable model with different configs"
	@echo "regress            runs regression tests"
	@echo "check_var          does various checks on distribution"
	@#echo "stable             makes stable distribution of last version"
	@echo "compiler_version   info on compiler"
	@echo "last_commit        name of last commit"
	@echo "dist               prepares distribution (Rules.make)"
	@echo "rules_save         copies back last saved Rules.make file"
	@echo "rules_dist         substitutes Rules.make with Rules.dist file"
	@echo "rules_new          copies Rules.make file to Rules.dist"
	@echo "rules_diff         difference between Rules.make and Rules.dist"
	@echo "advance_time       advances modification time of VERSION"
	@echo "make_executable    makes scripts executable"

test_compile:
	@femcheck/test_compile.sh

test_stable:
	@femcheck/test_stable.sh

check_var:
	@femcheck/check_var.sh

regress:
	if [ -d $(REGRESSDIR) ]; then cd $(REGRESSDIR); make regress; fi

revision:
	 $(FEMBIN)/revision_last.sh

rules_ggu_save:
	cp -f ./Rules.make arc/rules/Rules.ggu
	cp -f femcheck/Rules.dist ./Rules.make

rules_ggu_restore:
	cp -f arc/rules/Rules.ggu ./Rules.make

rules_dist:
	cp -f femcheck/Rules.dist ./Rules.make

rules_new:
	cp -f ./Rules.make femcheck/Rules.dist

rules_diff:
	@-diff femcheck/Rules.dist ./Rules.make || true

dist: cleandist
	mkdir -p arc/rules
	mv --backup=numbered ./Rules.make arc/rules/Rules.save
	cp -f femcheck/Rules.dist ./Rules.make
	make doc; make clean

stable:
	@stable/make_stable.sh $(FEMDIR)/arc/shyfem-$(VERSNAME).tar.gz

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

compiler_version:
	$(F77) $(FINFOFLAGS)
	$(CC) $(CINFOFLAGS)

last_commit:
	@gittags | tail -1

shyfemdir:
	@echo "shyfemdir: $(SHYFEMDIR)"

#---------------------------------------------------------------
# special targets for ggu
#---------------------------------------------------------------

nemon:
	fem3d/bin/nemunas_adjust.sh -nemunas

nemoff:
	fem3d/bin/nemunas_adjust.sh -original

rules_nemunas:
	cp -f arc/rules/Rules.nemunas ./Rules.make

rules_lagoon:
	cp -f arc/rules/Rules.lagoon ./Rules.make

rules_carbonium:
	cp -f arc/rules/Rules.carbonium ./Rules.make

git_nemunas:
	. arc/rules/nemunas-git.sh

#---------------------------------------------------------------
# check if routines are executable
#---------------------------------------------------------------

test_executable:
	@if [ ! -x fembin/make_executable.sh ]; then make make_executable; fi

make_executable:
	@echo "making programs executable"
	chmod +x fembin/*
	fembin/make_executable.sh

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

