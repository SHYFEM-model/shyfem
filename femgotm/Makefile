#$Id: Makefile,v 1.6 2005-07-06 14:07:17 kbk Exp $
#
# Makefile to build the the Turbulence library - libturbulence.a
#
# customized by Georg Umgiesser for SHYFEM (GGU)

FEMDIR = ..
FEMBIN = $(FEMDIR)/fembin

include ../Rules.make

LIB	= $(LIBDIR)/libturbulence$(buildtype).a

#------------------------------------------------------------

DIR     = femgotm
LIBGOTM  = gotm
LIB = $(DIRLIB)/lib$(LIBGOTM).a
LIBDIR = $(DIRLIB)
MODDIR = .

EXTRA_DOCSRC = 		\
tridiagonal.F90         \
diff_face.F90           \
fem_gotm_interface.F90	

EXTRA_MOD =		\
tridiagonal.o		\
eqstate.o		\
util.o

EXTRA_OBJ =		\
diff_face.o		\
fem_gotm_interface.o

EXTRA_SPECIAL =	cppdefs.h version.h gotmturb.nml Makefile.orig

#------------------------------------------------------------

DOCSRC	= \
$(EXTRA_DOCSRC)		\
turbulence.F90		\
production.F90          \
alpha_mnb.F90           \
r_ratio.F90             \
tkeeq.F90		\
q2over2eq.F90		\
kbeq.F90		\
lengthscaleeq.F90	\
dissipationeq.F90	\
genericeq.F90	        \
tkealgebraic.F90	\
kbalgebraic.F90	        \
algebraiclength.F90	\
epsbalgebraic.F90	\
variances.F90	        \
ispralength.F90		\
potentialml.F90	        \
cmue_a.F90              \
cmue_b.F90              \
cmue_c.F90              \
cmue_d.F90              \
cmue_ma.F90		\
cmue_sg.F90		\
cmue_rf.F90		\
compute_cpsi3.F90	\
compute_rist.F90	\
internal_wave.F90       \
fk_craig.F90            \
kpp.F90                 \
gotm_lib_version.F90

MOD   = \
$(EXTRA_MOD)		\
turbulence.o		\
kpp.o

OBJ   = \
$(EXTRA_OBJ)		\
gotm_lib_version.o	\
production.o            \
alpha_mnb.o             \
r_ratio.o               \
tkeeq.o			\
q2over2eq.o		\
tkealgebraic.o		\
kbeq.o	          	\
kbalgebraic.o		\
dissipationeq.o		\
genericeq.o		\
lengthscaleeq.o		\
algebraiclength.o	\
epsbalgebraic.o  	\
variances.o       	\
ispralength.o		\
potentialml.o		\
cmue_a.o		\
cmue_b.o		\
cmue_c.o		\
cmue_d.o		\
cmue_ma.o		\
cmue_rf.o		\
cmue_sg.o		\
compute_cpsi3.o		\
compute_rist.o		\
fk_craig.o		\
internal_wave.o

#---------------------------------------------------------------------

all: $(MOD) $(OBJ) $(LIB) nothing

$(LIB): $(MOD) $(OBJ)
	ar rvs $(LIB) *.o

$(OBJ): $(MOD)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/turbulence.tex 
	touch doc

clean:
	$(RM) -f $(LIB) $(MODDIR)/turbulence.{m,mod} $(MODDIR)/kpp.{m,mod} doc
	$(RM) -f *.mod *_genmod.f90
	$(RM) -f *.o

realclean: clean
	$(RM) -f *.o *~ *.stb

distclean: realclean

#------------------------------------------------------------

fem: all nothing

nothing:
	@true

list:
	pwd; ls

depend: nothing

cleanall: clean cleanmod realclean

cleanmod:
	$(RM) -f *.mod

.SUFFIXES:
.SUFFIXES: .F90

%.o: %.F90
	$(F95) $(FFLAGS) -c $<

#------------------------------------------------------------

#-----------------------------------------------------------------------
# Copyright by the GOTM-team under the GNU Public License - www.gnu.org
#-----------------------------------------------------------------------
