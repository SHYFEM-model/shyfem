#$Id: Makefile,v 1.6 2005-07-06 14:07:17 kbk Exp $
#
# Makefile to build the the Turbulence library - libturbulence.a
#

FEMDIR = ..
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
${LIB}(tridiagonal.o)	\
${LIB}(eqstate.o)	\
${LIB}(util.o)

EXTRA_OBJ =		\
${LIB}(diff_face.o)	\
${LIB}(fem_gotm_interface.o)

EXTRA_SPECIAL =	cppdefs.h version.h gotmturb.nml Makefile.orig

#------------------------------------------------------------

DOCSRC	= \
${EXTRA_DOCSRC}		\
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
${EXTRA_MOD}		\
${LIB}(turbulence.o)	\
${LIB}(kpp.o)

OBJ   = \
${EXTRA_OBJ}			\
${LIB}(gotm_lib_version.o)	\
${LIB}(production.o)            \
${LIB}(alpha_mnb.o)             \
${LIB}(r_ratio.o)               \
${LIB}(tkeeq.o)			\
${LIB}(q2over2eq.o)		\
${LIB}(tkealgebraic.o)		\
${LIB}(kbeq.o)	          	\
${LIB}(kbalgebraic.o)		\
${LIB}(dissipationeq.o)		\
${LIB}(genericeq.o)		\
${LIB}(lengthscaleeq.o)		\
${LIB}(algebraiclength.o)	\
${LIB}(epsbalgebraic.o)  	\
${LIB}(variances.o)       	\
${LIB}(ispralength.o)		\
${LIB}(potentialml.o)		\
${LIB}(cmue_a.o)		\
${LIB}(cmue_b.o)		\
${LIB}(cmue_c.o)		\
${LIB}(cmue_d.o)		\
${LIB}(cmue_ma.o)		\
${LIB}(cmue_rf.o)		\
${LIB}(cmue_sg.o)		\
${LIB}(compute_cpsi3.o)		\
${LIB}(compute_rist.o)		\
${LIB}(fk_craig.o)		\
${LIB}(internal_wave.o)


all: $(MOD) ${OBJ}
	$(MOVE_MODULES_COMMAND)

$(OBJ): $(MOD)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/turbulence.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/turbulence.{m,mod} $(MODDIR)/kpp.{m,mod} doc
	$(RM) -f *.mod *_genmod.f90


realclean: clean
	$(RM) -f *.o *~ *.stb

distclean: realclean

#------------------------------------------------------------

fem: all

list:
	pwd; ls

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
