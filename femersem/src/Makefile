## Makefile to build the BFM standalone
##
## uses makedepf90 to build dependencies "on the fly"
## http://www.helsinki.fi/~eedelman/makedepf90.html

FEMDIR = ../..
FEMBIN    = $(FEMDIR)/fembin

include ../../Rules.make

ifeq ($(FORTRAN_COMPILER),GNU_G77)
  FORTRAN_COMPILER=G77
endif
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  FORTRAN_COMPILER=GFORTRAN
endif
ifeq ($(FORTRAN_COMPILER),INTEL)
  FORTRAN_COMPILER=IFORT
endif
ifeq ($(FORTRAN_COMPILER),PORTLAND)
  FORTRAN_COMPILER=PGF90
endif

include ./Rules.make

FC = $(F77)
buildtype = 
LIBDIR = ../../femlib

EXE	= bfm_standalone$(buildtype)
MAINOBJ	= ./standalone/standalone_main.o 
LIB	= $(LIBDIR)/libbfm$(buildtype).a
BFMSRC = ./BFM

BFM_MOD = \
	${LIB}(${BFMSRC}/General/ModuleGlobalMem.o)			\
	${LIB}(${BFMSRC}/General/ModuleConstants.o)			\
	${LIB}(${BFMSRC}/General/ModuleGlobFun.o)			\
	${LIB}(${BFMSRC}/General/ModuleMem.o)				\
	${LIB}(${BFMSRC}/General/ModuleParam.o)				\
	${LIB}(${BFMSRC}/General/ModuleInterface.o)			\
	${LIB}(${BFMSRC}/Light/ModuleLightAdaptation.o)			\
	${LIB}(${BFMSRC}/Light/ModulePhotoAvailableRadiation.o)		\
	${LIB}(${BFMSRC}/Oxygen/ModuleWindOxReaeration_3.o)		\
	${LIB}(${BFMSRC}/PelB/ModuleMesoZoo.o)				\
	${LIB}(${BFMSRC}/PelB/ModuleMicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelBac.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelChem.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelGlobal.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePhyto.o)				\
	${LIB}(${BFMSRC}/PelBen/ModuleSettling.o)			\
	${LIB}(${BFMSRC}/PelBen/ModuleControlBenPartNutrientBuffers.o)	\
	${LIB}(${BFMSRC}/Ben/ModuleBenOrganism.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBenBac.o)				\
	${LIB}(${BFMSRC}/Ben/ModuleFilterFeeder.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBioturbation.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAmmonium.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAnoxic.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenDenitriDepth.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNitrate.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutConstants.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutType.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutInterface.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutVariables.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenOxygen.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenPhosphate.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenQ1Transport.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenSilica.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn1.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn2.o)			 

	
BFM_OBJ = \
	${LIB}(${BFMSRC}/General/AllocateMem.o)				\
	${LIB}(${BFMSRC}/General/Ecology.o)				\
	${LIB}(${BFMSRC}/General/InitBoxParams.o)			\
	${LIB}(${BFMSRC}/General/Initialize.o)				\
	${LIB}(${BFMSRC}/General/InitTransportStateTypes.o)		\
	${LIB}(${BFMSRC}/General/eTq.o)					\
	${LIB}(${BFMSRC}/Light/LightAdaptation.o)			\
	${LIB}(${BFMSRC}/Light/PhotoAvailableRadiation.o)		\
	${LIB}(${BFMSRC}/Oxygen/WindOxReaeration_3.o)			\
	${LIB}(${BFMSRC}/Oxygen/CalcOxygenSaturation_3.o)		\
	${LIB}(${BFMSRC}/Oxygen/CalcSchmidtNumberOx.o)		\
	${LIB}(${BFMSRC}/PelB/CalcChlorophylla.o)			\
	${LIB}(${BFMSRC}/PelB/CalcVerticalExtinction.o)			\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/MesoZoo.o)				\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/PelBac.o)					\
	${LIB}(${BFMSRC}/PelB/PelChem.o)				\
	${LIB}(${BFMSRC}/PelB/PelGlobal.o)				\
	${LIB}(${BFMSRC}/PelB/PelagicSystem.o)				\
	${LIB}(${BFMSRC}/PelB/Phyto.o)					\
	${LIB}(${BFMSRC}/Ben/ResetTotMassVar.o)				\
	${LIB}(${BFMSRC}/PelBen/BentoPelCoup.o)				\
	${LIB}(${BFMSRC}/PelBen/PelForcingForBen.o)			\
	${LIB}(${BFMSRC}/PelBen/Sedimentation.o)			\
	${LIB}(${BFMSRC}/PelBen/Settling.o)				\
	${LIB}(${BFMSRC}/PelBen/ControlBenPartNutrientBuffers.o)	\
	${LIB}(${BFMSRC}/Ben/BenBac.o)					\
	${LIB}(${BFMSRC}/Ben/BenOrganism.o)				\
	${LIB}(${BFMSRC}/Ben/BenthicSystem.o)				\
	${LIB}(${BFMSRC}/Ben/Bioturbation.o)				\
	${LIB}(${BFMSRC}/Ben/FilterFeeder.o)				\
	${LIB}(${BFMSRC}/Ben/ResetTotMassVar.o)				\
	${LIB}(${BFMSRC}/Ben/BenCheckMassConservation.o)		\
	${LIB}(${BFMSRC}/Bennut/BenAmmonium.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenAnoxic.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenDenitriDepth.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenNitrate.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenNitrogenShifting.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenOxygen.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenPhosphate.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenQ1Transport.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenSilica.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient2.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient3.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn1.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn2.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateFromSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateShift.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateTau.o) 			\
	${LIB}(${BFMSRC}/Bennut/CompleteSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/InitializeSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/DefineSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/GetInfoFromSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/PrintSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/bess_exp.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessi0.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessi1.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessk0.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessk1.o) 				\
	${LIB}(${BFMSRC}/Bennut/calculate_equation.o) 			\
	${LIB}(${BFMSRC}/Bennut/calculatelayer.o) 			\
	${LIB}(${BFMSRC}/Bennut/funcalc.o) 				\
	${LIB}(${BFMSRC}/Bennut/input_para.o) 				\
	${LIB}(${BFMSRC}/Bennut/kfind.o) 				\
	${LIB}(${BFMSRC}/Bennut/lubksb.o) 				\
	${LIB}(${BFMSRC}/Bennut/ludcmp.o) 				\
	${LIB}(${BFMSRC}/Bennut/manage_coeff.o) 			\
	${LIB}(${BFMSRC}/Bennut/noutput.o) 				\
	${LIB}(${BFMSRC}/Bennut/qgaus_exp.o) 				\
	${LIB}(${BFMSRC}/Bennut/re_store.o) 				\
	${LIB}(${BFMSRC}/Bennut/set_max_sing.o) 			\
	${LIB}(${BFMSRC}/Bennut/svbksb.o) 				\
	${LIB}(${BFMSRC}/Bennut/svdcmp.o) 				\
	${LIB}(${BFMSRC}/Bennut/transfer.o) 				 

APIBFM_MOD =	\
	${LIB}(./share/api_bfm.o)			\
	${LIB}(./share/netcdf_bfm.o)			\
	${LIB}(./share/string_functions.o)		\
	${LIB}(./share/init_var_bfm.o)			\
	${LIB}(./share/init_cnps.o)			\
	${LIB}(./share/calcmean_bfm.o)			\
	${LIB}(${BFMSRC}/General/set_var_info_bfm.o)	

STDALONE_MOD =	\
	${LIB}(./standalone/time.o)		\
	${LIB}(./standalone/standalone.o)	\
	${LIB}(./standalone/integration-efw.o)	\
	${LIB}(./standalone/integration-lf.o)	\
	${LIB}(./standalone/integration-RK2.o)	\
	${LIB}(./standalone/D3toD1.o)		\
	${LIB}(./standalone/D2toD1.o)		\
	${LIB}(./standalone/ResetFluxes.o)	\
	${LIB}(./standalone/GetDelta.o)		

FEM_OBJ = \
	${LIB}(./share/femtobfm.o)                      


fem: bfm
bfm: create_dirs all
all: ${BFM_MOD} ${APIBFM_MOD} ${STDALONE_MOD} ${FEM_OBJ} ${BFM_OBJ} 
	$(MOVE_MODULES_COMMAND)

create_dirs: ../lib/$(FORTRAN_COMPILER) ../modules/$(FORTRAN_COMPILER)
../lib/$(FORTRAN_COMPILER):
	mkdir -p ../lib/$(FORTRAN_COMPILER)
../modules/$(FORTRAN_COMPILER):
	mkdir -p ../modules/$(FORTRAN_COMPILER)

#all: ${BFM_MOD} ${APIBFM_MOD} ${STDALONE_MOD} ${BFM_OBJ} ${EXE}
#	$(MOVE_MODULES_COMMAND)

$(EXE): $(MAINOBJ) 
	$(FC) $(FFLAGS) -o ../bin/$@ $(MAINOBJ) $(LIB) -L$(NETCDFLIBDIR) -lnetcdf

$(BFM_MOD) : $(BFMSRC)/General/ModuleMem.F90

${BFMSRC}/General/ModuleMem.F90 : $(BFMSRC)/General/GlobalDefsBFM.model
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  -read ${BFMSRC}/General/GlobalDefsBFM.model \
		-from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem allocmem netcdfmem 
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  -read ${BFMSRC}/General/GlobalDefsBFM.model \
		-from ${BFMSRC}/proto -to ${BFMSRC}/include -actions headermem

distclean: clean

clean:
	$(RM) ../modules/$(FORTRAN_COMPILER)/*.{m,mod}
	$(RM) ../lib/$(FORTRAN_COMPILER)/lib*$(buildtype).a
	$(RM) ./share/*.o
	$(RM) ./standalone/*.o
	$(RM) ${BFMSRC}/*/*.o

	$(RM) -r ../modules/*/*.m
	$(RM) -r ../modules/*/*.mod
	$(RM) -r ../src/*.m
	$(RM) -r ../src/*.mod
	$(RM) ./gotm/*.o
	$(RM) ../../femlib/libbfm.a
	$(RM) ../lib/libbfm.a

depend:
	makedepf90 ${BFMSRC}/*/*.?90 ./share/*.?90 ./standalone/*.?90 > .depend

include .depend

