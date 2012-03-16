c
c alukas
c
c revision log :
c
c 14.03.2012    ggu     eliminated PAUSE statements (gfortran)
c
c  CONTENT:
c
c   SUBROUTINE ALUKAS       - main subroutine for WC kinetics
c   SUBROUTINE SMITH_ALUKAS - light limitation factor
c   FUNCTION ALIGHT_ALUKAS  - light extinction 
c   FUNCTION DO_SAT_ALUKAS  - DO saturation concentration
c   FUNCTION KAWIND_ALUKAS  - reaeration coefficient because of wind
c   FUNCTION KAHYDR_ALUKAS  - reaeration coefficient because of wind
c   FUNCTION FKPH_ALUKAS    - pH limiting factor(does it has sense?)
c   FUNCTION TSTROG_ALUKAS  - temperature limitation factor (Stroganov type (bell shaped))
C   SUBROUTINE cur_smith    - old version, not used
C   FUNCTION HFLUX          - 
C   FUNCTION HSWRAD         -
C   FUNCTION HLWRAD         -
C   FUNCTION HBRAD          -
C   FUNCTION HLEVAP         -
C   FUNCTION HTRANS         -
C    
C**********************************************************************
C                                                                     *
C    NAME    : ALUKAS                                                 *
C    PURPOSE : MAIN SUBROUTINE TO MANAGE THE KINETIC CALCULATIONS     *
C                                                                     *
C    ALUKAS IS A SUCCESSOR OF AQUABC THAT BECOMES THE NAME            *
C    OF THE MODEL CONSISTING INSIDE WC (ALUKAS) AND BS models         *
C                                                                     *
C                                                                     *
C    FIRST DEVELOPMENT : OCTOCER 2007, CORPI, KLAIPEDA, LITHUANIA     *
C    DEVELOPER         : ALI ERTURK                                   *
C    CONTRIBUTERS      : -                                            *
C    REVISIONS         :                                              *
C                                                                     *
C    COMPLETELY REVISED TO BE INDEPENDENT FORM ESTAS IN 29/06/2009    *
C                                                                     *
C    02/06/2010                                                       *
C    ADDED THREE STATE VARIABLES FOR SEDIMENT DISSOLVED ORGANIC       *
C    MATTER INTERACTION                                               *
C---------------------------------------------------------------------*
C    DOCUMENTATION                                                    *
C                                                                     *
C                                                                     *
C**********************************************************************

      SUBROUTINE ALUKAS
     *           (STATE_VARIABLES  , DERIVATIVES, NUM_VARS, 
     *            MODEL_CONSTANTS  , NUM_CONSTS , 
     *            DRIVING_FUNCTIONS, NUM_DRVS   , 
     *            PROCESSES        , NUM_PROCESSES,
     *            SAVED_OUTPUTS    , NUM_SAVED_OUTPUTS, 
     *            PSTIME,    CELLNO, LAYER, CALLED_BEFORE)

      IMPLICIT NONE


      INTEGER NUM_VARS
      INTEGER NUM_CONSTS
      INTEGER NUM_DRVS
      INTEGER NUM_PROCESSES
      INTEGER NUM_SAVED_OUTPUTS
      INTEGER CELLNO
      INTEGER LAYER
      INTEGER CALLED_BEFORE

      DOUBLE PRECISION PSTIME


      DOUBLE PRECISION PROCESSES
      DIMENSION PROCESSES(NUM_VARS, NUM_PROCESSES)

      DOUBLE PRECISION DERIVATIVES
      DIMENSION DERIVATIVES(NUM_VARS)

      DOUBLE PRECISION STATE_VARIABLES
      DIMENSION STATE_VARIABLES(NUM_VARS)

      DOUBLE PRECISION MODEL_CONSTANTS
      DIMENSION MODEL_CONSTANTS(NUM_CONSTS)

      DOUBLE PRECISION DRIVING_FUNCTIONS
      DIMENSION DRIVING_FUNCTIONS(NUM_DRVS)

      DOUBLE PRECISION SAVED_OUTPUTS
      DIMENSION SAVED_OUTPUTS(NUM_SAVED_OUTPUTS)


C     STATE VARIABLES
      DOUBLE PRECISION NH3
      DOUBLE PRECISION NOx
      DOUBLE PRECISION OPO4
      DOUBLE PRECISION PHY
      DOUBLE PRECISION EXLADDETC
      DOUBLE PRECISION DOO
      DOUBLE PRECISION EXLAPDETC
      DOUBLE PRECISION EXREDDETC
      DOUBLE PRECISION ZOO
      DOUBLE PRECISION ISI
      DOUBLE PRECISION EXREPDETC
      DOUBLE PRECISION PHY_2
      DOUBLE PRECISION PHY_3
      DOUBLE PRECISION INC
      DOUBLE PRECISION GPHYDDETC
      DOUBLE PRECISION GPHYPDETC
      DOUBLE PRECISION DPHYDDETC
      DOUBLE PRECISION DPHYPDETC
      DOUBLE PRECISION CPHYDDETC
      DOUBLE PRECISION CPHYPDETC
      DOUBLE PRECISION ZOOPDDETC
      DOUBLE PRECISION ZOOPPDETC
      DOUBLE PRECISION SED_DOC
      DOUBLE PRECISION SED_DON
      DOUBLE PRECISION SED_DOP
      

C     MODEL CONSTANTS
C     1=small 2=medium 3=large
      DOUBLE PRECISION WTYPE

C     NH4->NO3 Nitrification rate at 20C, day-1
      DOUBLE PRECISION kcnit

C     NH4->NO3 Nitrification rate temp. const
      DOUBLE PRECISION ktnit

C     NH4->NO3 Half sat. for nitrific., mgO2/l
      DOUBLE PRECISION knit      

C     NO3-> Denitrification rate at 20C, day-1day-1
      DOUBLE PRECISION kcdenit

C     NO3-> Denitrification rate temp. const.unitless
      DOUBLE PRECISION ktdenit

C     NO3-> Half sat. for denitrific., mgO2/L
      DOUBLE PRECISION kdenit    

C     PHYT   growth rate constant  day-1
      DOUBLE PRECISION k1c

C     PHYT_2 growth rate constant  day-1
      DOUBLE PRECISION k1c_2

C     PHYT_3 growth rate constant at 20C, day-1
      DOUBLE PRECISION k1c_3

C     PHY   Carbon half saturation, mgC/l
      DOUBLE PRECISION kc

C     PHY_2 Carbon half saturation, mgC/l
      DOUBLE PRECISION kc_2

C     PHY_3 Carbon half saturation, mgC/l
      DOUBLE PRECISION kc_3

C     PHY   Quantum yield const. mg C/mole photon
      DOUBLE PRECISION PHIMX

C     PHY_2 Quantum yield const. mg C/mole photon
      DOUBLE PRECISION PHIMX_2

C     PHY_3 Quantum yield const. mg C/mole photon, LGHTSW=2
      DOUBLE PRECISION PHIMX_3

C     Chloroph. extinction, mcg [Chla/l]/m
      DOUBLE PRECISION XKC

C     OPA Carbon to chlorophyl ratio
      DOUBLE PRECISION CCHL

C     Diatoms Carbon to chlorophyl ratio
      DOUBLE PRECISION CCHL_2

C     Cyanobacteria Carbon to chlorophyl ratio
      DOUBLE PRECISION CCHL_3

C     Nitrogen half saturation for OPA mgN/l
      DOUBLE PRECISION kn

C     Nitrogen half saturation for Diatoms mgN/l
      DOUBLE PRECISION kn_2      

C     Nitrogen half saturation for Cyanobacteria mgN/l
      DOUBLE PRECISION kn_3         

C     Phosphorus half saturation for OPA mgP/l
      DOUBLE PRECISION kp

C     Phosphorus half saturation for Diatoms mgP/l
      DOUBLE PRECISION kp_2   

C     Phosphorus half saturation for Cyanobacteria mgP/l
      DOUBLE PRECISION kp_3    

C     Silicon half saturation concentration for Diatoms
      DOUBLE PRECISION kmnsi

C     OPA respiration rate
      DOUBLE PRECISION k1rc

C     Diatoms respiration rate
      DOUBLE PRECISION k1rc_2    

C     Cyanobacteria respiration rate
      DOUBLE PRECISION k1rc_3    

C     OPA Phytoplankton death rate,day-1
      DOUBLE PRECISION k1d

C     Diatoms Phytoplankton death rate,day-1
      DOUBLE PRECISION k1d_2 

C     Cyanobacteria Phytoplankton death rate,day-1
      DOUBLE PRECISION k1d_3

C     PHYT Nutrient lim. option 1-min, 0-product
      DOUBLE PRECISION NUTLIM

C     mgP/mgC in OPA
      DOUBLE PRECISION pc

C     mgP/mgC in Diatoms
      DOUBLE PRECISION pc_2

C     mgP/mgC in Cyanobacteria
      DOUBLE PRECISION pc_3

C     mgN/mgC in OPA
      DOUBLE PRECISION nc

C     mgN/mgC in Diatoms
      DOUBLE PRECISION nc_2

C     mgN/mgC in Cyanobacteria
      DOUBLE PRECISION nc_3

C     OPA fraction of photorespiration in primary production
      DOUBLE PRECISION fpr

C     Diatoms fraction of photorespiration in primary production
      DOUBLE PRECISION fpr_2

C     Cyanobacteria fraction of photorespiration in primary production
      DOUBLE PRECISION fpr_3

C     OPA  excretion rate
      DOUBLE PRECISION kexcr

C     Diatoms excretion rate
      DOUBLE PRECISION kexcr_2 

C     Cyanobacteria excretipn rate
      DOUBLE PRECISION kexcr_3

C     PHY  excretion and respiration temperature coefficient
      DOUBLE PRECISION kert

C     PHY_2 excretion and respiration temperature coefficient
      DOUBLE PRECISION kert_2

C     PHY_3 excretion and respiration temperature coefficient
      DOUBLE PRECISION kert_3

C     PHY  death rate temperature coefficient
      DOUBLE PRECISION kdt

C     PHY_2 death rate temperature coefficient
      DOUBLE PRECISION kdt_2

C     PHY_3 death rate temperature coefficient
      DOUBLE PRECISION kdt_3

C     Ratio between N fixing and nonfixing growth rate
      DOUBLE PRECISION R_fix

C     Effectivity of switching to N fixation 
      DOUBLE PRECISION K_fix

C     Preference coefficient for greens
      DOUBLE PRECISION p_1

C     Preference coefficient for diatoms
      DOUBLE PRECISION p_2

C     Preference coefficient for cyanobacteria
      DOUBLE PRECISION p_3

C     Preference coefficient for detritus
      DOUBLE PRECISION p_4

C     Preference coefficient for greens based particulate detritus carbon
      DOUBLE PRECISION p_5

C     Preference coefficient for diatoms basesdparticulate detritus carbon
      DOUBLE PRECISION p_6

C     Preference coefficient for cyanobacteria based particulate 
c     detritus carbon
      DOUBLE PRECISION p_7


C     Minimum concentration of food for zooplankton
      DOUBLE PRECISION FOOD_min

C     Rearation rate at 20C, day-1,if k2!0 then use kawind or kahydra
      DOUBLE PRECISION K2

C     Light extinction coef. m-1
      DOUBLE PRECISION ke

C     min, max temp for pH mineralization factor 
      DOUBLE PRECISION pHminm
      DOUBLE PRECISION pHmaxm          

C     min, max temp for pH nitrification factor
      DOUBLE PRECISION pHminn  
      DOUBLE PRECISION pHmaxn     

C     min, max temp for pH denitrification factor
      DOUBLE PRECISION pHmindn
      DOUBLE PRECISION pHmaxdn

C     PHYT->ZOO Phytoplankton grazing rate by Zoo
      DOUBLE PRECISION kgrz

C     PHYT->ZOO half saturation constant for grazing
      DOUBLE PRECISION kpz

C     PHYT->ZOO zoo-phyto digestion efficiency
      DOUBLE PRECISION eff

C     ZOO death [with excrection]
      DOUBLE PRECISION kdz 

C     ZOO respiration rate
      DOUBLE PRECISION krz

C     ZOO excretion fraction in respiration
      DOUBLE PRECISION kexcrz

C     Silica to carbon ratio for diatoms
      DOUBLE PRECISION asc

C     ZOO growth Q10 factor
      DOUBLE PRECISION q10z

C     ZOO growth maximal temperature
      DOUBLE PRECISION tmaxz

C     ZOO growth optimal temperature
      DOUBLE PRECISION toptz

C     ZOO growth maximal acclimation temperature(delay)
      DOUBLE PRECISION tamaxz

C     ZOO growth maximal acclimation approaching rate
      DOUBLE PRECISION kamaxz

C     ZOO growth temperature below which no acclimation
      DOUBLE PRECISION taminz

C     ZOO respiration temperature coefficient
      DOUBLE PRECISION krtz

C     ZOO death temperature coefficient
      DOUBLE PRECISION kdtz

C     PHY growth Q10 factor
      DOUBLE PRECISION q10

C     PHY growth maximal temperature
      DOUBLE PRECISION tmax

C     PHY growth optimal temperature
      DOUBLE PRECISION topt

C     PHY growth maximal acclimation temperature(delay)
      DOUBLE PRECISION tamax

C     PHY growth maximal acclimation approaching rate
      DOUBLE PRECISION kamax

C     PHY growth temperature below which no acclimation
      DOUBLE PRECISION tamin

C     PHY_2 growth Q10 factor
      DOUBLE PRECISION q10_2

C     PHY_2 growth maximal temperature
      DOUBLE PRECISION tmax_2 

C     PHY_2 growth optimal temperature
      DOUBLE PRECISION topt_2

C     PHY_2 growth maximal acclimation temperature(delay)
      DOUBLE PRECISION tamax_2

C     PHY_2 growth maximal acclimation approaching rate
      DOUBLE PRECISION kamax_2

C     PHY_2 growth temperature below which no acclimation
      DOUBLE PRECISION tamin_2

C     PHY_3 growth Q10 factor
      DOUBLE PRECISION q10_3

C     PHY_3 growth maximal temperature
      DOUBLE PRECISION tmax_3

C     PHY_3 growth optimal temperature
      DOUBLE PRECISION topt_3

C     PHY_3 growth maximal acclimation temperature(delay)
      DOUBLE PRECISION tamax_3

C     PHY_3 growth maximal acclimation approaching rate
      DOUBLE PRECISION kamax_3

C     PHY_3 growth temperature below which no acclimation
      DOUBLE PRECISION tamin_3



C     Rate constant for ZOOPPDETC dissolution
      DOUBLE PRECISION C_DISS_ZOOPPDETC

C     Temperature correction for ZOOPPDETC dissolution
      DOUBLE PRECISION T_DISS_ZOOPPDETC

C     Rate constant for GPHYPDETC dissolution
      DOUBLE PRECISION C_DISS_GPHYPDETC

C     Temperature correction for GPHYPDETC dissolution
      DOUBLE PRECISION T_DISS_GPHYPDETC

C     Rate constant for DPHYPDETC dissolution
      DOUBLE PRECISION C_DISS_DPHYPDETC

C     Temperature correction for DPHYPDETC dissolution
      DOUBLE PRECISION T_DISS_DPHYPDETC

C     Rate constant for CPHYPDETC dissolution
      DOUBLE PRECISION C_DISS_CPHYPDETC

C     Temperature correction for CPHYPDETC dissolution
      DOUBLE PRECISION T_DISS_CPHYPDETC

C     Rate constant for EXLAPDETC dissolution
      DOUBLE PRECISION C_DISS_EXLAPDETC

C     Temperature correction for EXLAPDETC dissolution
      DOUBLE PRECISION T_DISS_EXLAPDETC

C     N:C ratio for EXLAPDETC
      DOUBLE PRECISION NC_EXLAPDETC

C     P:C ratio for EXLAPDETC
      DOUBLE PRECISION PC_EXLAPDETC

C     Si:C ratio for EXLAPDETC
      DOUBLE PRECISION SiC_EXLAPDETC



C     Rate constant for EXREPDETC dissolution
      DOUBLE PRECISION C_DISS_EXREPDETC

C     Temperature correction for EXREPDETC dissolution
      DOUBLE PRECISION T_DISS_EXREPDETC

C     N:C ratio for EXREPDETC
      DOUBLE PRECISION NC_EXREPDETC

C     P:C ratio for EXREPDETC
      DOUBLE PRECISION PC_EXREPDETC

C     Si:C ratio for EXREPDETC
      DOUBLE PRECISION SiC_EXREPDETC



C     Rate constant for ZOOPDDETC oxidation
      DOUBLE PRECISION C_OX_ZOOPDDETC

C     Temprature correction factor for ZOOPDDETC oxidation
      DOUBLE PRECISION T_OX_ZOOPDDETC

C     Minimum pH for ZOOPDDETC oxidation
      DOUBLE PRECISION pHminm_OX_ZOOPDDETC

C     Maximum pH for ZOOPDDETC oxidation
      DOUBLE PRECISION pHmaxm_OX_ZOOPDDETC

C     Half saturation for oxygen for ZOOPDDETC oxidation
      DOUBLE PRECISION k_OX_ZOOPDDETC

C     O:C ratio for ZOOPDDETC
      DOUBLE PRECISION OC_ZOOPDDETC



C     Rate constant for GPHYDDETC oxidation
      DOUBLE PRECISION C_OX_GPHYDDETC  

C     Temprature correction factor for GPHYDDETC oxidation
      DOUBLE PRECISION T_OX_GPHYDDETC

C     Minimum pH for GPHYDDETC oxidation
      DOUBLE PRECISION pHminm_OX_GPHYDDETC

C     Maximum pH for GPHYDDETC oxidation
      DOUBLE PRECISION pHmaxm_OX_GPHYDDETC

C     Half saturation for oxygen for GPHYDDETC oxidation
      DOUBLE PRECISION k_OX_GPHYDDETC

C     O:C ratio for GPHYDDETC
      DOUBLE PRECISION OC_GPHYDDETC



C     Rate constant for DPHYDDETC oxidation
      DOUBLE PRECISION C_OX_DPHYDDETC

C     Temprature correction factor for DPHYDDETC oxidation
      DOUBLE PRECISION T_OX_DPHYDDETC

C     Minimum pH for DPHYDDETC oxidation
      DOUBLE PRECISION pHminm_OX_DPHYDDETC

C     Maximum pH for DPHYDDETC oxidation
      DOUBLE PRECISION pHmaxm_OX_DPHYDDETC

C     Half saturation for oxygen for DPHYDDETC oxidation
      DOUBLE PRECISION k_OX_DPHYDDETC

C     O:C ratio for DPHYDDETC
      DOUBLE PRECISION OC_DPHYDDETC



C     Rate constant for CPHYDDETC oxidation
      DOUBLE PRECISION C_OX_CPHYDDETC

C     Temprature correction factor for CPHYDDETC oxidation
      DOUBLE PRECISION T_OX_CPHYDDETC

C     Minimum pH for CPHYDDETC oxidation
      DOUBLE PRECISION pHminm_OX_CPHYDDETC

C     Maximum pH for CPHYDDETC oxidation
      DOUBLE PRECISION pHmaxm_OX_CPHYDDETC

C     Half saturation for oxygen for CPHYDDETC oxidation
      DOUBLE PRECISION k_OX_CPHYDDETC

C     O:C ratio for CPHYDDETC
      DOUBLE PRECISION OC_CPHYDDETC



C     Rate constant for EXLADDETC oxidation
      DOUBLE PRECISION C_OX_EXLADDETC
      
C     Temprature correction factor for EXLADDETC oxidation
      DOUBLE PRECISION T_OX_EXLADDETC

C     Minimum pH for EXLADDETC oxidation
      DOUBLE PRECISION pHminm_OX_EXLADDETC

C     Maximum pH for EXLADDETC oxidation
      DOUBLE PRECISION pHmaxm_OX_EXLADDETC

C     Half saturation for oxygen for EXLADDETC oxidation
      DOUBLE PRECISION k_OX_EXLADDETC

C     O:C ratio for EXLADDETC
      DOUBLE PRECISION OC_EXLADDETC



C     Rate constant for EXREDDETC oxidation
      DOUBLE PRECISION C_OX_EXREDDETC

C     Temprature correction factor for EXREDDETC oxidation
      DOUBLE PRECISION T_OX_EXREDDETC

C     Minimum pH for EXREDDETC oxidation
      DOUBLE PRECISION pHminm_OX_EXREDDETC

C     Maximum pH for EXREDDETC oxidation
      DOUBLE PRECISION pHmaxm_OX_EXREDDETC

C     Half saturation for oxygen for EXREDDETC oxidation
      DOUBLE PRECISION k_OX_EXREDDETC

C     O:C ratio for EXREDDETC
      DOUBLE PRECISION OC_EXREDDETC






C     AUXILARY VARIABLES
      DOUBLE PRECISION TPHY
      DOUBLE PRECISION TCHLA
      DOUBLE PRECISION XRT
      DOUBLE PRECISION XRT_2
      DOUBLE PRECISION XRT_3
      DOUBLE PRECISION GP1_T
      DOUBLE PRECISION GP2_T
      DOUBLE PRECISION GP3_T      
      DOUBLE PRECISION LLIGHT
      DOUBLE PRECISION LLIGHT_2
      DOUBLE PRECISION LLIGHT_3
      DOUBLE PRECISION CCHLX0
      
      DOUBLE PRECISION X1
      DOUBLE PRECISION X1_2
      DOUBLE PRECISION X1_3
      DOUBLE PRECISION X1_3_fix 
      DOUBLE PRECISION X2
      DOUBLE PRECISION X2_2
      DOUBLE PRECISION X2_3
      DOUBLE PRECISION X3 
      DOUBLE PRECISION X4
      DOUBLE PRECISION X4_2
      DOUBLE PRECISION X4_3
      DOUBLE PRECISION LNUT
      DOUBLE PRECISION LNUT_2
      DOUBLE PRECISION LNUT_3
      DOUBLE PRECISION LNUT_3_fix 
      DOUBLE PRECISION GP1
      DOUBLE PRECISION GP2
      DOUBLE PRECISION GP3_1
      DOUBLE PRECISION GP3_2
      DOUBLE PRECISION GP3       
      DOUBLE PRECISION GPP
      DOUBLE PRECISION GPP_2
      DOUBLE PRECISION GPP_3
      DOUBLE PRECISION RRT
      DOUBLE PRECISION RRT_2
      DOUBLE PRECISION RRT_3
      
      DOUBLE PRECISION RESD
      DOUBLE PRECISION RESPH        
      DOUBLE PRECISION RES
      DOUBLE PRECISION RES_2
      DOUBLE PRECISION RES_3
      DOUBLE PRECISION EXCR
      DOUBLE PRECISION EXCR_2
      DOUBLE PRECISION EXCR_3
      
      DOUBLE PRECISION DR
      DOUBLE PRECISION DR_2
      DOUBLE PRECISION DR_3
      DOUBLE PRECISION DED
      DOUBLE PRECISION DED2
      DOUBLE PRECISION DED3 
      DOUBLE PRECISION DP1
      DOUBLE PRECISION DP2
      DOUBLE PRECISION DP3       
      DOUBLE PRECISION DEDP
      DOUBLE PRECISION DEDP_2
      DOUBLE PRECISION DEDP_3          
      DOUBLE PRECISION DPP
      DOUBLE PRECISION DPP_2
      DOUBLE PRECISION DPP_3
      
      DOUBLE PRECISION GRZ_1
      DOUBLE PRECISION GRZ_2
      DOUBLE PRECISION GRZ_3
      DOUBLE PRECISION GRZ_4
      DOUBLE PRECISION GRZ_5
      DOUBLE PRECISION GRZ_6
      DOUBLE PRECISION GRZ_7
      DOUBLE PRECISION GRZ
      
      DOUBLE PRECISION DZRT
      DOUBLE PRECISION DZ
      DOUBLE PRECISION RZRT
      DOUBLE PRECISION RZ 
      DOUBLE PRECISION EXZ
      DOUBLE PRECISION XZRT 
      DOUBLE PRECISION GZ 
      DOUBLE PRECISION DEFEC
      DOUBLE PRECISION ZSINK
      
      DOUBLE PRECISION kpHn
      DOUBLE PRECISION N1
      DOUBLE PRECISION PN      
      DOUBLE PRECISION PN_2      
      DOUBLE PRECISION PN_3      
      DOUBLE PRECISION NALG2
      
      DOUBLE PRECISION NOALG
      DOUBLE PRECISION kpHdn 
      DOUBLE PRECISION NIT1
      
      
      DOUBLE PRECISION OPALG2
      
      DOUBLE PRECISION OX
      
      DOUBLE PRECISION OSAT
      DOUBLE PRECISION K2WIND
      DOUBLE PRECISION KA
      DOUBLE PRECISION DO1
      DOUBLE PRECISION DO2
      DOUBLE PRECISION DO3
      DOUBLE PRECISION DO4
      DOUBLE PRECISION N2
      
      DOUBLE PRECISION UISI



      DOUBLE PRECISION Resp 
      DOUBLE PRECISION Kliq_co2
      DOUBLE PRECISION power10    
      DOUBLE PRECISION co2sat 
      DOUBLE PRECISION atmexch
      DOUBLE PRECISION phot

      DOUBLE PRECISION T0
      DOUBLE PRECISION FISI
      DOUBLE PRECISION SOD1
      DOUBLE PRECISION SODT
      DOUBLE PRECISION EPS
      DOUBLE PRECISION CCHLX
      DOUBLE PRECISION FOPO4
      DOUBLE PRECISION KPO4
      DOUBLE PRECISION CCHLX_2
      DOUBLE PRECISION CCHLX_3
      DOUBLE PRECISION K2HYDR
      DOUBLE PRECISION VEL
      DOUBLE PRECISION PH
      DOUBLE PRECISION DEPTH
      DOUBLE PRECISION SHORTR
      DOUBLE PRECISION FRADWL


C     C:N  ratio for zooplankton carbon (calculated)
      DOUBLE PRECISION NC_ZOOC

C     C:P  ratio for zooplankton carbon (calculated)
      DOUBLE PRECISION PC_ZOOC

C     C:Si ratio for zooplankton carbon (calculated)
      DOUBLE PRECISION SiC_ZOOC

C     O:C ratio for  zooplankton carbon (calculated)
      DOUBLE PRECISION OC_ZOOC



C     ZOOPPDETC generation rate by ZOOC death
      DOUBLE PRECISION DEATH_ZOOC_ZOOPPDETC

C     Dissolution rate of ZOOPPDETC to ZOOPDDETC
      DOUBLE PRECISION DISS_ZOOPPDETC



C     GPHYPDETC generation rate by PHYGC death
      DOUBLE PRECISION DEATH_PYHGC_GPHYPDETC !

C     Dissolution rate of GPHYPDETC to GPHYDDETC
      DOUBLE PRECISION DISS_GPHYPDETC



C     DPHYPDETC generation rate by PHYDC death
      DOUBLE PRECISION DEATH_PYHDC_DPHYPDETC

C     Dissolution rate of DPHYPDETC to DPHYDDETC
      DOUBLE PRECISION DISS_DPHYPDETC



C     CPHYPDETC generation rate by PHYCC death
      DOUBLE PRECISION DEATH_PYHCC_CPHYPDETC

C     Dissolution rate of CPHYPDETC to CPHYDDETC
      DOUBLE PRECISION DISS_CPHYPDETC



C     Dissolution rate of EXLAPDETC to EXLADDETC
      DOUBLE PRECISION DISS_EXLAPDETC



C     Dissolution rate of EXREPDETC to EXREDDETC
      DOUBLE PRECISION DISS_EXREPDETC



C     ZOOPDDETC oxidation rate
      DOUBLE PRECISION OX_ZOOPDDETC

C     pH limitation factor for ZOOPDDETC oxidation
      DOUBLE PRECISION kpH_OX_ZOOPDDETC

C     Inorganic carbon generation rate by  ZOOPDDETC oxidation
      DOUBLE PRECISION INC_OX_ZOOPDDETC

C     Ammonia nitrogen generation rate by ZOOPDDETC oxidation
      DOUBLE PRECISION N_OX_ZOOPDDETC

C     Phosphate phosphorus generation rate by ZOOPDDETC oxidation
      DOUBLE PRECISION P_OX_ZOOPDDETC

C     Avialable silicon generation rate by ZOOPDDETC oxidation
      DOUBLE PRECISION Si_OX_ZOOPDDETC

C     ZOOPDDETC generation because of excretion
      DOUBLE PRECISION EXCR_ZOOPDDETC



C     GPHYDDETC oxidation rate
      DOUBLE PRECISION OX_GPHYDDETC

C     pH limitation factor for GPHYDDETC oxidation
      DOUBLE PRECISION kpH_OX_GPHYDDETC

C     Inorganic carbon generation rate by GPHYDDETC oxidation
      DOUBLE PRECISION INC_OX_GPHYDDETC

C     Ammonia nitrogen generation rate by GPHYDDETC oxidation
      DOUBLE PRECISION N_OX_GPHYDDETC

C     Phosphate phosphorus generation rate by GPHYDDETC oxidation
      DOUBLE PRECISION P_OX_GPHYDDETC

C     GPHYDDETC generation because of excretion
      DOUBLE PRECISION EXCR_GPHYDDETC



C     DPHYDDETC oxidation rate
      DOUBLE PRECISION OX_DPHYDDETC

C     pH limitation factor for DPHYDDETC oxidation
      DOUBLE PRECISION kpH_OX_DPHYDDETC

C     Inorganic carbon generation rate by DPHYDDETC oxidation
      DOUBLE PRECISION INC_OX_DPHYDDETC

C     Ammonia nitrogen generation rate by DPHYDDETC oxidation
      DOUBLE PRECISION N_OX_DPHYDDETC

C     Phosphate phosphorus generation rate by DPHYDDETC oxidation
      DOUBLE PRECISION P_OX_DPHYDDETC

C     Avialable silicon generation rate by DPHYDDETC oxidation
      DOUBLE PRECISION Si_OX_DPHYDDETC

C     DPHYDDETC generation because of excretion
      DOUBLE PRECISION EXCR_DPHYDDETC



C     CPHYDDETC oxidation rate
      DOUBLE PRECISION OX_CPHYDDETC

C     pH limitation factor for CPHYDDETC oxidation
      DOUBLE PRECISION kpH_OX_CPHYDDETC

C     Inorganic carbon generation rate by CPHYDDETC oxidation
      DOUBLE PRECISION INC_OX_CPHYDDETC

C     Ammonia nitrogen generation rate by CPHYDDETC oxidation
      DOUBLE PRECISION N_OX_CPHYDDETC

C     Phosphate phosphorus generation rate by CPHYDDETC oxidation
      DOUBLE PRECISION P_OX_CPHYDDETC

C     CPHYDDETC generation because of excretion
      DOUBLE PRECISION EXCR_CPHYDDETC



C     EXLADDETC oxidation rate
      DOUBLE PRECISION OX_EXLADDETC

C     Sediment based DOC oxidation rate 
      DOUBLE PRECISION OX_SED_DOC

C     Sediment based DON oxidation rate 
      DOUBLE PRECISION OX_SED_DON

C     Sediment based DOP oxidation rate 
      DOUBLE PRECISION OX_SED_DOP
      
C     pH limitation factor for EXLADDETC oxidation
      DOUBLE PRECISION kpH_OX_EXLADDETC

C     Inorganic carbon generation rate by EXLADDETC oxidation
      DOUBLE PRECISION INC_OX_EXLADDETC

C     Ammonia nitrogen generation rate by EXLADDETC oxidation
      DOUBLE PRECISION N_OX_EXLADDETC

C     Phosphate phosphorus generation rate by EXLADDETC oxidation
      DOUBLE PRECISION P_OX_EXLADDETC

C     Avialable silicon generation rate by EXLADDETC oxidation
      DOUBLE PRECISION Si_OX_EXLADDETC

C     N:C ratio for EXLADDETC
      DOUBLE PRECISION NC_EXLADDETC

C     P:C ratio for EXLADDETC
      DOUBLE PRECISION PC_EXLADDETC

C     Si:C ratio for EXLADDETC
      DOUBLE PRECISION SiC_EXLADDETC




C     EXREDDETC oxidation rate
      DOUBLE PRECISION OX_EXREDDETC

C     pH limitation factor for EXREDDETC oxidation
      DOUBLE PRECISION kpH_OX_EXREDDETC

C     Inorganic carbon generation rate by EXREDDETC oxidation
      DOUBLE PRECISION INC_OX_EXREDDETC

C     Ammonia nitrogen generation rate by EXREDDETC oxidation
      DOUBLE PRECISION N_OX_EXREDDETC

C     Phosphate phosphorus generation rate by EXREDDETC oxidation
      DOUBLE PRECISION P_OX_EXREDDETC

C     Avialable silicon generation rate by EXREDDETC oxidation
      DOUBLE PRECISION Si_OX_EXREDDETC

C     N:C ratio for EXREDDETC
      DOUBLE PRECISION NC_EXREDDETC

C     P:C ratio for EXREDDETC
      DOUBLE PRECISION PC_EXREDDETC

C     Si:C ratio for EXREDDETC
      DOUBLE PRECISION SiC_EXREDDETC



C     Grazing ratio for greens
      DOUBLE PRECISION GRAZ_RAT_1

C     Grazing ratio for diatoms
      DOUBLE PRECISION GRAZ_RAT_2

C     Grazing ratio for cyanobacteria
      DOUBLE PRECISION GRAZ_RAT_3

C     Grazing ratio for external labile particulate detritus carbon
      DOUBLE PRECISION GRAZ_RAT_4

C     Grazing ratio for greens based particulate detritus carbon
      DOUBLE PRECISION GRAZ_RAT_5

C     Grazing ratio for diatoms based particulate detritus carbon 
      DOUBLE PRECISION GRAZ_RAT_6

C     Grazing ratio for cyanobacteria based particulate detritus carbon 
      DOUBLE PRECISION GRAZ_RAT_7


C     Rate of ammonium nitrogen generation by zooplankton respiration
      DOUBLE PRECISION N_RESP_ZOOC

C     Rate of ammonium nitrogen generation by green phytoplankton respiration
      DOUBLE PRECISION N_RESP_PHYGC

C     Rate of ammonium nitrogen generation by diatom respiration
      DOUBLE PRECISION N_RESP_PHYDC

C     Rate of ammonium nitrogen generation by cyanobacteria respiration
      DOUBLE PRECISION N_RESP_PHYCC



C     Rate of phosphate phosphorus generation by zooplankton respiration
      DOUBLE PRECISION P_RESP_ZOOC

C     Rate of phosphate phosphorus generation by green phytoplankton respiration
      DOUBLE PRECISION P_RESP_PHYGC

C     Rate of phosphate phosphorus generation by diatom respiration
      DOUBLE PRECISION P_RESP_PHYDC

C     Rate of phosphate phosphorus generation by cyanobacteria respiration
      DOUBLE PRECISION P_RESP_PHYCC


C     Rate of available silicon generation by zooplankton respiration
      DOUBLE PRECISION Si_RESP_ZOOC

C     Rate of phosphate phosphorus generation by diatom respiration
      DOUBLE PRECISION Si_RESP_PHYDC



C     Rate of oxygen consumption by zooplankton respiration
      DOUBLE PRECISION OX_RESP_ZOOC

C     Rate of oxygen consumption by green phytoplankton respiration
      DOUBLE PRECISION OX_RESP_PHYGC

C     Rate of oxygen consumption by diatom respiration
      DOUBLE PRECISION OX_RESP_PHYDC

C     Rate of oxygen consumption by cyanobacteria respiration
      DOUBLE PRECISION OX_RESP_PHYCC


      INTEGER RUNTIME_ERROR


C     ALTERNATIVES
C      INTEGER SEDSEG
C      INTEGER TOPSEG
      INTEGER GRTYPE


C     FORCING VARIABLES
      DOUBLE PRECISION LTEMP
      DOUBLE PRECISION STP20
      DOUBLE PRECISION FDAY
      DOUBLE PRECISION WIND
      DOUBLE PRECISION AIRTMP
      DOUBLE PRECISION SALIN
      DOUBLE PRECISION SOD1D
      DOUBLE PRECISION SODTA
      DOUBLE PRECISION BRATIO
      DOUBLE PRECISION XICECVR
      DOUBLE PRECISION SURF_ELEVATION


C     AUXILARY VARIABLES
      DOUBLE PRECISION IA
C      DOUBLE PRECISION PHYTCC
C      DIMENSION PHYTCC(MCELLS, MLAYER)

      DOUBLE PRECISION TOPSEG


C     FUNCTIONS CALLED
      DOUBLE PRECISION ALIGHT_ALUKAS
      DOUBLE PRECISION KAWIND_ALUKAS
      DOUBLE PRECISION KAHYDR_ALUKAS
      DOUBLE PRECISION DO_SAT_ALUKAS
      DOUBLE PRECISION TSTROG_ALUKAS
      DOUBLE PRECISION FKPH_ALUKAS


      INTEGER I


      NH3       = STATE_VARIABLES(1)
      NOx       = STATE_VARIABLES(2)
      OPO4      = STATE_VARIABLES(3)
      PHY       = STATE_VARIABLES(4)
      EXLADDETC = STATE_VARIABLES(5)
      DOO       = STATE_VARIABLES(6)
      EXLAPDETC = STATE_VARIABLES(7)
      EXREDDETC = STATE_VARIABLES(8)
      ZOO       = STATE_VARIABLES(9) 
      ISI       = STATE_VARIABLES(10)
      EXREPDETC = STATE_VARIABLES(11)
      PHY_2     = STATE_VARIABLES(12)
      PHY_3     = STATE_VARIABLES(13)
      INC       = STATE_VARIABLES(14)
      GPHYDDETC = STATE_VARIABLES(15)
      GPHYPDETC = STATE_VARIABLES(16)
      DPHYDDETC = STATE_VARIABLES(17)
      DPHYPDETC = STATE_VARIABLES(18)
      CPHYDDETC = STATE_VARIABLES(19)
      CPHYPDETC = STATE_VARIABLES(20)
      ZOOPDDETC = STATE_VARIABLES(21)
      ZOOPPDETC = STATE_VARIABLES(22)
      SED_DOC   = STATE_VARIABLES(23)
      SED_DON   = STATE_VARIABLES(24)
      SED_DOP   = STATE_VARIABLES(25)

      IF (CALLED_BEFORE.EQ.0) THEN
          CCHLX   = MODEL_CONSTANTS(46)
          CCHLX_2 = MODEL_CONSTANTS(246)
          CCHLX_3 = MODEL_CONSTANTS(346)
      ELSE
          CCHLX   = SAVED_OUTPUTS(1)
          CCHLX_2 = SAVED_OUTPUTS(2)
          CCHLX_3 = SAVED_OUTPUTS(3)
      END IF


      WTYPE               = MODEL_CONSTANTS(1)
      kcnit               = MODEL_CONSTANTS(11)
      ktnit               = MODEL_CONSTANTS(12)
      knit                = MODEL_CONSTANTS(13)
      kcdenit             = MODEL_CONSTANTS(21)
      ktdenit             = MODEL_CONSTANTS(22)
      kdenit              = MODEL_CONSTANTS(23)
      k1c                 = MODEL_CONSTANTS(41)
      k1c_2               = MODEL_CONSTANTS(241) 
      k1c_3               = MODEL_CONSTANTS(341)
      kc                  = MODEL_CONSTANTS(42)
      kc_2                = MODEL_CONSTANTS(242)
      kc_3                = MODEL_CONSTANTS(342)
      PHIMX               = MODEL_CONSTANTS(44)
      PHIMX_2             = MODEL_CONSTANTS(244)
      PHIMX_3             = MODEL_CONSTANTS(344)
      XKC                 = MODEL_CONSTANTS(45)
      kn                  = MODEL_CONSTANTS(48)
      kn_2                = MODEL_CONSTANTS(248)
      kn_3                = MODEL_CONSTANTS(348)
      kp                  = MODEL_CONSTANTS(49)
      kp_2                = MODEL_CONSTANTS(249)
      kp_3                = MODEL_CONSTANTS(349)
      kmnsi               = MODEL_CONSTANTS(51)
      k1rc                = MODEL_CONSTANTS(50)
      k1rc_2              = MODEL_CONSTANTS(250)
      k1rc_3              = MODEL_CONSTANTS(350)
      k1d                 = MODEL_CONSTANTS(52)
      k1d_2               = MODEL_CONSTANTS(252) 
      k1d_3               = MODEL_CONSTANTS(352)
      NUTLIM              = MODEL_CONSTANTS(54)
      pc                  = MODEL_CONSTANTS(57)
      pc_2                = MODEL_CONSTANTS(257)
      pc_3                = MODEL_CONSTANTS(357)
      nc                  = MODEL_CONSTANTS(58)
      nc_2                = MODEL_CONSTANTS(258)
      nc_3                = MODEL_CONSTANTS(358)
      fpr                 = MODEL_CONSTANTS(60)
      fpr_2               = MODEL_CONSTANTS(260)
      fpr_3               = MODEL_CONSTANTS(360)
      kexcr               = MODEL_CONSTANTS(61)
      kexcr_2             = MODEL_CONSTANTS(261)
      kexcr_3             = MODEL_CONSTANTS(361)
      kert                = MODEL_CONSTANTS(62)
      kert_2              = MODEL_CONSTANTS(262)
      kert_3              = MODEL_CONSTANTS(362)
      kdt                 = MODEL_CONSTANTS(63)
      kdt_2               = MODEL_CONSTANTS(263)
      kdt_3               = MODEL_CONSTANTS(363)
      R_fix               = MODEL_CONSTANTS(65)
      K_fix               = MODEL_CONSTANTS(66)
      p_1                 = MODEL_CONSTANTS(67)
      p_2                 = MODEL_CONSTANTS(267)
      p_3                 = MODEL_CONSTANTS(367)
      p_4                 = MODEL_CONSTANTS(467)
      p_5                 = MODEL_CONSTANTS(468)
      p_6                 = MODEL_CONSTANTS(469)
      p_7                 = MODEL_CONSTANTS(470)
      FOOD_min            = MODEL_CONSTANTS(68)
      K2                  = MODEL_CONSTANTS(82)
      ke                  = MODEL_CONSTANTS(83)
      pHminm              = MODEL_CONSTANTS(105)
      pHmaxm              = MODEL_CONSTANTS(106)
      pHminn              = MODEL_CONSTANTS(107)
      pHmaxn              = MODEL_CONSTANTS(108)
      pHmindn             = MODEL_CONSTANTS(109)
      pHmaxdn             = MODEL_CONSTANTS(110)
      kgrz                = MODEL_CONSTANTS(112)
      kpz                 = MODEL_CONSTANTS(113)
      eff                 = MODEL_CONSTANTS(114)
      kdz                 = MODEL_CONSTANTS(115)
      krz                 = MODEL_CONSTANTS(120)
      kexcrz              = MODEL_CONSTANTS(121)
      asc                 = MODEL_CONSTANTS(211)
      q10z                = MODEL_CONSTANTS(130)
      tmaxz               = MODEL_CONSTANTS(131)
      toptz               = MODEL_CONSTANTS(132)
      tamaxz              = MODEL_CONSTANTS(133)
      kamaxz              = MODEL_CONSTANTS(134)
      taminz              = MODEL_CONSTANTS(135)
      krtz                = MODEL_CONSTANTS(136)
      kdtz                = MODEL_CONSTANTS(137)
      q10                 = MODEL_CONSTANTS(180)
      tmax                = MODEL_CONSTANTS(181)
      topt                = MODEL_CONSTANTS(182)
      tamax               = MODEL_CONSTANTS(183)
      kamax               = MODEL_CONSTANTS(184)
      tamin               = MODEL_CONSTANTS(185) 
      q10_2               = MODEL_CONSTANTS(280)
      tmax_2              = MODEL_CONSTANTS(281)
      topt_2              = MODEL_CONSTANTS(282)
      tamax_2             = MODEL_CONSTANTS(283)
      kamax_2             = MODEL_CONSTANTS(284)
      tamin_2             = MODEL_CONSTANTS(285)
      q10_3               = MODEL_CONSTANTS(380)
      tmax_3              = MODEL_CONSTANTS(381) 
      topt_3              = MODEL_CONSTANTS(382)
      tamax_3             = MODEL_CONSTANTS(383)
      kamax_3             = MODEL_CONSTANTS(384)
      tamin_3             = MODEL_CONSTANTS(385)
      C_DISS_ZOOPPDETC    = MODEL_CONSTANTS(501)
      T_DISS_ZOOPPDETC    = MODEL_CONSTANTS(502)
      C_DISS_GPHYPDETC    = MODEL_CONSTANTS(503)
      T_DISS_GPHYPDETC    = MODEL_CONSTANTS(504)
      C_DISS_DPHYPDETC    = MODEL_CONSTANTS(505)
      T_DISS_DPHYPDETC    = MODEL_CONSTANTS(506)
      C_DISS_CPHYPDETC    = MODEL_CONSTANTS(507)
      T_DISS_CPHYPDETC    = MODEL_CONSTANTS(508)
      C_DISS_EXLAPDETC    = MODEL_CONSTANTS(509)
      T_DISS_EXLAPDETC    = MODEL_CONSTANTS(510)
      NC_EXLAPDETC        = MODEL_CONSTANTS(511)
      PC_EXLAPDETC        = MODEL_CONSTANTS(512)
      SiC_EXLAPDETC       = MODEL_CONSTANTS(513)
      C_DISS_EXREPDETC    = MODEL_CONSTANTS(514)
      T_DISS_EXREPDETC    = MODEL_CONSTANTS(515)
      NC_EXREPDETC        = MODEL_CONSTANTS(516)
      PC_EXREPDETC        = MODEL_CONSTANTS(517)
      SiC_EXREPDETC       = MODEL_CONSTANTS(518)
      C_OX_ZOOPDDETC      = MODEL_CONSTANTS(519)
      T_OX_ZOOPDDETC      = MODEL_CONSTANTS(520)
      pHminm_OX_ZOOPDDETC = MODEL_CONSTANTS(521)
      pHmaxm_OX_ZOOPDDETC = MODEL_CONSTANTS(522)
      k_OX_ZOOPDDETC      = MODEL_CONSTANTS(523)
      OC_ZOOPDDETC        = MODEL_CONSTANTS(524)
      C_OX_GPHYDDETC      = MODEL_CONSTANTS(525)
      T_OX_GPHYDDETC      = MODEL_CONSTANTS(526)
      pHminm_OX_GPHYDDETC = MODEL_CONSTANTS(527)
      pHmaxm_OX_GPHYDDETC = MODEL_CONSTANTS(528)
      k_OX_GPHYDDETC      = MODEL_CONSTANTS(529)
      OC_GPHYDDETC        = MODEL_CONSTANTS(530)
      C_OX_DPHYDDETC      = MODEL_CONSTANTS(531)
      T_OX_DPHYDDETC      = MODEL_CONSTANTS(532)
      pHminm_OX_DPHYDDETC = MODEL_CONSTANTS(533)
      pHmaxm_OX_DPHYDDETC = MODEL_CONSTANTS(534)
      k_OX_DPHYDDETC      = MODEL_CONSTANTS(535)
      OC_DPHYDDETC        = MODEL_CONSTANTS(536)
      C_OX_CPHYDDETC      = MODEL_CONSTANTS(537)
      T_OX_CPHYDDETC      = MODEL_CONSTANTS(538)
      pHminm_OX_CPHYDDETC = MODEL_CONSTANTS(539)
      pHmaxm_OX_CPHYDDETC = MODEL_CONSTANTS(540)
      k_OX_CPHYDDETC      = MODEL_CONSTANTS(541)
      OC_CPHYDDETC        = MODEL_CONSTANTS(542)
      C_OX_EXLADDETC      = MODEL_CONSTANTS(543)
      T_OX_EXLADDETC      = MODEL_CONSTANTS(544)
      pHminm_OX_EXLADDETC = MODEL_CONSTANTS(545)
      pHmaxm_OX_EXLADDETC = MODEL_CONSTANTS(546)
      k_OX_EXLADDETC      = MODEL_CONSTANTS(547)
      OC_EXLADDETC        = MODEL_CONSTANTS(548)
      C_OX_EXREDDETC      = MODEL_CONSTANTS(549)
      T_OX_EXREDDETC      = MODEL_CONSTANTS(550)
      pHminm_OX_EXREDDETC = MODEL_CONSTANTS(551)
      pHmaxm_OX_EXREDDETC = MODEL_CONSTANTS(552)
      k_OX_EXREDDETC      = MODEL_CONSTANTS(553)
      OC_EXREDDETC        = MODEL_CONSTANTS(554)

      SALIN          = DRIVING_FUNCTIONS(1)
      LTEMP          = DRIVING_FUNCTIONS(2)
      PH             = DRIVING_FUNCTIONS(3)
      BRATIO         = DRIVING_FUNCTIONS(4)
      IA             = DRIVING_FUNCTIONS(5)
      DEPTH          = DRIVING_FUNCTIONS(6)
      TOPSEG         = DRIVING_FUNCTIONS(7)
      WIND           = DRIVING_FUNCTIONS(8)
      XICECVR        = DRIVING_FUNCTIONS(9)
      SURF_ELEVATION = DRIVING_FUNCTIONS(10)
      SHORTR         = DRIVING_FUNCTIONS(11)
      FRADWL         = DRIVING_FUNCTIONS(12)

c      write(*,*) 'CELLNO         : ', CELLNO
c      write(*,*) 'LAYER          : ', LAYER
c      write(*,*) 'PSTIME         : ', PSTIME
c      write(*,*) 'SALIN          : ', SALIN
c      write(*,*) 'LTEMP          : ', LTEMP 
c      write(*,*) 'PH             : ', PH
c      write(*,*) 'BRATIO         : ', BRATIO
c      write(*,*) 'IA             : ', IA
c      write(*,*) 'DEPTH          : ', DEPTH
c      write(*,*) 'TOPSEG         : ', TOPSEG
c      write(*,*) 'WIND           : ', WIND
c      write(*,*) 'XICECVR        : ', XICECVR
c      write(*,*) 'SURF_ELEVATION : ', SURF_ELEVATION
c      write(*,*) 'SHORTR         : ', SHORTR
c      write(*,*) 'FRADWL         : ', FRADWL

c      pause

      STP20 = LTEMP - 20.0

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TIME SERIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FDAY  = FRADWL


 
      T0      = 20.0D0    ! Reference temperature
      FISI    = 0.6D0     ! fraction of dissolved nonbiogenic silica 
      sod1    = 2.        ! sediment oxygen demand
      sodt    = 1.08      ! temperature coeficient for sediment oxygen demand 
      EPS     = 1e-20     ! small number


      SOD1D   = 0.1  !Will come from time series
      SODTA   = 1.04 !Will be model constant


      FOPO4   = 1.0  !Will be model constant
      KPO4    = 1.0  !Will be model constant
      KE      = 1.0  !Will be model constant
      KA      = 0.0D0


      RUNTIME_ERROR = 0
      GRTYPE  = 1  !Will be handeled later


C     PROCESSES

C     DETRITUS CARBON STOCHIOMETRY
      NC_EXLADDETC  = NC_EXLAPDETC
      PC_EXLADDETC  = PC_EXLAPDETC
      SiC_EXLADDETC = SiC_EXLAPDETC

      NC_EXREDDETC  = NC_EXREPDETC
      PC_EXREDDETC  = PC_EXREPDETC
      SiC_EXREDDETC = SiC_EXREPDETC


C     ZOOPLANKTON STOCHIOMETRY (ESTIMATED)
C     Assumption : Zooplankton get the average stochiometry of its food

      GRAZ_RAT_1 = ((p_1 * (PHY - FOOD_min)) / 
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz)) 

      GRAZ_RAT_2 = ((p_2 * (PHY_2 - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz)) 

      GRAZ_RAT_3 = ((p_3 * (PHY_3 - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRAZ_RAT_4 = ((p_4 * (EXLAPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRAZ_RAT_5 = ((p_5 * (GPHYPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRAZ_RAT_6 = ((p_6 * (DPHYPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRAZ_RAT_7 = ((p_7 * (CPHYPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))


      NC_ZOOC  = (GRAZ_RAT_1 * NC)   + (GRAZ_RAT_2 * NC_2) + 
     *           (GRAZ_RAT_3 * NC_3) + (GRAZ_RAT_4 * NC_EXLAPDETC) +
     *           (GRAZ_RAT_5 * NC)   + (GRAZ_RAT_6 * NC_2) +
     *           (GRAZ_RAT_7 * NC_3)


      PC_ZOOC  = (GRAZ_RAT_1 * PC)   + (GRAZ_RAT_2 * PC_2) + 
     *           (GRAZ_RAT_3 * PC_3) + (GRAZ_RAT_4 * PC_EXLAPDETC) +
     *           (GRAZ_RAT_5 * PC)   + (GRAZ_RAT_6 * PC_2) +
     *           (GRAZ_RAT_7 * PC_3)

      SiC_ZOOC = (GRAZ_RAT_2 * ASC) + (GRAZ_RAT_4 * SiC_EXLAPDETC) +
     *           (GRAZ_RAT_6 * ASC) 

C     PHYTOPLANKTON CARBON
      TPHY = PHY + PHY_2 + PHY_3 !Total phyt. carbon


      IF ((CCHLX.eq.0.0D0).or.(CCHLX_2.eq.0.0D0).or.
     *    (CCHLX_3.eq.0.0D0))then

          WRITE(6,*) 'ERROR IN SUBROUTINE WCKIND'
          WRITE(6,*) '--------------------------'
          WRITE(6,*) 'Carbon to chlorophylla ratio should not be zero'
          WRITE(6,*) 'CCHLX    : ', CCHLX
          WRITE(6,*) 'CCHLX_2  : ', CCHLX_2
          WRITE(6,*) 'CCHLX_3  : ', CCHLX_3
          STOP

      END IF
       

      TCHLA = (PHY/CCHLX) + (PHY_2/CCHLX_2) + (PHY_3/CCHLX_3)
       

      XRT   = TSTROG_ALUKAS(LTEMP, q10, tmax, topt, tamax, kamax, tamin)
      GP1_T = k1c * XRT


      IF (IA.GT.1.0D-10) THEN

C          CALL SMITH_ALUKAS(FDAY,  IA,  TCHLA ,CCHLX, GP1_T,
C     *               DEPTH, KE,  XKC   ,PHIMX, 0, LLIGHT, CCHLX0)

          CALL CUR_SMITH(FDAY,  IA,  TCHLA ,CCHLX, GP1_T,
     *               DEPTH, KE,  XKC   ,PHIMX, LLIGHT, CCHLX0)

          CCHLX = CCHLX0

      ELSE

          LLIGHT = 0.0D0

      END IF

               
      SAVED_OUTPUTS(1) = CCHLX
                     
     
      XRT_2 = TSTROG_ALUKAS(LTEMP, q10_2, tmax_2, topt_2, 
     *               tamax_2, kamax_2, tamin_2)

      GP2_T = k1c_2 * XRT_2 


      IF (IA.GT.1.0D-10) THEN

C          CALL SMITH_ALUKAS(FDAY, IA, TCHLA, CCHLX_2, GP2_T, 
C     *               DEPTH, KE, XKC, PHIMX_2, 1, LLIGHT_2,CCHLX0) 

          CALL CUR_SMITH(FDAY, IA, TCHLA, CCHLX_2, GP2_T, 
     *               DEPTH, KE, XKC, PHIMX_2, LLIGHT_2,CCHLX0)

          CCHLX_2 = CCHLX0

      ELSE

          LLIGHT_2 = 0.0D0

      END IF


      SAVED_OUTPUTS(2) = CCHLX_2


      XRT_3 = TSTROG_ALUKAS(LTEMP, q10_3, tmax_3, topt_3, 
     *               tamax_3, kamax_3, tamin_3)

      GP3_T = k1c_3 * XRT_3 


      IF (IA.GT.1.0D-10) THEN

C          CALL SMITH_ALUKAS(FDAY, IA, TCHLA, CCHLX_3, GP3_T, 
C     *               DEPTH, KE, XKC, PHIMX_3, 1, LLIGHT_3, CCHLX0) 

          CALL CUR_SMITH(FDAY, IA, TCHLA, CCHLX_3, GP3_T, 
     *               DEPTH, KE, XKC, PHIMX_3, LLIGHT_3, CCHLX0)

          CCHLX_3 = CCHLX0

      ELSE

          LLIGHT_3 = 0.0D0

      END IF


      SAVED_OUTPUTS(3) = CCHLX_3


C     Nutrient limitation
      X1   = (NH3 + NOx) / (kn + NH3 + NOx) 
      X1_2 = (NH3 + NOx) / (kn_2 + NH3 + NOx) 
      X1_3 = (NH3 + NOx) / (kn_3 + NH3 + NOx)  
      X1_3_fix = R_fix * (K_fix / (K_fix + NH3 + NOx)) 
 
      X2   = OPO4 / ((kp   / fopo4) + OPO4) 
      X2_2 = OPO4 / ((kp_2 / fopo4) + OPO4) 
      X2_3 = OPO4 / ((kp_3 / fopo4) + OPO4) 

      X3   = ISI / ((kmnsi / fisi) + ISI) 

c      X4   = INC / (INC + kc) 
c      X4_2 = INC / (INC + kc_2) 
c      X4_3 = INC / (INC + kc_3) 

      X4   = 1.0D0 
      X4_2 = 1.0D0 
      X4_3 = 1.0D0

      LNUT       = DMIN1(X1, X2, X4) 
      LNUT_2     = DMIN1(X1_2, X2_2, X3, X4_2) 
      LNUT_3     = DMIN1(X2_3, X1_3, X4_3) 
      LNUT_3_fix = DMIN1(X2_3,X1_3_fix, X4_3) 


C     Relative primary production
      GP1 = GP1_T * DMIN1(LLIGHT, LNUT) 
      GP2 = GP2_T * DMIN1(LLIGHT_2, LNUT_2) 

      GP3_1 = GP3_T * DMIN1(LLIGHT_3, LNUT_3) 
      GP3_2 = GP3_T * DMIN1(LLIGHT_3, LNUT_3_fix) 
      GP3   = GP3_1 + GP3_2 
       
      GPP   = GP1 * PHY 
      GPP_2 = GP2 * PHY_2 
      GPP_3 = GP3 * PHY_3 


C     Respiration
      RRT    = kert   ** (LTEMP - 20.0D0) 
      RRT_2  = kert_2 ** (LTEMP - 20.0D0) 
      RRT_3  = kert_3 ** (LTEMP - 20.0D0) 

      RESD   = k1rc * RRT       !dark respiration
      RESPH  = fpr  * GP1       !photo respiration
      RES    = RESD + RESPH 

      RESD    = k1rc_2 * RRT_2  !dark respiration
      RESPH   = fpr_2 * GP2     !photo respiration
      RES_2   = RESD + RESPH 


      RESD  = k1rc_3 * RRT_3     !dark respiration
      RESPH = fpr_3 * GP3        !photo respiration
      RES_3 = RESD + RESPH 

C     Excretion, mortality
      EXCR   = kexcr   * RRT
      EXCR_2 = kexcr_2 * RRT_2
      EXCR_3 = kexcr_3 * RRT_3


      DR  = k1d * (kdt ** (LTEMP - 20.0D0)) 
      DED = EXCR + DR 
      DP1 = RES + DED 


      DR_2 = k1d_2 * (kdt_2 ** (LTEMP - 20.0D0)) 
      DED2 = EXCR_2 + DR_2 
      DP2  = RES_2 + DED2 


      DR_3 = k1d_3 * (kdt_3 ** (LTEMP - 20.0D0)) 
      DED3 = EXCR_3 + DR_3
      DP3  = RES_3 + DED3 

      DEDP   = DED  * PHY 
      DEDP_2 = DED2 * PHY_2 
      DEDP_3 = DED3 * PHY_3 

      DPP    = DP1 * PHY 
      DPP_2  = DP2 * PHY_2 
      DPP_3  = DP3 * PHY_3 



      XZRT = TSTROG_ALUKAS(LTEMP, q10z, tmaxz, toptz, 
     *              tamaxz, kamaxz, taminz)



      GRZ_1 = XZRT * kgrz * ZOO * ((p_1 * (PHY - FOOD_min)) / 
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz)) 

      GRZ_2 = XZRT * kgrz * ZOO * ((p_2 * (PHY_2 - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) +
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRZ_3 = XZRT * kgrz * ZOO * ((p_3 * (PHY_3 - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRZ_4 = XZRT * kgrz * ZOO * ((p_4 * (EXLAPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRZ_5 = XZRT * kgrz * ZOO * ((p_5 * (GPHYPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRZ_6 = XZRT * kgrz * ZOO * ((p_6 * (DPHYPDETC - FOOD_min)) /
     *            ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *             (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *             (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRZ_7 = XZRT * kgrz * ZOO * ((p_7 * (CPHYPDETC - FOOD_min)) /
     *        ((p_1 * PHY) + (p_2 * PHY_2) + (p_3 * PHY_3) + 
     *         (p_4 * EXLAPDETC) + (p_5 * GPHYPDETC) + 
     *         (p_6 * DPHYPDETC) + (p_7 * CPHYPDETC) + kpz))

      GRZ = GRZ_1 + GRZ_2 + GRZ_3 + GRZ_4 + GRZ_5 + GRZ_6 + GRZ_7 


C     THE TOTAL KINETIC DERIVATIVE FOR GREEN PHYTOPLANKTON
      DERIVATIVES(4)  = GPP - DPP - GRZ_1


C     INDIVIDUAL PROCESSES FOR GREEN PHYTOPLANKTON KINETICS

C     GROSS PRIMARY PRODUCTION
      PROCESSES(4, 1) = GPP

C     RESPIRATION
      PROCESSES(4, 2) = RES * PHY

C     EXCRETION
      PROCESSES(4, 3) = EXCR * PHY

C     NON PREDATORY DEATH
      PROCESSES(4, 4) = DR * PHY

C     DEATH BY ZOOPLANKTON GRAZING
      PROCESSES(4, 5) = GRZ_1


      IF ((.NOT.(DERIVATIVES(4).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(4).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 4th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'GPP   : ', GPP
          WRITE(*,*) 'DPP   : ', DPP
          WRITE(*,*) 'GRZ_1 : ', GRZ_1

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     THE TOTAL KINETIC DERIVATIVE FOR DIATOMS
      DERIVATIVES(12) = GPP_2 - DPP_2 - GRZ_2


C     INDIVIDUAL PROCESSES FOR GREEN PHYTOPLANKTON KINETICS

C     GROSS PRIMARY PRODUCTION
      PROCESSES(12, 1) = GPP_2

C     RESPIRATION
      PROCESSES(12, 2) = RES_2 * PHY_2

C     EXCRETION
      PROCESSES(12, 3) = EXCR_2 * PHY_2

C     NON PREDATORY DEATH
      PROCESSES(12, 4) = DR_2 * PHY_2

C     DEATH BY ZOOPLANKTON GRAZING
      PROCESSES(12, 5) = GRZ_2


      IF ((.NOT.(DERIVATIVES(12).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(12).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 12th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'GPP_2 : ', GPP_2
          WRITE(*,*) 'DPP_2 : ', DPP_2
          WRITE(*,*) 'GRZ_2 : ', GRZ_2

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     THE TOTAL KINETIC DERIVATIVE FOR CYANOBACTERIA
      DERIVATIVES(13) = GPP_3 - DPP_3 - GRZ_3


C     INDIVIDUAL PROCESSES FOR GREEN PHYTOPLANKTON KINETICS

C     GROSS PRIMARY PRODUCTION
      PROCESSES(13, 1) = GPP_3

C     RESPIRATION
      PROCESSES(13, 2) = RES_3 * PHY_3

C     EXCRETION
      PROCESSES(13, 3) = EXCR_3 * PHY_3

C     NON PREDATORY DEATH
      PROCESSES(13, 4) = DR_3 * PHY_3

C     DEATH BY ZOOPLANKTON GRAZING
      PROCESSES(13, 5) = GRZ_3


      IF ((.NOT.(DERIVATIVES(13).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(13).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 13th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'GPP_3 : ', GPP_3
          WRITE(*,*) 'DPP_3 : ', DPP_3
          WRITE(*,*) 'GRZ_3 : ', GRZ_3

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     FOR ECOPATH LINKAGE
C      PROD_1(CELLNO, LAYER) = GRZ_1 + (DR * PHY)
C      PROD_2(CELLNO, LAYER) = GRZ_2 + (DR_2 * PHY_2)
C      PROD_3(CELLNO, LAYER) = GRZ_3 + (DR_3 * PHY_3)

C     FOR ECOPATH LINKAGE
C      GRAZ_1(CELLNO, LAYER) = GRZ_1
C      GRAZ_2(CELLNO, LAYER) = GRZ_2
C      GRAZ_3(CELLNO, LAYER) = GRZ_3

C     FOR ECOPATH LINKAGE
C      DEAD_1(CELLNO, LAYER) = DR * PHY
C      DEAD_2(CELLNO, LAYER) = DR_2 * PHY_2
C      DEAD_3(CELLNO, LAYER) = DR_3 * PHY_3


C     ZOOPLANKTON CARBON
      DZRT  = kdtz ** (LTEMP - 20.0D0) 
      DZ    = kdz * DZRT * ZOO 
      RZRT  = krtz ** (LTEMP - 20.0D0) 
      RZ    = krz * RZRT * ZOO 
      EXZ   = kexcrz * RZ 


      GZ    = (EFF * GRZ) - RZ - DZ - EXZ 
      DEFEC = (1.0D0 - EFF) * GRZ 
      ZSINK = DEFEC + DZ + EXZ


C     TOTAL DERIVATIVE FOR ZOOPLANKTON CARBON
      DERIVATIVES(9) = GZ


C     INDIVIDUAL PROCESSES FOR ZOOPLANKTON KINETICS

C     GAIN OF BIOMASS IN GRAZING OF PHYTOPLANKTON
      PROCESSES(9, 1)  = (EFF * GRZ_1)

C     GAIN OF BIOMASS IN GRAZING OF DIATIOMS
      PROCESSES(9, 2)  = (EFF * GRZ_2)

C     GAIN OF BIOMASS IN GRAZING OF CYANOBACTERIA
      PROCESSES(9, 3)  = (EFF * GRZ_3)

C     GAIN OF BIOMASS IN GRAZING OF EXTERNAL LABILE PARTICULATE
C     PARTICULATE DETRITUS
      PROCESSES(9, 4)  = (EFF * GRZ_4)

C     GAIN OF BIOMASS IN GRAZING OF GREENS BASED PARTICULATE DETRITUS
      PROCESSES(9, 5)  = (EFF * GRZ_5)

C     GAIN OF BIOMASS IN GRAZING OF DIATIOMS BASED PARTICULATE DETRITUS
      PROCESSES(9, 6)  = (EFF * GRZ_6)

C     GAIN OF BIOMASS IN GRAZING OF CYANOBACTERIA BASED 
C     PARTICULATE DETRITUS
      PROCESSES(9, 7)  = (EFF * GRZ_7)



C     RESPITATION
      PROCESSES(9, 8)  = RZ

C     EXCRETION
      PROCESSES(9, 9)  = EXZ

C     NON PREDATORY DEATH
      PROCESSES(9, 10) = DZ


      IF ((.NOT.(DERIVATIVES(9).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(9).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 9th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'GZ    : ', GZ
          WRITE(*,*) 'EFF   : ', EFF
          WRITE(*,*) 'RZ    : ', RZ
          WRITE(*,*) 'DZ    : ', DZ
          WRITE(*,*) 'EXZ   : ', EXZ

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     *****************************************************************
C
C      P A R T I C U L A T E  D E T R I T U S  C O M P A R T M E N T S
C      ---------------------------------------------------------------
C
C     *****************************************************************

C     DISSOLUTION OF PARTICULATE DETRITUS COMPARTMENTS
C     Zooplankton particulate detritus carbon dissolution
      DISS_ZOOPPDETC = C_DISS_ZOOPPDETC * 
     *    (T_DISS_ZOOPPDETC ** (LTEMP - T0)) * ZOOPPDETC

C     Green phytoplankton particulate detritus carbon dissolution
      DISS_GPHYPDETC = C_DISS_GPHYPDETC * 
     *    (T_DISS_GPHYPDETC ** (LTEMP - T0)) * GPHYPDETC

C     Diatom particulate detritus carbon dissolution
      DISS_DPHYPDETC = C_DISS_DPHYPDETC * 
     *    (T_DISS_DPHYPDETC ** (LTEMP - T0)) * DPHYPDETC

C     Cyanobacteria particulate detritus carbon dissolution
      DISS_CPHYPDETC = C_DISS_DPHYPDETC * 
     *    (T_DISS_CPHYPDETC ** (LTEMP - T0)) * CPHYPDETC

C     External labile particulate detritus carbon dissolution
      DISS_EXLAPDETC = C_DISS_EXLAPDETC * 
     *    (T_DISS_EXLAPDETC ** (LTEMP - T0)) * EXLAPDETC

C     External refractory particulate detritus carbon dissolution
      DISS_EXREPDETC = C_DISS_EXREPDETC * 
     *    (T_DISS_EXREPDETC ** (LTEMP - T0)) * EXREPDETC


      IF (LTEMP.LT.7.0D0) THEN

          DISS_ZOOPPDETC = 0.0D0
          DISS_GPHYPDETC = 0.0D0
          DISS_DPHYPDETC = 0.0D0
          DISS_CPHYPDETC = 0.0D0
          DISS_EXLAPDETC = 0.0D0
          DISS_EXREPDETC = 0.0D0

      END IF


C     DEATH
C     Zooplnakton death
      DEATH_ZOOC_ZOOPPDETC = DZ

C     Green phytoplankton death
      DEATH_PYHGC_GPHYPDETC = DR * PHY

C     Diatom death
      DEATH_PYHDC_DPHYPDETC = DR_2 * PHY_2

C     Cyanobacteria death
      DEATH_PYHCC_CPHYPDETC = DR_3 * PHY_3


C     TOTAL DERIVATIVE FOR EXTERNAL LABILE PARTICULATE
C     DETRITUS CARBON KINETICS
      DERIVATIVES(7)  = -1.0D0 * DISS_EXLAPDETC - GRZ_4


C     INDIVIDUAL PROCESSES FOR EXTERNAL LABILE PARTICULATE
C     DETRITUS CARBON

C     DISSOLUTION
      PROCESSES(7, 1) = DISS_EXLAPDETC

C     GRAZING BY ZOOPLANKTON
      PROCESSES(7, 2) = GRZ_4


      IF ((.NOT.(DERIVATIVES(7).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(7).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 7th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_EXLAPDETC : ', DISS_EXLAPDETC
          WRITE(*,*) 'GRZ_4          : ', GRZ_4

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR EXTERNAL LABILE PARTICULATE
C     DETRITUS CARBON KINETICS
      DERIVATIVES(11) = -1.0D0 * DISS_EXREPDETC


C     INDIVIDUAL PROCESSES FOR EXTERNAL LABILE PARTICULATE
C     DETRITUS CARBON

C     DISSOLUTION
      PROCESSES(11, 1) = DISS_EXREPDETC


      IF ((.NOT.(DERIVATIVES(11).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(11).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 11th state variable ',
     *                   'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_EXLAPDETC : ', DISS_EXREPDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR GREEN PHYTOPLANKTON BASED
C     DETRITUS CARBON KINETICS
      DERIVATIVES(16) = DEATH_PYHGC_GPHYPDETC - DISS_GPHYPDETC - GRZ_5


C     INDIVIDUAL PROCESSES FOR  GREEN PHYTOPLANKTON BASED
C     DETRITUS CARBON KINETICS

C     DEATH OF GREEN PHYTOPLANKTON
      PROCESSES(16, 1) = DEATH_PYHGC_GPHYPDETC

C     DISSOLUTION
      PROCESSES(16, 2) = DISS_GPHYPDETC

C     GRAZINF BY ZOOPLANKTON
      PROCESSES(16, 3) = GRZ_5


      IF ((.NOT.(DERIVATIVES(16).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(16).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 16th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DEATH_PYHGC_GPHYPDETC : ', DEATH_PYHGC_GPHYPDETC
          WRITE(*,*) 'DISS_GPHYPDETC        : ', DISS_GPHYPDETC
          WRITE(*,*) 'GRZ_5                 : ', GRZ_5

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR DIATIOMS BASED
C     DETRITUS CARBON KINETICS
      DERIVATIVES(18) = DEATH_PYHDC_DPHYPDETC - DISS_DPHYPDETC - GRZ_6


C     INDIVIDUAL PROCESSES FOR DIATOMS BASED
C     DETRITUS CARBON KINETICS

C     DEATH OF DIATIOMS
      PROCESSES(18, 1) = DEATH_PYHDC_DPHYPDETC

C     DISSOLUTION
      PROCESSES(18, 2) = DISS_DPHYPDETC

C     GRAZING BY ZOOPLANKTON
      PROCESSES(18, 3) = GRZ_6


      IF ((.NOT.(DERIVATIVES(18).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(18).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 18th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DEATH_PYHDC_DPHYPDETC : ', DEATH_PYHDC_DPHYPDETC
          WRITE(*,*) 'DISS_DPHYPDETC        : ', DISS_DPHYPDETC
          WRITE(*,*) 'GRZ_6                 : ', GRZ_6

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR CYANOBACTERIA BASED
C     DETRITUS CARBON KINETICS
      DERIVATIVES(20) = DEATH_PYHCC_CPHYPDETC - DISS_CPHYPDETC - GRZ_7


C     INDIVIDUAL PROCESSES FOR CYANOBACTERIA BASED
C     DETRITUS CARBON KINETICS

C     DEATH OF CYANOBACTERIA
      PROCESSES(20, 1) = DEATH_PYHCC_CPHYPDETC

C     DISSOLUTION
      PROCESSES(20, 2) = DISS_CPHYPDETC

C     GRAZNG BY ZOOPLANKTON
      PROCESSES(20, 3) = GRZ_7


      IF ((.NOT.(DERIVATIVES(20).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(20).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 20th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'

          WRITE(*,*) 'DEATH_PYHCC_CPHYPDETC : ', DEATH_PYHCC_CPHYPDETC
          WRITE(*,*) 'DISS_CPHYPDETC        : ', DISS_CPHYPDETC
          WRITE(*,*) 'GRZ_7                 : ', GRZ_7

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR ZOOPLANKTON BASED
C     DETRITUS CARBON KINETICS
      DERIVATIVES(22) = DEATH_ZOOC_ZOOPPDETC  - DISS_ZOOPPDETC


C     INDIVIDUAL PROCESSES FOR ZOOPLANKTON BASED
C     DETRITUS CARBON KINETICS

C     DEATH OF ZOOPLANKTON
      PROCESSES(22, 1) = DEATH_ZOOC_ZOOPPDETC

C     DISSOLUTION
      PROCESSES(22, 2) = DISS_ZOOPPDETC


      IF ((.NOT.(DERIVATIVES(22).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(22).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 22th state variable ',
     *                   'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DEATH_ZOOC_ZOOPPDETC : ', DEATH_PYHCC_CPHYPDETC
          WRITE(*,*) 'DISS_ZOOPPDETC       : ', DISS_ZOOPPDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     *****************************************************************
C
C        D I S S O L V E D  D E T R I T U S  C O M P A R T M E N T S 
C        -----------------------------------------------------------
C
C     *****************************************************************


C     OXIDATION OF DISSOLVED DETRITUS COMPARTMENTS
C     pH limitation for zooplankton dissloved  detritus oxidation 
      kpH_OX_ZOOPDDETC = 
     *    FKPH_ALUKAS(pH, pHminm_OX_ZOOPDDETC, pHmaxm_OX_ZOOPDDETC)

C     Zooplankton dissolved detritus carbon oxidation rate
      OX_ZOOPDDETC = C_OX_ZOOPDDETC * 
     *    (T_OX_ZOOPDDETC ** (LTEMP - T0)) * kpH_OX_ZOOPDDETC * 
     *     ZOOPDDETC * (DOO / (k_OX_ZOOPDDETC + DOO))


C     pH limitation for green phytoplankton dissolved detritus 
C     carbon oxidation 
      kpH_OX_GPHYDDETC = 
     *    FKPH_ALUKAS(pH, pHminm_OX_GPHYDDETC, pHmaxm_OX_GPHYDDETC)

C     Green phytoplankton dissolved detritus oxidation rate
      OX_GPHYDDETC = C_OX_GPHYDDETC * 
     *    (T_OX_GPHYDDETC ** (LTEMP - T0)) * kpH_OX_GPHYDDETC * 
     *     GPHYDDETC * (DOO / (k_OX_GPHYDDETC + DOO))


C     pH limitation for diatom dissloved  detritus oxidation 
      kpH_OX_DPHYDDETC = 
     *    FKPH_ALUKAS(pH, pHminm_OX_DPHYDDETC, pHmaxm_OX_DPHYDDETC)

C     Diatom dissolved detritus carbon oxidation rate
      OX_DPHYDDETC = C_OX_DPHYDDETC * 
     *    (T_OX_DPHYDDETC ** (LTEMP - T0)) * kpH_OX_DPHYDDETC * 
     *     DPHYDDETC * (DOO / (k_OX_DPHYDDETC + DOO))


C     pH limitation for cyanobacteria dissloved  detritus oxidation 
      kpH_OX_CPHYDDETC = 
     *    FKPH_ALUKAS(pH, pHminm_OX_CPHYDDETC, pHmaxm_OX_CPHYDDETC)

C     Cyanobacteria dissolved detritus carbon oxidation rate
      OX_CPHYDDETC = C_OX_CPHYDDETC * 
     *    (T_OX_CPHYDDETC ** (LTEMP - T0)) * kpH_OX_CPHYDDETC * 
     *     CPHYDDETC * (DOO / (k_OX_CPHYDDETC + DOO))

C     pH limitation for external labile dissloved detritus oxidation 
      kpH_OX_EXLADDETC = 
     *    FKPH_ALUKAS(pH, pHminm_OX_EXLADDETC, pHmaxm_OX_EXLADDETC)

C     External labile dissolved detritus carbon oxidation rate
      OX_EXLADDETC = C_OX_EXLADDETC * 
     *    (T_OX_EXLADDETC ** (LTEMP - T0)) * kpH_OX_EXLADDETC * 
     *     EXLADDETC * (DOO / (k_OX_EXLADDETC + DOO))

C    Sediment based DOC
      OX_SED_DOC = C_OX_EXLADDETC * 
     *    (T_OX_EXLADDETC ** (LTEMP - T0)) * kpH_OX_EXLADDETC * 
     *     SED_DOC * (DOO / (k_OX_EXLADDETC + DOO))

C    Sediment based DON
      OX_SED_DON = C_OX_EXLADDETC * 
     *    (T_OX_EXLADDETC ** (LTEMP - T0)) * kpH_OX_EXLADDETC * 
     *     SED_DON * (DOO / (k_OX_EXLADDETC + DOO))

C    Sediment based DOP
      OX_SED_DOP = C_OX_EXLADDETC * 
     *    (T_OX_EXLADDETC ** (LTEMP - T0)) * kpH_OX_EXLADDETC * 
     *     SED_DOP * (DOO / (k_OX_EXLADDETC + DOO))

     
C     pH limitation for external refractory dissloved detritus oxidation
      kpH_OX_EXREDDETC = 
     *    FKPH_ALUKAS(pH, pHminm_OX_EXREDDETC, pHmaxm_OX_EXREDDETC)

C     External refractory dissolved detritus carbon oxidation rate
      OX_EXREDDETC = C_OX_EXREDDETC * 
     *    (T_OX_EXREDDETC ** (LTEMP - T0)) * kpH_OX_EXREDDETC * 
     *     EXREDDETC * (DOO / (k_OX_EXREDDETC + DOO))


C     CONSUMPTION OF DISSOLVED ORGANIC MATTER BY DENITRIFICATION
C     DENITRIFICATION IS ASSUMED TO OCCUR ONLY IN VERY LOW OXYGEN 
C     CONCENTRATION IN THE WATER COLUMN
C     Denitrification
C     NIT2 = ((5.0 * 32.0) / (4.0 * 14.0)) * NIT1


C     EXCRETION
C     Zooplankton
      EXCR_ZOOPDDETC = EXZ

C     Green Phytoplankton
      EXCR_GPHYDDETC = EXCR

C     Diatom
      EXCR_DPHYDDETC = EXCR_2

C     Cyanobacteria
      EXCR_CPHYDDETC = EXCR_3

C     DERIVATIVES


C     TOTAL DERIVATIVE FOR CEXTERNAL LABILE DISSOLVED 
C     DETRITUS CARBON
      DERIVATIVES(5)  = DISS_EXLAPDETC - OX_EXLADDETC

C     SED_DOC
      DERIVATIVES(23)  = (-1.0D0) * OX_SED_DOC
      PROCESSES(23, 1) = OX_SED_DOC
      
C     SED_DON
      DERIVATIVES(24) = (-1.0D0) * OX_SED_DON
      PROCESSES(24, 1) = OX_SED_DON
      
C     SED_DOP
      DERIVATIVES(25) = (-1.0D0) * OX_SED_DOP
      PROCESSES(25, 1) = OX_SED_DOP

C     INDIVIDUAL PROCESSES FOR EXTERNAL LABILE DISSOLVED 
C     DETRITUS CARBON

C     DISSOLUTION OF EXTERNAL LABILE DISSOLVED DETRITUS CARBON
      PROCESSES(5, 1) = DISS_EXLAPDETC

C     DEGRADATION OF EXTERNAL LABILE DISSOLVED DETRITUS CARBON
      PROCESSES(5, 2) = OX_EXLADDETC


      IF ((.NOT.(DERIVATIVES(5).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(5).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 5th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_EXLAPDETC : ', DISS_EXLAPDETC
          WRITE(*,*) 'OX_EXLADDETC   : ', OX_EXLADDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

          WRITE(*,*)


          IF ((.NOT.(OX_EXLADDETC.LE.0.0D0)).AND.
     *        (.NOT.(OX_EXLADDETC.GE.0.0D0))) THEN

              WRITE(*,*) 'OX_EXLADDETC causes the problem.'
              WRITE(*,*) 
              WRITE(*,*) 'RELATED VARIABLES'
              WRITE(*,*) '-----------------'
              WRITE(*,*) 'C_OX_EXLADDETC   : ', C_OX_EXLADDETC
              WRITE(*,*) 'T_OX_EXLADDETC   : ', T_OX_EXLADDETC
              WRITE(*,*) 'kpH_OX_EXLADDETC : ', kpH_OX_EXLADDETC
              WRITE(*,*) 'EXLADDETC        : ', k_OX_EXLADDETC
              WRITE(*,*) 'DOO              : ', DOO

          END IF


      END IF


C     TOTAL DERIVATIVE FOR EXTERNAL REFRACTORY DISSOLVED 
C     DETRITUS CARBON
      DERIVATIVES(8)  = DISS_EXREPDETC - OX_EXREDDETC


C     INDIVIDUAL PROCESSES FOR EXTERNAL REFRACTORY DISSOLVED 
C     DETRITUS CARBON

C     DISSOLUTION OF EXTERNAL REFRACTORY DISSOLVED DETRITUS CARBON
      PROCESSES(8, 1) = DISS_EXREPDETC

C     DEGRADATION OF EXTERNAL REFRACTORY DISSOLVED DETRITUS CARBON
      PROCESSES(8, 2) = OX_EXREDDETC


      IF ((.NOT.(DERIVATIVES(8).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(8).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 8th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_EXREPDETC : ', DISS_EXREPDETC
          WRITE(*,*) 'OX_EXREDDETC   : ', OX_EXREDDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF



C     TOTAL DERIVATIVE FOR GREEN PHYTOPLANKTON BASED DISSOLVED 
C     DETRITUS CARBON
      DERIVATIVES(15) = DISS_GPHYPDETC + EXCR_GPHYDDETC - OX_GPHYDDETC


C     INDIVIDUAL PROCESSES FOR GREEN PHYTOPLANKTON BASED DISSOLVED 
C     DETRITUS CARBON

C     DISSOLUTION OF GREEN PHYTOPLANKTON BASED DISSOLVED DETRITUS CARBON
      PROCESSES(15, 1) = DISS_GPHYPDETC

C     EXCRETION OF GREEN PHYTOPLANKTON
      PROCESSES(15, 2) = EXCR_GPHYDDETC

C     DEGRADATION OF GREEN PHYTOPLANKTON BASED DISSOLVED DETRITUS CARBON
      PROCESSES(15, 3) = OX_GPHYDDETC


      IF ((.NOT.(DERIVATIVES(15).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(15).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 15th state variable ',
     *                   'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_GPHYPDETC : ', DISS_GPHYPDETC
          WRITE(*,*) 'EXCR_GPHYDDETC : ', EXCR_GPHYDDETC
          WRITE(*,*) 'OX_GPHYDDETC   : ', OX_GPHYDDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR DIATIOMS BASED DISSOLVED 
C     DETRITUS CARBON
      DERIVATIVES(17) = DISS_DPHYPDETC + EXCR_DPHYDDETC - OX_DPHYDDETC


C     INDIVIDUAL PROCESSES FOR DIATIOMS BASED DISSOLVED 
C     DETRITUS CARBON

C     DISSOLUTION OF DIATIOMS BASED DISSOLVED DETRITUS CARBON
      PROCESSES(17, 1) = DISS_DPHYPDETC

C     EXCRETION OF DIATIOMS
      PROCESSES(17, 2) = EXCR_DPHYDDETC

C     DEGRADATION OF DIATIOMS BASED DISSOLVED DETRITUS CARBON
      PROCESSES(17, 3) = OX_DPHYDDETC


      IF ((.NOT.(DERIVATIVES(17).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(17).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 17th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_DPHYPDETC : ', DISS_DPHYPDETC
          WRITE(*,*) 'EXCR_DPHYDDETC : ', EXCR_DPHYDDETC
          WRITE(*,*) 'OX_DPHYDDETC   : ', OX_DPHYDDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF


C     TOTAL DERIVATIVE FOR CYANOBACTERIA BASED DISSOLVED 
C     DETRITUS CARBON
      DERIVATIVES(19) = DISS_CPHYPDETC + EXCR_CPHYDDETC - OX_CPHYDDETC


C     INDIVIDUAL PROCESSES FOR CYANOBACTERIA BASED DISSOLVED 
C     DETRITUS CARBON

C     DISSOLUTION OF CYANOBACTERIA BASED DISSOLVED DETRITUS CARBON
      PROCESSES(19, 1) = DISS_CPHYPDETC

C     EXCRETION OF CYANOBACTERIA
      PROCESSES(19, 2) = EXCR_CPHYDDETC

C     DEGRADATION OF CYANOBACTERIA BASED DISSOLVED DETRITUS CARBON
      PROCESSES(19, 3) = OX_CPHYDDETC


      IF ((.NOT.(DERIVATIVES(19).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(19).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 19th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_CPHYPDETC : ', DISS_CPHYPDETC
          WRITE(*,*) 'EXCR_CPHYDDETC : ', EXCR_CPHYDDETC
          WRITE(*,*) 'OX_CPHYDDETC   : ', OX_CPHYDDETC

          RUNTIME_ERROR = 1

C         Add additional error detection code to detect the exact
C         location of error.P
          !PAUSE

      END IF


C     TOTAL DERIVATIVE FOR ZOOPLANKTON BASED DISSOLVED 
C     DETRITUS CARBON
      DERIVATIVES(21) = DISS_ZOOPPDETC + EXCR_ZOOPDDETC - OX_ZOOPDDETC


C     INDIVIDUAL PROCESSES FOR ZOOPLANKTON BASED DISSOLVED 
C     DETRITUS CARBON

C     DISSOLUTION OF ZOOPLANKTON BASED DISSOLVED DETRITUS CARBON
      PROCESSES(21, 1) = DISS_ZOOPPDETC

C     EXCRETION OF ZOOPLANKTON
      PROCESSES(21, 2) = EXCR_ZOOPDDETC

C     DEGRADATION OF ZOOPLANKTON BASED DISSOLVED DETRITUS CARBON
      PROCESSES(21, 3) = OX_ZOOPDDETC


      IF ((.NOT.(DERIVATIVES(21).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(21).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 21st state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'DISS_ZOOPPDETC : ', DISS_ZOOPPDETC
          WRITE(*,*) 'EXCR_ZOOPDDETC : ', EXCR_ZOOPDDETC
          WRITE(*,*) 'OX_ZOOPDDETC   : ', OX_ZOOPDDETC

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

      END IF




C     FOR ECOPATH LINKAGE
C      DEAD_Z(CELLNO, LAYER) = DZ


C     AMMONIUM NITROGEN
      kpHn = FKPH_ALUKAS(pH, pHminn, pHmaxn) 


C     Nitrification      
	N1 = kcnit * (ktnit ** (LTEMP - t0)) * 
     *     (DOO / (knit+DOO)) * kpHn * NH3


C     Uptake preference factor for Greens
      PN = (NH3 * NOx) / ((kn + NH3) * (kn + NOx)) + 
     *      (kn * NH3) / ((NH3 + NOx) * (kn + NOx))
 

      IF (NH3.LT.EPS) THEN

          PN = 0.0D0

      END IF


C     Uptake preference factor for Diatoms
      PN_2 = (NH3 * NOx) / ((kn_2 + NH3) * (kn_2 + NOx)) + 
     *       (kn_2 * NH3) / ((NH3 + NOx) * (kn_2 + NOx)) 

      IF (NH3 .lt. eps) THEN

          PN_2 = 0.0D0

      END IF


C     Uptake preference factor for Cyanobacteria
      PN_3 = (NH3 * NOx) / ((kn_3 + NH3) * (kn_3 + NOx)) + 
     *       (kn_3 * NH3) / ((NH3 + NOx) * (kn_3 + NOx)) 
    
    
      IF (NH3 .lt. eps) THEN

          PN_3 = 0.0D0

      END IF


C     Uptake of Ammonium Nitrogen
      NALG2 = (nc * PN * GPP) + (nc_2 * PN_2 * GPP_2) + 
     *        (nc_3 * PN_3 * GP3_1 * PHY_3)


      IF (NH3.LT.eps) THEN

          NALG2 = 0.0D0

      END IF


C     Dissolved detritus oxidation
      N_OX_ZOOPDDETC = OX_ZOOPDDETC * NC_ZOOC
      N_OX_GPHYDDETC = OX_GPHYDDETC * NC
      N_OX_DPHYDDETC = OX_DPHYDDETC * NC_2
      N_OX_CPHYDDETC = OX_CPHYDDETC * NC_3
      N_OX_EXLADDETC = OX_EXLADDETC * NC_EXLADDETC
      N_OX_EXREDDETC = OX_EXREDDETC * NC_EXREDDETC


C     Respiration
      N_RESP_ZOOC  = RZ * NC_ZOOC
      N_RESP_PHYGC = RES * PHY * NC
      N_RESP_PHYDC = RES_2 * PHY_2 * NC_2
      N_RESP_PHYCC = RES_3 * PHY_3 * NC_3


C     DERIVATIVE


C     TOTAL DERIVATIVE FOR AMMONIA NITROGEN
      DERIVATIVES(1) = N_RESP_ZOOC + N_RESP_PHYGC + N_RESP_PHYDC + 
     *    N_RESP_PHYCC + N_OX_ZOOPDDETC + N_OX_GPHYDDETC +
     *    N_OX_DPHYDDETC + N_OX_CPHYDDETC + N_OX_EXLADDETC + 
     *    N_OX_EXREDDETC - NALG2 - N1 + OX_SED_DON


C     INDIVIDUAL PROCESSES FOR AMMONIA NITROGEN

C     ZOOPLANKTON RESPIRATION
      PROCESSES(1, 1) = N_RESP_ZOOC

C     GREEN PHYTOPLANKTON RESPIRATION
      PROCESSES(1, 2) = N_RESP_PHYGC

C     DIATOMS RESPIRATION
      PROCESSES(1, 3) = N_RESP_PHYDC

C     CYANOBACTERIA RESPIRATION
      PROCESSES(1, 4) = N_RESP_PHYCC

C     DEGRADATION OF ZOOPLANKTON BASED DISSOLVED DETRITUS
      PROCESSES(1, 5) = N_OX_ZOOPDDETC

C     DEGRADATION OF GREEN PHYTOPLANKTON BASED DISSOLVED DETRITUS
      PROCESSES(1, 6) = N_OX_GPHYDDETC

C     DEGRADATION OF DIATOMS BASED DISSOLVED DETRITUS
      PROCESSES(1, 7) = N_OX_DPHYDDETC

C     DEGRADATION OF CYANOBACTERIA BASED DISSOLVED DETRITUS
      PROCESSES(1, 8) = N_OX_CPHYDDETC

C     DEGRADATION OF EXTERNAL LABILE DISSOLVED DETRITUS
      PROCESSES(1, 9) = N_OX_EXLADDETC

C     DEGRADATION OF EXTERNAL REFRACTORY DISSOLVED DETRITUS
      PROCESSES(1, 10) = N_OX_EXREDDETC

C     NITROGEN UPTAKE OF GREEN PHYTOPLANKTON
      PROCESSES(1, 11) = nc * PN * GPP

C     NITROGEN UPTAKE OF DIATOMS
      PROCESSES(1, 12) = nc_2 * PN_2 * GPP_2

C     NITROGEN UPTAKE OF CYANOBACTERIA
      PROCESSES(1, 13) = nc_3 * PN_3 * GP3_1 * PHY_3

C     NITRIFICATION
      PROCESSES(1, 14) = N1

C     DEGRADATION OF SEDIMENT BASED DON 
      PROCESSES(1, 15) = OX_SED_DON


C     NITRATE NITROGEN

C     Uptake
      NOALG = ((1.0D0 - PN) * nc * GPP) + 
     *        ((1.0D0 - PN_2) * nc_2 * GPP_2) +
     *        ((1.0D0 - PN_3) * nc_3 * GP3_1 * PHY_3)


      IF (NOx .lt. eps) THEN
          NOALG = 0.0D0
      END IF


      kpHdn = FKPH_ALUKAS(pH, pHmindn, pHmaxdn) 


C     Denitrification
      NIT1 = kcdenit * (ktdenit ** (LTEMP - t0)) *
     *       (kdenit / (kdenit + DOO)) * kpHdn * NOx


C     TOTAL DERIVATIVE FOR NITRATE NITROGEN KINETICS
      DERIVATIVES(2)= N1-NOALG-NIT1


C     INDIVIDUAL PROCESSES FOR NITRATE NITROGEN

C     NITRIFICATION
      PROCESSES(2, 1) = N1

C     GREEN PHYTOPLANKTON UPTAKE
      PROCESSES(2, 2) = (1.0D0 - PN) * nc * GPP

C     DIATOMS UPTAKE
      PROCESSES(2, 3) = (1.0D0 - PN_2) * nc_2 * GPP_2

C     CYANOBACTERIA UPTAKE
      PROCESSES(2, 4) = (1.0D0 - PN_3) * nc_3 * GP3_1 * PHY_3

C     DENITRIFICATION
      PROCESSES(2, 5) = NIT1


C     PHOSPHATE PHOSPHORUS
C     Dissolved detritus oxidation
      P_OX_ZOOPDDETC = OX_ZOOPDDETC * PC_ZOOC
      P_OX_GPHYDDETC = OX_GPHYDDETC * PC
      P_OX_DPHYDDETC = OX_DPHYDDETC * PC_2
      P_OX_CPHYDDETC = OX_CPHYDDETC * PC_3
      P_OX_EXLADDETC = OX_EXLADDETC * PC_EXLADDETC
      P_OX_EXREDDETC = OX_EXREDDETC * PC_EXREDDETC


C     Respiration
      P_RESP_ZOOC  = RZ * PC_ZOOC
      P_RESP_PHYGC = RES * PHY * PC
      P_RESP_PHYDC = RES_2 * PHY_2 * PC_2
      P_RESP_PHYCC = RES_3 * PHY_3 * PC_3

C     Uptake
      OPALG2 = (pc * GPP) + (pc_2 * GPP_2) + (pc_3 * GPP_3)


C     TOTAL DERIVATIVE FOR PHOSPHATE PHOSPHORUS
      DERIVATIVES(3) = P_RESP_ZOOC + P_RESP_PHYGC + P_RESP_PHYDC + 
     *    P_RESP_PHYCC + P_OX_ZOOPDDETC + P_OX_GPHYDDETC +
     *    P_OX_DPHYDDETC + P_OX_CPHYDDETC + P_OX_EXLADDETC + 
     *    P_OX_EXREDDETC - OPALG2 + OX_SED_DOP


C     INDIVIDUAL PROCESSES FOR PHOSPHATE PHOSPHORUS KINETICS

C     ZOOPLANKTON RESPIRATION
      PROCESSES(3, 1) = P_RESP_ZOOC

C     GREEN PHYTOPLANKTON RESPIRATION
      PROCESSES(3, 2) = P_RESP_PHYGC

C     DIATOMS RESPIRATION
      PROCESSES(3, 3) = P_RESP_PHYDC

C     CYANOBACTERIA RESPIRATION
      PROCESSES(3, 4) = P_RESP_PHYCC

C     DEGRADATION OF ZOOPLANKTON BASED DISSOLVED DETRITUS
      PROCESSES(3, 5) = P_OX_ZOOPDDETC

C     DEGRADATION OF GREEN PHYTOPLANKTON BASED DISSOLVED DETRITUS
      PROCESSES(3, 6) = P_OX_GPHYDDETC

C     DEGRADATION OF DIATOMS BASED DISSOLVED DETRITUS
      PROCESSES(3, 7) = P_OX_DPHYDDETC

C     DEGRADATION OF CYANOBACTERIA BASED DISSOLVED DETRITUS
      PROCESSES(3, 8) = P_OX_CPHYDDETC

C     DEGRADATION OF EXTERNAL LABILE DISSOLVED DETRITUS
      PROCESSES(3, 9) = P_OX_EXLADDETC

C     DEGRADATION OF EXTERNAL REFRACTORY DISSOLVED DETRITUS
      PROCESSES(3, 10) = P_OX_EXREDDETC

C     PHOSPHORUS UPTAKE OF GREEN PHYTOPLANKTON
      PROCESSES(3, 11) = pc * GPP

C     PHOSPHORUS UPTAKE OF DIATOMS
      PROCESSES(3, 12) = pc_2 * GPP_2

C     PHOSPHORUS UPTAKE OF CYANOBACTERIA
      PROCESSES(3, 13) = pc_3 * GPP_3

C     DEGRADATION OF SEDIMENT BASED DISSOLVED ORGANIC PHOSPHORUS
      PROCESSES(3, 14) = OX_SED_DOP


C     DISSOLVED OXYGEN 

C     If surface element and not completly covered with ice
      IF ((TOPSEG.EQ.1.0D0) .AND. (XICECVR .GT. 0.0D0)) THEN

 
          IF (K2 .EQ. 0.0) THEN

              K2WIND = KAWIND_ALUKAS(WIND , LTEMP, AIRTMP,  DEPTH,
     +         WTYPE) 
              K2HYDR = KAHYDR_ALUKAS(DEPTH, VEL  , STP20)


              IF (K2WIND.GT.K2HYDR) THEN

                  KA = K2WIND * XICECVR

              ELSE

                  KA = K2HYDR * XICECVR

              END IF


          ELSE


              IF (K2 .GT. 0)THEN

                  KA = ((K2 * 1.028**STP20) * XICECVR)

              ELSE

                  KA = 0.0D0

              ENDIF


          END IF


      ELSE

          KA = 0.0D0

      END IF


C     Oxygen saturation concentration
      OSAT = DO_SAT_ALUKAS(LTEMP, SALIN, SURF_ELEVATION)

      if (OSAT.gt.25.0D0) then

         write(*,*) 'NODE  : ', CELLNO
         write(*,*) 'TIME  : ', PSTIME
         write(*,*) 'TEMP  : ', LTEMP
         write(*,*) 'SALIN : ', SALIN
         write(*,*) 'ELEV  : ', SURF_ELEVATION
         write(*,*) 'OSAT  : ', OSAT 
         write(*,*) 'KA    : ', KA

         RUNTIME_ERROR = 1
         !PAUSE

      end if

C     Aeration
      DO1 = KA * (OSAT - DOO)

C     Photosynthesis using NH3
      DO2 = (OC_GPHYDDETC * PN * GP1 * PHY) + 
     *      (OC_DPHYDDETC * PN_2 * GP2 * PHY_2) + 
     *      (OC_CPHYDDETC * PN_3 * GP3 * PHY_3)

C     Photosynthesis using NOx
      DO3 = (OC_GPHYDDETC * (1.0D0 - PN) * GP1 * PHY * 
     *       (1.0D0 + 32.0D0 * 1.5D0 * nc / 14.0D0)) +
     *      (OC_DPHYDDETC * (1.0D0 - PN_2) * GP2 * PHY_2 * 
     *       (1.0D0 + 32.0D0 * 1.5D0 * nc_2 / 14.0D0)) + 
     *      (OC_CPHYDDETC * (1.0D0 - PN_3) * GP3 * PHY_3 * 
     *       (1.0D0 + 32.0D0 * 1.5D0 * nc_3 / 14.0D0))

C     Respiration
      OX_RESP_ZOOC  = RZ * OC_ZOOPDDETC
      OX_RESP_PHYGC = RES * PHY * OC_GPHYDDETC
      OX_RESP_PHYDC = RES_2 * PHY_2 * OC_DPHYDDETC
      OX_RESP_PHYCC = RES_3 * PHY_3 * OC_CPHYDDETC


      DO4 = OX_RESP_ZOOC + OX_RESP_PHYGC + OX_RESP_PHYDC + 
     *      OX_RESP_PHYCC

C     Nitrification
      N2 = (64.0D0 / 14.0D0) * N1


C     Dissolved organic matter oxidation 
      OX = (OX_ZOOPDDETC * OC_ZOOPDDETC) + 
     *     (OX_GPHYDDETC * OC_GPHYDDETC) + 
     *     (OX_DPHYDDETC * OC_DPHYDDETC) + 
     *     (OX_CPHYDDETC * OC_CPHYDDETC) + 
     *     (OX_EXLADDETC * OC_EXLADDETC) + 
     *     (OX_EXREDDETC * OC_EXREDDETC) +
     *     (OX_SED_DOC   * OC_EXLADDETC)

 


C     TOTAL DERIVATIVE FOR DISSOLVED OXYGEN
      DERIVATIVES(6) = DO1 + DO2 + DO3 - DO4 - N2 - OX


C     INDIVIDUAL PROCESSES FOR DISSOLVED OXYGEN KINETICS

C     AERATION
      PROCESSES(6, 1) = DO1

C     PHOTOSNYTHESIS BY GREEN PHYTOPLANKTON - AMMONIA UPTAKE
      PROCESSES(6, 2) = OC_GPHYDDETC * PN * GP1 * PHY

C     PHOTOSNYTHESIS BY DIATOMS - AMMONIA UPTAKE
      PROCESSES(6, 3) = OC_DPHYDDETC * PN_2 * GP2 * PHY_2

C     PHOTOSNYTHESIS BY CYANOBACTERIA - AMMONIA UPTAKE
      PROCESSES(6, 4) = OC_CPHYDDETC * PN_3 * GP3 * PHY_3

C     PHOTOSNYTHESIS BY GREEN PHYTOPLANKTON - NITRATE UPTAKE
      PROCESSES(6, 5) = OC_GPHYDDETC * (1.0D0 - PN) *  
     *       GP1 * PHY * (1.0D0 + 32.0D0 * 1.5D0 * nc / 14.0D0)

C     PHOTOSNYTHESIS BY DIATOMS - NITRATE UPTAKE
      PROCESSES(6, 6) = OC_DPHYDDETC * (1.0D0 - PN_2) *  
     *       GP2 * PHY_2 * (1.0D0 + 32.0D0 * 1.5D0 * nc_2 / 14.0D0)

C     PHOTOSNYTHESIS BY CYANOBACTERIA - NITRATE UPTAKE
      PROCESSES(6, 7) = OC_CPHYDDETC * (1.0D0 - PN_3) *  
     *       GP3 * PHY_3 * (1.0D0 + 32.0D0 * 1.5D0 * nc_3 / 14.0D0)

C     ZOOPLANKTON RESPIRATION
      PROCESSES(6, 8) = OX_RESP_ZOOC

C     GREEN PHYTOPLANKTON RESPIRATION
      PROCESSES(6, 9) = OX_RESP_PHYGC

C     DIATOMS RESPIRATION
      PROCESSES(6, 10) = OX_RESP_PHYDC

C     CYANOBACTERIA RESPIRATION
      PROCESSES(6, 11) = OX_RESP_PHYCC

C     NITRIFICATION
      PROCESSES(6, 12) = N2

C     ZOOPLANKTON BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(6, 13) = OX_ZOOPDDETC * OC_ZOOPDDETC

C     GREEN PHYTOPLANKTON BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(6, 14) = OX_GPHYDDETC * OC_GPHYDDETC

C     DIATOM BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(6, 15) = OX_DPHYDDETC * OC_DPHYDDETC

C     CYANOBACTERIA BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(6, 16) = OX_CPHYDDETC * OC_CPHYDDETC 

C     EXTERNAL LABILE DISSOLVED DETRITUS DEGRADATION
      PROCESSES(6, 17) = OX_EXLADDETC * OC_EXLADDETC

C     EXTERNAL RERACTORY DISSOLVED DETRITUS DEGRADATION
      PROCESSES(6, 18) = OX_EXREDDETC * OC_EXREDDETC

C     DEGREDATION (OXIDATION) OF SEDIMENT BASED DISSOLVED ORGANIC CARBON
      PROCESSES(6, 19) = OX_SED_DOC * OC_EXLADDETC

C     AVAILABLE SILICA
C     Dissolved Detritus oxidation
      Si_OX_ZOOPDDETC = OX_ZOOPDDETC * SiC_ZOOC
      Si_OX_DPHYDDETC = OX_DPHYDDETC * asc
      Si_OX_EXLADDETC = OX_EXLADDETC * SiC_EXLADDETC 
      Si_OX_EXREDDETC = OX_EXREDDETC * SiC_EXREDDETC


C     Respiration
      Si_RESP_ZOOC  = RZ * SiC_ZOOC
      Si_RESP_PHYDC = RES_2 * asc


C     Uptake by diatoms
      UISI = asc * GP2 * PHY_2


C     TOTAL DERIVATIVE FOR AVAILABLE SILICON
      DERIVATIVES(10) = Si_OX_ZOOPDDETC + Si_OX_DPHYDDETC + 
     *             Si_OX_EXLADDETC + Si_OX_EXREDDETC +
     *             Si_RESP_ZOOC + Si_RESP_PHYDC - UISI


C     INDIVIDUAL PROCESSES FOR AVAILABLE SILICON

C     ZOOPLANKTON BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(10, 1) = Si_OX_ZOOPDDETC

C     DIATOMS BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(10, 2) = Si_OX_DPHYDDETC

C     EXTERNAL LABILE DISSOLVED DETRITUS DEGRADATION
      PROCESSES(10, 3) = Si_OX_EXLADDETC

C     EXTERNAL REFEACTORY DISSOLVED DETRITUS DEGRADATION
      PROCESSES(10, 4) = Si_OX_EXREDDETC

C     ZOOPLANKTON RESPIRATION
      PROCESSES(10, 5) = Si_RESP_ZOOC

C     DIATOMS RESPIRATION
      PROCESSES(10, 6) = Si_RESP_PHYDC

C     UPTAKE BY DIATOMS
      PROCESSES(10, 7) = UISI



C     INORGANIC CARBON (CO2)
      Resp = (RES * PHY) + (RES_2 * PHY_2) + (RES_3 * PHY_3) + RZ 

      Kliq_co2 = KA * ((32.0D0 / 44.0D0) ** 0.25D0) 

      power10 = (2385.73 / (LTEMP + 273.15)) - 14.0184 +
     *           (0.0152642 * (LTEMP + 273.15)) 

      co2sat =  0.00035 * 44.0D3 * (10.0D0 ** power10)
 
      atmexch = (12.0D0 / 44.0D0) * 
     *          (Kliq_co2 * (co2sat - (44.0D0 / 12.0D0) * INC)) 

      phot = GPP + GPP_2 + GPP_3 


C     TOTAL DERIVATIVE FOR INORGANIC CARBON
      DERIVATIVES(14) = Resp + OX_ZOOPDDETC + OX_GPHYDDETC + 
     *     OX_DPHYDDETC + OX_CPHYDDETC + OX_EXLADDETC + 
     *     OX_EXREDDETC + atmexch - phot + OX_SED_DOC


      IF ((.NOT.(DERIVATIVES(14).LE.0.0D0)).AND.
     *    (.NOT.(DERIVATIVES(14).GE.0.0D0))) THEN

          WRITE(*,*) 'ERROR IN FUNCTION WCKIND'
          WRITE(*,*) '------------------------'
          WRITE(*,*) 'PSTIME  : ', PSTIME
          WRITE(*,*) 'CELL NO : ', CELLNO
          WRITE(*,*) 'LAYER   : ', LAYER

          WRITE(*,*) 'Derivative for the 14th state variable ',
     *               'is not a number.'

          WRITE(*,*) 'RELATED VARIABLES'
          WRITE(*,*) '-----------------'
          WRITE(*,*) 'Resp : ', Resp
          WRITE(*,*) 'OX_ZOOPDDETC   : ', OX_ZOOPDDETC
          WRITE(*,*) 'OX_GPHYDDETC   : ', OX_GPHYDDETC
          WRITE(*,*) 'OX_DPHYDDETC   : ', OX_DPHYDDETC
          WRITE(*,*) 'OX_CPHYDDETC   : ', OX_CPHYDDETC
          WRITE(*,*) 'OX_EXLADDETC   : ', OX_EXLADDETC
          WRITE(*,*) 'OX_EXREDDETC   : ', OX_EXREDDETC
          WRITE(*,*) 'atmexch        : ', atmexch
          WRITE(*,*) 'phot           : ', phot

          RUNTIME_ERROR = 1
          !PAUSE

C         Add additional error detection code to detect the exact
C         location of error.

          IF ((.NOT.(atmexch.LE.0.0D0)).AND.
     *        (.NOT.(atmexch.GE.0.0D0))) THEN

              WRITE(*,*)
              WRITE(*,*)

              WRITE(*,*) 'PROCESS "atmech" casued not a number error'


              WRITE(*,*) 'RELATED VARIABLES'
              WRITE(*,*) '-----------------'
              WRITE(*,*) 'Kliq_co2 : ', Kliq_co2
              WRITE(*,*) 'co2sat   : ', co2sat
              WRITE(*,*) 'INC      : ', INC

          END IF


      END IF



C     INDIVIDUAL PROCESSES FOR INORGANIC CARBON KINETICS

C     ZOOPLANKTON RESPIRATION
      PROCESSES(14, 1) = RES * PHY

C     GREEN PHYTOPLANKTON RESPIRATION
      PROCESSES(14, 2) = RES_2 * PHY_2

C     DIATOMS RESPIRATION
      PROCESSES(14, 3) = RES_3 * PHY_3

C     CYANOBACTERIA RESPIRATION
      PROCESSES(14, 4) = RZ

C     ZOOPLANKTON BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(14, 5) = OX_ZOOPDDETC

C     GREEN PHYTOPLANKTON BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(14, 6) = OX_GPHYDDETC

C     DIATOMS DISSOLVED DETRITUS DEGRADATION
      PROCESSES(14, 7) = OX_DPHYDDETC

C     CYANOBACTERIA BASED DISSOLVED DETRITUS DEGRADATION
      PROCESSES(14, 8) = OX_DPHYDDETC

C     ATMOSPHERIC EXCHANGE
      PROCESSES(14, 9) = atmexch

C     GREEN PHYTOPLANKTON PHOTOSYNTHESIS
      PROCESSES(14, 10) = GPP

C     DIATOMS PHOTOSYNTHESIS
      PROCESSES(14, 11) = GPP_2

C     CYANOBACTERIA PHOTOSYNTHESIS
      PROCESSES(14, 12) = GPP_3

C     DEGRADATION (OXIDATION) OF SEDIMENT BASED DISSOLVED ORGANIC CARBON
      PROCESSES(14, 13) = OX_SED_DOC

      SAVED_OUTPUTS(4) = NC_ZOOC
      SAVED_OUTPUTS(5) = PC_ZOOC
      SAVED_OUTPUTS(6) = SiC_ZOOC


      IF (RUNTIME_ERROR.EQ.1) THEN

          WRITE(*,*) 'One or more kinetic derivative(s) caused ',
     *               'problems. Program halted.'

          STOP

      END IF


      END
c
c**********************************************************************************
c

C     USER SUBROUTINES FOR KINETIC MODEL
      SUBROUTINE SMITH_ALUKAS(PHOTO, Ia, TCHLA, CCHLXI, 
     *                 GITMAX, H, ke, XKC,  			!ggu
     *                 PHIMX, ERR_TOLARATE, LLIGHT, CCHLX) 

C         PHOTO  : Photo period
C         Ia     : Avaiable light at surface
C         TCHLA  : Chlorophyll-A
C         CCHLXI : Carbon to Chlorophyll ratio
C         GITMAX : 
C         LAYER  : NO
C         KE     : Background light extinction coefficient (1/m)
C         XKC    : Chlorophyll light extinction coefficient (1/m)
C         PHIMAX : Max. Quantum Yield
       
          IMPLICIT NONE

          DOUBLE PRECISION PHOTO
          DOUBLE PRECISION Ia
          DOUBLE PRECISION TCHLA
          DOUBLE PRECISION CCHLXI
          DOUBLE PRECISION GITMAX
          DOUBLE PRECISION H
          DOUBLE PRECISION ke
          DOUBLE PRECISION XKC
          DOUBLE PRECISION PHIMX
          INTEGER ERR_TOLARATE

C         LLIGHT : Calculated light limitation
C         CCHLX  : Calculated Carbon to Chlorophyll ratio 

          DOUBLE PRECISION LLIGHT
          DOUBLE PRECISION CCHLX


          DOUBLE PRECISION FDAY		! ggu
          DOUBLE PRECISION ITOT
          DOUBLE PRECISION CCHL1
          DOUBLE PRECISION KESHD
          DOUBLE PRECISION SKE
          DOUBLE PRECISION TEMP1
          DOUBLE PRECISION TEMP2
          DOUBLE PRECISION TEMP3
          DOUBLE PRECISION IMAX
          DOUBLE PRECISION SUM
          DOUBLE PRECISION DTDAY
          DOUBLE PRECISION I0
          DOUBLE PRECISION RLIGHT
          DOUBLE PRECISION IAV
          DOUBLE PRECISION IAVSG


          DOUBLE PRECISION PI
          INTEGER I

c          print *,'FDAY   : ', FDAY
c          print *,'IA     : ', IA
c          print *,'TCHLA  : ', TCHLA
c          print *,'CCHLX  : ', CCHLX
c          print *,'GP1_T  : ', GP1_T
c          print *,'DEPTH  : ', DEPTH
c          print *,'KE     : ', KE
c          print *,'XKC    : ', XKC
c          print *,'PHIMX  : ', PHIMX
c          print *,'LLIGHT : ', LLIGHT
c          print *,'CCHLX0 : ', CCHLX0

           print *,'PHOTO  : ', PHOTO
           print *,'IA     : ', IA
           print *,'TCHLA  : ', TCHLA
           print *,'CCHLXI : ', CCHLXI
           print *,'GITMAX : ', GITMAX
           print *,'H      : ', H
           print *,'KE     : ', KE
           print *,'XKC    : ', XKC
           print *,'PHIMX  : ', PHIMX
           print *,'LLIGHT : ', LLIGHT
           print *,'CCHLX  : ', CCHLX


          PI = 3.14159

          FDAY = PHOTO
          ITOT = Ia
          CCHL1 = CCHLXI
 
 
C         Dick Smith formulation integrated every day
          KESHD = XKC * 1000.0 * TCHLA
          SKE = ke
          SKE = SKE + KESHD
          TEMP1 = SKE * H
          TEMP2 = 0.083 * PHIMX * XKC / (GITMAX * CCHL1 * 2.718)

          print *, 'Location 1 in SMITH'

          IF ((.NOT.(TEMP1.LE.0.0D0)).AND.(.NOT.(TEMP1.GE.0.0D0))) THEN


              IF (ERR_TOLARATE.EQ.0) THEN

                  WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
                  WRITE(*,*) '-------------------------'
                  WRITE(*,*) 'TEMP1 IS NOT A NUMBER'
                  WRITE(*,*)
                  WRITE(*,*) 'VALUES OF RELATED VARIABLES'
                  WRITE(*,*) '---------------------------'
                  WRITE(*,*) 'SMITH_ALUKAS, XKC   : ', XKC
                  WRITE(*,*) 'SMITH_ALUKAS, TCHLA : ', TCHLA
                  WRITE(*,*) 'SMITH_ALUKAS, KESHD : ', KESHD
                  WRITE(*,*) 'SMITH_ALUKAS, KE    : ', KE
                  WRITE(*,*) 'SMITH_ALUKAS, SKE   : ', SKE
                  WRITE(*,*) 'SMITH_ALUKAS, TEMP1 : ', TEMP1
                  STOP

              ELSE

                  LLIGHT = 1.0D-3
                  CCHLX  = 3.0D1
                  GOTO 9999

              END IF


          END IF

          print *, 'Location 2 in SMITH'

          IF (TEMP1.GE.100.0) THEN


              IF (ERR_TOLARATE.EQ.0) THEN

                  WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
                  WRITE(*,*) '-------------------------'
                  WRITE(*,*) 'UNREALISTIC VALUE FOR TEMP1'
                  WRITE(*,*)
                  WRITE(*,*) 'VALUES OF RELATED VARIABLES'
                  WRITE(*,*) '---------------------------'
                  WRITE(*,*) 'SMITH_ALUKAS, XKC   : ', XKC
                  WRITE(*,*) 'SMITH_ALUKAS, TCHLA : ', TCHLA
                  WRITE(*,*) 'SMITH_ALUKAS, KESHD : ', KESHD
                  WRITE(*,*) 'SMITH_ALUKAS, KE    : ', KE
                  WRITE(*,*) 'SMITH_ALUKAS, SKE   : ', SKE
                  WRITE(*,*) 'SMITH_ALUKAS, TEMP1 : ', TEMP1
                  STOP

              ELSE

                  LLIGHT = 1.0D-3
                  CCHLX  = 3.0D1
                  GOTO 9999

              END IF


          END IF


          print *, 'Location 3 in SMITH'


          IF (TEMP1.EQ.0.0D0) THEN


              IF (ERR_TOLARATE.EQ.0) THEN

                  WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
                  WRITE(*,*) '-------------------------'
                  WRITE(*,*) 'TEMP1 SHOULD NOT BE ZERO'
                  WRITE(*,*)
                  WRITE(*,*) 'VALUES OF RELATED VARIABLES'
                  WRITE(*,*) '---------------------------'
                  WRITE(*,*) 'SMITH_ALUKAS, XKC   : ', XKC
                  WRITE(*,*) 'SMITH_ALUKAS, TCHLA : ', TCHLA
                  WRITE(*,*) 'SMITH_ALUKAS, KESHD : ', KESHD
                  WRITE(*,*) 'SMITH_ALUKAS, KE    : ', KE
                  WRITE(*,*) 'SMITH_ALUKAS, SKE   : ', SKE
                  WRITE(*,*) 'SMITH_ALUKAS, TEMP1 : ', TEMP1
                  STOP

              ELSE

                  LLIGHT = 1.0D-3
                  CCHLX  = 3.0D1
                  GOTO 9999

              END IF


          END IF


          print *, 'Location 4 in SMITH'


          IF ((.NOT.(TEMP2.LE.0.0D0)).AND.(.NOT.(TEMP2.GE.0.0D0))) THEN
    

              IF (ERR_TOLARATE.EQ.0) THEN

                  WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
                  WRITE(*,*) '-------------------------'
                  WRITE(*,*) 'TEMP2 IS NOT A NUMBER'
                  WRITE(*,*)
                  WRITE(*,*) 'VALUES OF RELATED VARIABLES'
                  WRITE(*,*) '---------------------------'
                  WRITE(*,*) 'SMITH_ALUKAS, CCHL1  : ', CCHL1
                  WRITE(*,*) 'SMITH_ALUKAS, GITMAX : ', GITMAX
                  WRITE(*,*) 'SMITH_ALUKAS, PHIMX  : ', PHIMX
                  STOP

              ELSE

                  LLIGHT = 1.0D-3
                  CCHLX  = 3.0D1
                  GOTO 9999

              END IF


          END IF

          print *, 'Location 5 in SMITH'


c          IF (TEMP2.LT.-1.0D0) THEN

c              WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
c              WRITE(*,*) '-------------------------'
c              WRITE(*,*) 'UNREALISTIC VALUE FOR TEMP2'
c              WRITE(*,*)
c              WRITE(*,*) 'VALUES OF RELATED VARIABLES'
c              WRITE(*,*) '---------------------------'
c              WRITE(*,*) 'SMITH_ALUKAS, CCHL1  : ', CCHL1
c              WRITE(*,*) 'SMITH_ALUKAS, GITMAX : ', GITMAX
c              WRITE(*,*) 'SMITH_ALUKAS, PHIMX  : ', PHIMX
c              STOP

c          END IF


          TEMP3 = EXP(-TEMP1)

          print *, 'Location 6 in SMITH'

          IF ((.NOT.(TEMP3.LE.0.0D0)).AND.(.NOT.(TEMP3.GE.0.0D0))) THEN


              IF (ERR_TOLARATE.EQ.0) THEN

                  WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
                  WRITE(*,*) '-------------------------'
                  WRITE(*,*) 'TEMP3 IS NOT A NUMBER'
                  WRITE(*,*)
                  WRITE(*,*) 'VALUES OF RELATED VARIABLES'
                  WRITE(*,*) '---------------------------'
                  WRITE(*,*) 'SMITH_ALUKAS, TEMP1 : ', TEMP1
                  WRITE(*,*)
                  WRITE(*,*) 'VALUES RELATED TO TEMP1'
                  WRITE(*,*) '---------------------------'
                  WRITE(*,*) 'SMITH_ALUKAS, XKC   : ', XKC
                  WRITE(*,*) 'SMITH_ALUKAS, TCHLA : ', TCHLA
                  WRITE(*,*) 'SMITH_ALUKAS, KESHD : ', KESHD
                  WRITE(*,*) 'SMITH_ALUKAS, KE    : ', KE
                  WRITE(*,*) 'SMITH_ALUKAS, SKE   : ', SKE
                  STOP

              ELSE

                  LLIGHT = 1.0D-3
                  CCHLX  = 3.0D1
                  GOTO 9999

              END IF


          END IF


          print *, 'Location 7 in SMITH'


          IMAX = PI * ITOT / (2.0D0 * FDAY)
          SUM = 0.0

          print *, 'Location 8 in SMITH'

          DO 1010 I = 1,25

              DTDAY = (I - 1) / 24.0


              IF (DTDAY.LE.FDAY) THEN 

                  I0 = IMAX * SIN(PI * DTDAY / FDAY)


c                  IF (((-TEMP2 * I0 * TEMP3).GE.100.0).OR.
c     *                ((-TEMP2 * I0 * TEMP3).LE.-100.0)) THEN

c                      WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
c                      WRITE(*,*) '-------------------------'

c                      WRITE(*,*) 'UNREALISTIC VALUE FOR ', 
c     *                           '(-TEMP2 * I0 * TEMP3)'

c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES OF RELATED VARIABLES'
c                      WRITE(*,*) '---------------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, I0    : ', I0
c                      WRITE(*,*) 'SMITH_ALUKAS, TEMP2 : ', TEMP2
c                      WRITE(*,*) 'SMITH_ALUKAS, TEMP3 : ', TEMP3
c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES RELATED TO I0'
c                      WRITE(*,*) '--------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, IMAX  : ', IMAX
c                      WRITE(*,*) 'SMITH_ALUKAS, DTDAY : ', DTDAY
c                      WRITE(*,*) 'SMITH_ALUKAS, FDAY  : ', FDAY
c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES RELATED TO TEMP2'
c                      WRITE(*,*) '-----------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, CCHL1  : ', CCHL1
c                      WRITE(*,*) 'SMITH_ALUKAS, GITMAX : ', GITMAX
c                      WRITE(*,*) 'SMITH_ALUKAS, PHIMX  : ', PHIMX
c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES RELATED TO TEMP3'
c                      WRITE(*,*) '---------------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, TEMP1 : ', TEMP1
c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES RELATED TO TEMP1'
c                      WRITE(*,*) '---------------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, XKC   : ', XKC
c                      WRITE(*,*) 'SMITH_ALUKAS, TCHLA : ', TCHLA
c                      WRITE(*,*) 'SMITH_ALUKAS, KESHD : ', KESHD
c                      WRITE(*,*) 'SMITH_ALUKAS, KE    : ', KE
c                      WRITE(*,*) 'SMITH_ALUKAS, SKE   : ', SKE
c                      STOP

c                  END IF


c                  IF (((-TEMP2 * I0).GE.100.0).OR.
c     *                ((-TEMP2 * I0).LE.-100.0)) THEN

c                      WRITE(*,*) 'ERROR IN SUBROUTINE SMITH_ALUKAS'
c                      WRITE(*,*) '-------------------------'

c                      WRITE(*,*) 'UNREALISTIC VALUE FOR ', 
c     *                           '(-TEMP2 * I0)'

c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES OF RELATED VARIABLES'
c                      WRITE(*,*) '---------------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, TEMP2 : ', TEMP2
c                      WRITE(*,*) 'SMITH_ALUKAS, I0    : ', I0
c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES RELATED TO I0'
c                      WRITE(*,*) '--------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, IMAX  : ', IMAX
c                      WRITE(*,*) 'SMITH_ALUKAS, DTDAY : ', DTDAY
c                      WRITE(*,*) 'SMITH_ALUKAS, FDAY  : ', FDAY
c                      WRITE(*,*)
c                      WRITE(*,*) 'VALUES RELATED TO TEMP2'
c                      WRITE(*,*) '-----------------------'
c                      WRITE(*,*) 'SMITH_ALUKAS, CCHL1  : ', CCHL1
c                      WRITE(*,*) 'SMITH_ALUKAS, GITMAX : ', GITMAX
c                      WRITE(*,*) 'SMITH_ALUKAS, PHIMX  : ', PHIMX
c                      STOP

c                  END IF

          print *, 'Location 9a in SMITH'
                  SUM = SUM + 2.7183 / TEMP1 * 
     *                  (EXP(-TEMP2 * I0 * TEMP3) - EXP(-TEMP2 * I0))

              ELSE

                  EXIT

          print *, 'Location 9b in SMITH'
              END IF


 1010     CONTINUE


          RLIGHT = SUM / 24.0
          LLIGHT = RLIGHT

C         Adapt carbon to chlorophyll ratio:
          IAV   = 0.9 * ITOT/FDAY
          IAVSG = IAV * (1.0 - TEMP3) / TEMP1
          CCHLX = 0.3 * 0.083 * PHIMX * XKC * IAVSG / (GITMAX * 2.718)

          print *, 'Location 10 in SMITH'


          IF (CCHLX.LT.1.5D1) THEN

              CCHLX = 1.5D1

          END IF



          IF (CCHLX.GT.6.0D1) THEN

              CCHLX = 6.0D1

          END IF

          print *, 'Location 11 in SMITH'

 9999     CONTINUE
 
      END

c********************************************************************
c********************************************************************

      DOUBLE PRECISION FUNCTION ALIGHT_ALUKAS(TLAYER, LAYER, LAYERS, 
     *                                 CELL  , CELLS, DEPTHS, PHYTCC,
     *                                 CCHLXI, KE   , XKC   , SURFL)

          IMPLICIT NONE

          INTEGER TLAYER
          INTEGER LAYER
          INTEGER LAYERS
          INTEGER CELL
          INTEGER CELLS
      
	    DOUBLE PRECISION DEPTHS(CELLS, LAYERS)
          DOUBLE PRECISION PHYTCC(CELLS, LAYERS)
          DOUBLE PRECISION CCHLXI
          DOUBLE PRECISION KE
          DOUBLE PRECISION XKC
          DOUBLE PRECISION SURFL

          DOUBLE PRECISION TCHLA
          DOUBLE PRECISION LIGHT
          DOUBLE PRECISION EXTIN
          INTEGER I

          LIGHT = SURFL


          DO 1010, I = TLAYER + 1, LAYER

              TCHLA = PHYTCC(CELL, I) / CCHLXI
              EXTIN = KE + (XKC * 1000.0 * TCHLA)
              LIGHT = LIGHT * EXP(-EXTIN * DEPTHS(CELL, I))

 1010     CONTINUE


          ALIGHT_ALUKAS = LIGHT

      END



C     Function, which returns saturation concentration of dissolved oxygen
      DOUBLE PRECISION FUNCTION DO_SAT_ALUKAS(T, S, H)

C         T : Water temperature (in Celcius)
C         S : Salinity (in ppt)
C         H : Water surface elevation (in m) 

          DOUBLE PRECISION T
          DOUBLE PRECISION S
          DOUBLE PRECISION H
 
C         T_KELVIN : Water temperature (in Kelvins)
C         H_FEET   : Water surface elevation (in feet)
          DOUBLE PRECISION T_KELVIN  
          DOUBLE PRECISION H_FEET

          DOUBLE PRECISION LN_CSF
          DOUBLE PRECISION LN_CSS
          DOUBLE PRECISION CSS
          DOUBLE PRECISION CSP

          DOUBLE PRECISION CS

C         P     : Atmospheric pressure at elevation H (in atm)
C         P0    : Standart pressure (in mmHg) 
C         PWV   : Partial pressure of water vapor (in atm)
C         THETA : A constant

          DOUBLE PRECISION P
          DOUBLE PRECISION P0
          DOUBLE PRECISION LN_PWV           
          DOUBLE PRECISION PWV
          DOUBLE PRECISION THETA

          T_KELVIN = T + 273.15
          H_FEET = H / 0.3048

C         Calculate the effect of temperature on dissolved oxygen saturation
          LN_CSF = -139.34411 + (157570.1 / T_KELVIN)
     *           - (66423080.0 / (T_KELVIN**2)) + 
     *               (12438000000.0 / (T_KELVIN**3))
     *           - (862194900000.0 / (T_KELVIN**4))

          !Calculate the effect of salinity on dissolved oxygen saturation
          LN_CSS = LN_CSF - S 
     *    * (0.017674 - (10.754 / T_KELVIN)  + (2140.7 / (T_KELVIN**2)))

          CSS = EXP(LN_CSS)
 
C         Calculate the effect of altitude on dissolved oxygen saturation
    
C         Calculate THETA
          THETA = 0.000975 - (0.00001426 * T) + (0.00000006436 * (T**2))

C         Set standard pressure to mean sea level
          P0 = 760

C         Calculate atmospheric pressure at altitude H
          P = (P0 - (0.02667 * H_FEET)) / 760.0

C         Calculate vapour pressure of water
          LN_PWV = 11.8571 - 
     *    (3840.7 / T_KELVIN) - (216961 / (T_KELVIN**2))

          PWV = EXP(LN_PWV)

C         Final calculation including elevation effect
          CSP = CSS  * P * (((1 - (PWV / P)) * (1 - (THETA * P)))
     *        / ((1 - PWV) * (1 - THETA)))

          CS = CSP

          DO_SAT_ALUKAS = CS

      END



      DOUBLE PRECISION FUNCTION KAWIND_ALUKAS(WS, TW, TA, DEPTH, WTYPE)

          IMPLICIT NONE
C         WS         wind speed, m/s
C         TW         water temperature C
C         TA         air temperature C
C         depth      defined depth(segmax) in geometrical
C         WTYPE      type of water body


          DOUBLE PRECISION WS
          DOUBLE PRECISION TW
          DOUBLE PRECISION TA
          DOUBLE PRECISION DEPTH
          DOUBLE PRECISION WTYPE

C         RK         reareation term calculated in kawind
          DOUBLE PRECISION RK
          DOUBLE PRECISION R0MIN
          INTEGER IWTYPE
          INTEGER N


C         UT     : SHEAR VELOCITY (CM/SEC)
C         UC     : CRITICAL SHEAR VELOCITY (CM/SEC) 
C         KARMAN : VONKARMAN'S CONSTANT
C         ZE     : EQUILIBRIUM ROUGHNESS (CM)
C         1/LAM  : A REYNOLD'S NUMBER
C         GAM    : NONDIMENSIONAL COEFFICIENT DEPENDENT ON WATER BODY SIZE.
          DOUBLE PRECISION UT
          DOUBLE PRECISION UC
          DOUBLE PRECISION KARMAN
          DOUBLE PRECISION ZE
          DOUBLE PRECISION LAM
          DOUBLE PRECISION GAM


C         DIFF : DIFFUSIVITY OF OXYGEN IN WATER (CM**2/SEC)
C         VW   : VISCOSITY OF WATER             (CM**2/SEC)
C         VA   : VISCOSITY OF AIR               (CM**2/SEC)
C         PW   : DENSITY OF WATER               (G/CM**3)
C         PA   : DENSITY OF AIR                 (G/CM**3)
          DOUBLE PRECISION DIFF
          DOUBLE PRECISION VW
          DOUBLE PRECISION VA
          DOUBLE PRECISION PA
          DOUBLE PRECISION PW

          DOUBLE PRECISION KA3
          DOUBLE PRECISION WH
          DOUBLE PRECISION SRCD
          DOUBLE PRECISION ERR
          DOUBLE PRECISION EF
          DOUBLE PRECISION F1
          DOUBLE PRECISION F2
          DOUBLE PRECISION FP1
          DOUBLE PRECISION FP2
          DOUBLE PRECISION FP3
          DOUBLE PRECISION FP4
          DOUBLE PRECISION SRCD2
          DOUBLE PRECISION CDDRAG
          DOUBLE PRECISION US
          DOUBLE PRECISION Z0
          DOUBLE PRECISION RK1
          DOUBLE PRECISION RK2
          DOUBLE PRECISION RK3
          DOUBLE PRECISION GAMU

  
          R0MIN = 1.e-15
          IWTYPE = WTYPE


          IF (IWTYPE.EQ.3) THEN

              UT  = 10.0
              UC  = 11.0
              ZE  = 0.35
              LAM = 3.0
              GAM = 5.0

          ELSE 


              IF (IWTYPE .EQ.1) THEN 

                  UT = 9.0
                  UC = 22.0
                  ZE = 0.25
                  LAM = 10.0
                  GAM = 10.0

              ELSE


                  IF (IWTYPE.EQ.2) THEN  

                      UT  = 10.0
                      UC  = 11.0
                      ZE  = 0.25
                      LAM = 3.0
                      GAM = 6.5

                  ELSE

                      WRITE(*,*) 
     *                'KAWIND_ALUKAS : WRONG VALUE FOR WATERBODY TYPE'

                      STOP

                  END IF


              END IF


          END IF


          DIFF = 4.58E-07 * TW + 1.2E-05
          VW = 0.0164 - 0.00024514*TW
          VA = 0.133 + 0.0009*TA
          PA = 0.00129 - 0.0000040*TA
          PW = 1.00
          WS = WS*100.0
          RK = 1.0

C         NEWTON RAPHSON METHOD TO CALCULATE THE SQUARE ROOT OF THE DRAG
C         COEFFICIENT
          N = 0

          KARMAN = 0.4
          KA3    = KARMAN**0.3333
          WH     = 1000.0

C         INITIAL GUESS FOR SQUARE ROOT OF THE DRAG COEFFICIENT
          SRCD = 0.04
          ERR = 1.0


          DO 1010, N = 1, 9

C             CALCULATE VALUE OF FUNCTION(F2) AND 
C             DERIVATIVE OF FUNCTION(FP)
              EF  = EXP( - SRCD*WS/UT)
              F1  = LOG((WH/ZE) + (WH*LAM/VA)*SRCD*WS*EF)
              F2  = F1 - KARMAN / SRCD
              FP1 = 1.0/((WH/ZE) + (LAM*WH/VA)*SRCD*WS*EF)
              FP2 = ((WH*LAM)/(VA*UT))*SRCD*(WS**2.0)*EF
              FP3 = (WH*(LAM / VA)) * WS * EF
              FP4 = FP1*(FP2 + FP3) + (KARMAN/(SRCD**2.0))

C             ANEW GUESS FOR SQUARE ROOT OF DRAG AND COMPARE TO
C             PREVIOUS GUESS AND LOOP BACK THROUGH N-R WITH NEW GUESS
C             IF APPROPRIATE
              SRCD2 = SRCD - F2/FP4
              ERR = ABS(SRCD - SRCD2)
              SRCD = SRCD2


              IF (ERR.LE.0.0005) THEN

                  EXIT

              END IF


 1010     CONTINUE   
       

          IF ((ERR.GT.0.005).AND.(N.EQ.9)) THEN 

              WRITE(*,*) 'KAWIND_ALUKAS : SOLUTION DID NOT CONVERGE'
              STOP

          END IF

  
          CDDRAG = SRCD**2.0
          US     = SRCD * WS
          Z0     = 1.0 / ((1.0 / ZE) + LAM * US * EXP(-US / UT) / VA)
          WS     = WS/100.0


          IF (WS.LT.6.0) THEN 

              RK1 = ((DIFF / VW)**0.666667) * SRCD * ((PA / PW)**0.5)
              RK = RK1 * KA3 * (WS/GAM)
              RK = RK * 3600.0 *24.0
              RK= RK / DEPTH

          END IF


          IF ((WS.GE.6.0).AND.(WS.LE.20.0)) THEN

              GAMU = GAM * US * (EXP(-(US / UC) + 1.0) / UC)

              RK1  = ((DIFF/VW)**0.6667) * KA3 * ((PA/PW)**0.5) * 
     *                (US / GAMU)

              RK2  = ((DIFF * US * PA * VA) / (KARMAN * Z0 * PW * VW))
     *               **0.5

              RK3  = (1.0 / RK1) + (1.0 / RK2)
              RK   = 1.0 / RK3
              RK   = RK * 3600.0 * (24.0 / 100.0)
              RK   = RK / DEPTH

          END IF


          IF (WS.GT.20.0) THEN

              RK = ((DIFF * PA * VA * US) / (KARMAN * ZE * PW * VW))
     *             **0.5

              RK = RK * 3600.0 * (24.0 / 100.0)
              RK = RK / DEPTH

          END IF


          KAWIND_ALUKAS = RK

      END
c
c*************************************************************************
c
      DOUBLE PRECISION FUNCTION KAHYDR_ALUKAS(H, VEL, STP20)

          IMPLICIT NONE

          DOUBLE PRECISION H
          DOUBLE PRECISION VEL
          DOUBLE PRECISION STP20


          DOUBLE PRECISION CFOREA
          DOUBLE PRECISION AVDEPE
          DOUBLE PRECISION AVVELE

          DOUBLE PRECISION REAK
          DOUBLE PRECISION EXPREV
          DOUBLE PRECISION EXPRED
          DOUBLE PRECISION K20
          DOUBLE PRECISION TRANDP
          DOUBLE PRECISION DIF

          CFOREA = 1.0
          AVDEPE = H
          AVVELE = ABS(VEL)


C         Calculate reaeration coefficient for free-flowing reach
C         Calculate reaeration coefficient as a power function of
C         average hydraulic depth and velocity; determine exponents
C         to depth and velocity terms and assign value to REAK


C         Use Owen's formulation for reaeration
          IF (AVDEPE.LE.0.61) THEN 

              REAK = 5.349
              EXPREV = 0.67
              EXPRED = - 1.85
              K20 = REAK*(AVVELE**EXPREV)*(AVDEPE**EXPRED)


              IF (K20.GT.24) THEN 

                  K20 = 24.0

              END IF


              KAHYDR_ALUKAS = K20*(1.028**STP20)
              RETURN

          END IF


C        Calculate transition depth; transition depth determines
C        which method of calculation is used given the current
C        velocity
         IF (AVVELE.LT.0.518) THEN

             TRANDP = 0.0

         ELSE

             TRANDP = 4.411*(AVVELE**2.9135)

         END IF


         DIF = AVDEPE - TRANDP


C        Use Churchill's formulation for reaeration
         IF (DIF.LE.0.0) THEN

             REAK   = 5.049
             EXPREV = 0.969
             EXPRED = - 1.673
             K20    = REAK * (AVVELE**EXPREV) * (AVDEPE**EXPRED)


             IF (K20.GT.24) THEN

                 K20 = 24.0

             END IF


             KAHYDR_ALUKAS = K20*(1.028**STP20)
             RETURN

          END IF
 

C         Use O'Connor-Dobbins formulation for reaeration
          REAK   = 3.93
          EXPREV = 0.5
          EXPRED = - 1.5
          K20 = REAK*(AVVELE**EXPREV)*(AVDEPE**EXPRED)

 
          IF (K20.GT.24.0) THEN

              K20 = 24.0

          END IF
 
 
          KAHYDR_ALUKAS = K20 * (1.028**STP20)

      END

c**************************************************************
c*************************************************************

      DOUBLE PRECISION FUNCTION FKPH_ALUKAS(pH, pHmin, pHmax)
C         pH correction factor
          IMPLICIT NONE

          DOUBLE PRECISION pH
          DOUBLE PRECISION pHmin
          DOUBLE PRECISION pHmax


          IF ((pH.LT.0.0).OR.(pH.GT.14.0)) THEN

              WRITE(*,*) 'PH : ', pH   
              STOP "FKPH_ALUKAS : pH value is not correct"

          END IF



          IF (pH.LT.pHmin) THEN

              FKPH_ALUKAS = EXP(pH-pHmin)

          ELSE


              IF (pH.GT.pHmax) THEN

                  FKPH_ALUKAS = EXP(pHmax - pH)

              ELSE

                  FKPH_ALUKAS = 1.0D0

              END IF


          END IF


      END

c*****************************************************************
c*****************************************************************

      DOUBLE PRECISION FUNCTION TSTROG_ALUKAS
     *        (Ts, Q10, TMax, TOpt, XM, KT, TRef)
 
C     Stroganov temperature correction for algal growth based on
C     assumption that organism is never fully acclimated, i.e. temperature change
C     is faster than acclimation time.
C     Acclimation is interpreted as temperature shift that it is equal to 0
C     below the temperature TRef and fast aproaches maximal value (XM) for
C     higher temperatures

C     rt=TSTROG_ALUKAS(Ts,Q10,TMax,TOpt,XM,KT,TRef)

          
C     Ts   - temperature
C     Q10  - Q10 factor (~2)
C     TMax - maximal temperature
C     TOpt - optimal temperature
C     XM   - maximal acclimation temperature(maximal shift)
C     KT   - exponential rate to aproach XM
C     Tref - temperature below which no aclimation take place
 

          DOUBLE PRECISION Ts
          DOUBLE PRECISION Q10
          DOUBLE PRECISION Tmax
          DOUBLE PRECISION TOpt
          DOUBLE PRECISION XM
          DOUBLE PRECISION KT
          DOUBLE PRECISION Tref

          DOUBLE PRECISION Acclimation
          DOUBLE PRECISION TMaxA
          DOUBLE PRECISION TOptA
          DOUBLE PRECISION WT
          DOUBLE PRECISION YT
          DOUBLE PRECISION VT
          DOUBLE PRECISION XT
          DOUBLE PRECISION T
          INTEGER SIGN

          T = Ts 


          IF (T.LE.TRef) THEN

              SIGN = -1

          ELSE

              SIGN = 1

          END IF
 
 
          Acclimation = sign * XM * (1.0 - EXP(-KT * ABS(T-TRef)))
          TMaxA = TMax + Acclimation
          TOptA = TOpt + Acclimation
          WT = LOG(Q10) * (TMaxA - TOptA)
          YT = WT + 2 * LOG(Q10)
          VT = (TMaxA - T) / (TMaxA - TOptA)
          XT = ((WT * (1.0 + SQRT(1.0 + 40.0 / YT)))**2.0) / 400.0


          IF (T.GE.TMaxA) THEN

              TSTROG_ALUKAS = 1e-10

          ELSE

              TSTROG_ALUKAS = (VT**XT) * EXP(XT*(1-VT))

          END IF         


      END



c
c********************************************************************      
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

      SUBROUTINE cur_smith(PHOTO,Ia,TCHLA,CCHLXI,GITMAX,H,ke,XKC,PHIMX, 
     *                 LLIGHT,CCHLX)   
C     Dick-Smith light limitation formulation(adapted from EUTRO)
     
C   Version to use in ALUKAS with dayly averaged light intensity
C   with  double precision variables
 
C Geting FDAY, ITOT

          double precision PHOTO
          double precision Ia
          double precision TCHLA
          double precision CCHLXI
          double precision GITMAX
          double precision H
          double precision ke
          double precision XKC
          double precision PHIMX

          double precision LLIGHT
          double precision CCHLX


          double precision FDAY
          double precision ITOT
          double precision CCHL1
          double precision KESHD
          double precision SKE
          double precision TEMP1
          double precision TEMP2
          double precision TEMP3
          double precision IMAX
          double precision SUM
          double precision DTDAY
          double precision I0
          double precision RLIGHT
          double precision IAV
          double precision IAVSG


          double precision PI
          INTEGER I          
          

          PI = 3.14159

          FDAY=PHOTO
          ITOT = Ia
c         CCHL1 = CCHLX(ISEG)      
          CCHL1 = CCHLXI
 
C         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C         %?IAVBOT=IAVBOTX(ISEG);
C         %PHYT=TPHY;
C         %GITMAX=k1c*rtmult(TEMP, t_1, t_2, t_3, t_4, k_1, 0.98, 0.98, k_4);
C         %TCHLA = PHYT/CCHL1;
C         %RLIGHT = RLGHTS (ISEG, 1)
C         %
C         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
C         %           Dick Smith formulation integrated every day
 
          KESHD = XKC * 1000.0 * TCHLA    
          SKE = ke
          SKE = SKE + KESHD
          TEMP1 = SKE * H
          
          if (GITMAX .lt. 1e-20 . or. CCHL1 .lt. 1e-20) then
           write(6,*) 'SMITH: TEMP2 is NaN ', 'GITMAX=', GITMAX,           
     *               'CCHL1=', CCHL1
           stop
          end if
                
          TEMP2 = 0.083 * PHIMX * XKC / (GITMAX * CCHL1 * 2.718)
          
           
          
          TEMP3 = EXP( - TEMP1)

C         %RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3
C         %IMAX = PI*ITOT*RLGHTS (ISEG, 2)/(2.*FDAY)
          IMAX = pi*ITOT/(2.*FDAY)
          SUM = 0.0


          DO 1010 I = 1,25

              DTDAY = (I - 1)/24.0
              IF (DTDAY.LE.FDAY) THEN
                  I0 = IMAX * SIN(PI * DTDAY/FDAY)
                  SUM = SUM + 2.7183 / TEMP1 * 
     *                  (EXP( -TEMP2 * I0 * TEMP3) - EXP( -TEMP2 * I0))
              ELSE
                  EXIT
              END IF
              
 1010     CONTINUE


          RLIGHT = SUM/24.0
          LLIGHT=RLIGHT

C         Adapt carbon to chlorophyll ratio:
          IAV=0.9 * ITOT/FDAY
          IAVSG=IAV*(1.0-TEMP3)/TEMP1
          CCHLX=0.3 * 0.083 * PHIMX * XKC * IAVSG / (GITMAX * 2.718)


          IF (CCHLX.LT.15.0) THEN
              CCHLX=15.0
c              write(6,*)PHIMX,XKC,IAVSG,GITMAX
          END IF


      END
      
c********************************************************************
c********************************************************************

!       DOUBLE PRECISION FUNCTION HFLUX(STIME, LAT, LONG, HEMIS, ELEV, 
!      *                                AITEMP, RELHUM, FCLOUD, WINDS, 
!      *                                G_REF, DUST_C, WTEMP, DENS, A, B)
! 
!           IMPLICIT NONE
! 
!           DOUBLE PRECISION STIME
!           DOUBLE PRECISION LAT
!           DOUBLE PRECISION LONG
!           INTEGER HEMIS
!           DOUBLE PRECISION ELEV
!           DOUBLE PRECISION AITEMP
!           DOUBLE PRECISION RELHUM
!           DOUBLE PRECISION FCLOUD
!           DOUBLE PRECISION WINDS
!           DOUBLE PRECISION G_REF
!           DOUBLE PRECISION DUST_C
!           DOUBLE PRECISION WTEMP
!           DOUBLE PRECISION DENS
!           DOUBLE PRECISION A
!           DOUBLE PRECISION B
! 
!           DOUBLE PRECISION HSWRAD
!           DOUBLE PRECISION SWRAD
! 
!           DOUBLE PRECISION HLWRAD
!           DOUBLE PRECISION LWRAD
! 
!           DOUBLE PRECISION HBRAD
!           DOUBLE PRECISION BRAD
! 
!           DOUBLE PRECISION HLEVAP
!           DOUBLE PRECISION LEVAP
! 
!           DOUBLE PRECISION HTRANS
!           DOUBLE PRECISION STRANS
! 
!           SWRAD = HSWRAD(STIME, LAT, LONG, HEMIS, ELEV, 
!      *                   AITEMP, RELHUM, FCLOUD, G_REF, DUST_C)
! 
!           LWRAD  = HLWRAD(AITEMP, FCLOUD)
!           BRAD   = HBRAD(WTEMP - 273.16)
! 
!           LEVAP  = HLEVAP(WINDS, (WTEMP - 273.16), 
!      *             AITEMP, RELHUM, DENS, A, B)
! 
!           STRANS = HTRANS(AITEMP, (WTEMP - 273.16), 
!      *             ELEV, WINDS, DENS, A, B)
! 
!           HFLUX = SWRAD + LWRAD + STRANS - BRAD - LEVAP 
! 
!       END
! 
! 
! 
!       DOUBLE PRECISION FUNCTION HSWRAD(STIME, LAT, LONG, HEMIS, ELEV, 
!      *                                 AITEMP, RELHUM, FCLOUD, G_REF, 
!      *                                 DUST_C)
! 
!           IMPLICIT NONE
! 
!           DOUBLE PRECISION STIME
!           DOUBLE PRECISION LAT
!           DOUBLE PRECISION LONG
!           INTEGER HEMIS
!           DOUBLE PRECISION ELEV
!           DOUBLE PRECISION AITEMP
!           DOUBLE PRECISION RELHUM
!           DOUBLE PRECISION FCLOUD
!           DOUBLE PRECISION G_REF
!           DOUBLE PRECISION DUST_C
! 
! 
! C         H_ZERO : Amount of solar radiation reaching the earths
! C                  outer atmospere
! C
! C         A_T    : Atmospheric transmission term
! C
! C         R_S    : Albedo (reflection coefficient) 
! C
! C         C_A    : Fraction of solar radiation not adsorbed by
! C                  the clouds
!           DOUBLE PRECISION H_ZERO
!           DOUBLE PRECISION A_T
!           DOUBLE PRECISION R_S
!           DOUBLE PRECISION C_A
! 
! 
! C         FOR CALCULATION OF H_ZERO
!           DOUBLE PRECISION H_SC
!           DOUBLE PRECISION PI
!           DOUBLE PRECISION R
!           DOUBLE PRECISION DELTA
!           DOUBLE PRECISION GAMMA
!           DOUBLE PRECISION THETA
!           DOUBLE PRECISION H_B
!           DOUBLE PRECISION H_E
! 
!           PARAMETER(H_SC = 1.44D3)
!           PARAMETER(PI   = 3.14159265)
! 
! C         FOR CALCULATION OF R AND DELTA
!           DOUBLE PRECISION DAY
! 
! C         FOR CALCULATION OF H_E AND H_B
!           DOUBLE PRECISION BRACK
!           DOUBLE PRECISION DEL_TS
!           DOUBLE PRECISION A
!           DOUBLE PRECISION B
!           DOUBLE PRECISION HOUR
! 
! C         FOR CALCULATION DEL_TS
!           DOUBLE PRECISION E_A
!           DOUBLE PRECISION L_SM
! 
! C         FOR CALCULATION OF GAMMA
!           DOUBLE PRECISION T_SS
!           DOUBLE PRECISION T_SU
! 
! 
! C         FOR CALCULATION OF A_T
!           DOUBLE PRECISION A_1
!           DOUBLE PRECISION A_2
!           DOUBLE PRECISION C_D
!           DOUBLE PRECISION R_G
! 
! C         FOR CALCULATION OF A_1 AND A_2
!           DOUBLE PRECISION P_WC
!           DOUBLE PRECISION OAMASS
! 
! C         FOR CALCULATION OF P_WC
!           DOUBLE PRECISION T_D
! 
! C         FOR CALCULATION OF T_D
!           DOUBLE PRECISION E_S
!           DOUBLE PRECISION R_H
! 
! 
! C         FOR CALCULATION OF OAMASS
!           DOUBLE PRECISION ALPHA
!           DOUBLE PRECISION Z
! 
! C         FOR CALCULATION OF ALPHA     
!           DOUBLE PRECISION ALPHA1
! 
! C         FOR CALCULATION OF ALPHA1
!           DOUBLE PRECISION OMEGA
! 
! 
!           DAY = INT(STIME)
! 
! 
!           IF (DAY.GT.3.65D2) THEN
! 
!               DAY = INT(STIME) - (3.65D2 * (INT(INT(STIME) / 365)))
! 
!           END IF
! 
! 
!           HOUR = (STIME - INT(STIME)) * 2.4D1
!           THETA = LAT
!           Z = ELEV
!           R_H = RELHUM
!           R_G = G_REF
!           C_D = DUST_C
! 
! 
! C         CALCULATION OF H_ZERO
!           R = 1.0D0 + (1.7D-2 * 
!      *                 COS(((2.0D0 * PI) / 3.65D2) * (1.86D2 - DAY)))
! 
!           DELTA = ((2.345D1 * PI) / 1.8D2) * 
!      *            COS(((2.0D0 * PI)/3.65D2)*(1.72D2 - DAY))
! 
! 
!           IF (HEMIS.LE.0) THEN
! 
!               E_A = -1.0D0
! 
!           ELSE
! 
!               E_A = 1.0D0
! 
!           END IF
! 
! 
!           L_SM = 1.5D1 * INT(LONG / 1.5D1)
!           DEL_TS = (E_A / 1.5D1) * (L_SM - LONG)
! 
! 
!           IF (HOUR.LE.1.2D1) THEN
! 
!               A = 1.0D0
! 
!           ELSE
! 
!               A = -1.0D0
! 
!           END IF
! 
! 
!           BRACK = (PI / 1.2D1) * ((HOUR - 1.0D0) - DEL_TS + (A * 1.2D1))
! 
! 
!           IF (BRACK.GT.(2 * PI)) THEN
! 
!               B = -1.0D0
! 
!           END IF
! 
! 
!           IF (BRACK.LT.0.0D0) THEN
! 
!               B = 1.0D0
! 
!           END IF
! 
! 
!           IF ((BRACK.GE.0.0D0).AND.(BRACK.LE.(2 * PI))) THEN
! 
!               B = 0.0D0
! 
!           END IF
! 
! 
!           H_B = BRACK + (B * 2 * PI)
! 
! 
!           BRACK = (PI / 1.2D1) * (HOUR - DEL_TS + (A * 1.2D1))
! 
! 
!           IF (BRACK.GT.(2 * PI)) THEN
! 
!               B = -1.0D0
! 
!           END IF
! 
! 
!           IF (BRACK.LT.0.0D0) THEN
! 
!               B = 1.0D0
! 
!           END IF
! 
! 
!           IF ((BRACK.GE.0.0D0).AND.(BRACK.LE.(2 * PI))) THEN
! 
!               B = 0.0D0
! 
!           END IF
! 
! 
!           H_E = BRACK + (B * 2 * PI)
! 
!           T_SS = ((1.2D1 / PI) * 
!      *           ACOS(-((SIN((PI * THETA) / 1.8D2) * SIN(DELTA)) / 
!      *                  (COS((PI * THETA)/1.8D2) * COS(DELTA))))) + 
!      *           DEL_TS + 1.2D1 
! 
!           T_SU = (-1.0D0 * T_SS) + (2.0D0 * DEL_TS) + 2.4D1
! 
! 
!           IF ((HOUR.LT.T_SU).OR.(HOUR.GT.T_SS)) THEN
! 
!               GAMMA = 0.0D0
! 
!           ELSE
! 
!               GAMMA = 1.0D0
! 
!           END IF
! 
! 
!           H_ZERO = (H_SC / (R ** 2.0D0)) * GAMMA *
!      *             ((SIN((PI * THETA) / 1.8D2) * SIN(DELTA)) +
!      *              ((1.2D1 / PI) * COS((PI * THETA) / 1.8D2) *
!      *               COS(DELTA) * (SIN(H_E) - SIN(H_B))))
! 
! 
! C         CALCULATION OF A_T
!           OMEGA = H_E
! 
!           ALPHA1 = ABS((SIN((PI * DELTA) / 1.8D2) * SIN(DELTA)) + 
!      *                 (COS((PI * DELTA) / 1.8D2) * COS(DELTA) * 
!      *                  COS(OMEGA)))
! 
!           ALPHA = ATAN(ALPHA1 / SQRT(1.0D0 - ALPHA1))
! 
!           OAMASS = (((2.88D2 - (6.5D-3 * Z)) / 2.88D2) ** 5.256D0) /
!      *             (SIN(ALPHA) + (1.5D-1 * ((((ALPHA * 1.8D2) / PI) + 
!      *                                         3.855D0) ** (-1.253D0))))
! 
!           E_S = 6.11D0 * (1.0D1 ** ((7.5D0 * AITEMP) / 
!      *                   (2.377D2 + AITEMP)))
! 
!           T_D = (2.377D2 * LOG10((E_S * R_H) / 6.11D2)) / 
!      *          (7.5D0 - LOG10((E_S * R_H) / 6.11D2))
! 
!           P_WC = 8.5D-1 * EXP(1.1D-1 + (6.14D-2 * T_D))
! 
!           A_1 = EXP(-(4.65D-1 + (1.34D-1 * P_WC)) * 
!      *               (1.29D-1 + (1.71D-1 * EXP(-8.8D-1 * OAMASS))) * 
!      *               OAMASS)
! 
!           A_2 = EXP(-(4.65D-1 + (1.34D-1 * P_WC)) * 
!      *               (1.79D-1 + (4.21D-1 * EXP(-7.21D-1 * OAMASS))) * 
!      *               OAMASS)
! 
!           A_T = (A_2 + 0.5D0 * (1.0D0 - A_1 - C_D)) / 
!      *          (1.0D0 - 0.5D0 * R_G * (1.0D0 - A_1 - C_D))
! 
! 
! C         CALCULATION OF R_S
! 
! C         OVERCAST
!           IF (FCLOUD.GT.9.0D-1) THEN
! 
!               A = 3.3D-1
!               B = -4.5D-1
! 
!           END IF
! 
! 
! C         BROKEN
!           IF ((FCLOUD.GT.5.0D-1).AND.(FCLOUD.LE.9.0D-1)) THEN
! 
!               A = 9.5D-1
!               B = -7.5D-1
! 
!           END IF
! 
! 
! C         SCATTERED
!           IF ((FCLOUD.GE.1.0D-1).AND.(FCLOUD.LE.5.0D-1)) THEN
! 
!               A = 2.2D0
!               B = -9.7D-1
! 
!           END IF
! 
! 
! C         CLEAR
!           IF (FCLOUD.LT.1.0D-1) THEN
! 
!               A = 1.18D0
!               B = -7.7D-1
! 
!           END IF
! 
! 
!           R_S = A * (((1.8D2 * ALPHA) / PI) ** B)
! 
! 
! C         CALCULATION OF C_A
!           C_A = 1.0D0 - (6.5D-1 * (FCLOUD ** 2.0D0))
! 
! 
! C         CALCULATION OF SOLAR WAVE RADIATION 
!           HSWRAD = H_ZERO * A_T * (1.0D0 - R_S) * C_A
! 
!       END
! 
! 
! 
!       DOUBLE PRECISION FUNCTION HLWRAD(AITEMP, FCLOUD)
! 
!           IMPLICIT NONE
! 
!           DOUBLE PRECISION AITEMP
!           DOUBLE PRECISION FCLOUD
! 
! 
!           DOUBLE PRECISION EPSI_A
!           DOUBLE PRECISION SIGMA
!           DOUBLE PRECISION ALPHA0
!           DOUBLE PRECISION REFLC
! 
!           PARAMETER(SIGMA  = 5.67D-8)
!           PARAMETER(ALPHA0 = 0.937D-5)
!           PARAMETER(REFLC  = 9.7D-1)
! 
!           EPSI_A = ALPHA0 * (1.0D0 + (1.7D-1 * FCLOUD)) * 
!      *             ((AITEMP + 273.16 ) ** 2.0D0)
! 
!           HLWRAD = REFLC * EPSI_A * SIGMA * ((AITEMP + 273.16) ** 4.0D0)
! 
!       END
! 
! 
! 
!       DOUBLE PRECISION FUNCTION HBRAD(WTEMP)
! 
!           IMPLICIT NONE
! 
!           DOUBLE PRECISION WTEMP
! 
! 
!           DOUBLE PRECISION EPSI_W
!           DOUBLE PRECISION SIGMA
! 
!           PARAMETER(SIGMA  = 5.67D-8)
!           PARAMETER(EPSI_W  = 9.7D-1)
! 
!           HBRAD = EPSI_W * SIGMA * ((WTEMP + 273.16) ** 4.0D0)
! 
!       END
! 
! 
! 
!       SUBROUTINE CSUN(STIME, LAT, LONG, HEMIS, FDAY, SRISE, SSET)
! 
!           IMPLICIT NONE
! 
!           DOUBLE PRECISION STIME
!           DOUBLE PRECISION LAT
!           DOUBLE PRECISION LONG
!           INTEGER HEMIS
!           DOUBLE PRECISION FDAY
!           DOUBLE PRECISION SRISE
!           DOUBLE PRECISION SSET
! 
! 
! C         FOR CALCULATION OF H_ZERO
!           DOUBLE PRECISION PI
!           DOUBLE PRECISION DELTA
!           DOUBLE PRECISION THETA
! 
!           PARAMETER(PI   = 3.14159265)
! 
! C         FOR CALCULATION OF R AND DELTA
!           DOUBLE PRECISION DAY
! 
! C         FOR CALCULATION OF H_E AND H_B
!           DOUBLE PRECISION DEL_TS
! C          DOUBLE PRECISION HOUR
! 
! C         FOR CALCULATION DEL_TS
!           DOUBLE PRECISION E_A
!           DOUBLE PRECISION L_SM
! 
! C         FOR CALCULATION OF GAMMA
!           DOUBLE PRECISION T_SS
!           DOUBLE PRECISION T_SU
! 
! 
!           DAY = INT(STIME)
! 
! 
!           IF (DAY.GT.3.65D2) THEN
! 
!               DAY = INT(STIME) - (3.65D2 * (INT(INT(STIME) / 365)))
! 
!           END IF
! 
! 
! c          HOUR = (STIME - INT(STIME)) * 2.4D1
!           THETA = LAT
! 
!           DELTA = ((2.345D1 * PI) / 1.8D2) * 
!      *            COS(((2.0D0 * PI)/3.65D2)*(1.72D2 - DAY))
! 
! 
!           IF (HEMIS.LE.0) THEN
! 
!               E_A = -1.0D0
! 
!           ELSE
! 
!               E_A = 1.0D0
! 
!           END IF
! 
! 
!           L_SM = 1.5D1 * INT(LONG / 1.5D1)
!           DEL_TS = (E_A / 1.5D1) * (L_SM - LONG)
! 
! 
!           T_SS = ((1.2D1 / PI) * 
!      *           ACOS(-((SIN((PI * THETA) / 1.8D2) * SIN(DELTA)) / 
!      *                  (COS((PI * THETA) / 1.8D2) * COS(DELTA))))) + 
!      *           DEL_TS + 1.2D1 
! 
!           T_SU = (-1.0D0 * T_SS) + (2.0D0 * DEL_TS) + 2.4D1
! 
!           FDAY  = (T_SS - T_SU) / 2.4D1
!           SRISE = T_SU
!           SSET  = T_SS
! 
!       END   
!       
!      
!       
! 
