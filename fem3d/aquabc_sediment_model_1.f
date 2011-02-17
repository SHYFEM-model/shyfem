
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Main routines for bottom sediment model No1 
c CONTENT:
c  SUBROUTINE SEDIMENT_MODEL_1   -  main routine
c  FUNCTION STRANGER(VALUE)      -  checks for NaNs (input type is double precision)
c  SUBROUTINE FLX_ALUKAS_TO_SED_MOD_1
c  FUNCTION SED_MOD_1_ALUKAS_MOLDI_C
c  SUBROUTINE SED_MOD_1_CVISC
c  SUBROUTINE FLX_SED_MOD_1_TO_ALUKAS 
c
c Produced by Ali Erturk 2010 June
c
c
c************************************************************************

      SUBROUTINE SEDIMENT_MODEL_1
     *           (INIT_SED_STATE_VARS  , SED_DEPTHS , SED_POROSITIES, 
     *            SED_DENSITIES        , PART_MIXING_COEFFS         ,
     *            SED_DIFFUSIONS       , SURF_MIXLEN, SED_BURRIALS  , 
     *            SURF_WATER_CONCS     , SED_TEMPS                  ,
     *            NUM_SED_VARS         , NUM_SED_LAYERS             ,      
     *            SED_MODEL_CONSTANTS  , NUM_SED_CONSTS             , 
     *            SED_DRIVING_FUNCTIONS, NUM_SED_DRIV               ,
     *            FLUXES_TO_SEDIMENTS  , NUM_FLUXES_TO_SEDIMENTS    ,
     *            NUM_FLUX_RECEIVING_SED_LAYERS,  
     *            CELLNO, LAYER, PSTIME, TIME_STEP                  , 
     *            FINAL_SED_STATE_VARS , 
     *            FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS  ,
     *            PROCESSES_sed        , NDIAGVAR_sed ,  
     *            SED_OUTPUTS          , NUM_SED_OUTPUTS)

      IMPLICIT NONE     

C     ARGUMENTS RELATED TO ARAY SIZES
      INTEGER NUM_SED_VARS
      INTEGER NUM_SED_LAYERS
      INTEGER NUM_SED_CONSTS
      INTEGER NUM_SED_DRIV
      INTEGER NUM_FLUXES_TO_SEDIMENTS
      INTEGER NUM_FLUXES_FROM_SEDIMENTS
      INTEGER NDIAGVAR_sed
      INTEGER NUM_FLUX_RECEIVING_SED_LAYERS 
      INTEGER NUM_SED_OUTPUTS
!      PARAMETER(NUM_SED_OUTPUTS = 14)
            
C     INPUT ARGUMENTS
      DOUBLE PRECISION INIT_SED_STATE_VARS(NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION SED_DEPTHS         (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_DENSITIES      (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_POROSITIES     (NUM_SED_LAYERS)
      DOUBLE PRECISION PART_MIXING_COEFFS (NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION SED_DIFFUSIONS     (NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION SURF_MIXLEN
      DOUBLE PRECISION SED_BURRIALS       (NUM_SED_LAYERS)
      DOUBLE PRECISION SURF_WATER_CONCS   (NUM_SED_VARS)
      DOUBLE PRECISION SED_TEMPS          (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_MODEL_CONSTANTS(NUM_SED_CONSTS)
      DOUBLE PRECISION PROCESSES_sed
     * (NUM_SED_LAYERS,NUM_SED_VARS, NDIAGVAR_sed)

      DOUBLE PRECISION SED_DRIVING_FUNCTIONS
     *                 (NUM_SED_LAYERS, NUM_SED_DRIV)

      DOUBLE PRECISION FLUXES_TO_SEDIMENTS  (NUM_FLUXES_TO_SEDIMENTS)
      
      
      INTEGER CELLNO
      INTEGER LAYER
      DOUBLE PRECISION PSTIME
      DOUBLE PRECISION TIME_STEP

C     OUTPUT ARGUMENTS
      DOUBLE PRECISION FINAL_SED_STATE_VARS 
     *                 (NUM_SED_LAYERS, NUM_SED_VARS)

      DOUBLE PRECISION FLUXES_FROM_SEDIMENTS(NUM_FLUXES_FROM_SEDIMENTS)
      DOUBLE PRECISION SED_OUTPUTS(NUM_SED_LAYERS, NUM_SED_OUTPUTS)
      
C     COUNTERS
      INTEGER I, J

C     SEDIMENT STATE VARIABLES
      DOUBLE PRECISION SED_NH4N(NUM_SED_LAYERS)
      DOUBLE PRECISION SED_NO3N(NUM_SED_LAYERS)
      DOUBLE PRECISION SED_DON (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_PON (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_PO4P(NUM_SED_LAYERS)
      DOUBLE PRECISION SED_DOP (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_POP (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_DOXY(NUM_SED_LAYERS)
      DOUBLE PRECISION SED_DOC (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_POC (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_DSi (NUM_SED_LAYERS)
      DOUBLE PRECISION SED_PSi (NUM_SED_LAYERS)

C     SEDIMENT MODEL COEFFICIENTS
C     Dissolution of particulate organic carbon
      DOUBLE PRECISION K_OXIC_DISS_POC   
      DOUBLE PRECISION K_ANOXIC_DISS_POC 
      DOUBLE PRECISION THETA_DISS_POC    
      DOUBLE PRECISION KHS_DISS_POC

C     Dissolution of particulate organic nitrogen
      DOUBLE PRECISION K_OXIC_DISS_PON   
      DOUBLE PRECISION K_ANOXIC_DISS_PON 
      DOUBLE PRECISION THETA_DISS_PON    
      DOUBLE PRECISION KHS_DISS_PON 

C     Dissolution of particulate organic phosphorus
      DOUBLE PRECISION K_OXIC_DISS_POP   
      DOUBLE PRECISION K_ANOXIC_DISS_POP 
      DOUBLE PRECISION THETA_DISS_POP    
      DOUBLE PRECISION KHS_DISS_POP

C     Dissolution of particulate silicon
      DOUBLE PRECISION K_OXIC_DISS_PSi   
      DOUBLE PRECISION K_ANOXIC_DISS_PSi 
      DOUBLE PRECISION THETA_DISS_PSi    
      DOUBLE PRECISION KHS_DISS_PSi

C     Mineralization of dissolved organic carbon
      DOUBLE PRECISION K_OXIC_MINER_DOC
      DOUBLE PRECISION K_ANOXIC_MINER_DOC
      DOUBLE PRECISION THETA_MINER_DOC
      DOUBLE PRECISION KHS_MINER_DOC 
      DOUBLE PRECISION O_TO_C

C     Mineralization of dissolved organic nitrogen
      DOUBLE PRECISION K_OXIC_MINER_DON
      DOUBLE PRECISION K_ANOXIC_MINER_DON
      DOUBLE PRECISION THETA_MINER_DON
      DOUBLE PRECISION KHS_MINER_DON
      
C     Mineralization of dissolved organic phosphorus
      DOUBLE PRECISION K_OXIC_MINER_DOP
      DOUBLE PRECISION K_ANOXIC_MINER_DOP
      DOUBLE PRECISION THETA_MINER_DOP
      DOUBLE PRECISION KHS_MINER_DOP
      
C     Nitrification
      DOUBLE PRECISION K_NITR
      DOUBLE PRECISION THETA_NITR
      DOUBLE PRECISION KHS_NITR_NH4N
      DOUBLE PRECISION KHS_NITR_DOXY

C     Denitrification
      DOUBLE PRECISION K_DENITR
      DOUBLE PRECISION THETA_DENITR
      DOUBLE PRECISION KHS_DENITR_NO3N
      DOUBLE PRECISION KHS_DENITR_DOC
      DOUBLE PRECISION KHS_DENITR_DOXY
      DOUBLE PRECISION DENITR_YIELD

C     Anoxia     
      DOUBLE PRECISION DOXY_AT_ANOXIA

C     Solid partition      
      DOUBLE PRECISION SOLID_PART_COEFF_NH4
      DOUBLE PRECISION SOLID_PART_COEFF_PO4
      

C     SEDIMENT KINETIC PROCESS RATES
      DOUBLE PRECISION R_DISS_POC
      DOUBLE PRECISION R_DISS_PON
      DOUBLE PRECISION R_DISS_POP
      DOUBLE PRECISION R_MINER_DOC
      DOUBLE PRECISION R_MINER_DON
      DOUBLE PRECISION R_MINER_DOP
      DOUBLE PRECISION R_NITR
      DOUBLE PRECISION R_DENITR
      DOUBLE PRECISION R_DISS_PSi

C     SEDIMENT TRANSPORT PROCESS RATES
      DOUBLE PRECISION NEIGHBOUR_CONC
      DOUBLE PRECISION UPPER_CONC_GRADIENT
      DOUBLE PRECISION SED_MIXLEN
      DOUBLE PRECISION SED_DIFFUSION_RATES(NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION PART_MIXING_RATES  (NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION SED_BURRIAL_RATES  (NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION UNIT_AREA_MASSES   (NUM_SED_LAYERS, NUM_SED_VARS)

C     DERIVS
      DOUBLE PRECISION DERIVS          (NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION TRANSPORT_DERIVS(NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION SETTLING_DERIVS (NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION KINETIC_DERIVS  (NUM_SED_LAYERS, NUM_SED_VARS)

      DOUBLE PRECISION TEMP
      DOUBLE PRECISION SETTLING_AFFECTED_DEPTH
      
      INTEGER STRANGER  !Function checking for strange values  
      INTEGER error     !Error indicator
      
      DOUBLE PRECISION DEOXYGENATION
      DOUBLE PRECISION DIFF_CORRECTION_FACTOR
      
      INTEGER TIME_LOOP
      INTEGER NUM_SUB_TIME_STEPS
      DOUBLE PRECISION INTERMED_RESULTS(NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION VOLUME_FRACTION
      INTEGER IN_WHICH_PHASE(NUM_SED_VARS)
      
      DOUBLE PRECISION SOLUTE_FRACTIONS(NUM_SED_LAYERS, NUM_SED_VARS)
      DOUBLE PRECISION WATER_DENSITY
      DOUBLE PRECISION SOLID_CONCS(NUM_SED_LAYERS)
      
      WATER_DENSITY = 1.0D0
      
      
      DO I = 1, NUM_SED_LAYERS

          SOLID_CONCS(I) = (SED_DENSITIES(I) - 
     *                     (WATER_DENSITY * SED_POROSITIES(I)))
     
C        Initialisation of solute fractions      
          DO J = 1, NUM_SED_VARS
              SOLUTE_FRACTIONS (I, J) = 1.0D0 ! fraction of variable in solute form. Not used for solids (Should be zero for solids)!
              PART_MIXING_RATES(I, J) = 0.0D0
          END DO

      END DO
      
C      IN_WHICH_PHASE
C      0 : SOLUTE
C      1 : SOLID
C      2 : ALL SEDIMENTS
  
      IN_WHICH_PHASE(1)  = 2 !2 
      IN_WHICH_PHASE(2)  = 0
      IN_WHICH_PHASE(3)  = 0
      IN_WHICH_PHASE(4)  = 1 !1
      IN_WHICH_PHASE(5)  = 2 !2
      IN_WHICH_PHASE(6)  = 0
      IN_WHICH_PHASE(7)  = 1 !1
      IN_WHICH_PHASE(8)  = 0
      IN_WHICH_PHASE(9)  = 0
      IN_WHICH_PHASE(10) = 1 !1
      IN_WHICH_PHASE(11) = 0
      IN_WHICH_PHASE(12) = 1 !1
      
      error = 0
       
      call set_3d_d_array(NUM_SED_LAYERS,NUM_SED_VARS,NDIAGVAR_sed
     +		,PROCESSES_sed,0.0D+0)

      !PROCESSES_sed(:,:,:) = 0.0D+0 !Fixme
       
      DO I = 1, NUM_SED_LAYERS
        DO J = 1, NUM_SED_VARS       
          if(STRANGER(INIT_SED_STATE_VARS(I,J)) .eq. 1) then
            print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j,
     *            'Cell ',CELLNO 
                print *, 'Initial state is NaN'
                print *, 'INITIAL(i,j)=',INIT_SED_STATE_VARS(I,J)
                error =1
               end if
         end do
      end do
      if (error .eq. 1) then
        print *, 'time= ', PSTIME
        stop
      endif
  
C     INITIALIZE SEDIMENT STATE VARIABLES
      DO I = 1, NUM_SED_LAYERS
          SED_NH4N(I) = INIT_SED_STATE_VARS(I, 1)
          SED_NO3N(I) = INIT_SED_STATE_VARS(I, 2)
          SED_DON (I) = INIT_SED_STATE_VARS(I, 3)
          SED_PON (I) = INIT_SED_STATE_VARS(I, 4)
          SED_PO4P(I) = INIT_SED_STATE_VARS(I, 5)
          SED_DOP (I) = INIT_SED_STATE_VARS(I, 6)
          SED_POP (I) = INIT_SED_STATE_VARS(I, 7)
          SED_DOXY(I) = INIT_SED_STATE_VARS(I, 8)
          SED_DOC (I) = INIT_SED_STATE_VARS(I, 9)
          SED_POC (I) = INIT_SED_STATE_VARS(I, 10)
          SED_DSi (I) = INIT_SED_STATE_VARS(I, 11)
          SED_PSi (I) = INIT_SED_STATE_VARS(I, 12)
      END DO


C     INITIALIZE SEDIMENT MODEL COEFFICIENTS
C     Dissolution of particulate organic carbon
      K_OXIC_DISS_POC    = SED_MODEL_CONSTANTS(1)
      K_ANOXIC_DISS_POC  = SED_MODEL_CONSTANTS(2)
      THETA_DISS_POC     = SED_MODEL_CONSTANTS(3)
      KHS_DISS_POC       = SED_MODEL_CONSTANTS(4)
      
C     Dissolution of particulate organic nitrogen
      K_OXIC_DISS_PON    = SED_MODEL_CONSTANTS(5)
      K_ANOXIC_DISS_PON  = SED_MODEL_CONSTANTS(6)
      THETA_DISS_PON     = SED_MODEL_CONSTANTS(7)
      KHS_DISS_PON       = SED_MODEL_CONSTANTS(8)
      
C     Dissolution of particulate organic phosphorus
      K_OXIC_DISS_POP    = SED_MODEL_CONSTANTS(9)
      K_ANOXIC_DISS_POP  = SED_MODEL_CONSTANTS(10)
      THETA_DISS_POP     = SED_MODEL_CONSTANTS(11)
      KHS_DISS_POP       = SED_MODEL_CONSTANTS(12)
      
C     Dissolution of particulate silicon
      K_OXIC_DISS_PSi    = SED_MODEL_CONSTANTS(13)
      K_ANOXIC_DISS_PSi  = SED_MODEL_CONSTANTS(14)
      THETA_DISS_PSi     = SED_MODEL_CONSTANTS(15)
      KHS_DISS_PSi       = SED_MODEL_CONSTANTS(16)
      
C     Mineralization of dissolved organic carbon
      K_OXIC_MINER_DOC   = SED_MODEL_CONSTANTS(17)
      K_ANOXIC_MINER_DOC = SED_MODEL_CONSTANTS(18)
      THETA_MINER_DOC    = SED_MODEL_CONSTANTS(19)
      KHS_MINER_DOC      = SED_MODEL_CONSTANTS(20)
      
C     Mineralization of dissolved organic nitrogen
      K_OXIC_MINER_DON   = SED_MODEL_CONSTANTS(21)
      K_ANOXIC_MINER_DON = SED_MODEL_CONSTANTS(22)
      THETA_MINER_DON    = SED_MODEL_CONSTANTS(23)
      KHS_MINER_DON      = SED_MODEL_CONSTANTS(24)
      
C     Mineralization of dissolved organic phosphorus
      K_OXIC_MINER_DOP   = SED_MODEL_CONSTANTS(25)
      K_ANOXIC_MINER_DOP = SED_MODEL_CONSTANTS(26)
      THETA_MINER_DOP    = SED_MODEL_CONSTANTS(27)
      KHS_MINER_DOP      = SED_MODEL_CONSTANTS(28)
      
      O_TO_C             = SED_MODEL_CONSTANTS(29)

C     Nitrification
      K_NITR             = SED_MODEL_CONSTANTS(30)
      THETA_NITR         = SED_MODEL_CONSTANTS(31)
      KHS_NITR_NH4N      = SED_MODEL_CONSTANTS(32)
      KHS_NITR_DOXY      = SED_MODEL_CONSTANTS(33)

C     Denitrification
      K_DENITR           = SED_MODEL_CONSTANTS(34)
      THETA_DENITR       = SED_MODEL_CONSTANTS(35)
      KHS_DENITR_NO3N    = SED_MODEL_CONSTANTS(36)
      KHS_DENITR_DOC     = SED_MODEL_CONSTANTS(37)
      KHS_DENITR_DOXY    = SED_MODEL_CONSTANTS(38)
      DENITR_YIELD       = SED_MODEL_CONSTANTS(39)

      DOXY_AT_ANOXIA     = SED_MODEL_CONSTANTS(40)
      
      SOLID_PART_COEFF_NH4 = SED_MODEL_CONSTANTS(41)
      SOLID_PART_COEFF_PO4 = SED_MODEL_CONSTANTS(42)
      
                        
!             print *,'NUM_SED_VARS= ', NUM_SED_VARS,
!      *        'NUM_SED_LAYERS= ', NUM_SED_LAYERS
!       print *, 'Init state(1,j): ', (INIT_SED_STATE_VARS(1,j),j=1,12)      
!       print *, 'Init state(2,j): ', (INIT_SED_STATE_VARS(2,j),j=1,12)
!       print *, 'Init state(3,j): ', (INIT_SED_STATE_VARS(3,j),j=1,12)
!       
!       print *, 'Constants:',(SED_MODEL_CONSTANTS(i), i = 1,33)     
!       
!       
!       return 


        
      DO I = 1, NUM_SED_LAYERS

          DO J = 1, NUM_SED_VARS
              INTERMED_RESULTS(I, J) = INIT_SED_STATE_VARS(I, J)
          END DO
         
      END DO

      NUM_SUB_TIME_STEPS = 1
      

 
      DO TIME_LOOP = 1, NUM_SUB_TIME_STEPS

      DO I = 1, NUM_SED_LAYERS

          DO J = 1, NUM_SED_VARS

              IF (IN_WHICH_PHASE(J).EQ.0) THEN
                  VOLUME_FRACTION = SED_POROSITIES(I)
              END IF
              
              IF (IN_WHICH_PHASE(J).EQ.1) THEN
                  VOLUME_FRACTION = 1.0D0 - SED_POROSITIES(I)
              END IF

              IF (IN_WHICH_PHASE(J).EQ.2) THEN
                  VOLUME_FRACTION = 1.0D0
              END IF
              
              UNIT_AREA_MASSES(I, J) = 
     *             INTERMED_RESULTS(I, J) * SED_DEPTHS(I) * 
     *             VOLUME_FRACTION                           !g/m3 * m =  g/m2
          END DO
         
      END DO
      
               
      DO I = 1, NUM_SED_LAYERS

C         Calculate solute fractions of NH4 and PO4      
          SOLUTE_FRACTIONS(I, 1) = (1.0D0 / SED_POROSITIES(I)) *
     *       (1.0D0 / (1.0D0 + 
     *                (SOLID_CONCS(I) * SOLID_PART_COEFF_NH4)))

          SOLUTE_FRACTIONS(I, 5) = (1.0D0 / SED_POROSITIES(I)) *
     *        (1.0D0 / (1.0D0 + 
     *                 (SOLID_CONCS(I) * SOLID_PART_COEFF_PO4)))

         
          DO J = 1, NUM_SED_VARS

              IF (IN_WHICH_PHASE(J).EQ.0) THEN
                  VOLUME_FRACTION = SED_POROSITIES(I)
              END IF
              
              IF (IN_WHICH_PHASE(J).EQ.1) THEN
                  VOLUME_FRACTION = 1.0D0 - SED_POROSITIES(I)
              END IF

              IF (IN_WHICH_PHASE(J).EQ.2) THEN
                  VOLUME_FRACTION = 1.0D0
              END IF
          
          
              IF ((IN_WHICH_PHASE(J) .EQ. 0).OR.
     *            (IN_WHICH_PHASE(J) .EQ. 2)) THEN
          
                  IF (I.EQ.1) THEN
                      NEIGHBOUR_CONC = SURF_WATER_CONCS(J)
                      SED_MIXLEN     = SURF_MIXLEN
                  ELSE
                      NEIGHBOUR_CONC = INTERMED_RESULTS(I - 1, J) * 
     *                                 SOLUTE_FRACTIONS(I - 1, J)
 
                      SED_MIXLEN     = 0.5D0 * 
     *                    (SED_DEPTHS(I - 1) + SED_DEPTHS(I))
                  END IF

                  UPPER_CONC_GRADIENT = 
     *                  (INTERMED_RESULTS(I, J) * 
     *                   SOLUTE_FRACTIONS(I, J)) - NEIGHBOUR_CONC

                  DIFF_CORRECTION_FACTOR = 1.0D0 / 
     *                (1.0D0 + (3.0D0 * (1.0D0 - SED_POROSITIES(I))))

                  SED_DIFFUSION_RATES(I, J) = DIFF_CORRECTION_FACTOR *
     *                (UPPER_CONC_GRADIENT * SED_DIFFUSIONS(I, J)) / 
     *                SED_MIXLEN

C                 CALCULATE FLUX FROM SEDIMENTS TO WATER COULUMN
                  IF (I.EQ.1) THEN
                      FLUXES_FROM_SEDIMENTS(J) = 
     *                    SED_DIFFUSION_RATES(I, J)
                  END IF

C                 CALCULATE PARTICLE MIXING
                  IF (IN_WHICH_PHASE(J).EQ.0) THEN
                      PART_MIXING_RATES(I, J) = 0.0D0
                  END IF

                  IF (IN_WHICH_PHASE(J).EQ.2) THEN

                      IF (I.EQ.1) THEN
                          PART_MIXING_RATES(I, J) = 0.0D0
                      ELSE
                          NEIGHBOUR_CONC = 
     *                        INTERMED_RESULTS(I - 1, J) * 
     *                        (1.0D0 - SOLUTE_FRACTIONS(I, J))

                          SED_MIXLEN = 0.5D0 * 
     *                        (SED_DEPTHS(I - J) + SED_DEPTHS(I))

                          UPPER_CONC_GRADIENT = 
     *                        (INTERMED_RESULTS(I, J) * 
     *                        (1.0D0 - SOLUTE_FRACTIONS(I, J))) - 
     *                        NEIGHBOUR_CONC

                          PART_MIXING_RATES(I, J) = 
     *                        (UPPER_CONC_GRADIENT * 
     *                         PART_MIXING_COEFFS(I, J)) / 
     *                         SED_MIXLEN
                      END IF

                  END IF

              ELSE
                  FLUXES_FROM_SEDIMENTS (J) = 0.0D0
                  SED_DIFFUSION_RATES(I, J) = 0.0D0
                 
                  IF (I.EQ.1) THEN
                       PART_MIXING_RATES(I, 1) = 0.0D0
                  ELSE
                      NEIGHBOUR_CONC = INTERMED_RESULTS(I - 1, J)
 
                      SED_MIXLEN     = 0.5D0 * 
     *                    (SED_DEPTHS(I - 1) + SED_DEPTHS(I))

                      UPPER_CONC_GRADIENT = 
     *                      (INTERMED_RESULTS(I, J) * 
     *                       SOLUTE_FRACTIONS(I, J)) - NEIGHBOUR_CONC

                      PART_MIXING_RATES(I, J) = 
     *                      (UPPER_CONC_GRADIENT * 
     *                       PART_MIXING_COEFFS(I, J)) / 
     *                       SED_MIXLEN
                  END IF

              END IF

                            
              SED_BURRIAL_RATES(I, J) = 
     *            INTERMED_RESULTS(I, J) * SED_BURRIALS(I) 

              DERIVS          (I, J) = 0.0D0
              TRANSPORT_DERIVS(I, J) = 0.0D0
              KINETIC_DERIVS  (I, J) = 0.0D0
              SETTLING_DERIVS (I, J) = 0.0D0

          END DO

      END DO


C     CALCULATE TRANSPORT DERIVS
      DO I = 1, (NUM_SED_LAYERS - 1)

          DO J = 1, NUM_SED_VARS

              IF (I.GT.1) THEN 

                  TRANSPORT_DERIVS(I, J) = 
     *                SED_DIFFUSION_RATES(I + 1, J) - 
     *                SED_DIFFUSION_RATES(I, J)     + 
     *                SED_BURRIAL_RATES  (I - 1, J) -
     *                SED_BURRIAL_RATES  (I, J)     +
     *                PART_MIXING_RATES  (I + 1, J) -
     *                PART_MIXING_RATES  (I, J)
                   ELSE
 
                  TRANSPORT_DERIVS(I, J) = 
     *                SED_DIFFUSION_RATES(I + 1, J) - 
     *                SED_DIFFUSION_RATES(I, J)     - 
     *                SED_BURRIAL_RATES  (I, J)     +
     *                PART_MIXING_RATES  (I + 1, J) -
     *                PART_MIXING_RATES  (I, J)

              END IF

          END DO

      END DO

      
      DO J = 1, NUM_SED_VARS

          IF (NUM_SED_LAYERS.GT.1) THEN 

              TRANSPORT_DERIVS(NUM_SED_LAYERS, J) = 
     *            ((-1.0D0) * SED_DIFFUSION_RATES(NUM_SED_LAYERS, J)) +
     *            SED_BURRIAL_RATES (NUM_SED_LAYERS - 1, J) -
     *            SED_BURRIAL_RATES (NUM_SED_LAYERS, J)
     
          ELSE
 
              TRANSPORT_DERIVS(NUM_SED_LAYERS, J) = 
     *            ((-1.0D0) * SED_DIFFUSION_RATES(NUM_SED_LAYERS, J)) + 
     *            SED_BURRIAL_RATES  (NUM_SED_LAYERS, J)

          END IF
  
      END DO

      
C     CALCULATE KINETIC DERIVS
      DO I = 1, NUM_SED_LAYERS
          
          TEMP = SED_TEMPS(I)

          IF (SED_DOXY(I).GE.DOXY_AT_ANOXIA) THEN

              R_DISS_POC  = K_OXIC_DISS_POC  * 
     *            (THETA_DISS_POC  ** (TEMP - 2.0D1)) * 
     *            (SED_POC(I) / (SED_POC(I) + KHS_DISS_POC)) * 
     *            SED_POC(I)
     
              R_DISS_PON  = K_OXIC_DISS_PON  * 
     *            (THETA_DISS_PON  ** (TEMP - 2.0D1)) *
     *            (SED_PON(I) / (SED_PON(I) + KHS_DISS_PON)) *
     *            SED_PON(I)

              R_DISS_POP  = K_OXIC_DISS_POP  * 
     *            (THETA_DISS_POP  ** (TEMP - 2.0D1)) * 
     *            (SED_POP(I) / (SED_POP(I) + KHS_DISS_POP)) *     
     *            SED_POP(I)

              R_DISS_PSi  = K_OXIC_DISS_PSi  * 
     *            (THETA_DISS_PSi  ** (TEMP - 2.0D1)) * 
     *            (SED_PSi(I) / (SED_PSi(I) + KHS_DISS_PSi)) *      
     *            SED_PSi(I)

              R_MINER_DOC = K_OXIC_MINER_DOC * 
     *            (THETA_MINER_DOC ** (TEMP - 2.0D1)) * 
     *            (SED_DOC(I) / (SED_DOC(I) + KHS_MINER_DOC)) *
     *            SED_DOC(I)

              DEOXYGENATION = O_TO_C * R_MINER_DOC
     
              R_MINER_DON = K_OXIC_MINER_DON * 
     *            (THETA_MINER_DON ** (TEMP - 2.0D1)) * 
     *            (SED_DON(I) / (SED_DON(I) + KHS_MINER_DON)) *
     *            SED_DON(I)

              R_MINER_DOP = K_OXIC_MINER_DOP * 
     *            (THETA_MINER_DOP ** (TEMP - 2.0D1)) * 
     *            (SED_DOP(I) / (SED_DOP(I) + KHS_MINER_DOP)) *
     *            SED_DOP(I)

          ELSE

              R_DISS_POC  = K_ANOXIC_DISS_POC  * 
     *            (THETA_DISS_POC  ** (TEMP - 2.0D1)) * 
     *            (SED_POC(I) / (SED_POC(I) + KHS_DISS_POC)) * 
     *            SED_POC(I)

              R_DISS_PON  = K_ANOXIC_DISS_PON  * 
     *            (THETA_DISS_PON  ** (TEMP - 2.0D1)) *
     *            (SED_PON(I) / (SED_PON(I) + KHS_DISS_PON)) *
     *            SED_PON(I)

              R_DISS_POP  = K_ANOXIC_DISS_POP  * 
     *            (THETA_DISS_POP  ** (TEMP - 2.0D1)) * 
     *            (SED_POP(I) / (SED_POP(I) + KHS_DISS_POP)) *     
     *            SED_POP(I)

              R_DISS_PSi  = K_ANOXIC_DISS_PSi  * 
     *            (THETA_DISS_PSi  ** (TEMP - 2.0D1)) * 
     *            (SED_PSi(I) / (SED_PSi(I) + KHS_DISS_PSi)) *      
     *            SED_PSi(I)

              R_MINER_DOC = K_ANOXIC_MINER_DOC * 
     *            (THETA_MINER_DOC ** (TEMP - 2.0D1)) * 
     *            (SED_DOC(I) / (SED_DOC(I) + KHS_MINER_DOC)) *
     *            SED_DOC(I)

              DEOXYGENATION = (O_TO_C * R_MINER_DOC) *
     *            (SED_DOXY(I) / 
     *                 (SED_DOXY(I) + (DOXY_AT_ANOXIA / 2.0D0)))

C              IF (SED_DOXY(I).LT.1.0D-2) THEN
C                  DEOXYGENATION = 0.0D0
C              ELSE
C                  DEOXYGENATION = O_TO_C * R_MINER_DOC
C              END IF

              R_MINER_DON = K_ANOXIC_MINER_DON * 
     *            (THETA_MINER_DON ** (TEMP - 2.0D1)) * 
     *            (SED_DON(I) / (SED_DON(I) + KHS_MINER_DON)) *
     *            SED_DON(I)

              R_MINER_DOP = K_ANOXIC_MINER_DOP * 
     *            (THETA_MINER_DOP ** (TEMP - 2.0D1)) * 
     *            (SED_DOP(I) / (SED_DOP(I) + KHS_MINER_DOP)) *
     *            SED_DOP(I)

          END IF
          
C         Nitrification
C          IF (SED_DOXY(I).LT.1.0D-2) THEN
C              R_NITR  = 0.0D0
C          ELSE
              R_NITR  = K_NITR * (THETA_NITR ** (TEMP - 2.0D1)) * 
     *                  (SED_DOXY(I) / (SED_DOXY(I) + KHS_NITR_DOXY)) *
     *                  (SED_NH4N(I) / (SED_NH4N(I) + KHS_NITR_NH4N)) *
     *                  SED_NH4N(I)
C          END IF

          R_DENITR = K_DENITR * (THETA_DENITR ** (TEMP - 2.0D1)) * 
     *             (KHS_DENITR_DOXY / (SED_DOXY(I) + KHS_DENITR_DOXY))*
     *             (SED_NO3N(I) / (SED_NO3N(I) + KHS_DENITR_NO3N)) *
     *             (SED_DOC (I) / (SED_DOC (I) + KHS_DENITR_DOC))  *
     *             SED_NO3N(I)

          R_DISS_POC = 
     *        R_DISS_POC * SED_DEPTHS(I) * (1.0D0 - SED_POROSITIES(I))
          
          R_DISS_PON = 
     *        R_DISS_PON * SED_DEPTHS(I) * (1.0D0 - SED_POROSITIES(I))

          R_DISS_POP = 
     *        R_DISS_POP * SED_DEPTHS(I) * (1.0D0 - SED_POROSITIES(I))

          R_DISS_PSi = 
     *        R_DISS_PSi * SED_DEPTHS(I) * (1.0D0 - SED_POROSITIES(I))

          R_MINER_DOC = R_MINER_DOC * SED_DEPTHS(I) * SED_POROSITIES(I)
          R_MINER_DON = R_MINER_DON * SED_DEPTHS(I) * SED_POROSITIES(I)
          R_MINER_DOP = R_MINER_DOP * SED_DEPTHS(I) * SED_POROSITIES(I)

          DEOXYGENATION = 
     *        DEOXYGENATION * SED_DEPTHS(I) * SED_POROSITIES(I)

          R_NITR   = R_NITR * SED_DEPTHS(I) * SED_POROSITIES(I)
          R_DENITR = R_DENITR * SED_DEPTHS(I) * SED_POROSITIES(I)

               
          KINETIC_DERIVS(I, 1) = R_MINER_DON - R_NITR
          KINETIC_DERIVS(I, 2) = R_NITR - R_DENITR
          KINETIC_DERIVS(I, 3) = R_DISS_PON - R_MINER_DON
          KINETIC_DERIVS(I, 4) = (-1.0D0) * R_DISS_PON
          KINETIC_DERIVS(I, 5) = R_MINER_DOP
          KINETIC_DERIVS(I, 6) = R_DISS_POP - R_MINER_DOP
          KINETIC_DERIVS(I, 7) = (-1.0D0) * R_DISS_POP
          KINETIC_DERIVS(I, 8) = ((-4.57D0) * R_NITR) - DEOXYGENATION

          KINETIC_DERIVS(I, 9) = R_DISS_POC - R_MINER_DOC -
     *        (R_DENITR / DENITR_YIELD)

          KINETIC_DERIVS(I, 10) = (-1.0D0) * R_DISS_POC
          KINETIC_DERIVS(I, 11) = R_DISS_PSi
          KINETIC_DERIVS(I, 12) = (-1.0D0) * R_DISS_PSi

      END DO

      SETTLING_AFFECTED_DEPTH = 0.0D0

      DO I = 1, NUM_FLUX_RECEIVING_SED_LAYERS

          SETTLING_AFFECTED_DEPTH = 
     *        SETTLING_AFFECTED_DEPTH + SED_DEPTHS(I)

      END DO

            
      DO I = 1, NUM_FLUX_RECEIVING_SED_LAYERS

          DO J = 1, NUM_SED_VARS
              SETTLING_DERIVS(I, J) = FLUXES_TO_SEDIMENTS(J) * 
     *            (SED_DEPTHS(I) / SETTLING_AFFECTED_DEPTH)
                       
          END DO

      END DO


      DO I = 1, NUM_SED_LAYERS

          DO J = 1, NUM_SED_VARS
          
              if(STRANGER(TRANSPORT_DERIVS(I, J)) .eq. 1) then
                print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j,
     *            'Cell ',CELLNO 
                print *, 'Transport derivative is NaN'
                print *, 'Deriv(i,j)=',TRANSPORT_DERIVS(I, J)
                error =1
               end if
                
              if(STRANGER(KINETIC_DERIVS(I, J)) .eq. 1) then
                print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j,
     *            'Cell ',CELLNO                  
                print *, 'Kinetic derivative is NaN'
                print *, 'Deriv(i,j)=',KINETIC_DERIVS(I, J)
                error =1
               end if
                
              if(STRANGER(SETTLING_DERIVS(I, J)) .eq. 1) then
                print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j,
     *            'Cell ',CELLNO  
                print *, 'Settling derivative is NaN'
                print *, 'Deriv(i,j)=',SETTLING_DERIVS(I, J)
                error =1
               end if 

                               
             DERIVS(I, J) = TRANSPORT_DERIVS(I, J) + 
     *                      KINETIC_DERIVS  (I, J) +
     *                      SETTLING_DERIVS (I, J)
     

              IF (IN_WHICH_PHASE(J).EQ.0) THEN
                  VOLUME_FRACTION = SED_POROSITIES(I)
              END IF
              
              IF (IN_WHICH_PHASE(J).EQ.1) THEN
                  VOLUME_FRACTION = 1.0D0 - SED_POROSITIES(I)
              END IF

              IF (IN_WHICH_PHASE(J).EQ.2) THEN
                  VOLUME_FRACTION = 1.0D0
              END IF

              
              UNIT_AREA_MASSES(I, J) = UNIT_AREA_MASSES(I, J) + 
     *             (DERIVS(I, J) * (TIME_STEP / NUM_SUB_TIME_STEPS))

              INTERMED_RESULTS(I, J) = UNIT_AREA_MASSES(I, J) / 
     *            (SED_DEPTHS(I) * VOLUME_FRACTION)
                            
          END DO

      END DO
      
      if (error .eq. 1) then
        print *, 'time= ', PSTIME 
        stop
      endif
      
      END DO !END OF SUB TIME STEPS : CAUTION, FLUXES ARE NOT AVERAGED YET
      
      
      DO I = 1, NUM_SED_LAYERS

          DO J = 1, NUM_SED_VARS
              FINAL_SED_STATE_VARS(I, J) = INTERMED_RESULTS(I, J)
              SED_OUTPUTS(I, J) = INTERMED_RESULTS(I, J)
          END DO
      
          SED_OUTPUTS(I, NUM_SED_VARS + 1) = 
     *        SOLUTE_FRACTIONS(I, 1) * INTERMED_RESULTS(I, 1)
          
          SED_OUTPUTS(I, NUM_SED_VARS + 2) = 
     *        SOLUTE_FRACTIONS(I, 5) * INTERMED_RESULTS(I, 5)
          
      END DO
      
!       if(CELLNO. eq. 2668 ) then
!       i=4
!        print *, '*****************************************'
!        print *, FINAL_SED_STATE_VARS(1, i),' ', SED_OUTPUTS(1, i)
!        print *, FINAL_SED_STATE_VARS(2, i),' ', SED_OUTPUTS(2, i)
!        print *, FINAL_SED_STATE_VARS(3, i),' ', SED_OUTPUTS(3, i)
!       end if
      
      
      
          
      END ! end of sediment routine
      
c************************************************************************ 
c************************************************************************
     
      INTEGER FUNCTION STRANGER(VALUE)
      
C cheks for NaN and Inf

      DOUBLE PRECISION VALUE, BIGNUMBER, RATIO
      
      BIGNUMBER=1.0D300
      STRANGER=0
      
      if (.not.(VALUE .lt. 0.D0).and..not.(VALUE .ge. 0.D0)) then
       STRANGER=1
      end if
      
      if ((VALUE .lt. 0.D0).and.(VALUE .ge. 0.D0)) then
       STRANGER=1
      end if
      
      RATIO = BIGNUMBER/VALUE
      
      if(RATIO .eq. 0.D0) then
       STRANGER=1
      endif
      return
      end
      
c************************************************************************
c************************************************************************
      
           SUBROUTINE FLX_ALUKAS_TO_SED_MOD_1
     *           (STATE_VARIABLES    , NUM_VARS   , 
     *            MODEL_CONSTANTS    , NUM_CONSTS , 
     *            DRIVING_FUNCTIONS  , NUM_DRIV   , 
     *            SETTLING_VELOCITIES, DISSOLVED_FRACTIONS,
     *            BOTTOM_FACTOR , CELLNO, LAYER, PSTIME,
     *            SETTLING_RATES, FLUXES, NUM_FLUXES,
     *            SEDIMENT_TYPE , FRACTION_OF_DEPOSITION, 
     *            NOT_DEPOSITED_FLUXES, NUM_NOT_DEPOSITED_FLUXES)

      IMPLICIT NONE

      INTEGER NUM_VARS
      INTEGER NUM_CONSTS
      INTEGER NUM_FLUXES
      INTEGER NUM_DRIV

      DOUBLE PRECISION STATE_VARIABLES
      DIMENSION STATE_VARIABLES(NUM_VARS)

      DOUBLE PRECISION MODEL_CONSTANTS
      DIMENSION MODEL_CONSTANTS(NUM_CONSTS)

      DOUBLE PRECISION SETTLING_VELOCITIES
      DIMENSION SETTLING_VELOCITIES(NUM_VARS)

      DOUBLE PRECISION DISSOLVED_FRACTIONS
      DIMENSION DISSOLVED_FRACTIONS(NUM_VARS)

      DOUBLE PRECISION SETTLING_RATES
      DIMENSION SETTLING_RATES(NUM_VARS)

      DOUBLE PRECISION FLUXES
      DIMENSION FLUXES(NUM_VARS)
      
      DOUBLE PRECISION DRIVING_FUNCTIONS
      DIMENSION DRIVING_FUNCTIONS(NUM_DRIV)

      INTEGER CELLNO
      INTEGER LAYER

      DOUBLE PRECISION PSTIME
      DOUBLE PRECISION BOTTOM_FACTOR

      INTEGER SEDIMENT_TYPE
      INTEGER NUM_NOT_DEPOSITED_FLUXES
      DOUBLE PRECISION FRACTION_OF_DEPOSITION(NUM_NOT_DEPOSITED_FLUXES)
      DOUBLE PRECISION NOT_DEPOSITED_FLUXES(NUM_NOT_DEPOSITED_FLUXES)
      
      DOUBLE PRECISION NH3_FLUX
      DOUBLE PRECISION NOx_FLUX
      DOUBLE PRECISION OPO4_FLUX
      DOUBLE PRECISION PHY_FLUX
      DOUBLE PRECISION EXLADDETC_FLUX
      DOUBLE PRECISION DOO_FLUX
      DOUBLE PRECISION EXLAPDETC_FLUX
      DOUBLE PRECISION EXREDDETC_FLUX
      DOUBLE PRECISION ZOO_FLUX
      DOUBLE PRECISION ISI_FLUX
      DOUBLE PRECISION EXREPDETC_FLUX
      DOUBLE PRECISION PHY_2_FLUX
      DOUBLE PRECISION PHY_3_FLUX
      DOUBLE PRECISION INC_FLUX
      DOUBLE PRECISION GPHYDDETC_FLUX
      DOUBLE PRECISION GPHYPDETC_FLUX
      DOUBLE PRECISION DPHYDDETC_FLUX
      DOUBLE PRECISION DPHYPDETC_FLUX
      DOUBLE PRECISION CPHYDDETC_FLUX
      DOUBLE PRECISION CPHYPDETC_FLUX
      DOUBLE PRECISION ZOOPDDETC_FLUX
      DOUBLE PRECISION ZOOPPDETC_FLUX
      DOUBLE PRECISION SED_DOC_FLUX
      DOUBLE PRECISION SED_DON_FLUX
      DOUBLE PRECISION SED_DOP_FLUX

      DOUBLE PRECISION SETTLING_FACTORS
      DIMENSION SETTLING_FACTORS(NUM_VARS)

      DOUBLE PRECISION NC_EXLADDETC
      DOUBLE PRECISION NC_EXREDDETC
      DOUBLE PRECISION NC
      DOUBLE PRECISION NC_2
      DOUBLE PRECISION NC_3
      DOUBLE PRECISION NC_ZOOC

      DOUBLE PRECISION PC_EXLADDETC
      DOUBLE PRECISION PC_EXREDDETC
      DOUBLE PRECISION PC
      DOUBLE PRECISION PC_2
      DOUBLE PRECISION PC_3
      DOUBLE PRECISION PC_ZOOC

      DOUBLE PRECISION SiC_EXLADDETC
      DOUBLE PRECISION SiC_EXREDDETC
      DOUBLE PRECISION SiC_2
      DOUBLE PRECISION SiC_ZOOC
      DOUBLE PRECISION ASC

      DOUBLE PRECISION OC_EXLADDETC
      DOUBLE PRECISION OC_EXREDDETC
      DOUBLE PRECISION OC_GPHYDDETC
      DOUBLE PRECISION OC_DPHYDDETC
      DOUBLE PRECISION OC_CPHYDDETC
      DOUBLE PRECISION OC_ZOOPDDETC
      
      INTEGER I


      NC_ZOOC  = DRIVING_FUNCTIONS(1)
      PC_ZOOC  = DRIVING_FUNCTIONS(2)
      SiC_ZOOC = DRIVING_FUNCTIONS(3)


      NC_EXLADDETC  = MODEL_CONSTANTS(511)
      NC_EXREDDETC  = MODEL_CONSTANTS(516)
      NC            = MODEL_CONSTANTS(58)
      NC_2          = MODEL_CONSTANTS(258)
      NC_3          = MODEL_CONSTANTS(358)

      PC_EXLADDETC  = MODEL_CONSTANTS(512)
      PC_EXREDDETC  = MODEL_CONSTANTS(517)
      PC            = MODEL_CONSTANTS(57)
      PC_2          = MODEL_CONSTANTS(257)
      PC_3          = MODEL_CONSTANTS(357)

      SiC_EXLADDETC = MODEL_CONSTANTS(513)
      SiC_EXREDDETC = MODEL_CONSTANTS(518)
      ASC           = MODEL_CONSTANTS(211)

      OC_EXLADDETC  = MODEL_CONSTANTS(548)
      OC_EXREDDETC  = MODEL_CONSTANTS(554)
      OC_GPHYDDETC  = MODEL_CONSTANTS(530)
      OC_DPHYDDETC  = MODEL_CONSTANTS(536)
      OC_CPHYDDETC  = MODEL_CONSTANTS(542)
      OC_ZOOPDDETC  = MODEL_CONSTANTS(524)


      DO I = 1, NUM_VARS

          SETTLING_FACTORS(I) = BOTTOM_FACTOR * 
     *                          (1.0D+0 - DISSOLVED_FRACTIONS(I))* 
     *                          SETTLING_VELOCITIES(I)
     
          IF (SETTLING_FACTORS(I).LT.0.0D0) THEN
              SETTLING_FACTORS(I) = 0.0D0
          END IF

          SETTLING_RATES(I)   = STATE_VARIABLES(I) * 
     *                          (1.0D+0 - DISSOLVED_FRACTIONS(I))*
     *                          SETTLING_VELOCITIES(I)
     
      END DO

      NH3_FLUX       = STATE_VARIABLES(1)  * SETTLING_FACTORS(1)
      NOx_FLUX       = STATE_VARIABLES(2)  * SETTLING_FACTORS(2)
      OPO4_FLUX      = STATE_VARIABLES(3)  * SETTLING_FACTORS(3)
      PHY_FLUX       = STATE_VARIABLES(4)  * SETTLING_FACTORS(4)
      EXLADDETC_FLUX = STATE_VARIABLES(5)  * SETTLING_FACTORS(5)
      DOO_FLUX       = STATE_VARIABLES(6)  * SETTLING_FACTORS(6)
      EXLAPDETC_FLUX = STATE_VARIABLES(7)  * SETTLING_FACTORS(7)
      EXREDDETC_FLUX = STATE_VARIABLES(8)  * SETTLING_FACTORS(8)
      ZOO_FLUX       = STATE_VARIABLES(9)  * SETTLING_FACTORS(9)
      ISI_FLUX       = STATE_VARIABLES(10) * SETTLING_FACTORS(10)
      EXREPDETC_FLUX = STATE_VARIABLES(11) * SETTLING_FACTORS(11)
      PHY_2_FLUX     = STATE_VARIABLES(12) * SETTLING_FACTORS(12)
      PHY_3_FLUX     = STATE_VARIABLES(13) * SETTLING_FACTORS(13)
      INC_FLUX       = STATE_VARIABLES(14) * SETTLING_FACTORS(14)
      GPHYDDETC_FLUX = STATE_VARIABLES(15) * SETTLING_FACTORS(15)
      GPHYPDETC_FLUX = STATE_VARIABLES(16) * SETTLING_FACTORS(16)
      DPHYDDETC_FLUX = STATE_VARIABLES(17) * SETTLING_FACTORS(17)
      DPHYPDETC_FLUX = STATE_VARIABLES(18) * SETTLING_FACTORS(18)
      CPHYDDETC_FLUX = STATE_VARIABLES(19) * SETTLING_FACTORS(19)
      CPHYPDETC_FLUX = STATE_VARIABLES(20) * SETTLING_FACTORS(20)
      ZOOPDDETC_FLUX = STATE_VARIABLES(21) * SETTLING_FACTORS(21)
      ZOOPPDETC_FLUX = STATE_VARIABLES(22) * SETTLING_FACTORS(22)
      SED_DOC_FLUX   = STATE_VARIABLES(23) * SETTLING_FACTORS(23)
      SED_DON_FLUX   = STATE_VARIABLES(24) * SETTLING_FACTORS(24)
      SED_DOP_FLUX   = STATE_VARIABLES(25) * SETTLING_FACTORS(25)

     
C     NOT DEPOSITED FLUXES      
      NOT_DEPOSITED_FLUXES(1)  =
     *     NH3_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(1))
      
      NOT_DEPOSITED_FLUXES(2)  = 
     *     NOx_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(2))
      
      NOT_DEPOSITED_FLUXES(3) =  
     *    OPO4_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(3))
      
      NOT_DEPOSITED_FLUXES(4)  = 
     *     PHY_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(4))
      
      NOT_DEPOSITED_FLUXES(5)  = 
     *    EXLADDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(5))
      
      NOT_DEPOSITED_FLUXES(6)  = 
     *     DOO_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(6))
      
      NOT_DEPOSITED_FLUXES(7)  = 
     *     EXLAPDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(7))
      
      NOT_DEPOSITED_FLUXES(8)  = 
     *     EXREDDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(8))
      
      NOT_DEPOSITED_FLUXES(9)  = 
     *     ZOO_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(9))
      
      NOT_DEPOSITED_FLUXES(10) = 
     *     ISI_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(10))
      
      NOT_DEPOSITED_FLUXES(11) = 
     *     EXREPDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(11))
      
      NOT_DEPOSITED_FLUXES(12) = 
     *     PHY_2_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(12))
      
      NOT_DEPOSITED_FLUXES(13) = 
     *     PHY_3_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(13))
      
      NOT_DEPOSITED_FLUXES(14) = 
     *     INC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(14))
      
      NOT_DEPOSITED_FLUXES(15) = 
     *     GPHYDDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(15))
      
      NOT_DEPOSITED_FLUXES(16) = 
     *    GPHYPDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(16))
      
      NOT_DEPOSITED_FLUXES(17) = 
     *     DPHYDDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(17))
      
      NOT_DEPOSITED_FLUXES(18) = 
     *    DPHYPDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(18))
      
      NOT_DEPOSITED_FLUXES(19) = 
     *     CPHYDDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(19))
      
      NOT_DEPOSITED_FLUXES(20) = 
     *     CPHYPDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(20))
      
      NOT_DEPOSITED_FLUXES(21) = 
     *     ZOOPDDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(21))
      
      NOT_DEPOSITED_FLUXES(22) = 
     *     ZOOPPDETC_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(22))

      NOT_DEPOSITED_FLUXES(23) = 
     *     SED_DOC_FLUX   * (1.0D0 - FRACTION_OF_DEPOSITION(23))

      NOT_DEPOSITED_FLUXES(24) = 
     *     SED_DON_FLUX   * (1.0D0 - FRACTION_OF_DEPOSITION(24))

      NOT_DEPOSITED_FLUXES(25) = 
     *     SED_DOP_FLUX   * (1.0D0 - FRACTION_OF_DEPOSITION(25))


C     AMMONIA FLUX
      FLUXES(1) = NH3_FLUX * FRACTION_OF_DEPOSITION(1)

C     NITRATE FLUX
      FLUXES(2) = NOx_FLUX * FRACTION_OF_DEPOSITION(2)

C     DISSOLVED ORGANIC NITROGEN FLUX
      FLUXES(3) = SED_DON_FLUX * FRACTION_OF_DEPOSITION(24)

C     PARTICULATE ORGANIC NITROGEN FLUX
      FLUXES(4) = 
     *   (EXLAPDETC_FLUX * NC_EXLADDETC * FRACTION_OF_DEPOSITION(7))  +
     *   (EXREPDETC_FLUX * NC_EXREDDETC * FRACTION_OF_DEPOSITION(11)) +
     *   (PHY_FLUX  * NC * FRACTION_OF_DEPOSITION(4))         +
     *   (GPHYPDETC_FLUX * NC * FRACTION_OF_DEPOSITION(16))   +
     *   (PHY_2_FLUX * NC_2 * FRACTION_OF_DEPOSITION(12))     +
     *   (DPHYPDETC_FLUX * NC_2 * FRACTION_OF_DEPOSITION(18)) +
     *   (PHY_3_FLUX * NC_3 * FRACTION_OF_DEPOSITION(13))     +
     *   (CPHYPDETC_FLUX * NC_3 * FRACTION_OF_DEPOSITION(20)) +
     *   (ZOO_FLUX * NC_ZOOC * FRACTION_OF_DEPOSITION(9))     +
     *   (ZOOPPDETC_FLUX * NC_ZOOC * FRACTION_OF_DEPOSITION(22))
      
C     PHOSPHATE FLUX
      FLUXES(5) = OPO4_FLUX * FRACTION_OF_DEPOSITION(3)

C     DISSOLVED ORGANIC PHOSPHORUS FLUX
      FLUXES(6) = SED_DON_FLUX * FRACTION_OF_DEPOSITION(25)

C     PARTICULATE ORGANIC PHOSPHORUS FLUX
      FLUXES(7) = 
     *   (EXLAPDETC_FLUX * PC_EXLADDETC * FRACTION_OF_DEPOSITION(7))  +
     *   (EXREPDETC_FLUX * PC_EXREDDETC * FRACTION_OF_DEPOSITION(11)) +
     *   (PHY_FLUX  * PC * FRACTION_OF_DEPOSITION(4))         +
     *   (GPHYPDETC_FLUX * PC * FRACTION_OF_DEPOSITION(16))   +
     *   (PHY_2_FLUX * PC_2 * FRACTION_OF_DEPOSITION(12))     +
     *   (DPHYPDETC_FLUX * PC_2 * FRACTION_OF_DEPOSITION(18)) +
     *   (PHY_3_FLUX * PC_3 * FRACTION_OF_DEPOSITION(13))     +
     *   (CPHYPDETC_FLUX * PC_3 * FRACTION_OF_DEPOSITION(20)) +
     *   (ZOO_FLUX * PC_ZOOC * FRACTION_OF_DEPOSITION(9))     +
     *   (ZOOPPDETC_FLUX * PC_ZOOC * FRACTION_OF_DEPOSITION(22))

     
C     DISSOLVED OXYGEN
      FLUXES(8) = 0.0D0

C     DISSOLVED ORGANIC CARBON
      FLUXES(9) = SED_DOC_FLUX * FRACTION_OF_DEPOSITION(23)

C     PARTICULATE ORGANIC CARBON FLUX
      FLUXES(10) = 
     *   (EXLAPDETC_FLUX * FRACTION_OF_DEPOSITION(7))  +
     *   (EXREPDETC_FLUX * FRACTION_OF_DEPOSITION(11)) +
     *   (PHY_FLUX * FRACTION_OF_DEPOSITION(4))        +
     *   (GPHYPDETC_FLUX * FRACTION_OF_DEPOSITION(16)) +
     *   (PHY_2_FLUX * FRACTION_OF_DEPOSITION(12))     +
     *   (DPHYPDETC_FLUX * FRACTION_OF_DEPOSITION(18)) +
     *   (PHY_3_FLUX * FRACTION_OF_DEPOSITION(13))     +
     *   (CPHYPDETC_FLUX * FRACTION_OF_DEPOSITION(20)) +
     *   (ZOO_FLUX * FRACTION_OF_DEPOSITION(9))        +
     *   (ZOOPPDETC_FLUX * FRACTION_OF_DEPOSITION(22))
      
C     DISSOLVED SILICON
      FLUXES(11) = ISI_FLUX * FRACTION_OF_DEPOSITION(10)

C     PARTICULATE  SILICON FLUX
      FLUXES(12) = 
     *   (EXLAPDETC_FLUX * SiC_EXLADDETC * FRACTION_OF_DEPOSITION(7))  +
     *   (EXREPDETC_FLUX * SiC_EXREDDETC * FRACTION_OF_DEPOSITION(11)) +
     *   (PHY_2_FLUX * ASC * FRACTION_OF_DEPOSITION(12))     +
     *   (DPHYPDETC_FLUX * ASC * FRACTION_OF_DEPOSITION(18)) +
     *   (ZOO_FLUX * SiC_ZOOC * FRACTION_OF_DEPOSITION(9))     +
     *   (ZOOPPDETC_FLUX * SiC_ZOOC * FRACTION_OF_DEPOSITION(22))

      
      END
      

      
      
c*************************************************************************
c*************************************************************************


C     THIS IS A USER PROGRAMMED FUNCTION. IT IS THE USERS RESPONSIBILITY
C     TO PROGRAM THE EQUATIONS WHICH GIVE THE MOLECULAR DIFFUSION COEFFICIENTS
C     OF STATE VARIABLES IN WATER AS A FUNCTION OF TEMPERATURE AND SALINITY

      DOUBLE PRECISION FUNCTION SED_MOD_1_ALUKAS_MOLDI_C
     *                          (SVARNO, T, SAL, TVAR)

      IMPLICIT NONE

C     SVARNO : Stave variable no
C     T      : Temperature in Celcisus
C     SAL    : Salinity (ppt)
C     TVAR   : Generic temporal variable

      INTEGER SVARNO
      DOUBLE PRECISION T
      DOUBLE PRECISION SAL
      DOUBLE PRECISION TVAR

      DOUBLE PRECISION TS
      DOUBLE PRECISION SS
      DOUBLE PRECISION V25
      DOUBLE PRECISION ONE
      DOUBLE PRECISION VTK
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION D
      DOUBLE PRECISION TK

      TK = T + 273.16

C     NH4N
C     Boudreau 1997, Springer-Verlag
      IF (SVARNO.EQ.1) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = (9.5D0 + 4.13D-1 * T) * 1.0D-6
      END IF

C     NO3N
C     Boudreau 1997, Springer-Verlag
      IF (SVARNO.EQ.2) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = (9.5D0 + 3.88D-1 * T) * 1.0D-6
      END IF

C     Dissolved organic nitrogen
      IF (SVARNO.EQ.3) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 1.0D-6
      END IF

C     Particulate organic nitrogen
      IF (SVARNO.EQ.4) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 0.0D0
      END IF

C     PO4P
C     Boudreau 1997, Springer-Verlag
      IF (SVARNO.EQ.5) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = (2.62D0 + 1.43D-1 * T) * 1.0D-6
      END IF

C     Dissolved organic phosphorus
      IF (SVARNO.EQ.6) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 1.0D-6
      END IF

C     Particulate organic phosphorus
      IF (SVARNO.EQ.7) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 0.0D0
      END IF

C     DOXY
C     Reference : Fossing et al., 2004 (NERI Technical report)
      IF (SVARNO.EQ.8) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = (1.17D1 + (3.44D-1 * T) + 
     *              (5.05D-3 * (T ** 2.0D0))) * 1.0D-6
      END IF

C     Dissolved organic carbon
      IF (SVARNO.EQ.9) THEN
           SED_MOD_1_ALUKAS_MOLDI_C = 1.0D-6
      END IF

C     Particulate organic carbon
      IF (SVARNO.EQ.10) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 0.0D0
      END IF

C     DSi
C     From Boudreau 1997, Springer-Verlag
C     Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
C     and 36.1 ppt S., Assume that this value can be scaled by
C     the Stokes-Einstein relationship to any other temperature.
      IF (SVARNO.EQ.10) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 1.0E-05
          TS = 25.0
          SS = 36.1
          CALL SED_MOD_1_CVISC(V25, SS , TS, 1.0D0)
          CALL SED_MOD_1_CVISC(VTK, SAL, T , 1.0D0)
 
          SED_MOD_1_ALUKAS_MOLDI_C = 
     *        SED_MOD_1_ALUKAS_MOLDI_C * V25 / 298.16 * TK / VTK
 
      END IF

C     Particulate silicon
      IF (SVARNO.EQ.12) THEN
          SED_MOD_1_ALUKAS_MOLDI_C = 0.0D0
      END IF

      SED_MOD_1_ALUKAS_MOLDI_C = SED_MOD_1_ALUKAS_MOLDI_C * 1.0D-4

      END
      
c****************************************************************************
c****************************************************************************


      SUBROUTINE SED_MOD_1_CVISC(V,S,T,P)

C     Calculates the shear viscosity of water using the equation
C     given by Kukulka et al. (1987).
C     Calculated viscosity is in centipoise.
C
C     Valid for 0<T<30 and 0<S<36.

      DOUBLE PRECISION V
      DOUBLE PRECISION S
      DOUBLE PRECISION T
      DOUBLE PRECISION P

      V =  1.7910 - T * (6.144D-02 - T*(1.4510D-03 - T*1.6826D-05))
     *  - 1.5290D-04 * P + 8.3885D-08 * P * P + 2.4727D-03 * S
     *  + (6.0574D-06*P - 2.6760D-09*P*P)*T + (T * (4.8429D-05
     *  - T * (4.7172D-06 - T * 7.5986D-08))) * S
      
      END
      
c************************************************************************
c************************************************************************


      SUBROUTINE FLX_SED_MOD_1_TO_ALUKAS
     *           (FLUXES_FROM_SEDIMENT, NUM_FLUXES_FROM_SEDIMENT, 
     *            FLUXES_TO_ALUKAS    , NUM_FLUXES_TO_ALUKAS)
     
      DOUBLE PRECISION FLUXES_FROM_SEDIMENT(NUM_FLUXES_FROM_SEDIMENT)
      DOUBLE PRECISION FLUXES_TO_ALUKAS    (NUM_FLUXES_TO_ALUKAS)

C     NH3 FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(1) = FLUXES_FROM_SEDIMENT(1)

C     NOx FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(2) = FLUXES_FROM_SEDIMENT(2)

C     OPO4 FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(3) = FLUXES_FROM_SEDIMENT(5)

C     PHY FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(4) = 0.0D0

C     EXLADDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(5) = 0.0D0

C     DOO FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(6) = FLUXES_FROM_SEDIMENT(8)

C     EXLAPDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(7) = 0.0D0

C     EXREDDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(8) = 0.0D0

C     ZOO FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(9) = 0.0D0

C     ISI FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(10) = FLUXES_FROM_SEDIMENT(12)

C     EXREPDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(11) = 0.0D0

C     PHY_2 FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(12) = 0.0D0

C     PHY_3 FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(13) = 0.0D0

C     INC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(14) = 0.0D0

C     GPHYDDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(15) = 0.0D0

C     GPHYPDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(16) = 0.0D0

C     DPHYDDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(17) = 0.0D0

C     DPHYPDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(18) = 0.0D0

C     CPHYDDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(19) = 0.0D0

C     CPHYPDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(20) = 0.0D0

C     ZOOPDDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(21) = 0.0D0

C     ZOOPPDETC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(22) = 0.0D0

C     SED_DOC FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(23) = FLUXES_FROM_SEDIMENT(9)

C     SED_DON FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(24) = FLUXES_FROM_SEDIMENT(3)

C     SED_DOP FLUX TO ALUKAS
      FLUXES_TO_ALUKAS(25) = FLUXES_FROM_SEDIMENT(6)
      
      END
      
c*****************************************************************
c*****************************************************************      
