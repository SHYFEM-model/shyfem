c
c aquabc programs
c
c revision log :
c
c 14.03.2012    ggu     changed BIGNUMBER because out of range (gfortran)
c
c CONTENTS:
c
c  subroutine AQUABC
c  subroutine cur_euler - calculates state variables from rates using Euler
c  FUNCTION STRANGERS(VALUE) 

c********************************************************************
c******************************************************************** 

     
      subroutine AQUABC(nnode,layer,layer_max, 
     *                  t,dt,
     *                  vol,volold,area,
     *                  depth,vel,TEMP,WIND,AIRTMP,sal,
     *                  pH,
     *                  ITOT,FDAY, 
     *                  loads,
     *                  par,
     *                  state,nstate,
     *                  INTER,NDIAGVAR,
     *                  ibsedim,
     *                  par_sed,
     *                  state_sed,nsstate, NOSLAY,
     *                  output_sed, nsoutput,     
     *                  INTER_sed, NDIAGVAR_sed)
     

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Aquatic biochemical cycling model
C This program serves as second level interface between SHYFEM 
C and WC, BS models
c
c Consists of WC module ALUKAS and sediment module (on development)
C Calls BS module when all WC layers are processed.
C Inputs:
c   nnode,layer,layer_max, 
c   t,dt
c   vol,volold,area
c   depth,vel,TEMP,WIND,AIRTMP,sal
c   pH
c   ITOT,FDAY 
c   loads
c   par - WC constans file
c   nstate
c   INTER,NDIAGVAR
c   ibsedim
c   par_sed
c   nsstate NOSLAY
c   NDIAGVAR_sed
c 
c Outputs: state     - array of values for WC state variables
C          INTER     - values of WC processes (components of derivatives) 
C          state_sed - array of values for BS state variables
C          INTER_sed - values of BS processes (components of derivatives)
C
C Notes:
C    Incomming time and step (t,dt) should be in days  
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      implicit none
 
        include 'param.h'

	integer nsdim,nssdim
	parameter ( nsdim = 25 , nssdim = 12 )
     
        integer ibsedim         !  0 - no bottom sediments processing
                                !  1 - b. sed. procesing
   
        integer nnode           !  node number
        integer LAYER           !  WC layer number
        integer LAYER_max       !  maximum number of WC layers for the node(reactor)
        real state(nsdim)             ! WC state variable [mg/L] == [g/m**3]
        integer nstate          !  number of state variables
        integer i,j
        real t 		            !current time [day]
        real dt                       !time step [day]
        real vol,volold               !volume [m**3]
        real area, depth              !depth of box [m]
        real VEL                      !velocity [m/s]
        real TEMP                     !water temperature [C]
        real WIND,AIRTMP
        real sal                      !salinity [psu] == [per mille]
        real pH
        real ITOT,FDAY		   !light intensity ly/day, photoperiod (fraction of day)

        real cold(nsdim)             !old state variable (for diagnostic purpose)
        real loads(nsdim)	          !loading for c [g/(m**3 day)]
        real rates(nsdim)            !WC state variables derivatives [g/m**3/day]        

        real par(800)

        integer NDIAGVAR
        real INTER(nsdim,NDIAGVAR)
        


                double precision STATE_VARIABLES(nsdim)
                double precision DERIVATIVES(nsdim)   !derivatives  [g/m3/day]
                double precision MODEL_CONSTANTS(554) 
                double precision DRIVING_FUNCTIONS(12) ! for AlUKAS
                double precision PROCESSES(nsdim,NDIAGVAR)                
                double precision SAVED_OUTPUTS(10)     ! outputs typical for the nodes saved for the next step   
                !double precision PSTIME
                double precision SURFL
                double precision IA

                double precision WCKIN_SAVE_ARRAY(nkndim, nlvdim, 10)
                save WCKIN_SAVE_ARRAY

                integer iprod
                parameter (iprod=nkndim*nlvdim)
                integer WCKIN_CALLED_BEFORE(nkndim, nlvdim)                
                data WCKIN_CALLED_BEFORE /iprod*0/
                save WCKIN_CALLED_BEFORE
                integer CALLED_BEFORE

C sediments variables
        real par_sed(100)       !  BS model parmeters 
        integer nsstate         !  number of botomm sediments state variables
        Integer NOSLAY          !  number of BS layers
        integer nsoutput        !  number of BS model outputs
!        parameter(nsoutput =14)  ! Temporary. fixme
        real state_sed(NOSLAY,nssdim)         ! BS state variable, g/m3 of pore water or solids or sediments       
        real output_sed(NOSLAY,nsoutput)       ! BS model outputs (states and and auxilaries
        integer NUM_SED_OUTPUTS  !nsoutput
        DOUBLE PRECISION SED_OUTPUTS(NOSLAY, nsoutput)
              

        DOUBLE PRECISION INIT_SED_STATE_VARS(NOSLAY, nssdim) ! to pass to BS subroutine,   g/m3 of pore water or solids
        DOUBLE PRECISION FINAL_SED_STATE_VARS(NOSLAY, nssdim)! received from BS subrotine, g/m3 of pore water or solids        


        integer NUM_FLUX_RECEIVING_SED_LAYERS
        parameter(NUM_FLUX_RECEIVING_SED_LAYERS = 1) !Number of first layers that receives material from WC        
        

        double precision DRIVING_FUNCTIONS_FOR_FLUXES(3) ! Fixme

 
        double precision SETTLING_VELOCITIES(nsdim)!g/m2/day
        double precision DISSOLVED_FRACTIONS(nsdim)                
        double precision SETTLING_RATES(nsdim)     !m/day
        
        integer NUM_FLUXES_TO                 
        double precision FLUXES_TO_SEDIMENTS(nssdim) !g/m2/day
        
        DOUBLE PRECISION SED_TEMPS(NOSLAY) ! BS Layers tempertures, C        
        double precision SALS              ! Salinity for Diff. coeff., psu
        double precision TS                ! Temperature for Diffusion coeff., C 
        
        INTEGER NUM_FLUXES_FROM_SEDIMENTS
        
       	DOUBLE PRECISION FLUXES_FROM_SEDIMENTS(nsdim) !g/m2/day
     *  	                                
       	
       	DOUBLE PRECISION SURF_WATER_CONCS(nssdim)  ! WC concentrations to calculate difusion between WC and 1-st BS layer, g/m3
        
       	INTEGER NUM_SED_DRIV
       	parameter(NUM_SED_DRIV = 1)
        DOUBLE PRECISION SED_DRIVING_FUNCTIONS
     *                 (NOSLAY, NUM_SED_DRIV)
     
        DOUBLE PRECISION SED_DIFFUSIONS(NOSLAY, nssdim)
     
        DOUBLE PRECISION SED_MOD_1_ALUKAS_MOLDI_C   ! Function name for Diff. coeff. calculation 
        DOUBLE PRECISION SED_DEPTHS (NOSLAY)        ! Thicknes of BS layers, m
        DOUBLE PRECISION SED_POROSITIES(NOSLAY)     ! Sediment porosities
        DOUBLE PRECISION SED_DENSITIES(NOSLAY) ! Sed. density , kg/m3
        DOUBLE PRECISION SURF_MIXLEN                ! Surface mixing length for exchange with WC, m
        DOUBLE PRECISION SED_BURRIALS(NOSLAY)       ! Burials rates (advection) for each of BS layers, m/day
        DOUBLE PRECISION FLUXES_TO_ALUKAS(nsdim)   ! Fluxes from BS to ALUKAS, g/m2/day        
        
        
        INTEGER NUM_SED_CONSTS
        parameter(NUM_SED_CONSTS = 100)
        DOUBLE PRECISION SED_MODEL_CONSTANTS(NUM_SED_CONSTS)
        
        INTEGER SEDIMENT_TYPE ! 1 - sandy, 2 - mudy
        
        INTEGER NUM_NOT_DEPOSITED_FLUXES               
        DOUBLE PRECISION FRACTION_OF_DEPOSITION(nsdim)
     *                       
        DOUBLE PRECISION NOT_DEPOSITED_FLUXES(nsdim) !g/m2/day 
        DOUBLE PRECISION PART_MIXING_COEFFS (NOSLAY, nssdim)                       
        
        DOUBLE PRECISION PSTIME, TIME_STEP ! days
       	
       	integer NDIAGVAR_sed 
        real INTER_sed(NOSLAY,nssdim,NDIAGVAR_sed)      				!Intermediate variables for BS (processes)
        DOUBLE PRECISION PROCESSES_sed(NOSLAY,nssdim, NDIAGVAR_sed)    !The same doble precision
        
     


c For debugging purposes
         integer ix, iy, iz, k, error
         real x,y,z
         
         integer STRANGERS
         
	if( nstate .ne. nsdim ) stop 'error stop aquabc: nstate'
	if( nsstate .ne. nssdim ) stop 'error stop aquabc: nsstate'
         
C     Asignment of some array bounds past to subloutines in order not to change their name         
         NUM_FLUXES_FROM_SEDIMENTS = nstate !nstate
         NUM_FLUXES_TO=nsstate              !nsstate Number of fluxes from WC to BS
         NUM_NOT_DEPOSITED_FLUXES = nstate  !nstate 
         NUM_SED_OUTPUTS = nsoutput         !nsoutput
         

         
C     Temporary initialisation. It is needed to asign to PROCESSES_sed and oposite later! Fixme

	    !INTER_sed(:,:,:)=0.
        call set_3d_r_array(NOSLAY,nssdim,NDIAGVAR_sed,INTER_sed,0.)     
       
C    Processing water column    
       
C     Converting single precision to double precision
        PSTIME = t
        do i=1,554
         MODEL_CONSTANTS(i)=par(i)
        enddo
        do i = 1, nstate
             STATE_VARIABLES(i) = state(i)
        end do

c     Driving functions management
       
      DRIVING_FUNCTIONS(1) = sal    ! salinity
      DRIVING_FUNCTIONS(2) = temp   ! w. temp.
      DRIVING_FUNCTIONS(3) = pH     ! pH
      DRIVING_FUNCTIONS(4) = 0.     ! area that is left form intersection of two layers (not important)

C     GENERATE IA. It comes as day average and used in light 
C     limitation by subroutine cur_smith in this version

      SURFL = ITOT / 2.0D0       ! Is ITOT short wave or includes infrared? Check method!
C      SURFL = (SURFL * 8.64D4 * 0.238846) / 1.0D4 ! for the case when ITOT is given in w/m2 to convert it to joules/s
                                                   ! multiplying by 86400s and convert to langlays multiplying 
                                                   ! by 0.238846/1.0D4. Should be used by smith_alukas
      IA = SURFL                 ! valid  for 2D only

c   Should be uncommented and fixed for 3D for use in alight cell is nnode
c      DO I = TOPLAY(CELLNO), NLAYER(CELLNO)
c              PHYTCC(CELLNO, I) = 
c     *            C_TRAC(CELLNO, I, 4) + C_TRAC(CELLNO, I, 12) + 
c     *            C_TRAC(CELLNO, I, 13)
c      END DO

c
c Should be uncommented and fixed for 3D
c      IA = ALIGHT_ALUKAS(TOPLAY(CELLNO), LAYER, NLAYER(CELLNO), 
c     *                   CELLNO, MCELLS, DEPTH, PHYTCC,
c     *                   CCHLX , KE    , XKC  , SURFL)

      DRIVING_FUNCTIONS(5) = SURFL  ! light to water surface
      DRIVING_FUNCTIONS(6) = depth  ! depth of layer
      DRIVING_FUNCTIONS(7) = 1.0    ! most upper layer number, fixme for 3d    
      DRIVING_FUNCTIONS(8) = WIND
      DRIVING_FUNCTIONS(9) = 1.0    ! ice free surface, 1-no ice
      DRIVING_FUNCTIONS(10) = 0.0   ! elevation (for pressure)
      DRIVING_FUNCTIONS(11) = 0.0   ! not used more 
      DRIVING_FUNCTIONS(12) = FDAY  ! day fraction of day light


 
c
c Should be uncommented and fixed for 3D
c      IA = ALIGHT_ALUKAS(TOPLAY(CELLNO), LAYER, NLAYER(CELLNO), 
c     *                   CELLNO, MCELLS, DEPTH, PHYTCC,
c     *                   CCHLX , KE    , XKC  , SURFL)

      DO I = 1, 10
          SAVED_OUTPUTS(I) = WCKIN_SAVE_ARRAY(nnode, LAYER, I)
      END DO

      DRIVING_FUNCTIONS_FOR_FLUXES(1) = SAVED_OUTPUTS(4)
      DRIVING_FUNCTIONS_FOR_FLUXES(2) = SAVED_OUTPUTS(5)
      DRIVING_FUNCTIONS_FOR_FLUXES(3) = SAVED_OUTPUTS(6)

      CALLED_BEFORE = WCKIN_CALLED_BEFORE(nnode, LAYER)
       

      CALL ALUKAS
     *           (STATE_VARIABLES  , DERIVATIVES, nstate, 
     *            MODEL_CONSTANTS  , 554, 
     *            DRIVING_FUNCTIONS, 12, 
     *            PROCESSES        , NDIAGVAR,
     *            SAVED_OUTPUTS    , 10, 
     *            PSTIME,    nnode , LAYER, CALLED_BEFORE)

        

         do i = 1, nstate
             do j = 1, NDIAGVAR
                 INTER(i,j) = PROCESSES(i,j) 
             end do 
         end do

         DO I = 1, 10
             WCKIN_SAVE_ARRAY(nnode, LAYER, I) = SAVED_OUTPUTS(I)
         END DO


         IF (WCKIN_CALLED_BEFORE(nnode, LAYER).EQ.0) THEN
             WCKIN_CALLED_BEFORE(nnode, LAYER) = 1
         END IF

         do i = 1, nstate
             rates(i) = DERIVATIVES(i)
         end do
         
C        Preparing fluxes fro WC to BS         


         SETTLING_VELOCITIES(1)  = 0.0D+0 !units: m/day
         SETTLING_VELOCITIES(2)  = 0.0D+0
         SETTLING_VELOCITIES(3)  = 0.0D+0
         SETTLING_VELOCITIES(4)  = 3.0D-1
         SETTLING_VELOCITIES(5)  = 0.0D+0
         SETTLING_VELOCITIES(6)  = 0.0D+0
         SETTLING_VELOCITIES(7)  = 1.0D-1
         SETTLING_VELOCITIES(8)  = 0.0D+0
         SETTLING_VELOCITIES(9)  = 0.0D+0
         SETTLING_VELOCITIES(10) = 0.0D+0
         SETTLING_VELOCITIES(11) = 1.0D-1
         SETTLING_VELOCITIES(12) = 3.0D-1
         SETTLING_VELOCITIES(13) = 3.0D-1
         SETTLING_VELOCITIES(14) = 0.0D+0
         SETTLING_VELOCITIES(15) = 0.0D+0
         SETTLING_VELOCITIES(16) = 1.0D-1
         SETTLING_VELOCITIES(17) = 0.0D+0
         SETTLING_VELOCITIES(18) = 1.0D-1
         SETTLING_VELOCITIES(19) = 0.0D+0
         SETTLING_VELOCITIES(20) = 1.0D-1
         SETTLING_VELOCITIES(21) = 0.0D+0
         SETTLING_VELOCITIES(22) = 1.0D-1
         SETTLING_VELOCITIES(23) = 0.0D+0
         SETTLING_VELOCITIES(24) = 0.0D+0
         SETTLING_VELOCITIES(25) = 0.0D+0

         
         DISSOLVED_FRACTIONS(1)  = 1.0D+0
         DISSOLVED_FRACTIONS(2)  = 1.0D+0
         DISSOLVED_FRACTIONS(3)  = 1.0D+0
         DISSOLVED_FRACTIONS(4)  = 0.0D+0
         DISSOLVED_FRACTIONS(5)  = 1.0D+0
         DISSOLVED_FRACTIONS(6)  = 1.0D+0
         DISSOLVED_FRACTIONS(7)  = 0.0D+0
         DISSOLVED_FRACTIONS(8)  = 1.0D+0
         DISSOLVED_FRACTIONS(9)  = 0.0D+0
         DISSOLVED_FRACTIONS(10) = 1.0D+0
         DISSOLVED_FRACTIONS(11) = 0.0D+0
         DISSOLVED_FRACTIONS(12) = 0.0D+0
         DISSOLVED_FRACTIONS(13) = 0.0D+0
         DISSOLVED_FRACTIONS(14) = 1.0D+0
         DISSOLVED_FRACTIONS(15) = 1.0D+0
         DISSOLVED_FRACTIONS(16) = 0.0D+0
         DISSOLVED_FRACTIONS(17) = 1.0D+0
         DISSOLVED_FRACTIONS(18) = 0.0D+0
         DISSOLVED_FRACTIONS(19) = 1.0D+0
         DISSOLVED_FRACTIONS(20) = 0.0D+0
         DISSOLVED_FRACTIONS(21) = 1.0D+0
         DISSOLVED_FRACTIONS(22) = 0.0D+0
         DISSOLVED_FRACTIONS(23) = 1.0D+0
         DISSOLVED_FRACTIONS(24) = 1.0D+0
         DISSOLVED_FRACTIONS(25) = 1.0D+0


C Fluxes to sediments by settling     
        CALL FLX_ALUKAS_TO_SED_MOD_1
     *           (STATE_VARIABLES               , nstate, 
     *            MODEL_CONSTANTS               , 554, 
     *            DRIVING_FUNCTIONS_FOR_FLUXES  , 3, 
     *            SETTLING_VELOCITIES, DISSOLVED_FRACTIONS, !inputs
     *            1.0D+0 , nnode, LAYER, PSTIME,
     *            SETTLING_RATES, FLUXES_TO_SEDIMENTS, NUM_FLUXES_TO,
     *            SEDIMENT_TYPE , FRACTION_OF_DEPOSITION, 
     *            NOT_DEPOSITED_FLUXES, NUM_NOT_DEPOSITED_FLUXES)
     

C settling rates comes in g/m2/day
         do i=1,nstate 
           rates(i)=rates(i) - settling_rates(i)/depth !for 3d  settling from upper layer should be added. Fixme! 
         enddo

C     Bottom sediments processing
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(ibsedim.eq.1.and.LAYER.eq.LAYER_max) then 
       
       
!        DO I = 1, NOSLAY
!         DO J = 1, nsstate      
!           if(STRANGERS(state_sed(I,J)) .eq. 1) then
!             print *, 'aquabc: Layer ', i, 'Variable ',j,
!      *            'Cell ',nnode 
!                 print *, 'Initial state is NaN'
!                 print *, 'INITIAL(i,j)=',state_sed(I,J)
!                 error =1
!                end if
!          end do
!         end do 
!



c ************************************************************************
c Temporary assignments, later should be read from files. fixme 
 
C        Fraction of deposited material (sandy sediments are assumed temporary)   
         SEDIMENT_TYPE =1
          
         if (SEDIMENT_TYPE .eq. 1) then
         
           !FRACTION_OF_DEPOSITION(:) = 0.3
           call set_1d_d_array(nsdim,FRACTION_OF_DEPOSITION,0.3D+0)     
           
          else if (SEDIMENT_TYPE .eq. 2) then
          
           !FRACTION_OF_DEPOSITION(:) = 1.
           call set_1d_d_array(nsdim,FRACTION_OF_DEPOSITION,1.D+0)     
           
          else
           print *, 'aquabc: incorrect sediment type'
           stop
          end if 

         !print *, 'layer= ',layer,'layermax= ',layer_max 
         
c       Thickness of sediment layers in meters (It might be spatialy variable)           
           SED_DEPTHS (1) = 0.02D+0
           SED_DEPTHS (2) = 0.03D+0
           SED_DEPTHS (3) = 0.05D+0
           
C       Porosities of BS layers (later should come from se. transport) Fixme!          
           SED_POROSITIES(1) = 0.40D+0
           SED_POROSITIES(2) = 0.30D+0
           SED_POROSITIES(3) = 0.25D+0
           
C       Densities of BS layers(later should come from se. transport) Fixme!      
           SED_DENSITIES(1) = 1.75D+0
           SED_DENSITIES(2) = 1.75D+0
           SED_DENSITIES(3) = 1.75D+0 
                              
C       Surface mixing length for exchange with WC(later should be in parameters) Fixme!  
           SURF_MIXLEN = 0.10  !m 0.1
           
C       Burial rate of BS state variables (around 1cm/year)                                
           SED_BURRIALS(1) =  2.7397D-5 !2.7397D-5    !m/day 
           SED_BURRIALS(2) =  2.7397D-5 !2.7397D-5
           SED_BURRIALS(3) =  2.7397D-5 !2.7397D-5 
           
!            PART_MIXING_COEFFS (1, :) = 2.64D-5 ! 0.d+0    !particle mixing diffusion coeff.  Made zero for solutes in sediment routine!
!            PART_MIXING_COEFFS (2, :) = 2.64D-5 ! 0.d+0
!            PART_MIXING_COEFFS (3, :) = 2.64D-5 ! 0.d+0 

           call set_2d_d_array(noslay,nssdim,PART_MIXING_COEFFS,2.64D-5)     

c ********************************************************************
       
C           Calculation of diffusion coefficients       
             do i=1,NOSLAY
              TS = TEMP
              SED_TEMPS(i) = TEMP ! sed. temperatures for sediments subroutine. Assumed equal to water column temp.
              SALS = sal          ! sed. salinity. Assumed equal to water column salinity 
              do j=1,nsstate
               SED_DIFFUSIONS(i,j) = SED_MOD_1_ALUKAS_MOLDI_C
     *                                   (j, TS, SALS, 0.0D+0)*86400.
              enddo
             enddo
             
C      Initial state of sediment variables to pass to BS model 
           do i=1,NOSLAY
            do j=1,nsstate                    
             INIT_SED_STATE_VARS(i,j) = state_sed(i,j)
            enddo
           enddo
            
C       WC concentrations corresponding to BS variables
                
           SURF_WATER_CONCS(1)  =  state(1)
           SURF_WATER_CONCS(2)  =  state(2)
           SURF_WATER_CONCS(3)  =  state(24)
           SURF_WATER_CONCS(4)  =  0.0D+0
           SURF_WATER_CONCS(5)  =  state(3)
           SURF_WATER_CONCS(6)  =  state(25)
           SURF_WATER_CONCS(7)  =  0.0D+0
           SURF_WATER_CONCS(8)  =  state(6)
           SURF_WATER_CONCS(9)  =  state(23)
           SURF_WATER_CONCS(10) =  0.0D+0
           SURF_WATER_CONCS(11) =  state(10)
           SURF_WATER_CONCS(12) =  0.0D+0
           
C 
           
C       BS constants           
           do i=1,NUM_SED_CONSTS
             SED_MODEL_CONSTANTS(i) = par_sed(i)
           enddo
           
C       Sediment drivig functions not used yet           
C     Call to the main sediment routine                      
        !print *, 'Calling sediment routine'
         
         PSTIME = t
         TIME_STEP = dt  
       call SEDIMENT_MODEL_1
     *           (INIT_SED_STATE_VARS  , SED_DEPTHS , SED_POROSITIES,
     *            SED_DENSITIES        ,PART_MIXING_COEFFS          ,  
     *            SED_DIFFUSIONS       , SURF_MIXLEN, SED_BURRIALS  , 
     *            SURF_WATER_CONCS     , SED_TEMPS                  ,
     *            nsstate              , NOSLAY                     , 
     *            SED_MODEL_CONSTANTS  , NUM_SED_CONSTS             , 
     *            SED_DRIVING_FUNCTIONS, NUM_SED_DRIV               ,
     *            FLUXES_TO_SEDIMENTS  , NUM_FLUXES_TO              ,
     *            1                     ,  
     *            nnode, LAYER, PSTIME, TIME_STEP                   ,                         
     *            FINAL_SED_STATE_VARS, 
     *            FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS  ,
     *            PROCESSES_sed        , NDIAGVAR_sed ,
     *            SED_OUTPUTS          , NUM_SED_OUTPUTS)
     
    ! NUM_FLUX_RECEIVING_SED_LAYERS changed by constant
     
C      Final state of BS variables and outputs to pass to to FEM interface
           do i=1,NOSLAY
            do j=1,nsstate                               
             state_sed(i,j) = FINAL_SED_STATE_VARS(i,j)
            enddo
            do j=1,nsoutput
              output_sed(i,j)= SED_OUTPUTS(i,j)
!             output_sed(i,j)= FINAL_SED_STATE_VARS(i,j)
            end do 
           enddo
                    

         
       DO I = 1, NOSLAY
        DO J = 1, nsstate      
          if(STRANGERS(state_sed(I,J)) .eq. 1) then
            print *, 'aquabc: Layer ', i, 'Variable ',j,
     *            'Cell ',nnode 
                print *, 'Initial state is NaN'
                print *, 'FINAL(i,j)=',state_sed(I,J)
                error =1
               end if
         end do
        end do

C     Updating WC model state variables rates
               
            
         CALL FLX_SED_MOD_1_TO_ALUKAS         
     *           (FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS, 
     *            FLUXES_TO_ALUKAS    , nstate)
     
     
!       print *, 't= ', t, 'nnode', nnode
!       print *, ' Fluxes from: ', (FLUXES_FROM_SEDIMENTS(j),j=1,12)


  
        
         do i=1,nstate 
           rates(i)=rates(i) + (FLUXES_TO_ALUKAS(i)/depth)
     *                       + (NOT_DEPOSITED_FLUXES(i) / depth)
         enddo
         
       end if
c     End of BS processing
ccccccccccccccccccccccccccc

c      Updating derivatives by loading
ccccccccccccccccccccccccccccccccccccccc
c       call load0d(rates,loads,vol) !Atention, check cds updates(*Vol!) fixme for alukas         

c      Updating state variables by new derivatives
cccccccccccccccccccccccccccccccccccccccccccccccccc

        INTER(23,2) = NOT_DEPOSITED_FLUXES(23)
        INTER(23,3) = rates(23)
        INTER(23,4) = (FLUXES_TO_ALUKAS(23)/depth)
        INTER(23,5) = (NOT_DEPOSITED_FLUXES(i) / depth)
        
        call cur_euler(nstate,dt,vol,volold,rates, state)


C debug
c         if(nnode.eq.2668)then
c          print *,'t=',t,'dt=',dt
c 
c          if(mod(int((t-730.)/dt),288).eq.0) then                 
c           print *,'nnode=', nnode,'t=',t,'dt=',dt
c           do i=1,nstate
c            print *,'i=',i, 'state=',state(i),'settl.rate=',
c     *                                  settling_rates(i)
c           enddo
c          endif
c         endif        
c end debug

       if( error .eq. 1) then
         print *, 'time= ', t 
         stop
       endif
       
        return
        end

c**************************************************************************
c**************************************************************************


	subroutine cur_euler(nstate,dt,vol,volold,cds,c)

c new c is computed and replaced by new values


	implicit none

	integer nstate	!number of state variables
	real dt			!time step [day]
	real vol		    !volume [m**3]
	real c(nstate)	!state variable [mg/L] == [g/m**3]	
	real cds(nstate)	!source term [g/day]
	real massnew
	integer i
	real volold	    !volume [m**3]
	real mass		    !mass [g]
	real mder		    !derivative [g/day]
		
c      volold=vol
      
	do i=1,nstate          
	  mass    = c(i) * volold
	  mder    = cds(i) * volold
	  massnew = mass + c(i)*(vol-volold) + dt*mder
	  c(i)    = massnew/vol
c       c(i) = c(i) + dt*cds(i)
	  if(c(i).lt.0.) c(i)=1.e-20
	end do	
      
	end
	

c********************************************************************
c********************************************************************
     
      INTEGER FUNCTION STRANGERS(VALUE)
      
C Cheks for NaN and Inf
C Input is single precision!

      REAL VALUE, BIGNUMBER, RATIO
      
      BIGNUMBER=1.0E30
      STRANGER=0
      
      if (.not.(VALUE .lt. 0.).and..not.(VALUE .ge. 0.)) then
       STRANGER=1
      end if
      
      if ((VALUE .lt. 0.).and.(VALUE .ge. 0.D0)) then
       STRANGER=1
      end if
      
      RATIO = BIGNUMBER/VALUE
      
      if(RATIO .eq. 0.) then
       STRANGERS=1
      endif
      return
      end
      
c************************************************************************
