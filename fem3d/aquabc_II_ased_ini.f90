!==================================================================
 module aquabc_II_sed_ini
!==================================================================
!  Initialises some sediments properties and constants also
!  
!  Initialization is done using routine sed_properties_ini(bsedim)
!  called in aquabc_II_fem_interface initialisation part
!  Recalculates BS state variables
!

 
	  implicit none  
	  
	  integer, parameter, private :: noslay = 6
	  integer, parameter, private :: nstate = 24
	  
	   ! Variables used in BS model
	   double precision, save :: SED_DENSITIES      (noslay) 
       double precision, save :: SED_POROSITIES     (noslay)
       double precision, save :: SED_DEPTHS         (noslay)
       double precision, save :: SED_BURRIALS       (noslay)
       double precision, save :: PART_MIXING_COEFFS (noslay)
       double precision, save :: SURF_MIXLEN
       double precision, save :: ADVECTIVE_VELOCITY
       ! Variables used in WC and BS model     
       double precision, save :: SETTLING_VELOCITIES(nstate)
       double precision, save :: DISSOLVED_FRACTIONS(nstate)  
      
      
      real, save :: WATER_DENSITY     
      


!==================================================================
 contains
!==================================================================

 subroutine sed_properties_ini(bsedim)
       include 'aquabc_II.h'

!  Assignment of sediment properties
!  Should be read from files or obtained from 
!  sediment transport in future. Fixme
!  bsedim - .true. if BS are processed
!  Called:
!        by aquabc_II_fem_interface
!
      logical bsedim

      if(nstate .ne. 24) then
       print *, 'SED_PROPERTIES_INI:' 
       print *, 'Number of WC variables is not equal 24'
       stop 
      end if
      
! Settling velocities for WC variables
      SETTLING_VELOCITIES(1)  = 0.0D+0 !units: m/day ! AMMONIUM NITROGEN
      SETTLING_VELOCITIES(2)  = 0.0D+0               ! NITRATE NITROGEN
      SETTLING_VELOCITIES(3)  = 0.0D+0               ! ORTHOPHOSPHATE PHOSPHORUS
      SETTLING_VELOCITIES(4)  = 0.0D+0               ! DISSOLVED OXYGEN
      SETTLING_VELOCITIES(5)  = 1.0D-2               ! Nitrification bacteria carbon
      SETTLING_VELOCITIES(6)  = 1.0D-2               ! Aerobic heterotropic bacteria carbon
      SETTLING_VELOCITIES(7)  = 1.0D-2               ! Denitrification bacteria carbon
      SETTLING_VELOCITIES(8)  = 3.0D-1               ! DIATOMS CARBON
      SETTLING_VELOCITIES(9)  = 2.5D-1               ! ZOOPLANKTON CARBON
      SETTLING_VELOCITIES(10) = 0.0D+0               ! ZOOPLANKTON NITROGEN
      SETTLING_VELOCITIES(11) = 0.0D+0               ! ZOOPLANKTON PHOSPHORUS
      SETTLING_VELOCITIES(12) = 1.0D-1               ! DETRITUS PARTICULATE ORG. CARBON
      SETTLING_VELOCITIES(13) = 1.0D-1               ! DETRITUS PARTICULATE ORG. NITROGEN
      SETTLING_VELOCITIES(14) = 1.0D-1               ! DETRITUS PARTICULATE ORG. PHOSPHORUS
      SETTLING_VELOCITIES(15) = 0.0D+0               ! DISSOLVED ORGANIC CARBON
      SETTLING_VELOCITIES(16) = 0.0D+0               ! DISSOLVED ORGANIC NITROGEN
      SETTLING_VELOCITIES(17) = 0.0D+0               ! DISSOLVED ORGANIC PHOSPHORUS
      SETTLING_VELOCITIES(18) = 2.0D-1               ! CYANOBACTERIA CARBON
      SETTLING_VELOCITIES(19) = 2.0D-1               ! OTHER PHYTOPLANKTON CARBON
      SETTLING_VELOCITIES(20) = 0.0D+0               ! DISSOLOVED SILICA
      SETTLING_VELOCITIES(21) = 2.0D-1               ! BIOGENIC SILICA
      SETTLING_VELOCITIES(22) = 2.0D-1               ! Nitrogen fixing cyano-bateria
      SETTLING_VELOCITIES(23) = 0.0D+0               ! Dissolved inorganic carbon
      SETTLING_VELOCITIES(24) = 0.0D+0               ! Alkalinity
!        SETTLING_VELOCITIES(25) = 0.0D+0

! Dissolved fractions for WC variables
      DISSOLVED_FRACTIONS(1)  = 1.0D+0         ! AMMONIUM NITROGEN
      DISSOLVED_FRACTIONS(2)  = 1.0D+0         ! NITRATE NITROGEN
      DISSOLVED_FRACTIONS(3)  = 1.D+0          ! ORTHOPHOSPHATE PHOSPHORUS
      DISSOLVED_FRACTIONS(4)  = 1.0D+0         ! DISSOLVED OXYGEN
      DISSOLVED_FRACTIONS(5)  = 0.0D+0         ! Nitrification bacteria carbon
      DISSOLVED_FRACTIONS(6)  = 0.0D+0         ! Aerobic heterotropic bacteria carbon
      DISSOLVED_FRACTIONS(7)  = 0.0D+0         ! Denitrification bacteria carbon
      DISSOLVED_FRACTIONS(8)  = 0.0D+0         ! DIATOMS CARBON
      DISSOLVED_FRACTIONS(9)  = 0.0D+0         ! ZOOPLANKTON CARBON
      DISSOLVED_FRACTIONS(10) = 0.0D+0         ! ZOOPLANKTON CARBON
      DISSOLVED_FRACTIONS(11) = 0.0D+0         ! ZOOPLANKTON PHOSPHORUS
      DISSOLVED_FRACTIONS(12) = 0.0D+0         ! DETRITUS PARTICULATE ORG. CARBON
      DISSOLVED_FRACTIONS(13) = 0.0D+0         ! DETRITUS PARTICULATE ORG. NITROGEN
      DISSOLVED_FRACTIONS(14) = 0.0D+0         ! DETRITUS PARTICULATE ORG. PHOSPHORUS
      DISSOLVED_FRACTIONS(15) = 1.0D+0         ! DISSOLVED ORGANIC CARBON
      DISSOLVED_FRACTIONS(16) = 1.0D+0         ! DISSOLVED ORGANIC NITROGEN
      DISSOLVED_FRACTIONS(17) = 1.0D+0         ! DISSOLVED ORGANIC PHOSPHORUS
      DISSOLVED_FRACTIONS(18) = 0.0D+0         ! CYANOBACTERIA CARBON
      DISSOLVED_FRACTIONS(19) = 0.0D+0         ! OTHER PHYTOPLANKTON CARBON
      DISSOLVED_FRACTIONS(20) = 1.0D+0         ! DISSOLOVED SILICA
      DISSOLVED_FRACTIONS(21) = 0.0D+0         ! BIOGENIC SILICA
      DISSOLVED_FRACTIONS(22) = 0.0D+0         ! Nitrogen fixing cyano-bateria
      DISSOLVED_FRACTIONS(23) = 1.0D+0         ! Inorganic carbon (dissolved in this version)
      DISSOLVED_FRACTIONS(24) = 1.0D+0         ! Alcalinitity
!         DISSOLVED_FRACTIONS(25) = 1.0D+0

      if(.not. bsedim) return

      if(noslay .ne. 6) then
       print *, 'SED_PROPERTIES_INI:' 
       print *, 'Number of BS variables is not equal 6'
       stop 
      end if
                
!     Porosities of BS layers       
      SED_POROSITIES(1) = 0.40
      SED_POROSITIES(2) = 0.40
      SED_POROSITIES(3) = 0.40
      SED_POROSITIES(4) = 0.40
      SED_POROSITIES(5) = 0.30
      SED_POROSITIES(6) = 0.25
        
!     Densities for BS layers      
      SED_DENSITIES(1) = 1.75
      SED_DENSITIES(2) = 1.75
      SED_DENSITIES(3) = 1.75
      SED_DENSITIES(4) = 1.75
      SED_DENSITIES(5) = 1.75
      SED_DENSITIES(6) = 1.75
      
! Thickness of sediment layers in meters (It might be spatialy variable)      
      SED_DEPTHS(1) = 0.005D+0
      SED_DEPTHS(2) = 0.005D+0
      SED_DEPTHS(3) = 0.005D+0
      SED_DEPTHS(4) = 0.005D+0
      SED_DEPTHS(5) = 0.03D+0
      SED_DEPTHS(6) = 0.05D+0
      
!     Burial rate of BS state variables (around 1cm/year) !m/day      
      SED_BURRIALS(1) =  2.7397D-5
      SED_BURRIALS(2) =  2.7397D-5
      SED_BURRIALS(3) =  2.7397D-5
      SED_BURRIALS(4) =  2.7397D-5
      SED_BURRIALS(5) =  2.7397D-5
      SED_BURRIALS(6) =  2.7397D-5
      
!     Particle mixing diffusion coeff.  
!     Made zero for solutes in sediment routine!     
      PART_MIXING_COEFFS(1) = 2.64D-5
      PART_MIXING_COEFFS(2) = 2.64D-5
      PART_MIXING_COEFFS(3) = 2.64D-5
      PART_MIXING_COEFFS(4) = 2.64D-5
      PART_MIXING_COEFFS(5) = 2.64D-5
      PART_MIXING_COEFFS(6) = 2.64D-5
      
!     Surface mixing length for exchange with WC     
      SURF_MIXLEN = 0.10
      
!     advective velocity in BS  
      ADVECTIVE_VELOCITY = 0
!     
      WATER_DENSITY = 1.
      
 end subroutine sed_properties_ini      
      
!**********************************************
 subroutine sed_recalc_ini(es, sed_output)
  
! Recalculates initial NH4 and PO4 as  concentrations in porewater and solids
       
        use para_aqua
        
        implicit none
        
        include 'aquabc_II.h'
        include 'nbasin.h'
        
        real es(noslay,nkn,nsstate)  
        real sed_output(NOSLAY,nkn,nsoutput)
        real SOLUTE_FRACTIONS_NH4N(NOSLAY)
        real SOLUTE_FRACTIONS_PO4P(NOSLAY)
        real SOLID_CONCS(noslay)
        
        double precision SOLID_PART_COEFF_NH4
        double precision SOLID_PART_COEFF_PO4
        
        integer i,j,k



          call para_get_value('SOLID_PART_COEFF_NH4',SOLID_PART_COEFF_NH4)  !41
          call para_get_value('SOLID_PART_COEFF_PO4',SOLID_PART_COEFF_PO4)  !42
          
          do k=1,nkn
           do i=1,NOSLAY

            SOLID_CONCS(I) = (SED_DENSITIES(I) - &
                          (WATER_DENSITY * SED_POROSITIES(I)))

            SOLUTE_FRACTIONS_NH4N(i) =  &
              (1.0/SED_POROSITIES(i)) * &
              (1.0 / (1.0 + (SOLID_CONCS(I) * SOLID_PART_COEFF_NH4)))

            SOLUTE_FRACTIONS_PO4P(I) =  &
              (1.0/SED_POROSITIES(i)) * &
              (1.0 / (1.0 + (SOLID_CONCS(I) * SOLID_PART_COEFF_PO4)))

            sed_output(i,k,nsstate+1) = es(i,k,1)
            es(i,k,1) = es(i,k,1)/ &
                              SOLUTE_FRACTIONS_NH4N(i)
                              
            sed_output(i,k,nsstate+2) = es(i,k,5)
            
            es(i,k,5) = es(i,k,5)/ &
                              SOLUTE_FRACTIONS_PO4P(i)
            do j = 1, nsstate
              sed_output(i,k,j) = es(i,k,j)
            end do

           end do !noslay

          end do  !nkn

 end subroutine sed_recalc_ini

!************************************************

 subroutine sed_recalc_final(es)
 
! Recalculates NH4 and PO4 as solute concentrations in porewater
 
        use para_aqua
        
        implicit none
        
        include 'aquabc_II.h'
        include 'nbasin.h'
        
        real es(noslay,nkn,nsstate)  
        real sed_output(NOSLAY,nkn,nsoutput)
      
        
        double precision SOLUTE_FRACTION_NH4N
        double precision SOLUTE_FRACTION_PO4P
        double precision SOLID_CONCS(noslay)
        
        double precision SOLID_PART_COEFF_NH4
        double precision SOLID_PART_COEFF_PO4
        
        integer i,j,k


           call para_get_value('SOLID_PART_COEFF_NH4',SOLID_PART_COEFF_NH4)  !41
           call para_get_value('SOLID_PART_COEFF_PO4',SOLID_PART_COEFF_PO4)  !4

!          Use this loop if initial NH4 and PO4 have to be considered as
!          solute concentrations on porewater
            do i = 1, nkn
               do j = 1, noslay
               SOLID_CONCS(j) = (SED_DENSITIES(i) - &
                                       (WATER_DENSITY * SED_POROSITIES(j)))
                   SOLUTE_FRACTION_NH4N =      &
                     (1./ SED_POROSITIES(j)) * &
                     (1.0 / (1.0 + (SOLID_CONCS(j) * SOLID_PART_COEFF_NH4)))

                  SOLUTE_FRACTION_PO4P =      &
                     (1./ SED_POROSITIES(j))* &
                     (1.0 / (1.0 + (SOLID_CONCS(j) * SOLID_PART_COEFF_PO4)))

!                   write(*,*) j, SOLUTE_FRACTION_NH4N,
!     +              SOLUTE_FRACTION_PO4P,SED_POROSITIES(j)

!                   write(*,*) 'NH4 before : ' ,j, es(j,i,1)
                   es(j,i,1) = es(j,i,1) * &
                       SOLUTE_FRACTION_NH4N
!                         write(*,*) 'Factor : ',
!      +                   (SOLUTE_FRACTION_NH4N / SED_POROSITIES(j))
!                         write(*,*) 'NH4 after : ' ,j, es(j,i,1)

                   es(j,i,5) = es(j,i,5) *  &
                       SOLUTE_FRACTION_PO4P
               end do
            end do
            
 end subroutine sed_recalc_final

!==================================================================
 end module aquabc_II_sed_ini
!==================================================================
