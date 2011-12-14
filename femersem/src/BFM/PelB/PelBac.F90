#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBac
!
! DESCRIPTION
!   This process describes the dynamics of bacterioplankton
!    
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelBacDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: B1c, R6c, B1n, R6n, &
  ! B1p, R6p, R1c, R1n, R1p, R2c, O2o, N6r, N4n, N1p, N3n, R7c
  ! The following Pelagic 1-d global boxvars are modified : flPTN6r
  ! The following Pelagic 1-d global boxvars are used: ETW, qnB1c, qpB1c, &
  ! eO2mO2, qpR6c, qnR6c
  ! The following 0-d global box parametes are used: p_pe_R1c, p_pe_R1n, &
  ! p_pe_R1p, p_qro
  ! The following global constants are used: RLEN
  ! The following constants are used: MW_C, ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: B1c, R6c, B1n, R6n, B1p, R6p, R1c, R1n, R1p, R2c, O2o, N6r, &
    N4n, N1p, N3n, R7c
#endif
  use mem, ONLY: ppB1c, ppR6c, ppB1n, ppR6n, ppB1p, ppR6p, ppR1c, &
    ppR1n, ppR1p, ppR2c, ppO2o, ppN6r, ppN4n, ppN1p, ppN3n, ppR7c, flPTN6r, &
    ETW, qnB1c, qpB1c, eO2mO2, qpR6c, qnR6c, NO_BOXES, iiBen, iiPel, flux_vector
  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro
  use mem_PelBac


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector, insw_vector



!  
!
! !AUTHORS
!   Original version by J.W. Baretta 
!    Giovanni Coppini (UNIBO), Hanneke Baretta-Bekker, Marcello Vichi (INGV) 
!    Piet Ruardij (NIOZ) 
!
!
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)  :: runn
  real(RLEN),dimension(NO_BOXES)  :: runp
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: eO2
  real(RLEN),dimension(NO_BOXES)  :: r
  real(RLEN),dimension(NO_BOXES)  :: flB1N6r
  real(RLEN),dimension(NO_BOXES)  :: rrc
  real(RLEN),dimension(NO_BOXES)  :: rd
  real(RLEN),dimension(NO_BOXES)  :: ruR1c
  real(RLEN),dimension(NO_BOXES)  :: ruR1n
  real(RLEN),dimension(NO_BOXES)  :: ruR1p
  real(RLEN),dimension(NO_BOXES)  :: ruR2c
  real(RLEN),dimension(NO_BOXES)  :: ruR6c
  real(RLEN),dimension(NO_BOXES)  :: ruR6p
  real(RLEN),dimension(NO_BOXES)  :: ruR6n
  real(RLEN),dimension(NO_BOXES)  :: cqun3
  real(RLEN),dimension(NO_BOXES)  :: rump
  real(RLEN),dimension(NO_BOXES)  :: rumn
  real(RLEN),dimension(NO_BOXES)  :: rumn3
  real(RLEN),dimension(NO_BOXES)  :: rumn4
  real(RLEN),dimension(NO_BOXES)  :: ren
  real(RLEN),dimension(NO_BOXES)  :: rep
  real(RLEN),dimension(NO_BOXES)  :: reR2c
  real(RLEN),dimension(NO_BOXES)  :: reR7c
  real(RLEN),dimension(NO_BOXES)  :: rut
  real(RLEN),dimension(NO_BOXES)  :: rum
  real(RLEN),dimension(NO_BOXES)  :: run
  real(RLEN),dimension(NO_BOXES)  :: sun
  real(RLEN),dimension(NO_BOXES)  :: rug
  real(RLEN),dimension(NO_BOXES)  :: suR1
  real(RLEN),dimension(NO_BOXES)  :: suR1n
  real(RLEN),dimension(NO_BOXES)  :: suR1p
  real(RLEN),dimension(NO_BOXES)  :: suR2
  real(RLEN),dimension(NO_BOXES)  :: cuR6
  real(RLEN),dimension(NO_BOXES)  :: cuR1
  real(RLEN),dimension(NO_BOXES)  :: iN1p
  real(RLEN),dimension(NO_BOXES)  :: iNIn
  real(RLEN),dimension(NO_BOXES)  :: iN
  real(RLEN),dimension(NO_BOXES)  :: qpR1c
  real(RLEN),dimension(NO_BOXES)  :: qnR1c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW(:),  p_q10)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: old definition
  !   2. density dependent mortality due to virus infection
  !
  !   It is assumed the mortality is distributed in the same way over
  !   LOC (R1) and detritus (R6) s for phytoplankton and microzooplankton.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rd  =  ( p_sd* et+( p_sd2* B1c(:)))* B1c(:)

  call flux_vector( iiPel, ppB1c,ppR6c, rd*( 1.0D+00- p_pe_R1c) )
  call flux_vector( iiPel, ppB1n,ppR6n, rd* qnB1c(:)*( 1.0D+00- p_pe_R1n) )
  call flux_vector( iiPel, ppB1p,ppR6p, rd* qpB1c(:)*( 1.0D+00- p_pe_R1p) )

  call flux_vector( iiPel, ppB1c,ppR1c, rd* p_pe_R1c )
  call flux_vector( iiPel, ppB1n,ppR1n, rd* qnB1c(:)* p_pe_R1n )
  call flux_vector( iiPel, ppB1p,ppR1p, rd* qpB1c(:)* p_pe_R1p )

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate quota in R1c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  qpR1c  =   R1p(:)/ R1c(:)
  qnR1c  =   R1n(:)/ R1c(:)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version)

    case ( 1 )  !LUCA

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rum  =   p_sum* et* B1c(:)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! No correction of food avilabilities:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      cuR1  =   1.0D+00
      cuR6  =   1.0D+00

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! oxygen environment:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      eO2  =   min(  1.0D+00,  eO2mO2(:))



    case ( 2 )  ! BFM option

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      iN1p  =   min(  1.0D+00,  max(  0.0D+00,   qpB1c(:)/ p_qpc))  !Phosphorus
      iNIn  =   min(  1.0D+00,  max(  0.0D+00,   qnB1c(:)/ p_qnc))  !Nitrogen

      iN  =   min(  iN1p,  iNIn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  =   p_sum* iN* et* B1c(:)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of food avilabilities dependent on internal quota
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      cuR1  =   min(  qpR1c/ p_qpc,  qnR1c/ p_qnc)
      cuR1  =   min(  1.0D+00,  cuR1)

      cuR6  =   min(  qpR6c(:)/ p_qpc,  qnR6c(:)/ p_qnc)
      cuR6  =   min(  1.0D+00,  cuR6)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! oxygen environment:
      ! To provide a faster switching between the two metabolic pathways the
      ! oxygen dependence eO2 has been changed from the standard
      !     eO2 = MM(O2.o, p_chdo) to the cubic one written below.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      eO2  =   MM_vector(  (O2o(:))**(3.0D+00),  (p_chdo)**(3.0D+00))




  end select


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate amount for R1, R6, and R2 and total amount of substrate avilable
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ruR1c  =   p_suR1* cuR1* R1c(:)
  ruR6c  =   p_suR6* cuR6* R6c(:)
  ruR2c  =   p_suR2* R2c(:)
  rut  =   ruR6c+ ruR2c+ ruR1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rug  =   min(  rum,  rut)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR6c  =   rug* ruR6c/ rut
  ruR2c  =   rug* ruR2c/ rut
  ruR1c  =   rug* ruR1c/ rut

  call flux_vector( iiPel, ppR6c,ppB1c, ruR6c )
  call flux_vector( iiPel, ppR2c,ppB1c, ruR2c )
  call flux_vector( iiPel, ppR1c,ppB1c, ruR1c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruR6n  =   qnR6c(:)* ruR6c
  ruR1n  =   qnR1c* ruR1c

  call flux_vector( iiPel, ppR6n,ppB1n, ruR6n )
  call flux_vector( iiPel, ppR1n,ppB1n, ruR1n )

  ruR6p  =   qpR6c(:)* ruR6c
  ruR1p  =   qpR1c* ruR1c

  call flux_vector( iiPel, ppR6p,ppB1p, ruR6p )
  call flux_vector( iiPel, ppR1p,ppB1p, ruR1p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration calculation + flux
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =  ( p_pu_ra+ p_pu_ra_o*( 1.0D+00- eO2))* rug+ p_srs* B1c(:)* et
  call flux_vector( iiPel, ppB1c,ppB1c,-( rrc) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic bacteria are a wide functional group comprising both aerobic and
  ! anaerobic bacteria. At (very) low Oxygen concentrations bacteria use
  ! N6.r as electron acceptor in the respiration process. However, if N3.n
  ! is still present, the rate of consumption of N6.r is converted in N3.n
  ! consumption (see ChemicalProcesses.p).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiPel, ppO2o,ppO2o,-( eO2* rrc/ MW_C) )
  flB1N6r  =  ( 1.0D+00- eO2)* rrc/ MW_C* p_qro
  call flux_vector( iiPel, ppN6r,ppN6r, flB1N6r )
  flPTN6r(:)  =   flPTN6r(:)+ flB1N6r

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  run  =   rug- rrc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version)
    case ( 1 )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Only uptake of ammonium possible
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ren  =  ( qnB1c(:)- p_qnc)* B1c(:)* ONE_PER_DAY
      call flux_vector( iiPel, ppB1n,ppN4n, ren* insw_vector( ren) )
      call flux_vector(iiPel, ppN4n,ppB1n,- ren* insw_vector( - ren)* N4n(:)/( &
        1.0D+00+ N4n(:)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rep  =  ( qpB1c(:)- p_qpc)* B1c(:)* ONE_PER_DAY
      call flux_vector( iiPel, ppB1p,ppN1p, rep* insw_vector( rep) )
      call flux_vector(iiPel, ppN1p,ppB1p,- rep* insw_vector( - rep)* N1p(:)/( &
        0.5D+00+ N1p(:)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Activity exrecetion (defined as reR7c) + stress excetion (defined as &
      ! reR2c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      reR7c  =   p_pu_ea_R7* run

      r  =   max(  1.0D+00- qpB1c(:)/ p_qpc,  1.0D+00- qnB1c(:)/ p_qnc)
      reR2c  =   ONE_PER_DAY* r* insw_vector(  r)* B1c(:)

      run  =   run- reR7c- reR2c



    case ( 2 )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Nutrient uptake
      !
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cqun3  =   p_lN4/( p_lN4+ N4n(:))
      rumn3  =   p_qun* N3n(:)* B1c(:)* cqun3  ! max pot. uptake of N3
      rumn4  =   p_qun* N4n(:)* B1c(:)  ! max pot. uptake of N4
      rumn  =   rumn3+ rumn4

      rump  =   p_qup* N1p(:)* B1c(:)  ! max pot. uptake

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Only stress excetion (defined as reR7c) , no other excretion (reR2c=0)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r  =   min(  run, ( ruR6n+ ruR1n+ rumn)/ p_qlnc)
      reR7c  =   run- min(  r, ( ruR6p+ ruR1p+ rump)/ p_qlpc)
      reR7c  =   max(  0.0D+00,  reR7c)

      reR2c  =   0.0D+00
      run  =   run- reR7c

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! insw: No excretion if net. growth <0
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ren  =   run*(( ruR6n+ ruR1n)/ run- p_qnc)* insw_vector(  run)
      ! excess of nutrients : ren > 0
      r  =   insw_vector(  ren)
      call flux_vector( iiPel, ppB1n,ppN4n, r* ren )

      ! shortage of nutrients : ren < 0 --> Nutrient uptake
      runn  =   min( - ren,  rumn)*( 1.0D+00- r)
      call flux_vector(iiPel, ppN4n,ppB1n, runn* rumn4/( 1.0D-80+ rumn))
      call flux_vector(iiPel, ppN3n,ppB1n, runn* rumn3/( 1.0D-80+ rumn))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! insw: No excretion if net. growth <0
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rep  =   run*(( ruR6p+ ruR1p)/ run- p_qpc)* insw_vector(  run)
      ! excess of nutrients : rep > 0
      r  =   insw_vector(  rep)
      call flux_vector( iiPel, ppB1p,ppN1p, rep* r )

      ! shortage of nutrients : rep < 0 --> Nutrient uptake
      runp  =   min( - rep,  rump)*( 1.0D+00- r)
      call flux_vector( iiPel, ppN1p,ppB1p, runp )



  end select


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion fluxes + correction net prod.:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiPel, ppB1c,ppR2c, reR2c )
  call flux_vector( iiPel, ppB1c,ppR7c, reR7c )



  ! Compute section of PelagicBacteria



  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
